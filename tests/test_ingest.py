"""Tests for ingest.py barcode classification, truncation detection, and backward compatibility.

Tests cover:
- Classification mode: barcode start/end classification, target assignment,
  confidence scores, edit distances
- Truncation detection: 5-level classification when construct TOML is provided
- Backward compatibility: single-target mode without sample sheet
"""

from __future__ import annotations

import sqlite3
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from barcodes import BARCODES

# ---------------------------------------------------------------------------
# Constants for synthetic read construction
# ---------------------------------------------------------------------------

NB01 = "CACAAAGACACCGACAACTTTCTT"
NB02 = "ACAGACGACTACAAACGGAATCGA"
NB05 = "AAGGTTACACAAACCCTGGACAAG"
NB10 = "GAGAGGACAAAGGTTTCAACGCTT"

FLANK_F = "AAGGTTAA"
FLANK_R = "CAGCACCT"
REV_FLANK_F = "AGGTGCTG"
REV_FLANK_R = "TTAACCTTAGCAAT"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def rc(seq: str) -> str:
    """Simple reverse complement for test data construction."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))


def make_read_sequence(bc_up: str, bc_down: str, target_seq: str) -> str:
    """Build a full construct sequence: flank-barcode-flank-target-flank-rc(barcode)-flank."""
    return FLANK_F + bc_up + FLANK_R + target_seq + REV_FLANK_F + rc(bc_down) + REV_FLANK_R


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def test_env(tmp_path: Path):
    """Create a complete test environment with BAM, sample sheet, refs, summary, and DB.

    Reads:
      read_0000: nb05/nb10 -> target_A
      read_0001: nb01/nb02 -> target_B
      read_0002: nb05/nb10 -> target_A

    Targets:
      target_A: 100bp (ATCGATCGATCG... repeated)
      target_B: 100bp (GCTAGCTAGCTA... repeated)
    """
    target_a_seq = "ATCGATCGATCGATCGATCG" * 5  # 100bp
    target_b_seq = "GCTAGCTAGCTAGCTAGCTA" * 5  # 100bp

    reads = [
        ("read_0000", NB05, NB10, target_a_seq),
        ("read_0001", NB01, NB02, target_b_seq),
        ("read_0002", NB05, NB10, target_a_seq),
    ]

    # --- BAM ---
    bam_name = "FAL12345_20260129_IF_sup_v5.2.0_trim0_0.bam"
    bam_path = tmp_path / bam_name

    header = {"HD": {"VN": "1.6", "SO": "unknown"}}
    save_verbosity = pysam.set_verbosity(0)
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_out:
        for read_name, bc_up, bc_down, tgt_seq in reads:
            seq = make_read_sequence(bc_up, bc_down, tgt_seq)
            a = pysam.AlignedSegment()
            a.query_name = read_name
            a.query_sequence = seq
            a.flag = 4  # unmapped
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            bam_out.write(a)
    pysam.set_verbosity(save_verbosity)

    # --- Sample sheet ---
    ss_path = tmp_path / "sample_sheet.csv"
    ss_path.write_text(
        "flow_cell_id,kit,barcode,alias,type\n"
        "FC001,SQK-NBD114.96,barcode05--barcode10,target_A,test_sample\n"
        "FC001,SQK-NBD114.96,barcode01--barcode02,target_B,test_sample\n"
    )

    # --- Reference FASTAs ---
    ref_dir = tmp_path / "refs"
    ref_dir.mkdir()
    (ref_dir / "target_A.fasta").write_text(f">target_A\n{target_a_seq}\n")
    (ref_dir / "target_B.fasta").write_text(f">target_B\n{target_b_seq}\n")

    # --- Summary TSV ---
    summary_path = tmp_path / "summary.tsv"
    lines = ["read_id\tend_reason"]
    for read_name, _, _, _ in reads:
        lines.append(f"{read_name}\tsignal_positive")
    summary_path.write_text("\n".join(lines) + "\n")

    # --- Create DB via mkdb.py ---
    result = subprocess.run(
        [
            sys.executable, "bin/mkdb.py",
            "-e", "FAL12345_20260129_IF",
            "-o", str(tmp_path),
        ],
        capture_output=True, text=True,
        cwd="/tmp/ont-sma-seq",
    )
    assert result.returncode == 0, f"mkdb.py failed: {result.stderr}"
    db_path = tmp_path / "SMA_FAL12345_20260129_IF.db"

    # --- Run ingest.py in classification mode ---
    output_bam = tmp_path / "tagged.bam"
    ingest_result = subprocess.run(
        [
            sys.executable, "bin/ingest.py",
            "-e", "FAL12345_20260129_IF",
            "-b", str(bam_path),
            "-s", str(summary_path),
            "-d", str(db_path),
            "-o", str(output_bam),
            "-ss", str(ss_path),
            "-rd", str(ref_dir),
        ],
        capture_output=True, text=True,
        cwd="/tmp/ont-sma-seq",
        env={**__import__("os").environ, "PYTHONPATH": "bin"},
    )

    return {
        "tmp_path": tmp_path,
        "bam_path": bam_path,
        "ss_path": ss_path,
        "ref_dir": ref_dir,
        "summary_path": summary_path,
        "db_path": db_path,
        "output_bam": output_bam,
        "ingest_result": ingest_result,
    }


# ---------------------------------------------------------------------------
# TestIngestWithClassification
# ---------------------------------------------------------------------------


class TestIngestWithClassification:
    """Validate ingest.py classification mode with barcode demultiplexing."""

    def test_ingest_succeeds(self, test_env):
        """ingest.py should exit with returncode 0 in classification mode."""
        r = test_env["ingest_result"]
        assert r.returncode == 0, (
            f"ingest.py failed (rc={r.returncode}).\n"
            f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
        )

    def test_all_reads_ingested(self, test_env):
        """All 3 reads should be in the Reads table."""
        assert test_env["ingest_result"].returncode == 0
        conn = sqlite3.connect(test_env["db_path"])
        count = conn.execute("SELECT COUNT(*) FROM Reads").fetchone()[0]
        conn.close()
        assert count == 3

    def test_barcode_start_classified(self, test_env):
        """Start barcodes should be correctly classified for each read."""
        assert test_env["ingest_result"].returncode == 0
        conn = sqlite3.connect(test_env["db_path"])
        rows = conn.execute(
            "SELECT read_id, bc_start_id FROM Reads ORDER BY read_id"
        ).fetchall()
        conn.close()

        expected = {
            "read_0000": "nb05",
            "read_0001": "nb01",
            "read_0002": "nb05",
        }
        for read_id, bc_start_id in rows:
            assert bc_start_id == expected[read_id], (
                f"{read_id}: expected bc_start_id={expected[read_id]}, got {bc_start_id}"
            )

    def test_barcode_end_classified(self, test_env):
        """End barcodes should be correctly classified for each read."""
        assert test_env["ingest_result"].returncode == 0
        conn = sqlite3.connect(test_env["db_path"])
        rows = conn.execute(
            "SELECT read_id, bc_end_id FROM Reads ORDER BY read_id"
        ).fetchall()
        conn.close()

        expected = {
            "read_0000": "nb10",
            "read_0001": "nb02",
            "read_0002": "nb10",
        }
        for read_id, bc_end_id in rows:
            assert bc_end_id == expected[read_id], (
                f"{read_id}: expected bc_end_id={expected[read_id]}, got {bc_end_id}"
            )

    def test_target_assigned_from_barcode_pair(self, test_env):
        """Target IDs should be assigned based on barcode pair -> sample sheet alias."""
        assert test_env["ingest_result"].returncode == 0
        conn = sqlite3.connect(test_env["db_path"])
        rows = conn.execute(
            "SELECT read_id, tgt_id FROM Reads ORDER BY read_id"
        ).fetchall()
        conn.close()

        expected = {
            "read_0000": "target_A",
            "read_0001": "target_B",
            "read_0002": "target_A",
        }
        for read_id, tgt_id in rows:
            assert tgt_id == expected[read_id], (
                f"{read_id}: expected tgt_id={expected[read_id]}, got {tgt_id}"
            )

    def test_confidence_is_high(self, test_env):
        """All barcode confidences should be >= 0.9 for perfect barcode matches."""
        assert test_env["ingest_result"].returncode == 0
        conn = sqlite3.connect(test_env["db_path"])
        rows = conn.execute(
            "SELECT read_id, bc_start_conf, bc_end_conf FROM Reads"
        ).fetchall()
        conn.close()

        for read_id, start_conf, end_conf in rows:
            assert start_conf >= 0.9, (
                f"{read_id}: bc_start_conf={start_conf} < 0.9"
            )
            assert end_conf >= 0.9, (
                f"{read_id}: bc_end_conf={end_conf} < 0.9"
            )

    def test_ed_calculated(self, test_env):
        """Edit distance should be calculated (not None) for all reads."""
        assert test_env["ingest_result"].returncode == 0
        conn = sqlite3.connect(test_env["db_path"])
        rows = conn.execute(
            "SELECT read_id, ed FROM Reads"
        ).fetchall()
        conn.close()

        for read_id, ed in rows:
            assert ed is not None, f"{read_id}: ed is None"


# ---------------------------------------------------------------------------
# TestIngestBackwardCompatibility
# ---------------------------------------------------------------------------


class TestIngestBackwardCompatibility:
    """Validate that ingest.py still works in single-target mode (no sample sheet)."""

    def test_single_target_mode(self, tmp_path: Path):
        """Run ingest without -ss/-rd, just -r. tgt_id should be set, bc_start_id should be None."""
        target_seq = "ATCGATCGATCGATCGATCG" * 5  # 100bp

        # --- Single read BAM ---
        bam_name = "FAL12345_20260129_IF_sup_v5.2.0_trim0_0.bam"
        bam_path = tmp_path / bam_name

        header = {"HD": {"VN": "1.6", "SO": "unknown"}}
        save_verbosity = pysam.set_verbosity(0)
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_out:
            a = pysam.AlignedSegment()
            a.query_name = "read_single"
            a.query_sequence = target_seq
            a.flag = 4
            a.query_qualities = pysam.qualitystring_to_array("I" * len(target_seq))
            bam_out.write(a)
        pysam.set_verbosity(save_verbosity)

        # --- Reference FASTA ---
        ref_path = tmp_path / "target.fasta"
        ref_path.write_text(f">my_target\n{target_seq}\n")

        # --- Summary TSV ---
        summary_path = tmp_path / "summary.tsv"
        summary_path.write_text("read_id\tend_reason\nread_single\tsignal_positive\n")

        # --- DB ---
        result = subprocess.run(
            [
                sys.executable, "bin/mkdb.py",
                "-e", "FAL12345_20260129_IF",
                "-o", str(tmp_path),
            ],
            capture_output=True, text=True,
            cwd="/tmp/ont-sma-seq",
        )
        assert result.returncode == 0, f"mkdb.py failed: {result.stderr}"
        db_path = tmp_path / "SMA_FAL12345_20260129_IF.db"

        # --- Run ingest in single-target mode ---
        output_bam = tmp_path / "tagged.bam"
        ingest_result = subprocess.run(
            [
                sys.executable, "bin/ingest.py",
                "-e", "FAL12345_20260129_IF",
                "-b", str(bam_path),
                "-s", str(summary_path),
                "-r", str(ref_path),
                "-d", str(db_path),
                "-o", str(output_bam),
            ],
            capture_output=True, text=True,
            cwd="/tmp/ont-sma-seq",
            env={**__import__("os").environ, "PYTHONPATH": "bin"},
        )
        assert ingest_result.returncode == 0, (
            f"ingest.py single-target failed (rc={ingest_result.returncode}).\n"
            f"stdout:\n{ingest_result.stdout}\nstderr:\n{ingest_result.stderr}"
        )

        # --- Verify ---
        conn = sqlite3.connect(db_path)
        row = conn.execute(
            "SELECT tgt_id, bc_start_id FROM Reads WHERE read_id = 'read_single'"
        ).fetchone()
        conn.close()

        assert row is not None, "read_single not found in DB"
        tgt_id, bc_start_id = row
        assert tgt_id == "my_target", f"Expected tgt_id='my_target', got '{tgt_id}'"
        assert bc_start_id is None, f"Expected bc_start_id=None, got '{bc_start_id}'"


# ---------------------------------------------------------------------------
# TestTruncationClassification
# ---------------------------------------------------------------------------


def _make_construct_toml(path: Path) -> None:
    """Write a minimal construct TOML for truncation testing."""
    path.write_text(
        '[arrangement]\n'
        'name = "test_construct"\n'
        'kit = "SQK-NBD114.96"\n'
        'mask1_front = "AAGGTTAA"\n'
        'mask1_rear = "CAGCACCT"\n'
        'mask2_front = "CAGCACCT"\n'
        'mask2_rear = "TTAACCTT"\n'
        'barcode1_pattern = "NB%02i"\n'
        'barcode2_pattern = "NB%02i"\n'
        'first_index = 1\n'
        'last_index = 96\n'
        '\n'
        '[scoring]\n'
        '\n'
        '[sma]\n'
        'mode = "dual_independent"\n'
        '\n'
        '[sma.confidence]\n'
        'full_length_threshold = 0.75\n'
        'start_barcode_min = 0.6\n'
        '\n'
        '[sma.truncation]\n'
        'min_target_length = 20\n'
        '\n'
        '[[sma.targets]]\n'
        'barcode1 = "NB05"\n'
        'barcode2 = "NB10"\n'
        'alias = "target_A"\n'
        'reference = "target_A.fasta"\n'
    )


@pytest.fixture()
def construct_env(tmp_path: Path):
    """Create a test environment with construct TOML for truncation testing.

    Reads:
      read_full:   full construct (bc1 + flank + target + rc_flank + rc_bc2) -> full_length
      read_trunc:  bc1 + flank + target only (no rc_bc2)                   -> bc1_target
      read_bc1:    bc1 + short junk (too short for target)                 -> bc1_only

    Uses poly-G target to ensure the truncated read has very low end-barcode
    confidence (< 0.3), placing it in the bc1_target bucket rather than
    bc1_target_bc2.
    """
    target_a_seq = "G" * 100  # 100bp poly-G for clean separation

    # Full-length read: flank_f + bc1 + flank_r + target + rev_flank_f + rc(bc2) + rev_flank_r
    full_seq = FLANK_F + NB05 + FLANK_R + target_a_seq + REV_FLANK_F + rc(NB10) + REV_FLANK_R

    # Truncated read: has bc1 + flank + target, but NO rc_bc2 region
    trunc_seq = FLANK_F + NB05 + FLANK_R + target_a_seq

    # BC1-only read: has bc1 but too short for any target content
    bc1_seq = FLANK_F + NB05 + FLANK_R + "AAAA"

    reads = [
        ("read_full", full_seq),
        ("read_trunc", trunc_seq),
        ("read_bc1", bc1_seq),
    ]

    # --- BAM ---
    bam_name = "FAL12345_20260129_IF_sup_v5.2.0_trim0_0.bam"
    bam_path = tmp_path / bam_name

    header = {"HD": {"VN": "1.6", "SO": "unknown"}}
    save_verbosity = pysam.set_verbosity(0)
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_out:
        for read_name, seq in reads:
            a = pysam.AlignedSegment()
            a.query_name = read_name
            a.query_sequence = seq
            a.flag = 4
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            bam_out.write(a)
    pysam.set_verbosity(save_verbosity)

    # --- Construct TOML ---
    construct_path = tmp_path / "construct.toml"
    _make_construct_toml(construct_path)

    # --- Reference FASTA ---
    ref_dir = tmp_path / "refs"
    ref_dir.mkdir()
    (ref_dir / "target_A.fasta").write_text(f">target_A\n{target_a_seq}\n")

    # --- Sample sheet (nb05/nb10 -> target_A) ---
    ss_path = tmp_path / "sample_sheet.csv"
    ss_path.write_text(
        "flow_cell_id,kit,barcode,alias,type\n"
        "FC001,SQK-NBD114.96,barcode05--barcode10,target_A,test_sample\n"
    )

    # --- Summary TSV ---
    summary_path = tmp_path / "summary.tsv"
    lines = ["read_id\tend_reason"]
    for read_name, _ in reads:
        lines.append(f"{read_name}\tsignal_positive")
    summary_path.write_text("\n".join(lines) + "\n")

    # --- Create DB via mkdb.py ---
    result = subprocess.run(
        [
            sys.executable, "bin/mkdb.py",
            "-e", "FAL12345_20260129_IF",
            "-o", str(tmp_path),
        ],
        capture_output=True, text=True,
        cwd="/tmp/ont-sma-seq",
    )
    assert result.returncode == 0, f"mkdb.py failed: {result.stderr}"
    db_path = tmp_path / "SMA_FAL12345_20260129_IF.db"

    # --- Run ingest.py with -c construct TOML ---
    output_bam = tmp_path / "tagged.bam"
    ingest_result = subprocess.run(
        [
            sys.executable, "bin/ingest.py",
            "-e", "FAL12345_20260129_IF",
            "-b", str(bam_path),
            "-s", str(summary_path),
            "-d", str(db_path),
            "-o", str(output_bam),
            "-ss", str(ss_path),
            "-rd", str(ref_dir),
            "-c", str(construct_path),
        ],
        capture_output=True, text=True,
        cwd="/tmp/ont-sma-seq",
        env={**__import__("os").environ, "PYTHONPATH": "bin"},
    )

    return {
        "tmp_path": tmp_path,
        "db_path": db_path,
        "ingest_result": ingest_result,
    }


class TestTruncationClassification:
    """Validate 5-level truncation classification when construct TOML is provided."""

    def test_ingest_with_construct_succeeds(self, construct_env):
        """ingest.py should succeed with -c construct TOML flag."""
        r = construct_env["ingest_result"]
        assert r.returncode == 0, (
            f"ingest.py failed (rc={r.returncode}).\n"
            f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
        )

    def test_full_length_classified(self, construct_env):
        """Full-length read (bc1+target+rc_bc2) should get trunc_level='full_length'."""
        r = construct_env["ingest_result"]
        assert r.returncode == 0, f"ingest.py failed:\n{r.stdout}\n{r.stderr}"

        conn = sqlite3.connect(construct_env["db_path"])
        row = conn.execute(
            "SELECT trunc_level FROM Reads WHERE read_id = 'read_full'"
        ).fetchone()
        conn.close()

        assert row is not None, "read_full not found in DB"
        assert row[0] == "full_length", (
            f"Expected trunc_level='full_length', got '{row[0]}'"
        )

    def test_bc1_target_classified(self, construct_env):
        """Truncated read (bc1+target, no bc2) should get trunc_level='bc1_target'."""
        r = construct_env["ingest_result"]
        assert r.returncode == 0, f"ingest.py failed:\n{r.stdout}\n{r.stderr}"

        conn = sqlite3.connect(construct_env["db_path"])
        row = conn.execute(
            "SELECT trunc_level FROM Reads WHERE read_id = 'read_trunc'"
        ).fetchone()
        conn.close()

        assert row is not None, "read_trunc not found in DB"
        assert row[0] == "bc1_target", (
            f"Expected trunc_level='bc1_target', got '{row[0]}'"
        )

    def test_bc1_only_classified(self, construct_env):
        """Very short read (bc1 + tiny junk) should get trunc_level='bc1_only'."""
        r = construct_env["ingest_result"]
        assert r.returncode == 0, f"ingest.py failed:\n{r.stdout}\n{r.stderr}"

        conn = sqlite3.connect(construct_env["db_path"])
        row = conn.execute(
            "SELECT trunc_level FROM Reads WHERE read_id = 'read_bc1'"
        ).fetchone()
        conn.close()

        assert row is not None, "read_bc1 not found in DB"
        assert row[0] == "bc1_only", (
            f"Expected trunc_level='bc1_only', got '{row[0]}'"
        )

    def test_backward_compat_no_construct(self, tmp_path: Path):
        """Without -c flag, trunc_level should be NULL for all reads."""
        target_seq = "ATCGATCGATCGATCGATCG" * 5

        # --- BAM ---
        bam_name = "FAL12345_20260129_IF_sup_v5.2.0_trim0_0.bam"
        bam_path = tmp_path / bam_name

        header = {"HD": {"VN": "1.6", "SO": "unknown"}}
        save_verbosity = pysam.set_verbosity(0)
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_out:
            seq = make_read_sequence(NB05, NB10, target_seq)
            a = pysam.AlignedSegment()
            a.query_name = "read_noconst"
            a.query_sequence = seq
            a.flag = 4
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            bam_out.write(a)
        pysam.set_verbosity(save_verbosity)

        # --- Sample sheet ---
        ss_path = tmp_path / "sample_sheet.csv"
        ss_path.write_text(
            "flow_cell_id,kit,barcode,alias,type\n"
            "FC001,SQK-NBD114.96,barcode05--barcode10,target_A,test_sample\n"
        )

        # --- Reference FASTA ---
        ref_dir = tmp_path / "refs"
        ref_dir.mkdir()
        (ref_dir / "target_A.fasta").write_text(f">target_A\n{target_seq}\n")

        # --- Summary TSV ---
        summary_path = tmp_path / "summary.tsv"
        summary_path.write_text("read_id\tend_reason\nread_noconst\tsignal_positive\n")

        # --- DB ---
        result = subprocess.run(
            [
                sys.executable, "bin/mkdb.py",
                "-e", "FAL12345_20260129_IF",
                "-o", str(tmp_path),
            ],
            capture_output=True, text=True,
            cwd="/tmp/ont-sma-seq",
        )
        assert result.returncode == 0, f"mkdb.py failed: {result.stderr}"
        db_path = tmp_path / "SMA_FAL12345_20260129_IF.db"

        # --- Run ingest WITHOUT -c ---
        output_bam = tmp_path / "tagged.bam"
        ingest_result = subprocess.run(
            [
                sys.executable, "bin/ingest.py",
                "-e", "FAL12345_20260129_IF",
                "-b", str(bam_path),
                "-s", str(summary_path),
                "-d", str(db_path),
                "-o", str(output_bam),
                "-ss", str(ss_path),
                "-rd", str(ref_dir),
            ],
            capture_output=True, text=True,
            cwd="/tmp/ont-sma-seq",
            env={**__import__("os").environ, "PYTHONPATH": "bin"},
        )
        assert ingest_result.returncode == 0, (
            f"ingest.py failed:\n{ingest_result.stdout}\n{ingest_result.stderr}"
        )

        conn = sqlite3.connect(db_path)
        row = conn.execute(
            "SELECT trunc_level FROM Reads WHERE read_id = 'read_noconst'"
        ).fetchone()
        conn.close()

        assert row is not None, "read_noconst not found in DB"
        assert row[0] is None, (
            f"Expected trunc_level=NULL without -c flag, got '{row[0]}'"
        )
