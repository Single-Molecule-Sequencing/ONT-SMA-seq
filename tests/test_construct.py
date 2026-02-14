"""Tests for the construct module.

Tests cover:
- TestParseConstructToml: parsing dual mode, start_only mode, arrangement flanks,
  scoring defaults, confidence defaults, truncation defaults, barcode_pair_to_alias
  mapping, used_barcodes set, flank_sequences
- TestConstructValidation: rejects missing arrangement, missing sma, invalid mode,
  dual target without barcode2, invalid barcode ID, non-DNA flanks
"""

from __future__ import annotations

from pathlib import Path

import pytest

from construct import ConstructConfig, ValidationError, parse_construct_toml


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

DUAL_TOML = """\
[arrangement]
name = "SMA_CYP2D6_dual"
kit = "SQK-NBD114-96"
mask1_front = "AAGGTTAA"
mask1_rear  = "CAGCACCT"
mask2_front = "AGGTGCTG"
mask2_rear  = "TTAACCTTAGCAAT"
barcode1_pattern = "NB%02i"
barcode2_pattern = "NB%02i"
first_index = 1
last_index = 96

[scoring]
max_barcode_penalty = 11
min_barcode_penalty_dist = 3
front_barcode_window = 100
rear_barcode_window = 100
min_flank_score = 0.5

[sma]
mode = "dual_independent"

[sma.confidence]
full_length_threshold = 0.75
start_barcode_min = 0.6
flank_max_error_rate = 0.5

[sma.truncation]
auto_generate_refs = true
min_target_length = 20

[[sma.targets]]
barcode1 = "NB05"
barcode2 = "NB10"
alias = "CYP2D6_v04_fwd"
reference = "references/CYP2D6_v04_fwd.fasta"

[[sma.targets]]
barcode1 = "NB10"
barcode2 = "NB05"
alias = "CYP2D6_v04_rev"
reference = "references/CYP2D6_v04_rev.fasta"
"""

START_ONLY_TOML = """\
[arrangement]
name = "SMA_simple"
kit = "SQK-NBD114-96"
mask1_front = "AAGGTTAA"
mask1_rear  = "CAGCACCT"
barcode1_pattern = "NB%02i"
first_index = 1
last_index = 96

[sma]
mode = "start_only"

[[sma.targets]]
barcode1 = "NB01"
alias = "GeneX_fwd"
reference = "references/GeneX_fwd.fasta"
"""

MINIMAL_DUAL_TOML = """\
[arrangement]
name = "SMA_minimal"
kit = "SQK-NBD114-96"
mask1_front = "AAGGTTAA"
mask1_rear  = "CAGCACCT"
mask2_front = "AGGTGCTG"
mask2_rear  = "TTAACCTT"
barcode1_pattern = "NB%02i"
barcode2_pattern = "NB%02i"
first_index = 1
last_index = 96

[sma]
mode = "dual_independent"

[[sma.targets]]
barcode1 = "NB01"
barcode2 = "NB02"
alias = "target_A"
reference = "references/target_A.fasta"
"""


@pytest.fixture()
def dual_toml_path(tmp_path: Path) -> Path:
    """Write a full dual-mode construct TOML and return the path."""
    p = tmp_path / "dual_construct.toml"
    p.write_text(DUAL_TOML)
    return p


@pytest.fixture()
def start_only_toml_path(tmp_path: Path) -> Path:
    """Write a start_only-mode construct TOML and return the path."""
    p = tmp_path / "start_only_construct.toml"
    p.write_text(START_ONLY_TOML)
    return p


@pytest.fixture()
def minimal_dual_toml_path(tmp_path: Path) -> Path:
    """Write a minimal dual-mode construct TOML and return the path."""
    p = tmp_path / "minimal_dual.toml"
    p.write_text(MINIMAL_DUAL_TOML)
    return p


# ---------------------------------------------------------------------------
# TestParseConstructToml
# ---------------------------------------------------------------------------


class TestParseConstructToml:
    """Validate parse_construct_toml and ConstructConfig for valid inputs."""

    def test_parses_dual_mode(self, dual_toml_path: Path):
        """Dual-mode TOML must parse with mode='dual_independent'."""
        cfg = parse_construct_toml(dual_toml_path)
        assert cfg.mode == "dual_independent"

    def test_parses_start_only_mode(self, start_only_toml_path: Path):
        """Start-only TOML must parse with mode='start_only'."""
        cfg = parse_construct_toml(start_only_toml_path)
        assert cfg.mode == "start_only"

    def test_arrangement_flanks(self, dual_toml_path: Path):
        """Arrangement section flanking sequences must be parsed."""
        cfg = parse_construct_toml(dual_toml_path)
        assert cfg.arrangement.mask1_front == "AAGGTTAA"
        assert cfg.arrangement.mask1_rear == "CAGCACCT"
        assert cfg.arrangement.mask2_front == "AGGTGCTG"
        assert cfg.arrangement.mask2_rear == "TTAACCTTAGCAAT"

    def test_scoring_defaults(self, start_only_toml_path: Path):
        """When [scoring] is omitted, defaults must be populated."""
        cfg = parse_construct_toml(start_only_toml_path)
        assert cfg.scoring.max_barcode_penalty == 11
        assert cfg.scoring.min_barcode_penalty_dist == 3
        assert cfg.scoring.front_barcode_window == 100
        assert cfg.scoring.rear_barcode_window == 100
        assert cfg.scoring.min_flank_score == pytest.approx(0.5)

    def test_confidence_defaults(self, start_only_toml_path: Path):
        """When [sma.confidence] is omitted, defaults must be populated."""
        cfg = parse_construct_toml(start_only_toml_path)
        assert cfg.confidence.full_length_threshold == pytest.approx(0.75)
        assert cfg.confidence.start_barcode_min == pytest.approx(0.6)
        assert cfg.confidence.flank_max_error_rate == pytest.approx(0.5)

    def test_truncation_defaults(self, start_only_toml_path: Path):
        """When [sma.truncation] is omitted, defaults must be populated."""
        cfg = parse_construct_toml(start_only_toml_path)
        assert cfg.truncation.auto_generate_refs is True
        assert cfg.truncation.min_target_length == 20

    def test_barcode_pair_to_alias_dual(self, dual_toml_path: Path):
        """barcode_pair_to_alias() must return dict compatible with sample_sheet."""
        cfg = parse_construct_toml(dual_toml_path)
        mapping = cfg.barcode_pair_to_alias()
        assert mapping[("nb05", "nb10")] == "CYP2D6_v04_fwd"
        assert mapping[("nb10", "nb05")] == "CYP2D6_v04_rev"
        assert len(mapping) == 2

    def test_barcode_pair_to_alias_start_only(self, start_only_toml_path: Path):
        """start_only mode: barcode_pair_to_alias() uses barcode1 for both elements."""
        cfg = parse_construct_toml(start_only_toml_path)
        mapping = cfg.barcode_pair_to_alias()
        # In start_only mode, we use (barcode1, barcode1) as the key
        assert mapping[("nb01", "nb01")] == "GeneX_fwd"

    def test_used_barcode_ids(self, dual_toml_path: Path):
        """used_barcode_ids() must return all barcode IDs referenced in targets."""
        cfg = parse_construct_toml(dual_toml_path)
        used = cfg.used_barcode_ids()
        assert used == {"nb05", "nb10"}

    def test_flank_properties(self, dual_toml_path: Path):
        """flank_front and flank_rear must return mask1_rear and mask2_front."""
        cfg = parse_construct_toml(dual_toml_path)
        # flank_front = the trailing flank after barcode1 = mask1_rear
        assert cfg.flank_front == "CAGCACCT"
        # flank_rear = the leading flank before RC(barcode2) = mask2_front
        assert cfg.flank_rear == "AGGTGCTG"


# ---------------------------------------------------------------------------
# TestConstructValidation
# ---------------------------------------------------------------------------


class TestConstructValidation:
    """Validate that parse_construct_toml rejects invalid inputs."""

    def test_missing_arrangement(self, tmp_path: Path):
        """TOML without [arrangement] must raise ValidationError."""
        p = tmp_path / "no_arrangement.toml"
        p.write_text(
            '[sma]\nmode = "dual_independent"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\nbarcode2 = "NB02"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="arrangement"):
            parse_construct_toml(p)

    def test_missing_sma(self, tmp_path: Path):
        """TOML without [sma] must raise ValidationError."""
        p = tmp_path / "no_sma.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'barcode1_pattern = "NB%02i"\nfirst_index = 1\nlast_index = 96\n'
        )
        with pytest.raises(ValidationError, match="sma"):
            parse_construct_toml(p)

    def test_invalid_mode(self, tmp_path: Path):
        """mode must be 'dual_independent' or 'start_only'."""
        p = tmp_path / "bad_mode.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'barcode1_pattern = "NB%02i"\nfirst_index = 1\nlast_index = 96\n'
            '[sma]\nmode = "single_end"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="mode"):
            parse_construct_toml(p)

    def test_dual_target_without_barcode2(self, tmp_path: Path):
        """In dual_independent mode, targets must have barcode2."""
        p = tmp_path / "missing_bc2.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'mask2_front = "AGGTGCTG"\nmask2_rear = "TTAACCTT"\n'
            'barcode1_pattern = "NB%02i"\nbarcode2_pattern = "NB%02i"\n'
            "first_index = 1\nlast_index = 96\n"
            '[sma]\nmode = "dual_independent"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="barcode2"):
            parse_construct_toml(p)

    def test_invalid_barcode_id(self, tmp_path: Path):
        """Barcode ID outside NB01-NB96 must raise ValidationError."""
        p = tmp_path / "bad_bc.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'barcode1_pattern = "NB%02i"\nfirst_index = 1\nlast_index = 96\n'
            '[sma]\nmode = "start_only"\n'
            '[[sma.targets]]\nbarcode1 = "NB99"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="NB99"):
            parse_construct_toml(p)

    def test_non_dna_flanks(self, tmp_path: Path):
        """Flanking sequences with non-DNA characters must raise ValidationError."""
        p = tmp_path / "bad_flank.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAXGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'barcode1_pattern = "NB%02i"\nfirst_index = 1\nlast_index = 96\n'
            '[sma]\nmode = "start_only"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="DNA"):
            parse_construct_toml(p)

    def test_dual_mode_without_mask2_fails(self, tmp_path: Path):
        """dual_independent mode without mask2_front/mask2_rear must fail."""
        p = tmp_path / "dual_no_mask2.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'barcode1_pattern = "NB%02i"\nbarcode2_pattern = "NB%02i"\n'
            "first_index = 1\nlast_index = 96\n"
            '[sma]\nmode = "dual_independent"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\nbarcode2 = "NB02"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="mask2"):
            parse_construct_toml(p)

    def test_dual_mode_without_barcode2_pattern_fails(self, tmp_path: Path):
        """dual_independent mode without barcode2_pattern must fail."""
        p = tmp_path / "dual_no_bc2_pattern.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'mask2_front = "AGGTGCTG"\nmask2_rear = "TTAACCTT"\n'
            'barcode1_pattern = "NB%02i"\n'
            "first_index = 1\nlast_index = 96\n"
            '[sma]\nmode = "dual_independent"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\nbarcode2 = "NB02"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="barcode2_pattern"):
            parse_construct_toml(p)

    def test_start_only_with_mask2_fails(self, tmp_path: Path):
        """start_only mode with mask2_front or mask2_rear must fail."""
        p = tmp_path / "start_with_mask2.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'mask2_front = "AGGTGCTG"\nmask2_rear = "TTAACCTT"\n'
            'barcode1_pattern = "NB%02i"\nfirst_index = 1\nlast_index = 96\n'
            '[sma]\nmode = "start_only"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="should not have mask2"):
            parse_construct_toml(p)


# ---------------------------------------------------------------------------
# TestConstructPipelineIntegration
# ---------------------------------------------------------------------------


class TestConstructPipelineIntegration:
    """End-to-end integration test for TOML → mkrefs → ingest → DB verification."""

    def test_full_pipeline_with_construct(self, tmp_path: Path):
        """Test complete pipeline: TOML → mkrefs → synthetic reads → ingest → verify truncation levels."""
        import os
        import sqlite3
        import subprocess
        import sys

        # Import barcode sequences and reverse_complement
        sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
        from barcodes import BARCODES, reverse_complement

        # Get actual barcode sequences from barcodes.py
        NB05 = BARCODES["nb05"]  # "AAGGTTACACAAACCCTGGACAAG"
        NB10 = BARCODES["nb10"]  # "GAGAGGACAAAGGTTTCAACGCTT"

        # Define flanking sequences
        FLANK_F = "AAGGTTAA"
        FLANK_R = "CAGCACCT"
        REV_FLANK_F = "AGGTGCTG"
        REV_FLANK_R = "TTAACCTT"

        # Use poly-G target to avoid false barcode matches
        TARGET = "G" * 100  # 100bp poly-G

        # ---------------------------------------------------------------------------
        # 1. Write construct TOML
        # ---------------------------------------------------------------------------
        construct_toml = tmp_path / "construct.toml"
        construct_toml.write_text(
            f"""\
[arrangement]
name = "test_integration"
kit = "SQK-NBD114-96"
mask1_front = "{FLANK_F}"
mask1_rear = "{FLANK_R}"
mask2_front = "{REV_FLANK_F}"
mask2_rear = "{REV_FLANK_R}"
barcode1_pattern = "NB%02i"
barcode2_pattern = "NB%02i"
first_index = 1
last_index = 96

[sma]
mode = "dual_independent"

[sma.truncation]
auto_generate_refs = false
min_target_length = 20

[[sma.targets]]
barcode1 = "NB05"
barcode2 = "NB10"
alias = "target_A"
reference = "references/target_A.fasta"
"""
        )

        # ---------------------------------------------------------------------------
        # 2. Write reference FASTA for target_A
        # ---------------------------------------------------------------------------
        ref_dir = tmp_path / "references"
        ref_dir.mkdir()
        ref_fasta = ref_dir / "target_A.fasta"
        ref_fasta.write_text(f">target_A\n{TARGET}\n")

        # ---------------------------------------------------------------------------
        # 3. Run mkrefs.py to generate truncated references
        # ---------------------------------------------------------------------------
        bin_dir = Path(__file__).parent.parent / "bin"
        mkrefs_script = bin_dir / "mkrefs.py"

        result = subprocess.run(
            [
                sys.executable,
                str(mkrefs_script),
                "-c",
                str(construct_toml),
                "-o",
                str(ref_dir),
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"mkrefs.py failed: {result.stderr}"

        # ---------------------------------------------------------------------------
        # 4. Assert truncated/manifest.tsv exists
        # ---------------------------------------------------------------------------
        truncated_dir = ref_dir / "truncated"
        manifest = truncated_dir / "manifest.tsv"
        assert manifest.exists(), "truncated/manifest.tsv not created"

        # ---------------------------------------------------------------------------
        # 5. Build synthetic reads
        # ---------------------------------------------------------------------------
        # Full-length read: FLANK_F + NB05 + FLANK_R + TARGET + REV_FLANK_F + RC(NB10) + REV_FLANK_R
        full_seq = (
            FLANK_F + NB05 + FLANK_R + TARGET + REV_FLANK_F + reverse_complement(NB10) + REV_FLANK_R
        )

        # Truncated read: FLANK_F + NB05 + FLANK_R + TARGET (no bc2 or reverse flanks)
        trunc_seq = FLANK_F + NB05 + FLANK_R + TARGET

        # ---------------------------------------------------------------------------
        # 6. Write BAM and summary TSV
        # ---------------------------------------------------------------------------
        import pysam

        # BAM filename must follow format: {exp_id}_{tier}_v{ver}_{trim}_{mods}.bam
        bam_path = tmp_path / "TEST_INTEGRATION_001_hac_v5.0.0_trim0_0.bam"
        header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "target_A", "LN": len(TARGET)}]}

        with pysam.AlignmentFile(bam_path, "wb", header=header) as bam_out:
            # Full-length read
            read1 = pysam.AlignedSegment()
            read1.query_name = "read_full"
            read1.query_sequence = full_seq
            read1.reference_id = -1  # unmapped
            read1.reference_start = -1
            read1.flag = 4  # unmapped
            read1.mapping_quality = 0
            read1.query_qualities = [30] * len(full_seq)
            bam_out.write(read1)

            # Truncated read
            read2 = pysam.AlignedSegment()
            read2.query_name = "read_trunc"
            read2.query_sequence = trunc_seq
            read2.reference_id = -1
            read2.reference_start = -1
            read2.flag = 4
            read2.mapping_quality = 0
            read2.query_qualities = [30] * len(trunc_seq)
            bam_out.write(read2)

        # Index the BAM
        pysam.index(str(bam_path))

        # Write summary TSV with required columns including end_reason
        summary_tsv = tmp_path / "summary.tsv"
        summary_tsv.write_text(
            "read_id\tsequence_length_template\tmean_qscore_template\tend_reason\n"
            f"read_full\t{len(full_seq)}\t15.0\tsignal_positive\n"
            f"read_trunc\t{len(trunc_seq)}\t15.0\tdata_service_unblock_mux_change\n"
        )

        # ---------------------------------------------------------------------------
        # 7. Write sample sheet CSV
        # ---------------------------------------------------------------------------
        sample_sheet = tmp_path / "sample_sheet.csv"
        sample_sheet.write_text(
            "barcode,alias\n"
            "barcode05--barcode10,target_A\n"
        )

        # ---------------------------------------------------------------------------
        # 8. Run mkdb.py to create database
        # ---------------------------------------------------------------------------
        mkdb_script = bin_dir / "mkdb.py"

        result = subprocess.run(
            [
                sys.executable,
                str(mkdb_script),
                "-e",
                "TEST_INTEGRATION_001",
                "-o",
                str(tmp_path),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"mkdb.py failed: {result.stderr}"

        # mkdb creates SMA_<expid>.db
        db_path = tmp_path / "SMA_TEST_INTEGRATION_001.db"
        assert db_path.exists(), "Database not created by mkdb.py"

        # ---------------------------------------------------------------------------
        # 9. Run ingest.py with construct TOML, sample sheet, and reference directory
        # ---------------------------------------------------------------------------
        ingest_script = bin_dir / "ingest.py"
        output_bam = tmp_path / "tagged.bam"

        result = subprocess.run(
            [
                sys.executable,
                str(ingest_script),
                "-e",
                "TEST_INTEGRATION_001",
                "-b",
                str(bam_path),
                "-s",
                str(summary_tsv),
                "-d",
                str(db_path),
                "-o",
                str(output_bam),
                "-ss",
                str(sample_sheet),
                "-rd",
                str(ref_dir),
                "-c",
                str(construct_toml),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"ingest.py failed: {result.stderr}"

        # ---------------------------------------------------------------------------
        # 10. Open DB and verify truncation levels
        # ---------------------------------------------------------------------------
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Verify full-length read
        cursor.execute(
            "SELECT trunc_level FROM Reads WHERE read_id = ?", ("read_full",)
        )
        full_result = cursor.fetchone()
        assert full_result is not None, "Full-length read not found in database"
        assert (
            full_result[0] == "full_length"
        ), f"Expected trunc_level='full_length', got '{full_result[0]}'"

        # Verify truncated read
        cursor.execute(
            "SELECT trunc_level FROM Reads WHERE read_id = ?", ("read_trunc",)
        )
        trunc_result = cursor.fetchone()
        assert trunc_result is not None, "Truncated read not found in database"
        # Should be bc1_target or bc1_target_bc2 (but bc2 is missing, so likely bc1_target)
        assert trunc_result[0] in (
            "bc1_target",
            "bc1_target_bc2",
        ), f"Expected truncated level, got '{trunc_result[0]}'"

        conn.close()
