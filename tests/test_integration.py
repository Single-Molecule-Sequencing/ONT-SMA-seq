#!/usr/bin/env python3
"""End-to-end integration test for the full SMA-seq pipeline with classification."""

import csv
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

NB05 = "AAGGTTACACAAACCCTGGACAAG"
NB10 = "GAGAGGACAAAGGTTTCAACGCTT"
FLANK_F = "AAGGTTAA"
FLANK_R = "CAGCACCT"
REV_FLANK_F = "AGGTGCTG"
REV_FLANK_R = "TTAACCTTAGCAAT"


def rc(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))


@pytest.fixture
def pipeline_env(tmp_path):
    """Full pipeline test environment."""
    target_seq = "ATCGATCGATCGATCGATCG" * 5

    # Build 5 reads: 3 with nb05--nb10, 2 with slightly mutated barcodes (still classifiable)
    reads = []
    for i in range(3):
        seq = FLANK_F + NB05 + FLANK_R + target_seq + REV_FLANK_F + rc(NB10) + REV_FLANK_R
        reads.append((f"read_{i:04d}", seq))
    # 2 reads with slightly mutated barcodes (still classifiable)
    for i in range(3, 5):
        mutated = list(NB05)
        mutated[0] = "G"
        seq = FLANK_F + "".join(mutated) + FLANK_R + target_seq + REV_FLANK_F + rc(NB10) + REV_FLANK_R
        reads.append((f"read_{i:04d}", seq))

    # Write BAM
    bam_path = tmp_path / "FAL99999_20260214_TEST_sup_v5.2.0_trim0_0.bam"
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for read_id, seq in reads:
            a = pysam.AlignedSegment()
            a.query_name = read_id
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.flag = 4
            bam.write(a)

    # Write pod5 summary (simulate extractMeta output)
    summary = tmp_path / "summary.tsv"
    lines = ["read_id\tend_reason\n"]
    for read_id, _ in reads:
        lines.append(f"{read_id}\tsignal_positive\n")
    summary.write_text("".join(lines))

    # Sample sheet
    ss = tmp_path / "sample_sheet.csv"
    rows = [{"flow_cell_id": "FAL99999", "kit": "SQK-NBD114-96",
             "sample_id": "20260214_TEST", "experiment_id": "exp_test",
             "barcode": "barcode05--barcode10", "alias": "test_target"}]
    with open(ss, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    # Reference
    ref_dir = tmp_path / "refs"
    ref_dir.mkdir()
    (ref_dir / "test_target.fasta").write_text(f">test_target\n{target_seq}\n")

    out_dir = tmp_path / "Output"
    return {
        "tmp": tmp_path, "bam": bam_path, "summary": summary,
        "ss": ss, "ref_dir": ref_dir, "out_dir": out_dir,
    }


class TestFullPipeline:

    def test_mkdb_then_ingest(self, pipeline_env):
        env = pipeline_env

        # Step 1: mkdb
        r1 = subprocess.run(
            [sys.executable, "bin/mkdb.py",
             "-e", "FAL99999_20260214_TEST", "-o", str(env["out_dir"])],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
        )
        assert r1.returncode == 0, r1.stderr

        db = env["out_dir"] / "SMA_FAL99999_20260214_TEST.db"
        assert db.exists()

        # Step 2: ingest with classification
        r2 = subprocess.run(
            [sys.executable, "bin/ingest.py",
             "-e", "FAL99999_20260214_TEST",
             "-b", str(env["bam"]),
             "-s", str(env["summary"]),
             "-d", str(db),
             "-o", str(env["tmp"] / "tagged.bam"),
             "-ss", str(env["ss"]),
             "-rd", str(env["ref_dir"])],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
            env={**os.environ, "PYTHONPATH": "bin"},
        )
        assert r2.returncode == 0, f"{r2.stderr}\n{r2.stdout}"

        # Verify all 5 reads ingested
        conn = sqlite3.connect(db)
        count = conn.execute("SELECT COUNT(*) FROM Reads").fetchone()[0]
        assert count == 5

        # Verify all classified as nb05 start
        starts = conn.execute(
            "SELECT DISTINCT bc_start_id FROM Reads"
        ).fetchall()
        assert starts == [("nb05",)]

        # Verify all classified as nb10 end
        ends = conn.execute(
            "SELECT DISTINCT bc_end_id FROM Reads"
        ).fetchall()
        assert ends == [("nb10",)]

        # Verify target assigned
        targets = conn.execute(
            "SELECT DISTINCT tgt_id FROM Reads"
        ).fetchall()
        assert targets == [("test_target",)]

        # Verify confidence for mutated reads is lower but still decent
        confs = conn.execute(
            "SELECT read_id, bc_start_conf FROM Reads ORDER BY read_id"
        ).fetchall()
        # First 3: perfect match -> conf ~1.0
        for read_id, conf in confs[:3]:
            assert conf >= 0.95, f"{read_id} conf={conf}"
        # Last 2: 1 mutation -> conf ~0.96
        for read_id, conf in confs[3:]:
            assert conf >= 0.9, f"{read_id} conf={conf}"

        conn.close()
