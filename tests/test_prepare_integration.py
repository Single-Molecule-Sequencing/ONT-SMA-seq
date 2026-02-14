"""Integration test for the full prepare pipeline."""
from __future__ import annotations

import csv
import sqlite3
from pathlib import Path

import pysam
import pytest


def _make_fake_experiment(tmp_path: Path) -> tuple[Path, Path]:
    """Create a minimal fake experiment directory with BAMs and a reference."""
    exp_dir = tmp_path / "experiment"
    run_dir = exp_dir / "no_sample" / "20260101_1200_MD-100_FC001_abc123"
    pod5_dir = run_dir / "pod5_pass"
    bam_dir = run_dir / "bam_pass" / "unclassified"
    pod5_dir.mkdir(parents=True)
    bam_dir.mkdir(parents=True)

    # Fake POD5 (just a marker file)
    (pod5_dir / "chunk.pod5").write_bytes(b"")

    # Create real BAMs
    seq_a = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32bp
    seq_b = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
    bam_path = bam_dir / "chunk1.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
        for i, seq in enumerate([seq_a, seq_b, seq_a, seq_b]):
            a = pysam.AlignedSegment(header)
            a.query_name = f"read_{i}"
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.flag = 4
            a.set_tag("er", "signal_positive", "Z")
            f.write(a)

    # Reference FASTA
    ref_path = tmp_path / "targets.fa"
    ref_path.write_text(">targetA\nACGTACGTACGTACGTACGTACGTACGTACGT\n>targetB\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n")

    return exp_dir, ref_path


def test_full_pipeline_produces_outputs(tmp_path):
    """Test that prepare produces all expected output files."""
    from prepare import Manifest, init_database, insert_run_metadata, collect_bam_files
    from prepare import symlink_pod5s
    from align import parse_fasta, process_bam
    from qc import run_qc
    from sma_merge.models import RunInfo, RunGroup

    exp_dir, ref_path = _make_fake_experiment(tmp_path)
    outdir = tmp_path / "output"
    outdir.mkdir()

    # Simulate discovery (we can't use real POD5 discovery in tests)
    run_dir = list((exp_dir / "no_sample").iterdir())[0]
    run = RunInfo(
        run_dir=run_dir, flow_cell_id="FC001", device_id="MD-100",
        protocol_group_id="test", basecall_model="sup@v5.2.0",
        sample_id="", run_id="20260101_1200_MD-100_FC001_abc123",
        sample_rate=5000, pod5_dir=run_dir / "pod5_pass",
        pod5_count=1, mod_base_models="",
    )
    group = RunGroup(flow_cell_id="FC001", runs=[run],
                     basecall_model="sup@v5.2.0", is_consistent=True)

    # Stage 3: Merge BAMs (REUSE path — just merge existing BAMs)
    bams = collect_bam_files([run])
    assert len(bams) == 1

    merged_bam = outdir / "merged.bam"
    pysam.merge("-f", str(merged_bam), *[str(b) for b in bams])

    # Stage 4: Init DB
    db_path = outdir / "SMA_TEST.db"
    init_database(db_path, "TEST_EXP", "FC001", "S001", "alias")

    # Stage 5a: Align
    refs = parse_fasta(ref_path)
    align_tsv = outdir / "alignments.tsv"
    class_tsv = outdir / "classification.tsv"
    n = process_bam(merged_bam, refs, align_tsv, class_tsv)
    assert n == 4

    # Check classification
    with open(class_tsv) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == 4
    # Reads 0,2 should match targetA, reads 1,3 should match targetB
    for row in rows:
        if row["read_id"] in ("read_0", "read_2"):
            assert row["assigned_ref"] == "targetA"
        else:
            assert row["assigned_ref"] == "targetB"

    # Stage 5b: QC
    qc_dir = outdir / "qc"
    plots = run_qc(align_tsv, class_tsv, qc_dir)
    assert len(plots) >= 10
    for p in plots:
        assert p.exists()

    # Check DB
    conn = sqlite3.connect(str(db_path))
    c = conn.cursor()
    c.execute("SELECT COUNT(*) FROM Mods")
    assert c.fetchone()[0] == 10
    conn.close()
