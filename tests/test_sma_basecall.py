"""Tests for sma_basecall.py BAM merging and basecalling."""

from pathlib import Path
import json
import subprocess
import sys

import pytest
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))

from sma_basecall import (
    collect_bam_files,
    build_dorado_command,
    subsample_read_ids,
    merge_sequencing_summaries,
    parse_sequencing_summary,
)


# ---------------------------------------------------------------------------
# Helper to create minimal valid BAMs
# ---------------------------------------------------------------------------


def _make_bam(path: Path, reads: list[tuple[str, str]] | None = None) -> None:
    """Create a minimal valid BAM file."""
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
    })
    with pysam.AlignmentFile(str(path), "wb", header=header) as f:
        if reads:
            for read_id, seq in reads:
                seg = pysam.AlignedSegment(header)
                seg.query_name = read_id
                seg.query_sequence = seq
                seg.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
                seg.flag = 4  # unmapped
                f.write(seg)


# ---------------------------------------------------------------------------
# collect_bam_files tests
# ---------------------------------------------------------------------------


class TestCollectBamFiles:
    def test_collects_from_demuxed_dirs(self, tmp_path):
        """Collect BAMs from demuxed subdirectories."""
        bam_dir = tmp_path / "bam_pass"
        bam_dir.mkdir()
        for sub in ["unclassified", "V04_2", "V04_4"]:
            d = bam_dir / sub
            d.mkdir()
            _make_bam(d / "reads.bam")
        bams = collect_bam_files([str(bam_dir)])
        assert len(bams) == 3

    def test_collects_from_flat_dir(self, tmp_path):
        """Collect BAMs from flat directory."""
        bam_dir = tmp_path / "bam_pass"
        bam_dir.mkdir()
        for i in range(3):
            _make_bam(bam_dir / f"reads_{i}.bam")
        bams = collect_bam_files([str(bam_dir)])
        assert len(bams) == 3

    def test_skips_empty_bams(self, tmp_path):
        """Empty BAM files should be skipped."""
        bam_dir = tmp_path / "bam_pass"
        bam_dir.mkdir()
        (bam_dir / "empty.bam").write_bytes(b"")
        _make_bam(bam_dir / "real.bam")
        bams = collect_bam_files([str(bam_dir)])
        assert len(bams) == 1

    def test_multiple_dirs(self, tmp_path):
        """Collect from multiple directories."""
        d1 = tmp_path / "run1" / "bam_pass"
        d2 = tmp_path / "run2" / "bam_pass"
        d1.mkdir(parents=True)
        d2.mkdir(parents=True)
        _make_bam(d1 / "r1.bam")
        _make_bam(d2 / "r2.bam")
        bams = collect_bam_files([str(d1), str(d2)])
        assert len(bams) == 2

    def test_nonexistent_dir(self):
        """Non-existent directory returns empty."""
        bams = collect_bam_files(["/no/such/dir"])
        assert bams == []


# ---------------------------------------------------------------------------
# build_dorado_command tests
# ---------------------------------------------------------------------------


class TestBuildDoradoCommand:
    def test_basic_command(self):
        cmd = build_dorado_command(
            model="sup",
            pod5_dir="/data/pod5",
            dorado_bin="dorado",
        )
        assert cmd == [
            "dorado", "basecaller", "sup", "/data/pod5",
            "--no-trim", "--emit-moves",
        ]

    def test_custom_model(self):
        cmd = build_dorado_command(
            model="dna_r10.4.1_e8.2_400bps_hac@v5.0.0",
            pod5_dir="/data/pod5",
            dorado_bin="/home/user/dorado/bin/dorado",
        )
        assert cmd[0] == "/home/user/dorado/bin/dorado"
        assert "dna_r10.4.1_e8.2_400bps_hac@v5.0.0" in cmd

    def test_subsample_adds_read_ids(self, tmp_path):
        ids_file = tmp_path / "ids.txt"
        ids_file.write_text("read1\nread2\n")
        cmd = build_dorado_command(
            model="sup", pod5_dir="/data/pod5",
            dorado_bin="dorado",
            read_ids_file=str(ids_file),
        )
        assert "--read-ids" in cmd
        assert str(ids_file) in cmd


# ---------------------------------------------------------------------------
# subsample_read_ids tests
# ---------------------------------------------------------------------------


class TestSubsampleReadIds:
    def test_subsamples_correct_count(self):
        all_ids = [f"read_{i}" for i in range(1000)]
        sampled = subsample_read_ids(all_ids, n=100, seed=42)
        assert len(sampled) == 100
        assert len(set(sampled)) == 100  # unique

    def test_returns_all_when_n_exceeds_total(self):
        all_ids = [f"read_{i}" for i in range(50)]
        sampled = subsample_read_ids(all_ids, n=100, seed=42)
        assert len(sampled) == 50

    def test_deterministic_with_seed(self):
        all_ids = [f"read_{i}" for i in range(1000)]
        s1 = subsample_read_ids(all_ids, n=100, seed=42)
        s2 = subsample_read_ids(all_ids, n=100, seed=42)
        assert s1 == s2


# ---------------------------------------------------------------------------
# merge_sequencing_summaries tests
# ---------------------------------------------------------------------------


class TestMergeSequencingSummaries:
    def test_merges_two_files(self, tmp_path):
        header = "read_id\trun_id\tend_reason\tduration\n"
        (tmp_path / "s1.txt").write_text(header + "r1\trun1\tsignal_positive\t1.5\n")
        (tmp_path / "s2.txt").write_text(header + "r2\trun2\tunblock\t0.3\n")
        out = tmp_path / "merged.tsv"
        n = merge_sequencing_summaries(
            [str(tmp_path / "s1.txt"), str(tmp_path / "s2.txt")],
            str(out),
        )
        assert n == 2
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 3  # header + 2 data lines
        assert "r1" in lines[1]
        assert "r2" in lines[2]

    def test_single_file_copies(self, tmp_path):
        header = "read_id\trun_id\tend_reason\n"
        (tmp_path / "s1.txt").write_text(header + "r1\trun1\tsignal_positive\n")
        out = tmp_path / "merged.tsv"
        merge_sequencing_summaries([str(tmp_path / "s1.txt")], str(out))
        assert out.exists()
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 2


# ---------------------------------------------------------------------------
# parse_sequencing_summary tests
# ---------------------------------------------------------------------------


class TestParseSequencingSummary:
    def test_parses_end_reason_and_duration(self, tmp_path):
        ss = tmp_path / "summary.tsv"
        ss.write_text(
            "read_id\tend_reason\tduration\tmean_qscore_template\n"
            "r1\tsignal_positive\t1.5\t14.2\n"
            "r2\tdata_service_unblock_mux_change\t0.3\t10.1\n"
        )
        result = parse_sequencing_summary(str(ss))
        assert len(result) == 2
        assert result["r1"]["end_reason"] == "signal_positive"
        assert result["r1"]["duration"] == pytest.approx(1.5)
        assert result["r1"]["mean_qscore"] == pytest.approx(14.2)
        assert result["r2"]["end_reason"] == "data_service_unblock_mux_change"


# ---------------------------------------------------------------------------
# CLI integration test (--from-bam)
# ---------------------------------------------------------------------------


class TestBasecallCLI:
    def test_from_bam_merges_and_indexes(self, tmp_path):
        """Integration: --from-bam collects, merges, sorts, indexes."""
        bam_dir = tmp_path / "run1" / "bam_pass"
        bam_dir.mkdir(parents=True)

        # Create 2 small BAMs with reads
        _make_bam(bam_dir / "reads_0.bam", [("read_0", "ACGT" * 25)])
        _make_bam(bam_dir / "reads_1.bam", [("read_1", "TGCA" * 30)])

        manifest = {
            "experiment_id": "test_exp",
            "experiment_dir": str(tmp_path),
            "flow_cells": [
                {"runs": [{"run_dir": str(tmp_path / "run1"), "bam_dir": "bam_pass"}]}
            ],
            "action": {
                "pod5_sources": [],
                "sequencing_summary_paths": [],
            },
        }
        manifest_path = tmp_path / "manifest.json"
        manifest_path.write_text(json.dumps(manifest))

        out_dir = tmp_path / "output"
        result = subprocess.run(
            [sys.executable, "bin/sma_basecall.py",
             str(manifest_path), "--from-bam", "-o", str(out_dir)],
            capture_output=True, text=True,
            cwd="/tmp/ont-sma-seq",
        )
        assert result.returncode == 0, f"stderr: {result.stderr}\nstdout: {result.stdout}"
        merged_bam = out_dir / "test_exp_merged.bam"
        assert merged_bam.exists()
        assert (out_dir / "test_exp_merged.bam.bai").exists()

        # Verify read count
        with pysam.AlignmentFile(str(merged_bam), check_sq=False) as f:
            reads = list(f)
        assert len(reads) == 2

        # Verify provenance
        prov = json.loads((out_dir / "test_exp_basecall.json").read_text())
        assert prov["mode"] == "from_bam"
        assert prov["total_reads"] == 2
