"""Integration test for the full calibrate pipeline."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pysam
import pytest

from calibrate.discover import discover_runs
from calibrate.merge import group_runs, find_bam_files, merge_bams
from calibrate.signal import load_sequencing_summary


@pytest.fixture
def full_experiment(tmp_path: Path) -> Path:
    """Create a complete minimal experiment with BAMs and summaries."""
    root = tmp_path / "experiment" / "sample" / "20251228_2219_MD_FBD69411_abc123"
    root.mkdir(parents=True)

    # Final summary
    (root / "final_summary_FBD69411_abc123_xyz.txt").write_text(
        "flow_cell_id=FBD69411\ndevice_id=MD-100098\n"
        "sample_id=sample\nexperiment_id=exp1\n"
        "started=2025-12-28T22:19:00Z\nprotocol_run_id=proto-1\n"
    )

    # Sequencing summary
    (root / "sequencing_summary_FBD69411_abc123_xyz.txt").write_text(
        "read_id\tduration\tend_reason\tsequence_length_template\t"
        "mean_qscore_template\tbarcode_arrangement\n"
        "read-001\t1.5\tsignal_positive\t500\t12.3\tbarcode05\n"
        "read-002\t0.8\tsignal_positive\t300\t10.1\tbarcode10\n"
    )

    # BAM with 2 unmapped reads
    bam_dir = root / "bam_pass"
    bam_dir.mkdir()
    bam_path = bam_dir / "reads_0.bam"
    header = {"HD": {"VN": "1.6", "SO": "unknown"}}
    save = pysam.set_verbosity(0)
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_out:
        for rid, seq in [("read-001", "ACGT" * 125), ("read-002", "ACGT" * 75)]:
            a = pysam.AlignedSegment()
            a.query_name = rid
            a.query_sequence = seq
            a.flag = 4
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            bam_out.write(a)
    pysam.set_verbosity(save)

    return tmp_path / "experiment"


@pytest.fixture
def two_run_experiment(tmp_path: Path) -> Path:
    """Create experiment with two runs on same flow cell (restart scenario)."""
    for run_name, started, reads in [
        ("20251228_2219_MD_FBD69411_abc123", "2025-12-28T22:19:00Z",
         [("read-001", "ACGT" * 125), ("read-002", "ACGT" * 75)]),
        ("20251229_0800_MD_FBD69411_def456", "2025-12-29T08:00:00Z",
         [("read-003", "ACGT" * 100), ("read-004", "ACGT" * 50)]),
    ]:
        run_dir = tmp_path / "experiment" / "sample" / run_name
        run_dir.mkdir(parents=True)

        parts = run_name.split("_")
        fc_id = parts[3]
        short_id = parts[4]

        (run_dir / f"final_summary_{fc_id}_{short_id}_xyz.txt").write_text(
            f"flow_cell_id={fc_id}\ndevice_id=MD-100098\n"
            f"sample_id=sample\nexperiment_id=exp1\n"
            f"started={started}\nprotocol_run_id={short_id}\n"
        )

        bam_dir = run_dir / "bam_pass"
        bam_dir.mkdir()
        bam_path = bam_dir / "reads_0.bam"
        header = {"HD": {"VN": "1.6", "SO": "unknown"}}
        save = pysam.set_verbosity(0)
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_out:
            for rid, seq in reads:
                a = pysam.AlignedSegment()
                a.query_name = rid
                a.query_sequence = seq
                a.flag = 4
                a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
                bam_out.write(a)
        pysam.set_verbosity(save)

    return tmp_path / "experiment"


class TestFullPipeline:
    def test_discover_finds_experiment(self, full_experiment: Path):
        runs = discover_runs(full_experiment)
        assert len(runs) == 1
        assert runs[0].flow_cell_id == "FBD69411"

    def test_groups_into_single_merge_group(self, full_experiment: Path):
        runs = discover_runs(full_experiment)
        groups = group_runs(runs)
        assert len(groups) == 1

    def test_finds_bam_files(self, full_experiment: Path):
        runs = discover_runs(full_experiment)
        bams = find_bam_files(runs[0].run_dir)
        assert len(bams) == 1

    def test_loads_sequencing_summary(self, full_experiment: Path):
        runs = discover_runs(full_experiment)
        assert runs[0].sequencing_summary_path is not None
        data = load_sequencing_summary(runs[0].sequencing_summary_path)
        assert len(data) == 2
        assert data["read-001"]["duration"] == pytest.approx(1.5)
        assert data["read-001"]["end_reason"] == "signal_positive"

    def test_merge_creates_merged_bam(self, full_experiment: Path, tmp_path: Path):
        runs = discover_runs(full_experiment)
        groups = group_runs(runs)
        merged_path = tmp_path / "merged.bam"
        provenance = merge_bams(groups[0], merged_path)
        assert merged_path.exists()
        assert provenance["num_bams"] == 1
        # Verify read count in merged BAM
        save = pysam.set_verbosity(0)
        with pysam.AlignmentFile(str(merged_path), "rb", check_sq=False) as bam:
            read_count = sum(1 for _ in bam)
        pysam.set_verbosity(save)
        assert read_count == 2


class TestTwoRunMerge:
    def test_discovers_two_runs(self, two_run_experiment: Path):
        runs = discover_runs(two_run_experiment)
        assert len(runs) == 2

    def test_groups_into_one_merge_group(self, two_run_experiment: Path):
        """Two runs on same flow cell within 24h should merge."""
        runs = discover_runs(two_run_experiment)
        groups = group_runs(runs)
        assert len(groups) == 1
        assert len(groups[0].runs) == 2

    def test_merge_combines_all_reads(self, two_run_experiment: Path, tmp_path: Path):
        runs = discover_runs(two_run_experiment)
        groups = group_runs(runs)
        merged_path = tmp_path / "merged.bam"
        provenance = merge_bams(groups[0], merged_path)
        assert provenance["num_bams"] == 2
        save = pysam.set_verbosity(0)
        with pysam.AlignmentFile(str(merged_path), "rb", check_sq=False) as bam:
            read_count = sum(1 for _ in bam)
        pysam.set_verbosity(save)
        assert read_count == 4


class TestVisualizationImport:
    def test_app_imports(self):
        """Calibration viz app must import cleanly."""
        from calibrate_viz.app import app
        assert app is not None

    def test_distributions_import(self):
        """Distribution computation must import cleanly."""
        from calibrate_viz.distributions import compute_kde
        assert compute_kde is not None

    def test_confusion_import(self):
        from calibrate_viz.confusion import compute_confusion_matrix
        assert compute_confusion_matrix is not None

    def test_separation_import(self):
        from calibrate_viz.separation import compute_pairwise_distances
        assert compute_pairwise_distances is not None

    def test_thresholds_import(self):
        from calibrate_viz.thresholds import compute_roc
        assert compute_roc is not None

    def test_comparison_import(self):
        from calibrate_viz.comparison import compute_experiment_summary
        assert compute_experiment_summary is not None
