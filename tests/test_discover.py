"""Tests for MinKNOW output directory discovery."""
from __future__ import annotations
from pathlib import Path
import pytest
from calibrate.discover import discover_runs, parse_final_summary, RunInfo


@pytest.fixture
def minknow_tree(tmp_path: Path) -> Path:
    """Create a minimal MinKNOW output directory tree with two runs."""
    root = tmp_path / "experiment"
    root.mkdir()

    # Run 1
    run1 = root / "no_sample_id" / "20251228_2219_MD-100098_FBD69411_34fa833d"
    run1.mkdir(parents=True)
    (run1 / "final_summary_FBD69411_34fa833d_abcd1234.txt").write_text(
        "protocol_run_id=run-uuid-1\nacquisition_run_id=acq-uuid-1\n"
        "started=2025-12-28T22:19:00Z\nacquisition_stopped=2025-12-29T03:20:00Z\n"
        "flow_cell_id=FBD69411\ndevice_id=MD-100098\n"
        "sample_id=no_sample_id\nexperiment_id=my_experiment\n"
    )
    seq_summary = run1 / "sequencing_summary_FBD69411_34fa833d_abcd1234.txt"
    seq_summary.write_text(
        "read_id\tduration\tend_reason\tsequence_length_template\tmean_qscore_template\n"
        "read-001\t1.5\tsignal_positive\t500\t12.3\n"
    )
    (run1 / "bam_pass").mkdir()

    # Run 2 — same flow cell, later start
    run2 = root / "no_sample_id" / "20251229_1055_MD-100098_FBD69411_5b5c57a9"
    run2.mkdir(parents=True)
    (run2 / "final_summary_FBD69411_5b5c57a9_efgh5678.txt").write_text(
        "protocol_run_id=run-uuid-2\nacquisition_run_id=acq-uuid-2\n"
        "started=2025-12-29T10:55:00Z\nacquisition_stopped=2025-12-29T14:00:00Z\n"
        "flow_cell_id=FBD69411\ndevice_id=MD-100098\n"
        "sample_id=no_sample_id\nexperiment_id=my_experiment\n"
    )
    (run2 / "bam_pass").mkdir()

    return root


class TestParsingFinalSummary:
    def test_parses_key_fields(self, minknow_tree):
        fs_path = next(minknow_tree.rglob("final_summary_FBD69411_34fa833d_*.txt"))
        info = parse_final_summary(fs_path)
        assert info.flow_cell_id == "FBD69411"
        assert info.device_id == "MD-100098"
        assert info.sample_id == "no_sample_id"
        assert info.experiment_id == "my_experiment"

    def test_extracts_start_time(self, minknow_tree):
        fs_path = next(minknow_tree.rglob("final_summary_FBD69411_34fa833d_*.txt"))
        info = parse_final_summary(fs_path)
        assert info.started == "2025-12-28T22:19:00Z"

    def test_finds_sequencing_summary(self, minknow_tree):
        fs_path = next(minknow_tree.rglob("final_summary_FBD69411_34fa833d_*.txt"))
        info = parse_final_summary(fs_path)
        assert info.sequencing_summary_path is not None
        assert info.sequencing_summary_path.name.startswith("sequencing_summary_")

    def test_no_sample_sheet_is_none(self, minknow_tree):
        fs_path = next(minknow_tree.rglob("final_summary_FBD69411_34fa833d_*.txt"))
        info = parse_final_summary(fs_path)
        assert info.sample_sheet_path is None


class TestDiscoverRuns:
    def test_finds_two_runs(self, minknow_tree):
        runs = discover_runs(minknow_tree)
        assert len(runs) == 2

    def test_groups_by_flow_cell(self, minknow_tree):
        runs = discover_runs(minknow_tree)
        flow_cells = {r.flow_cell_id for r in runs}
        assert flow_cells == {"FBD69411"}

    def test_returns_sorted_by_start_time(self, minknow_tree):
        runs = discover_runs(minknow_tree)
        assert runs[0].started < runs[1].started

    def test_returns_run_dirs(self, minknow_tree):
        runs = discover_runs(minknow_tree)
        for r in runs:
            assert r.run_dir.is_dir()
            assert (r.run_dir / "bam_pass").is_dir()
