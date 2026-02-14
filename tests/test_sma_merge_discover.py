"""Tests for MinKNOW run discovery."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


def _make_mock_reader(flow_cell_id="FBD66244", device_id="MD-101527",
                      protocol_group="test", model="dna_r10.4.1_e8.2_400bps_hac@v5.2.0",
                      sample_id="", run_id="abc123", sample_rate=5000):
    """Create a mock pod5 Reader that yields one read with given metadata."""
    mock_read = MagicMock()
    mock_read.run_info.context_tags = {
        "basecall_model_simplex": model,
        "basecall_models_modified": "",
    }
    mock_read.run_info.tracking_id = {
        "flow_cell_id": flow_cell_id,
        "device_id": device_id,
        "protocol_group_id": protocol_group,
        "sample_id": sample_id,
        "run_id": run_id,
    }
    mock_read.run_info.sample_rate = sample_rate

    reader = MagicMock()
    reader.__enter__ = MagicMock(return_value=reader)
    reader.__exit__ = MagicMock(return_value=False)
    reader.reads.side_effect = lambda: iter([mock_read])
    return reader


def _make_run_dir(root: Path, name: str = "20251230_1709_MD-101527_FBD66244_abc123"):
    """Create a MinKNOW-style directory with a dummy pod5 file."""
    run_dir = root / "no_sample_id" / name
    pod5_dir = run_dir / "pod5_pass" / "mixed"
    pod5_dir.mkdir(parents=True)
    (pod5_dir / "chunk_0.pod5").touch()
    return run_dir


class TestDiscoverRuns:

    def test_finds_single_run(self, tmp_path):
        _make_run_dir(tmp_path)
        mock_reader = _make_mock_reader()

        with patch("sma_merge.discover.p5.Reader", return_value=mock_reader):
            from sma_merge.discover import discover_runs
            runs = discover_runs(tmp_path)

        assert len(runs) == 1
        assert runs[0].flow_cell_id == "FBD66244"
        assert runs[0].device_id == "MD-101527"
        assert runs[0].basecall_model == "dna_r10.4.1_e8.2_400bps_hac@v5.2.0"
        assert runs[0].pod5_count == 1

    def test_finds_multiple_runs(self, tmp_path):
        _make_run_dir(tmp_path, "20251230_run1_MD_FBD66244_aaa")
        _make_run_dir(tmp_path, "20251231_run2_MD_FBD66244_bbb")
        mock_reader = _make_mock_reader()

        with patch("sma_merge.discover.p5.Reader", return_value=mock_reader):
            from sma_merge.discover import discover_runs
            runs = discover_runs(tmp_path)

        assert len(runs) == 2

    def test_skips_dirs_without_pod5(self, tmp_path):
        run_dir = tmp_path / "no_sample_id" / "run1"
        (run_dir / "pod5_pass" / "mixed").mkdir(parents=True)

        from sma_merge.discover import discover_runs
        runs = discover_runs(tmp_path)
        assert len(runs) == 0

    def test_counts_pod5_files(self, tmp_path):
        run_dir = tmp_path / "no_sample_id" / "run1"
        pod5_dir = run_dir / "pod5_pass" / "mixed"
        pod5_dir.mkdir(parents=True)
        for i in range(5):
            (pod5_dir / f"chunk_{i}.pod5").touch()

        mock_reader = _make_mock_reader()
        with patch("sma_merge.discover.p5.Reader", return_value=mock_reader):
            from sma_merge.discover import discover_runs
            runs = discover_runs(tmp_path)

        assert runs[0].pod5_count == 5


class TestFormatDiscovery:

    def test_format_table(self):
        from sma_merge.discover import format_discovery_table
        from sma_merge.models import RunInfo

        runs = [
            RunInfo(
                run_dir=Path("/tmp/run1"),
                flow_cell_id="FBD66244",
                device_id="MD-101527",
                protocol_group_id="exp_one_nick",
                basecall_model="hac@v5.2.0",
                sample_id="",
                run_id="abc",
                sample_rate=5000,
                pod5_dir=Path("/tmp/run1/pod5_pass"),
                pod5_count=437,
                mod_base_models="",
            ),
        ]
        table = format_discovery_table(runs)
        assert "FBD66244" in table
        assert "MD-101527" in table
        assert "437" in table
