# tests/test_sma_merge_integration.py
"""Integration tests for sma-merge with real experiment data.

These tests require real data at /tmp/one_nick_data/.
Skip if data is not available.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

ONE_NICK_DATA = Path("/tmp/one_nick_data")
HAS_REAL_DATA = ONE_NICK_DATA.exists() and any(ONE_NICK_DATA.rglob("pod5_pass"))

skip_no_data = pytest.mark.skipif(
    not HAS_REAL_DATA,
    reason="Real experiment data not available at /tmp/one_nick_data/",
)


@skip_no_data
class TestDiscoverRealData:

    def test_discovers_one_nick_run(self):
        from sma_merge.discover import discover_runs
        runs = discover_runs(ONE_NICK_DATA)
        assert len(runs) >= 1
        run = runs[0]
        assert run.flow_cell_id == "FBD66244"
        assert run.device_id == "MD-101527"
        assert "hac" in run.basecall_model
        assert run.pod5_count > 0
        assert run.sample_rate == 5000

    def test_validates_one_nick(self):
        from sma_merge.discover import discover_runs
        from sma_merge.validate import validate_runs
        runs = discover_runs(ONE_NICK_DATA)
        groups = validate_runs(runs)
        assert len(groups) == 1
        assert groups[0].is_consistent


@skip_no_data
class TestSubsampleRealData:

    def test_subsample_50_reads(self, tmp_path):
        from sma_merge.discover import discover_runs
        from sma_merge.subsample import subsample_pod5

        runs = discover_runs(ONE_NICK_DATA)
        pod5_dir = runs[0].pod5_dir

        output = tmp_path / "sub50.pod5"
        selected = subsample_pod5(pod5_dir, output, n_reads=50, seed=42)

        assert len(selected) == 50
        assert output.exists()
        assert output.stat().st_size > 0

    def test_end_reason_lookup_from_subsampled(self, tmp_path):
        from sma_merge.discover import discover_runs
        from sma_merge.subsample import subsample_pod5
        from sma_merge.tag import build_pod5_lookup

        runs = discover_runs(ONE_NICK_DATA)
        pod5_dir = runs[0].pod5_dir

        output = tmp_path / "sub50.pod5"
        selected = subsample_pod5(pod5_dir, output, n_reads=50, seed=42)

        lookup = build_pod5_lookup(output)
        assert len(lookup) == 50

        valid_reasons = {
            "unknown", "mux_change", "unblock_mux_change",
            "data_service_unblock_mux_change", "signal_positive",
            "signal_negative", "api_request", "device_data_error",
            "analysis_config_change", "paused",
        }
        for rid, (er, ns) in lookup.items():
            assert er in valid_reasons, f"Unexpected end_reason: {er}"
            assert ns > 0, f"Signal length should be > 0, got {ns}"


@skip_no_data
class TestEndToEndSubsample:

    def test_cli_subsample_pipeline(self, tmp_path, capsys):
        """Full pipeline: discover -> validate -> subsample -> build lookup."""
        from sma_merge.cli import cmd_discover

        cmd_discover(ONE_NICK_DATA)
        captured = capsys.readouterr()
        assert "FBD66244" in captured.out
        assert "OK" in captured.out

        from sma_merge.discover import discover_runs
        from sma_merge.subsample import subsample_pod5
        from sma_merge.tag import build_pod5_lookup

        runs = discover_runs(ONE_NICK_DATA)
        pod5_dir = runs[0].pod5_dir

        sub_pod5 = tmp_path / "sub20.pod5"
        selected = subsample_pod5(pod5_dir, sub_pod5, n_reads=20, seed=123)
        assert len(selected) == 20

        lookup = build_pod5_lookup(sub_pod5)
        assert len(lookup) == 20
        for rid in selected:
            assert rid in lookup

        reasons = [er for er, ns in lookup.values()]
        assert len(set(reasons)) >= 1
