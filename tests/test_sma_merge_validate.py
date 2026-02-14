"""Tests for run validation and grouping."""
from __future__ import annotations

from pathlib import Path

import pytest

from sma_merge.models import RunInfo


def _run(flow_cell="FBD66244", device="MD-101527", protocol="exp1",
         model="hac@v5.2.0", sample_id="", run_id="abc"):
    return RunInfo(
        run_dir=Path("/tmp/run"),
        flow_cell_id=flow_cell,
        device_id=device,
        protocol_group_id=protocol,
        basecall_model=model,
        sample_id=sample_id,
        run_id=run_id,
        sample_rate=5000,
        pod5_dir=Path("/tmp/run/pod5_pass"),
        pod5_count=10,
        mod_base_models="",
    )


class TestValidateRuns:

    def test_single_run_is_consistent(self):
        from sma_merge.validate import validate_runs
        groups = validate_runs([_run()])
        assert len(groups) == 1
        assert groups[0].is_consistent
        assert groups[0].flow_cell_id == "FBD66244"

    def test_same_flowcell_grouped(self):
        from sma_merge.validate import validate_runs
        runs = [_run(run_id="aaa"), _run(run_id="bbb")]
        groups = validate_runs(runs)
        assert len(groups) == 1
        assert len(groups[0].runs) == 2

    def test_different_flowcells_separate_groups(self):
        from sma_merge.validate import validate_runs
        runs = [_run(flow_cell="FC1"), _run(flow_cell="FC2")]
        groups = validate_runs(runs)
        assert len(groups) == 2

    def test_inconsistent_model_flagged(self):
        from sma_merge.validate import validate_runs
        runs = [_run(model="hac@v5.2.0"), _run(model="sup@v5.2.0")]
        groups = validate_runs(runs)
        assert len(groups) == 1
        assert not groups[0].is_consistent
        assert any("model" in issue.lower() for issue in groups[0].issues)

    def test_inconsistent_protocol_flagged(self):
        from sma_merge.validate import validate_runs
        runs = [_run(protocol="exp1"), _run(protocol="exp2")]
        groups = validate_runs(runs)
        assert not groups[0].is_consistent
        assert any("protocol" in issue.lower() for issue in groups[0].issues)

    def test_empty_sample_id_warned(self):
        from sma_merge.validate import validate_runs
        groups = validate_runs([_run(sample_id="")])
        assert any("sample_id" in issue.lower() for issue in groups[0].issues)

    def test_nonempty_sample_id_no_warning(self):
        from sma_merge.validate import validate_runs
        groups = validate_runs([_run(sample_id="my_sample")])
        assert not any("sample_id" in issue.lower() for issue in groups[0].issues)


class TestFormatValidation:

    def test_format_includes_status(self):
        from sma_merge.validate import validate_runs, format_validation
        groups = validate_runs([_run()])
        output = format_validation(groups)
        assert "PASS" in output or "OK" in output
