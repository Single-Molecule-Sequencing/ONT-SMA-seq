"""Tests for importing config from MinKNOW output directory."""

from __future__ import annotations

from pathlib import Path

import pytest

from viz.config_store import ConfigStore


@pytest.fixture()
def minknow_output(tmp_path: Path) -> Path:
    """Create a fake MinKNOW output directory structure."""
    run_dir = tmp_path / "data" / "experiment_group" / "sample1" / "20260214_1200_MN12345_FAL12345_abcdef01"
    run_dir.mkdir(parents=True)

    (run_dir / "sample_sheet_FAL12345_20260214_1200_abcdef01.csv").write_text(
        "protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,flow_cell_product_code,kit,barcode,alias,type\n"
        "abc-def,MN12345,FAL12345,sample1,exp_group,FLO-PRO114,SQK-NBD114-96,barcode05--barcode10,target_A,test_sample\n"
    )

    (run_dir / "final_summary_FAL12345_abcdef01_12345678.txt").write_text(
        "protocol_run_id=abc-def\nflow_cell_id=FAL12345\nsample_id=sample1\n"
    )

    return run_dir


class TestMinKNOWImport:
    def test_import_from_minknow_dir(self, minknow_output: Path, tmp_path: Path):
        exp_dir = tmp_path / "imported_experiment"
        exp_dir.mkdir()
        store = ConfigStore(exp_dir)
        store.import_from_minknow(minknow_output)

        entries = store.read_sample_sheet()
        assert len(entries) == 1
        assert entries[0].alias == "target_A"

    def test_import_sets_description(self, minknow_output: Path, tmp_path: Path):
        exp_dir = tmp_path / "imported_experiment"
        exp_dir.mkdir()
        store = ConfigStore(exp_dir)
        store.import_from_minknow(minknow_output)
        assert "sample1" in store.experiment_config.description

    def test_import_no_sample_sheet(self, tmp_path: Path):
        """Import from dir without sample sheet doesn't crash."""
        empty_minknow = tmp_path / "empty_minknow"
        empty_minknow.mkdir()
        exp_dir = tmp_path / "imported"
        exp_dir.mkdir()
        store = ConfigStore(exp_dir)
        store.import_from_minknow(empty_minknow)
        entries = store.read_sample_sheet()
        assert len(entries) == 0
