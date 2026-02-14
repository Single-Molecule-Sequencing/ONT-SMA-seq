"""Tests for SMA-seq config file generation."""
from __future__ import annotations

from pathlib import Path

import pytest

from calibrate.config_gen import generate_construct_toml, generate_sample_sheet


class TestGenerateConstructToml:
    def test_writes_valid_toml(self, tmp_path: Path):
        targets = [
            {"barcode1": "NB05", "barcode2": "NB10", "alias": "V04_2_fwd"},
            {"barcode1": "NB10", "barcode2": "NB05", "alias": "V04_2_rev"},
        ]
        out = tmp_path / "construct.toml"
        generate_construct_toml(
            targets=targets,
            kit="SQK-NBD114-96",
            mode="dual_independent",
            output_path=out,
            mask1_front="AAGGTTAA",
            mask1_rear="CAGCACCT",
            mask2_front="AGGTGCTG",
            mask2_rear="TTAACCTT",
        )
        assert out.exists()
        import tomllib
        with out.open("rb") as f:
            data = tomllib.load(f)
        assert data["arrangement"]["kit"] == "SQK-NBD114-96"
        assert data["sma"]["mode"] == "dual_independent"
        assert len(data["sma"]["targets"]) == 2

    def test_start_only_mode_omits_mask2(self, tmp_path: Path):
        targets = [{"barcode1": "NB05", "alias": "V04_2_fwd"}]
        out = tmp_path / "construct.toml"
        generate_construct_toml(
            targets=targets,
            kit="SQK-NBD114-96",
            mode="start_only",
            output_path=out,
            mask1_front="AAGGTTAA",
            mask1_rear="CAGCACCT",
        )
        import tomllib
        with out.open("rb") as f:
            data = tomllib.load(f)
        assert "mask2_front" not in data["arrangement"]


class TestGenerateSampleSheet:
    def test_writes_csv(self, tmp_path: Path):
        entries = [
            {"barcode": "barcode05--barcode10", "alias": "V04_2_fwd", "type": "test_sample"},
        ]
        out = tmp_path / "sample_sheet.csv"
        generate_sample_sheet(
            entries=entries,
            flow_cell_id="FBD69411",
            kit="SQK-NBD114-96",
            experiment_id="my_exp",
            output_path=out,
        )
        assert out.exists()
        lines = out.read_text().strip().splitlines()
        assert len(lines) == 2  # header + 1 entry
        assert "barcode05--barcode10" in lines[1]
