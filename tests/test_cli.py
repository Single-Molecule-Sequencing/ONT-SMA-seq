"""Tests for the calibrate CLI entry point."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


class TestCLIDiscoverSubcommand:
    def test_discover_prints_runs(self, tmp_path: Path):
        """discover subcommand should list found runs."""
        run_dir = tmp_path / "sample" / "20251228_2219_MD_FBD69411_abc123"
        run_dir.mkdir(parents=True)
        (run_dir / "final_summary_FBD69411_abc123_xyz.txt").write_text(
            "flow_cell_id=FBD69411\n"
            "device_id=MD-100098\n"
            "sample_id=sample\n"
            "experiment_id=exp1\n"
            "started=2025-12-28T22:19:00Z\n"
            "protocol_run_id=proto-1\n"
        )
        (run_dir / "bam_pass").mkdir()

        result = subprocess.run(
            [sys.executable, "-m", "calibrate.cli", "discover", str(tmp_path)],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq/bin",
        )
        assert result.returncode == 0
        assert "FBD69411" in result.stdout

    def test_validate_subcommand(self, tmp_path: Path):
        """validate subcommand should report validation status."""
        run_dir = tmp_path / "sample" / "20251228_2219_MD_FBD69411_abc123"
        run_dir.mkdir(parents=True)
        (run_dir / "final_summary_FBD69411_abc123_xyz.txt").write_text(
            "flow_cell_id=FBD69411\n"
            "device_id=MD-100098\n"
            "sample_id=sample\n"
            "experiment_id=exp1\n"
            "started=2025-12-28T22:19:00Z\n"
            "protocol_run_id=proto-1\n"
        )
        (run_dir / "bam_pass").mkdir()

        result = subprocess.run(
            [sys.executable, "-m", "calibrate.cli", "validate", str(tmp_path)],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq/bin",
        )
        assert result.returncode == 0
