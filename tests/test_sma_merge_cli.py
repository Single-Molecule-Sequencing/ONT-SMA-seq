"""Tests for sma-merge CLI."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest


class TestCliParsing:

    def test_discover_subcommand(self):
        from sma_merge.cli import build_parser
        parser = build_parser()
        args = parser.parse_args(["discover", "/tmp/experiment"])
        assert args.command == "discover"
        assert args.path == Path("/tmp/experiment")

    def test_merge_subcommand(self):
        from sma_merge.cli import build_parser
        parser = build_parser()
        args = parser.parse_args(["merge", "/tmp/experiment", "-o", "/tmp/out"])
        assert args.command == "merge"
        assert args.path == Path("/tmp/experiment")
        assert args.output == Path("/tmp/out")

    def test_subsample_subcommand(self):
        from sma_merge.cli import build_parser
        parser = build_parser()
        args = parser.parse_args(["subsample", "/tmp/exp", "-n", "5000", "-o", "/tmp/out"])
        assert args.command == "subsample"
        assert args.num_reads == 5000
        assert args.seed == 42  # default

    def test_subsample_custom_seed(self):
        from sma_merge.cli import build_parser
        parser = build_parser()
        args = parser.parse_args(["subsample", "/tmp/exp", "-n", "100", "-o", "/tmp/out", "--seed", "99"])
        assert args.seed == 99

    def test_merge_custom_device(self):
        from sma_merge.cli import build_parser
        parser = build_parser()
        args = parser.parse_args(["merge", "/tmp/exp", "-o", "/tmp/out", "--device", "cuda:0"])
        assert args.device == "cuda:0"

    def test_dorado_path_flag(self):
        from sma_merge.cli import build_parser
        parser = build_parser()
        args = parser.parse_args(["merge", "/tmp/exp", "-o", "/tmp/out", "--dorado", "/opt/dorado"])
        assert args.dorado == Path("/opt/dorado")


class TestDiscoverCommand:

    def test_discover_prints_table(self, tmp_path, capsys):
        from sma_merge.models import RunInfo

        mock_runs = [
            RunInfo(
                run_dir=tmp_path,
                flow_cell_id="FBD66244",
                device_id="MD-101527",
                protocol_group_id="test",
                basecall_model="hac@v5.2.0",
                sample_id="",
                run_id="abc",
                sample_rate=5000,
                pod5_dir=tmp_path / "pod5_pass",
                pod5_count=10,
                mod_base_models="",
            ),
        ]

        with patch("sma_merge.cli.discover_runs", return_value=mock_runs), \
             patch("sma_merge.cli.validate_runs") as mock_val:
            mock_val.return_value = []
            from sma_merge.cli import cmd_discover
            cmd_discover(tmp_path)

        captured = capsys.readouterr()
        assert "FBD66244" in captured.out
