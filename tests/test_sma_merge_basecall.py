"""Tests for POD5 merge and dorado basecalling wrappers."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest


class TestMergePod5s:

    def test_calls_pod5_merge(self, tmp_path):
        dir1 = tmp_path / "run1" / "pod5_pass"
        dir1.mkdir(parents=True)
        (dir1 / "a.pod5").touch()
        (dir1 / "b.pod5").touch()

        dir2 = tmp_path / "run2" / "pod5_pass"
        dir2.mkdir(parents=True)
        (dir2 / "c.pod5").touch()

        output = tmp_path / "merged.pod5"

        with patch("sma_merge.basecall.subprocess.run") as mock_run:
            from sma_merge.basecall import merge_pod5s
            merge_pod5s([dir1, dir2], output)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "pod5"
        assert cmd[1] == "merge"
        assert "--output" in cmd
        assert str(output) in cmd
        pod5_args = [a for a in cmd if a.endswith(".pod5") and a != str(output)]
        assert len(pod5_args) == 3

    def test_raises_on_no_pod5_files(self, tmp_path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()

        with pytest.raises(ValueError, match="No POD5 files"):
            from sma_merge.basecall import merge_pod5s
            merge_pod5s([empty_dir], tmp_path / "out.pod5")


class TestFindDorado:

    def test_returns_string(self):
        from sma_merge.basecall import find_dorado
        result = find_dorado()
        assert isinstance(result, str)


class TestBasecall:

    def test_calls_dorado_with_correct_flags(self, tmp_path):
        pod5 = tmp_path / "input.pod5"
        pod5.touch()
        output = tmp_path / "output.bam"

        with patch("sma_merge.basecall.subprocess.run") as mock_run, \
             patch("builtins.open", MagicMock()):
            from sma_merge.basecall import basecall
            basecall(
                pod5_path=pod5,
                output_bam=output,
                model="dna_r10.4.1_e8.2_400bps_hac@v5.2.0",
            )

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "--no-trim" in cmd
        assert "--kit-name" not in cmd
        assert "--emit-moves" in cmd
        assert "dna_r10.4.1_e8.2_400bps_hac@v5.2.0" in cmd

    def test_custom_dorado_path(self, tmp_path):
        pod5 = tmp_path / "input.pod5"
        pod5.touch()
        output = tmp_path / "output.bam"

        with patch("sma_merge.basecall.subprocess.run") as mock_run, \
             patch("builtins.open", MagicMock()):
            from sma_merge.basecall import basecall
            basecall(
                pod5_path=pod5,
                output_bam=output,
                model="hac@v5.2.0",
                dorado_path="/custom/dorado",
            )

        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "/custom/dorado"
