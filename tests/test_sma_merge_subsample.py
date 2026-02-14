"""Tests for POD5 read ID collection and subsampling."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch, call
from uuid import UUID

import pytest


def _mock_pod5_reader(read_ids: list[str]):
    """Create a mock pod5 Reader that yields reads with given IDs."""
    reads = []
    for rid in read_ids:
        mock_read = MagicMock()
        mock_read.read_id = UUID(rid) if len(rid) == 36 else rid
        reads.append(mock_read)

    reader = MagicMock()
    reader.__enter__ = MagicMock(return_value=reader)
    reader.__exit__ = MagicMock(return_value=False)
    reader.reads.return_value = iter(reads)
    return reader


TEST_UUIDS = [f"00000000-0000-0000-0000-{i:012d}" for i in range(20)]


class TestCollectReadIds:

    def test_collects_from_single_file(self, tmp_path):
        pod5_dir = tmp_path / "pod5_pass"
        pod5_dir.mkdir()
        (pod5_dir / "chunk_0.pod5").touch()

        reader = _mock_pod5_reader(TEST_UUIDS[:5])
        with patch("sma_merge.subsample.p5.Reader", return_value=reader):
            from sma_merge.subsample import collect_read_ids
            ids = collect_read_ids(pod5_dir)

        assert len(ids) == 5
        assert ids[0] == TEST_UUIDS[0]

    def test_collects_from_subdirectories(self, tmp_path):
        pod5_dir = tmp_path / "pod5_pass"
        (pod5_dir / "mixed").mkdir(parents=True)
        (pod5_dir / "mixed" / "chunk_0.pod5").touch()
        (pod5_dir / "mixed" / "chunk_1.pod5").touch()

        reader1 = _mock_pod5_reader(TEST_UUIDS[:3])
        reader2 = _mock_pod5_reader(TEST_UUIDS[3:6])
        with patch("sma_merge.subsample.p5.Reader", side_effect=[reader1, reader2]):
            from sma_merge.subsample import collect_read_ids
            ids = collect_read_ids(pod5_dir)

        assert len(ids) == 6


class TestSubsample:

    def test_subsample_selects_n_reads(self, tmp_path):
        pod5_dir = tmp_path / "pod5_pass"
        pod5_dir.mkdir()
        (pod5_dir / "chunk_0.pod5").touch()
        output = tmp_path / "out.pod5"

        reader = _mock_pod5_reader(TEST_UUIDS[:10])
        with patch("sma_merge.subsample.p5.Reader", return_value=reader), \
             patch("sma_merge.subsample.subprocess.run") as mock_run:
            from sma_merge.subsample import subsample_pod5
            selected = subsample_pod5(pod5_dir, output, n_reads=5, seed=42)

        assert len(selected) == 5
        # With seed=42, results are deterministic
        selected2 = None
        reader2 = _mock_pod5_reader(TEST_UUIDS[:10])
        with patch("sma_merge.subsample.p5.Reader", return_value=reader2), \
             patch("sma_merge.subsample.subprocess.run"):
            selected2 = subsample_pod5(pod5_dir, output, n_reads=5, seed=42)
        assert selected == selected2

    def test_subsample_caps_at_total(self, tmp_path):
        pod5_dir = tmp_path / "pod5_pass"
        pod5_dir.mkdir()
        (pod5_dir / "chunk_0.pod5").touch()
        output = tmp_path / "out.pod5"

        reader = _mock_pod5_reader(TEST_UUIDS[:3])
        with patch("sma_merge.subsample.p5.Reader", return_value=reader), \
             patch("sma_merge.subsample.subprocess.run"):
            from sma_merge.subsample import subsample_pod5
            selected = subsample_pod5(pod5_dir, output, n_reads=100, seed=42)

        assert len(selected) == 3

    def test_subsample_calls_pod5_filter(self, tmp_path):
        pod5_dir = tmp_path / "pod5_pass"
        pod5_dir.mkdir()
        (pod5_dir / "chunk_0.pod5").touch()
        output = tmp_path / "out.pod5"

        reader = _mock_pod5_reader(TEST_UUIDS[:5])
        with patch("sma_merge.subsample.p5.Reader", return_value=reader), \
             patch("sma_merge.subsample.subprocess.run") as mock_run:
            from sma_merge.subsample import subsample_pod5
            subsample_pod5(pod5_dir, output, n_reads=3, seed=42)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "pod5"
        assert cmd[1] == "filter"
        assert "--output" in cmd
        assert "--ids" in cmd
