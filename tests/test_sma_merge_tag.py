"""Tests for POD5 end_reason extraction and BAM tagging."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pysam
import pytest


def _make_ubam(bam_path: Path, read_ids: list[str], seq_len: int = 100):
    """Create a minimal unaligned BAM for testing."""
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "unsorted"}})
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as af:
        for rid in read_ids:
            seg = pysam.AlignedSegment(header)
            seg.query_name = rid
            seg.query_sequence = "A" * seq_len
            seg.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
            seg.flag = 4  # unmapped
            af.write(seg)


TEST_UUIDS = [f"00000000-0000-0000-0000-{i:012d}" for i in range(5)]


class TestBuildPod5Lookup:

    def test_builds_lookup_from_mock(self, tmp_path):
        pod5_file = tmp_path / "test.pod5"
        pod5_file.touch()

        mock_reads = []
        for i, rid in enumerate(TEST_UUIDS[:3]):
            r = MagicMock()
            r.read_id = rid
            r.end_reason.reason.name = ["SIGNAL_POSITIVE", "MUX_CHANGE", "DATA_SERVICE_UNBLOCK_MUX_CHANGE"][i]
            r.num_samples = [5000, 3000, 8000][i]
            mock_reads.append(r)

        reader = MagicMock()
        reader.__enter__ = MagicMock(return_value=reader)
        reader.__exit__ = MagicMock(return_value=False)
        reader.reads.return_value = iter(mock_reads)

        with patch("sma_merge.tag.p5.Reader", return_value=reader):
            from sma_merge.tag import build_pod5_lookup
            lookup = build_pod5_lookup(pod5_file)

        assert len(lookup) == 3
        assert lookup[TEST_UUIDS[0]] == ("signal_positive", 5000)
        assert lookup[TEST_UUIDS[1]] == ("mux_change", 3000)
        assert lookup[TEST_UUIDS[2]] == ("data_service_unblock_mux_change", 8000)


class TestTagBam:

    def test_adds_er_and_sl_tags(self, tmp_path):
        input_bam = tmp_path / "input.bam"
        output_bam = tmp_path / "output.bam"
        _make_ubam(input_bam, TEST_UUIDS[:3])

        lookup = {
            TEST_UUIDS[0]: ("signal_positive", 5000),
            TEST_UUIDS[1]: ("mux_change", 3000),
            TEST_UUIDS[2]: ("data_service_unblock_mux_change", 8000),
        }

        from sma_merge.tag import tag_bam
        tagged = tag_bam(input_bam, output_bam, lookup)

        assert tagged == 3

        with pysam.AlignmentFile(str(output_bam), check_sq=False) as af:
            reads = list(af)
        assert len(reads) == 3
        assert reads[0].get_tag("er") == "signal_positive"
        assert reads[0].get_tag("sl") == 5000
        assert reads[2].get_tag("er") == "data_service_unblock_mux_change"

    def test_handles_missing_reads_in_lookup(self, tmp_path):
        input_bam = tmp_path / "input.bam"
        output_bam = tmp_path / "output.bam"
        _make_ubam(input_bam, TEST_UUIDS[:3])

        lookup = {TEST_UUIDS[0]: ("signal_positive", 5000)}

        from sma_merge.tag import tag_bam
        tagged = tag_bam(input_bam, output_bam, lookup)

        assert tagged == 1

        with pysam.AlignmentFile(str(output_bam), check_sq=False) as af:
            reads = list(af)
        assert len(reads) == 3
        assert reads[0].has_tag("er")
        assert not reads[1].has_tag("er")
