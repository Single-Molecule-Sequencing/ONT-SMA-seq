"""Tests for alignment engine."""
from __future__ import annotations

import pytest


class TestParseCigar:
    """Test CIGAR string parsing into (op, length) tuples."""

    def test_simple_match(self):
        from align import parse_cigar
        assert parse_cigar("10=") == [("=", 10)]

    def test_mixed_ops(self):
        from align import parse_cigar
        result = parse_cigar("5=1X3=2I4=1D2=")
        assert result == [("=", 5), ("X", 1), ("=", 3), ("I", 2), ("=", 4), ("D", 1), ("=", 2)]

    def test_empty_string(self):
        from align import parse_cigar
        assert parse_cigar("") == []


class TestCigarMetrics:
    """Test whole-alignment metrics computed from CIGAR."""

    def test_perfect_match(self):
        from align import compute_cigar_metrics
        m = compute_cigar_metrics("10=", ref_len=10, read_len=10)
        assert m["matches"] == 10
        assert m["identity"] == 1.0
        assert m["ref_coverage"] == 1.0
        assert m["max_ins"] == 0
        assert m["max_del"] == 0
        assert m["n_sig_indels"] == 0

    def test_with_mismatch(self):
        from align import compute_cigar_metrics
        m = compute_cigar_metrics("4=1X5=", ref_len=10, read_len=10)
        assert m["matches"] == 9
        assert m["mismatches"] == 1
        assert m["identity"] == pytest.approx(0.9)

    def test_with_insertion(self):
        from align import compute_cigar_metrics
        # 5 match, 3 insert, 5 match = read_len 13, ref consumed 10
        m = compute_cigar_metrics("5=3I5=", ref_len=10, read_len=13)
        assert m["max_ins"] == 3
        assert m["n_sig_indels"] == 0  # 3 < 5 threshold

    def test_significant_indel(self):
        from align import compute_cigar_metrics
        m = compute_cigar_metrics("5=6I5=", ref_len=10, read_len=16)
        assert m["max_ins"] == 6
        assert m["n_sig_indels"] == 1

    def test_with_deletion(self):
        from align import compute_cigar_metrics
        # 5 match, 2 del, 3 match = read_len 8, ref consumed 10
        m = compute_cigar_metrics("5=2D3=", ref_len=10, read_len=8)
        assert m["max_del"] == 2
        assert m["ref_coverage"] == pytest.approx(0.8)  # 8 ref bases covered by alignment matches


class TestSegmentedMetrics:
    """Test metrics computed within 5'/mid/3' reference segments."""

    def test_perfect_match_all_segments_identical(self):
        from align import compute_segmented_metrics
        # 30bp perfect match, ref_len=30 -> segments at 10, 20
        segs = compute_segmented_metrics("30=", ref_len=30)
        for seg in ["five_prime", "middle", "three_prime"]:
            assert segs[f"{seg}_identity"] == 1.0
            assert segs[f"{seg}_contiguity"] == 10
            assert segs[f"{seg}_gap_count"] == 0

    def test_mismatch_in_five_prime(self):
        from align import compute_segmented_metrics
        # First 10bp of ref: 3=1X6= -> 5' segment has identity 9/10
        segs = compute_segmented_metrics("3=1X6=10=10=", ref_len=30)
        assert segs["five_prime_identity"] == pytest.approx(0.9)
        assert segs["middle_identity"] == 1.0
        assert segs["three_prime_identity"] == 1.0

    def test_insertion_in_middle(self):
        from align import compute_segmented_metrics
        # ref positions 0-9 perfect, 10-19 has a 3bp insertion, 20-29 perfect
        segs = compute_segmented_metrics("10=5=3I5=10=", ref_len=30)
        assert segs["five_prime_gap_count"] == 0
        assert segs["middle_gap_count"] == 1
        assert segs["three_prime_gap_count"] == 0

    def test_contiguity_broken_by_deletion(self):
        from align import compute_segmented_metrics
        # 5' segment: 4=2D4= -> contiguity = 4 (longest run without gap)
        segs = compute_segmented_metrics("4=2D4=10=10=", ref_len=30)
        assert segs["five_prime_contiguity"] == 4
