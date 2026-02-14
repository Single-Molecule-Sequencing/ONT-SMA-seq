"""Tests for barcode separation analysis."""
from __future__ import annotations

import pytest

from calibrate_viz.separation import (
    compute_pairwise_distances,
    compute_separation_metrics,
)


class TestPairwiseDistances:
    def test_returns_symmetric_matrix(self):
        barcodes = {"nb05": "CACAAAGACACCGACAACTTTCTT", "nb10": "GAGAGGACAAAGGTTTCAACGCTT"}
        matrix = compute_pairwise_distances(barcodes)
        assert matrix["nb05"]["nb10"] == matrix["nb10"]["nb05"]

    def test_diagonal_is_zero(self):
        barcodes = {"nb05": "CACAAAGACACCGACAACTTTCTT", "nb10": "GAGAGGACAAAGGTTTCAACGCTT"}
        matrix = compute_pairwise_distances(barcodes)
        assert matrix["nb05"]["nb05"] == 0
        assert matrix["nb10"]["nb10"] == 0

    def test_different_barcodes_have_nonzero_distance(self):
        barcodes = {"nb05": "CACAAAGACACCGACAACTTTCTT", "nb10": "GAGAGGACAAAGGTTTCAACGCTT"}
        matrix = compute_pairwise_distances(barcodes)
        assert matrix["nb05"]["nb10"] > 0

    def test_three_barcodes(self):
        barcodes = {
            "nb05": "CACAAAGACACCGACAACTTTCTT",
            "nb10": "GAGAGGACAAAGGTTTCAACGCTT",
            "nb15": "ATGAACACTTCGGATTCTATGCTT",
        }
        matrix = compute_pairwise_distances(barcodes)
        assert len(matrix) == 3
        for bc in barcodes:
            assert len(matrix[bc]) == 3


class TestSeparationMetrics:
    def test_returns_per_barcode_metrics(self):
        read_eds = {
            "nb05": {"true_match": [2, 3, 1, 2], "next_best": [15, 16, 14, 15]},
            "nb10": {"true_match": [1, 2, 3, 1], "next_best": [14, 13, 15, 14]},
        }
        metrics = compute_separation_metrics(read_eds)
        assert "nb05" in metrics
        assert metrics["nb05"]["mean_true_ed"] < metrics["nb05"]["mean_next_best_ed"]
        assert metrics["nb05"]["separation_gap"] > 10

    def test_error_rate_zero_for_well_separated(self):
        read_eds = {
            "nb05": {"true_match": [1, 2, 1, 2], "next_best": [15, 16, 14, 15]},
        }
        metrics = compute_separation_metrics(read_eds)
        assert metrics["nb05"]["estimated_error_rate"] == 0.0

    def test_high_error_rate_for_overlapping(self):
        read_eds = {
            "nb05": {"true_match": [10, 12, 11, 13], "next_best": [11, 10, 12, 11]},
        }
        metrics = compute_separation_metrics(read_eds)
        assert metrics["nb05"]["estimated_error_rate"] > 0.5
