"""Tests for threshold optimization and ROC computation."""
from __future__ import annotations

import numpy as np
import pytest

from calibrate_viz.thresholds import (
    compute_roc,
    recommend_threshold,
    export_thresholds_to_toml,
)


class TestComputeROC:
    def test_returns_fpr_tpr_thresholds(self):
        true_confs = np.array([0.9, 0.85, 0.8, 0.95, 0.88])
        false_confs = np.array([0.3, 0.2, 0.4, 0.15, 0.35])
        fpr, tpr, thresholds = compute_roc(true_confs, false_confs)
        assert len(fpr) == len(tpr) == len(thresholds)
        # FPR should be sorted ascending
        assert np.all(np.diff(fpr) >= -1e-10)

    def test_perfect_separation_gives_auc_near_1(self):
        true_confs = np.array([0.9, 0.95, 0.85, 0.92])
        false_confs = np.array([0.1, 0.05, 0.15, 0.08])
        fpr, tpr, _ = compute_roc(true_confs, false_confs, n_points=500)
        auc = float(np.trapz(tpr, fpr))
        assert auc > 0.95

    def test_random_gives_auc_near_half(self):
        rng = np.random.default_rng(42)
        true_confs = rng.uniform(0, 1, 100)
        false_confs = rng.uniform(0, 1, 100)
        fpr, tpr, _ = compute_roc(true_confs, false_confs)
        auc = float(np.trapz(tpr, fpr))
        assert 0.3 < auc < 0.7


class TestRecommendThreshold:
    def test_returns_threshold_and_metrics(self):
        true_confs = np.array([0.9, 0.85, 0.8, 0.95])
        false_confs = np.array([0.3, 0.2, 0.4, 0.15])
        result = recommend_threshold(true_confs, false_confs, target_fpr=0.05)
        assert "threshold" in result
        assert "sensitivity" in result
        assert "specificity" in result
        assert "f1" in result
        assert result["threshold"] > 0.4

    def test_high_sensitivity_for_separated_data(self):
        true_confs = np.array([0.9, 0.85, 0.8, 0.95, 0.88])
        false_confs = np.array([0.1, 0.05, 0.15, 0.08, 0.12])
        result = recommend_threshold(true_confs, false_confs, target_fpr=0.1)
        assert result["sensitivity"] > 0.8


class TestExportToToml:
    def test_writes_valid_toml(self, tmp_path):
        out = tmp_path / "construct.toml"
        export_thresholds_to_toml(
            output_path=out,
            start_barcode_min=0.65,
            full_length_threshold=0.78,
            flank_max_error_rate=0.5,
        )
        assert out.exists()
        import tomllib
        with out.open("rb") as f:
            data = tomllib.load(f)
        assert data["sma"]["confidence"]["start_barcode_min"] == 0.65
        assert data["sma"]["confidence"]["full_length_threshold"] == 0.78
