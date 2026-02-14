"""Tests for QC visualization module."""
from __future__ import annotations

import csv
from pathlib import Path

import pytest


def _make_classification_tsv(tmp_path: Path, rows: list[dict]) -> Path:
    """Write a classification.tsv for testing."""
    tsv = tmp_path / "classification.tsv"
    fields = ["read_id", "read_len", "assigned_ref", "best_ned", "second_ned", "margin", "confidence_flag", "end_reason"]
    with open(tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return tsv


def _make_alignments_tsv(tmp_path: Path, rows: list[dict]) -> Path:
    """Write an alignments.tsv for testing."""
    from align import METRIC_COLUMNS
    tsv = tmp_path / "alignments.tsv"
    with open(tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=METRIC_COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return tsv


class TestLoadData:
    """Test TSV loading utilities."""

    def test_load_classification(self, tmp_path):
        from qc import load_classification
        tsv = _make_classification_tsv(tmp_path, [
            {"read_id": "r1", "read_len": "100", "assigned_ref": "t1",
             "best_ned": "0.05", "second_ned": "0.4", "margin": "0.35",
             "confidence_flag": "HIGH", "end_reason": "signal_positive"},
        ])
        df = load_classification(tsv)
        assert len(df) == 1
        assert df["margin"].dtype == float

    def test_load_alignments(self, tmp_path):
        from qc import load_alignments
        from align import METRIC_COLUMNS
        row = {col: "0" for col in METRIC_COLUMNS}
        row["read_id"] = "r1"
        row["ref_id"] = "t1"
        tsv = _make_alignments_tsv(tmp_path, [row])
        df = load_alignments(tsv)
        assert len(df) == 1


def _sample_classification_df(n=100):
    """Generate a sample classification DataFrame for plotting tests."""
    import pandas as pd
    import numpy as np
    rng = np.random.default_rng(42)
    targets = ["t1", "t2"]
    ers = ["signal_positive", "data_service_unblock_mux_change"]
    return pd.DataFrame({
        "read_id": [f"r{i}" for i in range(n)],
        "read_len": rng.integers(100, 300, n),
        "assigned_ref": rng.choice(targets, n),
        "best_ned": rng.uniform(0.01, 0.3, n),
        "second_ned": rng.uniform(0.3, 0.8, n),
        "margin": rng.uniform(0.05, 0.5, n),
        "confidence_flag": rng.choice(["HIGH", "LOW", "POOR"], n, p=[0.8, 0.15, 0.05]),
        "end_reason": rng.choice(ers, n),
    })


def _sample_alignments_df(n_reads=50, n_refs=2):
    """Generate a sample alignments DataFrame."""
    import pandas as pd
    import numpy as np
    rng = np.random.default_rng(42)
    rows = []
    for i in range(n_reads):
        for j in range(n_refs):
            rows.append({
                "read_id": f"r{i}", "ref_id": f"t{j}",
                "read_len": 200, "ref_len": 200,
                "ed": int(rng.integers(1, 50)),
                "ned": float(rng.uniform(0.01, 0.3)),
                "identity": float(rng.uniform(0.7, 1.0)),
                "ref_coverage": float(rng.uniform(0.8, 1.0)),
                "read_to_ref_ratio": float(rng.uniform(0.8, 1.2)),
                "seg5_identity": float(rng.uniform(0.7, 1.0)),
                "seg5_contiguity": int(rng.integers(10, 60)),
                "segM_identity": float(rng.uniform(0.7, 1.0)),
                "segM_contiguity": int(rng.integers(10, 60)),
                "seg3_identity": float(rng.uniform(0.7, 1.0)),
                "seg3_contiguity": int(rng.integers(10, 60)),
                "max_ins": int(rng.integers(0, 10)),
                "max_del": int(rng.integers(0, 10)),
                "n_sig_indels": int(rng.integers(0, 3)),
                "five_prime_offset": int(rng.integers(0, 5)),
                "five_prime_identity_20": float(rng.uniform(0.7, 1.0)),
                "three_prime_offset": int(rng.integers(0, 5)),
                "three_prime_identity_20": float(rng.uniform(0.7, 1.0)),
                "rank": j + 1,
                "margin": float(rng.uniform(0.05, 0.5)),
            })
    return pd.DataFrame(rows)


class TestMarginKDE:
    """Test margin KDE plot generation."""

    def test_creates_png(self, tmp_path):
        from qc import plot_margin_kde
        tsv = _make_classification_tsv(tmp_path, [
            {"read_id": f"r{i}", "read_len": "100", "assigned_ref": "t1",
             "best_ned": "0.05", "second_ned": "0.4", "margin": f"{0.1 + i * 0.01}",
             "confidence_flag": "HIGH", "end_reason": "signal_positive"}
            for i in range(50)
        ])
        from qc import load_classification
        df = load_classification(tsv)
        out = tmp_path / "qc"
        out.mkdir()
        plot_margin_kde(df, out)
        assert (out / "margin_kde.png").exists()


class TestNedByTarget:
    def test_creates_png(self, tmp_path):
        from qc import plot_ned_by_target
        df_align = _sample_alignments_df()
        df_class = _sample_classification_df()
        out = tmp_path / "qc"
        out.mkdir()
        plot_ned_by_target(df_align, df_class, out)
        assert (out / "ned_by_target.png").exists()


class TestLengthVsNed:
    def test_creates_png(self, tmp_path):
        from qc import plot_length_vs_ned
        df = _sample_classification_df()
        out = tmp_path / "qc"
        out.mkdir()
        plot_length_vs_ned(df, out)
        assert (out / "length_vs_ned_scatter.png").exists()


class TestSegmentedViolins:
    def test_creates_png(self, tmp_path):
        from qc import plot_segmented_violins
        df_align = _sample_alignments_df()
        df_class = _sample_classification_df(50)
        out = tmp_path / "qc"
        out.mkdir()
        plot_segmented_violins(df_align, df_class, out)
        assert (out / "segmented_identity_violins.png").exists()


class TestFivePrimeQuality:
    def test_creates_png(self, tmp_path):
        from qc import plot_five_prime_quality
        df_align = _sample_alignments_df()
        df_class = _sample_classification_df(50)
        out = tmp_path / "qc"
        out.mkdir()
        plot_five_prime_quality(df_align, df_class, out)
        assert (out / "five_prime_quality.png").exists()


class TestIndelProfile:
    def test_creates_png(self, tmp_path):
        from qc import plot_indel_profile
        df_align = _sample_alignments_df()
        df_class = _sample_classification_df(50)
        out = tmp_path / "qc"
        out.mkdir()
        plot_indel_profile(df_align, df_class, out)
        assert (out / "indel_profile.png").exists()


class TestThresholdSweep:
    def test_creates_png(self, tmp_path):
        from qc import plot_threshold_sweep
        df = _sample_classification_df()
        out = tmp_path / "qc"
        out.mkdir()
        plot_threshold_sweep(df, out)
        assert (out / "threshold_sweep.png").exists()


class TestClassificationHeatmap:
    def test_creates_png(self, tmp_path):
        from qc import plot_classification_heatmap
        df_align = _sample_alignments_df()
        df_class = _sample_classification_df(50)
        out = tmp_path / "qc"
        out.mkdir()
        plot_classification_heatmap(df_align, df_class, out)
        assert (out / "classification_heatmap.png").exists()


class TestPerReadProfiles:
    def test_creates_png(self, tmp_path):
        from qc import plot_per_read_profiles
        df_align = _sample_alignments_df()
        out = tmp_path / "qc"
        out.mkdir()
        plot_per_read_profiles(df_align, out, n_sample=10)
        assert (out / "per_read_profiles_sample.png").exists()


class TestEndReasonFacets:
    def test_creates_faceted_plots(self, tmp_path):
        from qc import plot_end_reason_facets
        df_align = _sample_alignments_df()
        df_class = _sample_classification_df(50)
        out = tmp_path / "qc"
        out.mkdir()
        plot_end_reason_facets(df_align, df_class, out)
        assert (out / "margin_by_endreason.png").exists()


class TestRunQC:
    def test_generates_all_plots(self, tmp_path):
        from qc import run_qc

        df_class = _sample_classification_df(100)
        df_align = _sample_alignments_df(100, 2)

        class_tsv = tmp_path / "classification.tsv"
        align_tsv = tmp_path / "alignments.tsv"
        df_class.to_csv(class_tsv, sep="\t", index=False)
        df_align.to_csv(align_tsv, sep="\t", index=False)

        qc_dir = tmp_path / "qc"
        run_qc(align_tsv, class_tsv, qc_dir)

        expected = [
            "classification_heatmap.png",
            "margin_kde.png",
            "ned_by_target.png",
            "segmented_identity_violins.png",
            "length_vs_ned_scatter.png",
            "five_prime_quality.png",
            "indel_profile.png",
            "threshold_sweep.png",
            "per_read_profiles_sample.png",
            "margin_by_endreason.png",
        ]
        for name in expected:
            assert (qc_dir / name).exists(), f"Missing: {name}"
