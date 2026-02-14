"""Tests for cross-experiment comparison."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from calibrate_viz.comparison import (
    compute_experiment_summary,
    compute_overlay_distributions,
)


def _make_db(tmp_path: Path, name: str, reads: list[tuple]) -> Path:
    db_path = tmp_path / f"{name}.db"
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE Reads (
            read_id TEXT PRIMARY KEY, readlen INTEGER,
            signal_duration_s REAL, mean_qscore REAL, ER TEXT,
            bc_start_id TEXT, bc_start_conf REAL,
            bc_end_id TEXT, bc_end_conf REAL,
            tgt_id TEXT, trunc_level TEXT, ed INTEGER
        )
    """)
    conn.executemany(
        "INSERT INTO Reads VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        reads,
    )
    conn.commit()
    conn.close()
    return db_path


class TestExperimentSummary:
    def test_returns_key_metrics(self, tmp_path: Path):
        db = _make_db(tmp_path, "exp1", [
            ("r1", 500, 1.5, 12.0, "signal_positive", "nb05", 0.9, "nb10", 0.8, "t1", "full_length", 5),
            ("r2", 300, 0.9, 10.0, "signal_positive", "nb05", 0.7, None, 0.1, "t1", "bc1_target", 20),
        ])
        summary = compute_experiment_summary(db)
        assert summary["total_reads"] == 2
        assert summary["mean_readlen"] == 400.0
        assert "trunc_proportions" in summary
        assert "end_reason_proportions" in summary

    def test_trunc_proportions_sum_to_one(self, tmp_path: Path):
        db = _make_db(tmp_path, "exp2", [
            ("r1", 500, 1.5, 12.0, "sp", "nb05", 0.9, "nb10", 0.8, "t1", "full_length", 5),
            ("r2", 300, 0.9, 10.0, "sp", "nb05", 0.7, None, 0.1, "t1", "bc1_target", 20),
            ("r3", 100, 0.3, 8.0, "sp", "nb05", 0.3, None, 0.0, "t1", "adapter_only", 50),
        ])
        summary = compute_experiment_summary(db)
        total = sum(summary["trunc_proportions"].values())
        assert abs(total - 1.0) < 0.01

    def test_empty_db(self, tmp_path: Path):
        db = _make_db(tmp_path, "empty", [])
        summary = compute_experiment_summary(db)
        assert summary["total_reads"] == 0

    def test_mean_qscore(self, tmp_path: Path):
        db = _make_db(tmp_path, "qsc", [
            ("r1", 500, 1.5, 12.0, "sp", "nb05", 0.9, "nb10", 0.8, "t1", "fl", 5),
            ("r2", 300, 0.9, 8.0, "sp", "nb05", 0.7, None, 0.1, "t1", "bt", 20),
        ])
        summary = compute_experiment_summary(db)
        assert summary["mean_qscore"] == 10.0

    def test_barcode_and_target_counts(self, tmp_path: Path):
        db = _make_db(tmp_path, "counts", [
            ("r1", 500, 1.5, 12.0, "sp", "nb05", 0.9, "nb10", 0.8, "t1", "fl", 5),
            ("r2", 300, 0.9, 10.0, "sp", "nb10", 0.7, None, 0.1, "t2", "bt", 20),
            ("r3", 400, 1.2, 11.0, "sp", "nb05", 0.8, "nb10", 0.7, "t1", "fl", 10),
        ])
        summary = compute_experiment_summary(db)
        assert summary["barcode_count"] == 2
        assert summary["target_count"] == 2

    def test_end_reason_proportions(self, tmp_path: Path):
        db = _make_db(tmp_path, "er", [
            ("r1", 500, 1.5, 12.0, "signal_positive", "nb05", 0.9, "nb10", 0.8, "t1", "fl", 5),
            ("r2", 300, 0.9, 10.0, "signal_positive", "nb05", 0.7, None, 0.1, "t1", "bt", 20),
            ("r3", 100, 0.3, 8.0, "mux_change", "nb05", 0.3, None, 0.0, "t1", "ao", 50),
        ])
        summary = compute_experiment_summary(db)
        er = summary["end_reason_proportions"]
        assert "signal_positive" in er
        assert "mux_change" in er
        assert abs(er["signal_positive"] - 2 / 3) < 0.01
        assert abs(er["mux_change"] - 1 / 3) < 0.01


class TestOverlayDistributions:
    def test_overlay_returns_per_db_kde(self, tmp_path: Path):
        db1 = _make_db(tmp_path, "exp1", [
            ("r1", 500, 1.5, 12.0, "sp", "nb05", 0.9, "nb10", 0.8, "t1", "fl", 5),
            ("r2", 480, 1.4, 11.5, "sp", "nb05", 0.85, "nb10", 0.7, "t1", "fl", 8),
            ("r3", 510, 1.6, 12.5, "sp", "nb05", 0.95, "nb10", 0.9, "t1", "fl", 3),
        ])
        db2 = _make_db(tmp_path, "exp2", [
            ("r1", 300, 0.9, 10.0, "sp", "nb05", 0.7, "nb10", 0.5, "t1", "bt", 20),
            ("r2", 290, 0.85, 9.5, "sp", "nb05", 0.65, "nb10", 0.4, "t1", "bt", 22),
            ("r3", 310, 0.95, 10.5, "sp", "nb05", 0.75, "nb10", 0.6, "t1", "bt", 18),
        ])
        result = compute_overlay_distributions([db1, db2], column="readlen")
        assert len(result) == 2
        assert result[0]["db_name"] == "exp1"
        assert result[1]["db_name"] == "exp2"
        assert len(result[0]["x"]) > 0
        assert len(result[1]["x"]) > 0

    def test_overlay_single_db(self, tmp_path: Path):
        db = _make_db(tmp_path, "single", [
            ("r1", 500, 1.5, 12.0, "sp", "nb05", 0.9, "nb10", 0.8, "t1", "fl", 5),
            ("r2", 480, 1.4, 11.5, "sp", "nb05", 0.85, "nb10", 0.7, "t1", "fl", 8),
            ("r3", 510, 1.6, 12.5, "sp", "nb05", 0.95, "nb10", 0.9, "t1", "fl", 3),
        ])
        result = compute_overlay_distributions([db], column="readlen")
        assert len(result) == 1
        assert result[0]["db_name"] == "single"

    def test_overlay_empty_db(self, tmp_path: Path):
        db = _make_db(tmp_path, "empty", [])
        result = compute_overlay_distributions([db], column="readlen")
        assert len(result) == 1
        assert result[0]["x"] == []
        assert result[0]["count"] == 0

    def test_overlay_different_column(self, tmp_path: Path):
        db = _make_db(tmp_path, "sig", [
            ("r1", 500, 1.5, 12.0, "sp", "nb05", 0.9, "nb10", 0.8, "t1", "fl", 5),
            ("r2", 480, 1.4, 11.5, "sp", "nb05", 0.85, "nb10", 0.7, "t1", "fl", 8),
            ("r3", 510, 1.6, 12.5, "sp", "nb05", 0.95, "nb10", 0.9, "t1", "fl", 3),
        ])
        result = compute_overlay_distributions([db], column="signal_duration_s")
        assert len(result) == 1
        assert len(result[0]["x"]) > 0
