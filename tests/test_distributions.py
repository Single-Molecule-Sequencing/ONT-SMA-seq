"""Tests for KDE distribution computation."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
import pytest

from calibrate_viz.distributions import (
    compute_kde,
    load_distribution_data,
    find_kde_peaks,
    compute_grouped_kde,
)


@pytest.fixture
def sample_db(tmp_path: Path) -> Path:
    """Create a minimal SQLite DB with reads for distribution testing."""
    db_path = tmp_path / "test.db"
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE Reads (
            read_id TEXT PRIMARY KEY,
            readlen INTEGER,
            signal_duration_s REAL,
            mean_qscore REAL,
            ER TEXT,
            bc_start_id TEXT,
            bc_end_id TEXT,
            bc_start_conf REAL,
            bc_end_conf REAL,
            tgt_id TEXT,
            trunc_level TEXT
        )
    """)
    reads = [
        ("r1", 500, 1.5, 12.0, "signal_positive", "nb05", "nb10", 0.9, 0.8, "V04_2_fwd", "full_length"),
        ("r2", 480, 1.4, 11.5, "signal_positive", "nb05", "nb10", 0.85, 0.75, "V04_2_fwd", "full_length"),
        ("r3", 510, 1.6, 12.5, "signal_positive", "nb05", "nb10", 0.95, 0.9, "V04_2_fwd", "full_length"),
        ("r4", 200, 0.7, 9.0, "data_service_unblock_mux_change", "nb05", None, 0.7, 0.1, "V04_2_fwd", "bc1_target"),
        ("r5", 300, 0.9, 10.0, "signal_positive", "nb10", "nb05", 0.88, 0.82, "V04_2_rev", "full_length"),
        ("r6", 490, 1.45, 11.8, "signal_positive", "nb05", "nb10", 0.87, 0.77, "V04_2_fwd", "full_length"),
        ("r7", 520, 1.65, 12.7, "signal_positive", "nb05", "nb10", 0.92, 0.85, "V04_2_fwd", "full_length"),
        ("r8", 310, 0.95, 10.2, "signal_positive", "nb10", "nb05", 0.86, 0.80, "V04_2_rev", "full_length"),
        ("r9", 495, 1.48, 11.9, "signal_positive", "nb05", "nb10", 0.89, 0.78, "V04_2_fwd", "full_length"),
        ("r10", 505, 1.55, 12.2, "signal_positive", "nb05", "nb10", 0.91, 0.83, "V04_2_fwd", "full_length"),
    ]
    conn.executemany(
        "INSERT INTO Reads VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", reads,
    )
    conn.commit()
    conn.close()
    return db_path


class TestComputeKDE:
    def test_returns_x_and_y_arrays(self):
        data = np.array([100, 200, 300, 400, 500] * 10)
        x, y = compute_kde(data)
        assert len(x) == len(y)
        assert len(x) == 512

    def test_y_values_are_nonnegative(self):
        data = np.array([100, 200, 300, 400, 500] * 10)
        _, y = compute_kde(data)
        assert np.all(y >= 0)

    def test_empty_for_single_value(self):
        data = np.array([42.0])
        x, y = compute_kde(data)
        assert len(x) == 0

    def test_custom_n_points(self):
        data = np.array([1, 2, 3, 4, 5] * 10)
        x, _ = compute_kde(data, n_points=256)
        assert len(x) == 256


class TestLoadDistributionData:
    def test_loads_readlen(self, sample_db):
        data = load_distribution_data(sample_db, column="readlen")
        assert len(data) == 10

    def test_loads_signal_duration(self, sample_db):
        data = load_distribution_data(sample_db, column="signal_duration_s")
        assert len(data) == 10

    def test_groups_by_barcode(self, sample_db):
        groups = load_distribution_data(sample_db, column="readlen", group_by="bc_start_id")
        assert "nb05" in groups
        assert "nb10" in groups
        assert len(groups["nb05"]) == 8
        assert len(groups["nb10"]) == 2

    def test_groups_by_end_reason(self, sample_db):
        groups = load_distribution_data(sample_db, column="readlen", group_by="ER")
        assert "signal_positive" in groups
        assert "data_service_unblock_mux_change" in groups

    def test_groups_by_target(self, sample_db):
        groups = load_distribution_data(sample_db, column="readlen", group_by="tgt_id")
        assert "V04_2_fwd" in groups
        assert "V04_2_rev" in groups


class TestFindKDEPeaks:
    def test_finds_single_peak(self):
        x = np.linspace(0, 10, 500)
        y = np.exp(-((x - 5) ** 2) / 0.5)
        peaks = find_kde_peaks(x, y)
        assert len(peaks) >= 1
        assert any(abs(p - 5.0) < 0.5 for p in peaks)

    def test_finds_two_peaks(self):
        x = np.linspace(0, 20, 500)
        y = np.exp(-((x - 5) ** 2) / 0.5) + np.exp(-((x - 15) ** 2) / 0.5)
        peaks = find_kde_peaks(x, y)
        assert len(peaks) >= 2

    def test_returns_empty_for_short_data(self):
        x = np.array([1, 2, 3])
        y = np.array([0.1, 0.5, 0.1])
        peaks = find_kde_peaks(x, y)
        assert peaks == []


class TestComputeGroupedKDE:
    def test_ungrouped_returns_all(self, sample_db):
        result = compute_grouped_kde(sample_db, column="readlen")
        assert "all" in result["groups"]
        assert len(result["groups"]["all"]["x"]) == 512

    def test_grouped_returns_per_barcode(self, sample_db):
        result = compute_grouped_kde(sample_db, column="readlen", group_by="bc_start_id")
        assert "nb05" in result["groups"]
        assert "nb10" in result["groups"]
        assert "nb05" in result["peaks"]
