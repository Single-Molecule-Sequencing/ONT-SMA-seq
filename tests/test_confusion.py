"""Tests for barcode confusion matrix computation."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from calibrate_viz.confusion import (
    compute_confusion_matrix,
    compute_threshold_impact,
    get_affected_reads,
)


@pytest.fixture
def classified_db(tmp_path: Path) -> Path:
    """DB with reads that have known correct/incorrect barcode assignments."""
    db_path = tmp_path / "test.db"
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE Reads (
            read_id TEXT PRIMARY KEY,
            readlen INTEGER,
            ed INTEGER,
            bc_start_id TEXT,
            bc_start_conf REAL,
            bc_end_id TEXT,
            bc_end_conf REAL,
            tgt_id TEXT,
            trunc_level TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE Target (
            tgt_id TEXT PRIMARY KEY,
            tgt_reflen INTEGER
        )
    """)
    conn.execute("INSERT INTO Target VALUES ('V04_2_fwd', 500)")
    conn.execute("INSERT INTO Target VALUES ('V04_2_rev', 500)")
    # Correctly classified reads (low ed relative to target length)
    conn.execute("INSERT INTO Reads VALUES ('r1', 500, 5, 'nb05', 0.9, 'nb10', 0.8, 'V04_2_fwd', 'full_length')")
    conn.execute("INSERT INTO Reads VALUES ('r2', 490, 8, 'nb05', 0.85, 'nb10', 0.7, 'V04_2_fwd', 'full_length')")
    # Misclassified read (high ed = wrong barcode assignment likely)
    conn.execute("INSERT INTO Reads VALUES ('r3', 480, 200, 'nb05', 0.4, 'nb10', 0.3, 'V04_2_fwd', 'bc1_target')")
    # Different target
    conn.execute("INSERT INTO Reads VALUES ('r4', 500, 3, 'nb10', 0.92, 'nb05', 0.85, 'V04_2_rev', 'full_length')")
    conn.commit()
    conn.close()
    return db_path


class TestComputeConfusionMatrix:
    def test_returns_matrix_dict(self, classified_db: Path):
        matrix = compute_confusion_matrix(classified_db, ed_threshold=0.1)
        assert isinstance(matrix, dict)
        assert "labels" in matrix
        assert "counts" in matrix

    def test_correct_assignments_on_diagonal(self, classified_db: Path):
        matrix = compute_confusion_matrix(classified_db, ed_threshold=0.1)
        total_diagonal = sum(
            matrix["counts"][i][i] for i in range(len(matrix["labels"]))
        )
        assert total_diagonal >= 3


class TestComputeThresholdImpact:
    def test_returns_counts_per_trunc_level(self, classified_db: Path):
        result = compute_threshold_impact(
            classified_db,
            start_barcode_min=0.6,
            full_length_threshold=0.75,
        )
        assert isinstance(result, dict)
        # Should have some classification levels
        total = sum(result.values())
        assert total == 4  # 4 reads total

    def test_high_threshold_reduces_full_length(self, classified_db: Path):
        low = compute_threshold_impact(classified_db, 0.5, 0.7)
        high = compute_threshold_impact(classified_db, 0.5, 0.95)
        # Higher full_length_threshold means fewer full_length reads
        low_fl = low.get("full_length", 0)
        high_fl = high.get("full_length", 0)
        assert high_fl <= low_fl


class TestGetAffectedReads:
    def test_finds_changed_reads(self, classified_db: Path):
        affected = get_affected_reads(
            classified_db,
            old_start_min=0.6, old_fl_thresh=0.75,
            new_start_min=0.6, new_fl_thresh=0.95,
        )
        # Some reads should change classification
        assert isinstance(affected, list)
        for r in affected:
            assert "read_id" in r
            assert r["old_level"] != r["new_level"]
