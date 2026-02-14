#!/usr/bin/env python3
"""Tests for mkdb.py database schema."""

import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.fixture
def db_path(tmp_path):
    """Run mkdb.py and return the created DB path."""
    result = subprocess.run(
        [sys.executable, "bin/mkdb.py", "-e", "FAL12345_20260129_IF", "-o", str(tmp_path)],
        capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
    )
    assert result.returncode == 0, f"mkdb.py failed: {result.stderr}"
    return tmp_path / "SMA_FAL12345_20260129_IF.db"


class TestReadsSchema:

    def test_db_created(self, db_path):
        assert db_path.exists()

    def test_barcode_columns_exist(self, db_path):
        conn = sqlite3.connect(db_path)
        cursor = conn.execute("PRAGMA table_info(Reads)")
        columns = {row[1] for row in cursor.fetchall()}
        conn.close()

        expected_new = {
            "bc_start_id", "bc_start_ed", "bc_start_conf",
            "bc_end_id", "bc_end_ed", "bc_end_conf",
        }
        assert expected_new.issubset(columns), (
            f"Missing columns: {expected_new - columns}"
        )

    def test_barcode_column_types(self, db_path):
        conn = sqlite3.connect(db_path)
        cursor = conn.execute("PRAGMA table_info(Reads)")
        col_types = {row[1]: row[2] for row in cursor.fetchall()}
        conn.close()

        assert col_types["bc_start_id"] == "TEXT"
        assert col_types["bc_start_ed"] == "INTEGER"
        assert col_types["bc_start_conf"] == "REAL"
        assert col_types["bc_end_id"] == "TEXT"
        assert col_types["bc_end_ed"] == "INTEGER"
        assert col_types["bc_end_conf"] == "REAL"

    def test_original_columns_still_present(self, db_path):
        conn = sqlite3.connect(db_path)
        cursor = conn.execute("PRAGMA table_info(Reads)")
        columns = {row[1] for row in cursor.fetchall()}
        conn.close()

        original = {
            "uniq_id", "exp_id", "tgt_id", "read_id", "readseq",
            "readlen", "model_tier", "model_ver", "trim", "mod_bitflag",
            "ed", "q_bc", "q_ld", "ER",
        }
        assert original.issubset(columns)
