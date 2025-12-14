"""Tests for mkdb.py - Database Initialization Script"""

import sqlite3
import tempfile
from pathlib import Path
import sys
import os

# Add parent directory to path to import mkdb
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import mkdb


def test_create_tables():
    """Test database table creation."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / 'test.db'
        conn = sqlite3.connect(str(db_path))
        
        mkdb.create_tables(conn)
        
        # Verify tables exist
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = {row[0] for row in cursor.fetchall()}
        
        assert 'Exp' in tables
        assert 'Refseq' in tables
        assert 'Mods' in tables
        assert 'Reads' in tables
        
        conn.close()


def test_populate_mods_table():
    """Test populating the Mods table with modification bitflags."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / 'test.db'
        conn = sqlite3.connect(str(db_path))
        
        mkdb.create_tables(conn)
        mkdb.populate_mods_table(conn)
        
        # Verify modifications are populated
        cursor = conn.cursor()
        cursor.execute("SELECT mod_bitflag, mods FROM Mods ORDER BY mod_bitflag")
        mods = cursor.fetchall()
        
        assert len(mods) > 0
        assert (0, 'non') in mods
        assert (1, '6mA') in mods
        assert (2, '5mCG_5hmCG') in mods
        assert (4, '5mC_5hmC') in mods
        assert (8, '4mC_5mC') in mods
        assert (16, '5mC') in mods
        
        conn.close()


def test_main_creates_database():
    """Test main function creates database correctly."""
    with tempfile.TemporaryDirectory() as tmpdir:
        exp_id = 'test_exp'
        
        # Mock sys.argv
        sys.argv = ['mkdb.py', exp_id, '--output-dir', tmpdir]
        
        result = mkdb.main()
        
        assert result == 0
        
        db_path = Path(tmpdir) / f'SMA_{exp_id}.db'
        assert db_path.exists()
        
        # Verify database contents
        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()
        
        # Check tables exist
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = {row[0] for row in cursor.fetchall()}
        assert len(tables) == 4
        
        # Check Mods table is populated
        cursor.execute("SELECT COUNT(*) FROM Mods")
        count = cursor.fetchone()[0]
        assert count > 0
        
        conn.close()


def test_main_overwrites_existing_database():
    """Test main function overwrites existing database."""
    with tempfile.TemporaryDirectory() as tmpdir:
        exp_id = 'test_exp'
        db_path = Path(tmpdir) / f'SMA_{exp_id}.db'
        
        # Create an existing database
        db_path.touch()
        
        # Mock sys.argv
        sys.argv = ['mkdb.py', exp_id, '--output-dir', tmpdir]
        
        result = mkdb.main()
        
        assert result == 0
        assert db_path.exists()
        
        # Verify it's a valid database (not just an empty file)
        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = cursor.fetchall()
        assert len(tables) > 0
        conn.close()


if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-v'])
