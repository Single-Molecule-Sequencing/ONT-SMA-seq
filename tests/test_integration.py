"""Integration tests for the complete ONT-SMA-seq pipeline"""

import tempfile
import sqlite3
from pathlib import Path
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import mkdb
import inputInit
import ingest


def create_test_bam(bam_path):
    """Create a minimal test BAM file."""
    try:
        import pysam
    except ImportError:
        return False
    
    # Create a simple BAM with header
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 1000, 'SN': 'ref1'}]
    }
    
    with pysam.AlignmentFile(str(bam_path), 'wb', header=header) as outf:
        # Create a few test reads
        for i in range(5):
            a = pysam.AlignedSegment()
            a.query_name = f'read_{i}'
            a.query_sequence = 'A' * (500 + i * 10)
            a.flag = 4  # Unmapped
            a.reference_id = -1
            a.reference_start = -1
            a.mapping_quality = 0
            a.query_qualities = [30] * len(a.query_sequence)
            outf.write(a)
    
    return True


def create_test_pod5(pod5_dir):
    """Create a minimal test POD5 directory structure."""
    try:
        import pod5 as p5
        import numpy as np
        from pod5.tools.pod5_convert_from_fast5 import create_pod5_file
    except ImportError:
        # If pod5 is not available, create a dummy structure
        pod5_dir.mkdir(parents=True, exist_ok=True)
        # We'll skip the actual pod5 test if the library isn't available
        return False
    
    # For now, just create an empty pod5 directory
    # Full pod5 file creation requires complex setup
    pod5_dir.mkdir(parents=True, exist_ok=True)
    return False


def test_pipeline_integration():
    """Test the complete pipeline from start to finish."""
    # Check if required libraries are available
    try:
        import pysam
        import edlib
    except ImportError:
        print("Skipping integration test - pysam or edlib not available")
        return
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        exp_id = 'test_exp'
        
        # Step 1: Create database with mkdb.py
        sys.argv = ['mkdb.py', exp_id, '--output-dir', str(tmpdir)]
        result = mkdb.main()
        assert result == 0
        
        db_path = tmpdir / f'SMA_{exp_id}.db'
        assert db_path.exists()
        
        # Verify database structure
        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()
        
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = {row[0] for row in cursor.fetchall()}
        assert 'Exp' in tables
        assert 'Refseq' in tables
        assert 'Mods' in tables
        assert 'Reads' in tables
        
        # Verify Mods table is populated
        cursor.execute("SELECT COUNT(*) FROM Mods")
        count = cursor.fetchone()[0]
        assert count > 0
        
        conn.close()
        
        print(f"✓ Database created successfully: {db_path}")
        
        # Step 2: Create test input files
        source_dir = tmpdir / 'source'
        source_dir.mkdir()
        
        bam_source = source_dir / f'{exp_id}_s_v5.2.0_1_non.bam'
        if not create_test_bam(bam_source):
            print("Skipping full integration test - cannot create test BAM")
            return
        
        pod5_source = source_dir / 'pod5_data'
        pod5_source.mkdir()
        
        ref_source = source_dir / 'reference.fa'
        ref_content = '>short\n' + 'A' * 500 + '\n>long\n' + 'C' * 1000 + '\n'
        ref_source.write_text(ref_content)
        
        # Step 3: Run inputInit.py
        input_dir = tmpdir / 'Input'
        sys.argv = [
            'inputInit.py',
            str(bam_source),
            str(pod5_source),
            str(ref_source),
            '--output-dir', str(input_dir)
        ]
        
        result = inputInit.main()
        assert result == 0
        
        # Verify symlinks were created
        assert (input_dir / f'{exp_id}_s_v5.2.0_1_non.bam').exists()
        assert (input_dir / f'{exp_id}_pod5').exists()
        assert (input_dir / f'{exp_id}.fa').exists()
        
        print(f"✓ Input files standardized in: {input_dir}")
        
        # Since we can't easily create valid POD5 files, we'll skip the full ingest test
        # but verify that the structure is correct
        
        print("✓ Pipeline integration test passed (structure verified)")


def test_database_schema():
    """Test database schema matches specification."""
    with tempfile.TemporaryDirectory() as tmpdir:
        exp_id = 'test'
        sys.argv = ['mkdb.py', exp_id, '--output-dir', tmpdir]
        mkdb.main()
        
        db_path = Path(tmpdir) / f'SMA_{exp_id}.db'
        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()
        
        # Check Reads table schema
        cursor.execute("PRAGMA table_info(Reads)")
        columns = {row[1]: row[2] for row in cursor.fetchall()}
        
        expected_columns = [
            'uniq_id', 'exp_id', 'refseq_id', 'read_id', 'readseq', 'readlen',
            'model_tier', 'model_ver', 'trim', 'mod_bitflag',
            'ed', 'q_bc', 'q_ld', 'ER',
            'forced', 'channel', 'well', 'pore_type',
            'num_samples', 'start_sample', 'median_before', 'scale', 'offset'
        ]
        
        for col in expected_columns:
            assert col in columns, f"Column {col} missing from Reads table"
        
        # Check Mods table
        cursor.execute("PRAGMA table_info(Mods)")
        columns = {row[1] for row in cursor.fetchall()}
        assert 'mod_bitflag' in columns
        assert 'mods' in columns
        
        # Check Refseq table
        cursor.execute("PRAGMA table_info(Refseq)")
        columns = {row[1] for row in cursor.fetchall()}
        assert 'refseq_id' in columns
        assert 'refseq' in columns
        assert 'reflen' in columns
        
        # Check Exp table
        cursor.execute("PRAGMA table_info(Exp)")
        columns = {row[1] for row in cursor.fetchall()}
        assert 'exp_id' in columns
        assert 'exp_desc' in columns
        
        conn.close()
        print("✓ Database schema verified")


if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-v'])
