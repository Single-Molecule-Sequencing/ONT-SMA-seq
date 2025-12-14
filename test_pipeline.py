#!/usr/bin/env python3
"""
Comprehensive test suite for the ONT-SMA-seq pipeline.

Tests all three main scripts (mkdb.py, inputInit.py, ingest.py) and
verifies the complete pipeline workflow.
"""

import os
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path


class TestResult:
    """Simple test result tracker."""
    
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.tests = []
    
    def record(self, name, passed, message=""):
        """Record a test result."""
        self.tests.append((name, passed, message))
        if passed:
            self.passed += 1
            print(f"✓ {name}")
        else:
            self.failed += 1
            print(f"✗ {name}: {message}")
    
    def summary(self):
        """Print test summary."""
        total = self.passed + self.failed
        print(f"\n{'='*70}")
        print(f"Test Results: {self.passed}/{total} passed")
        if self.failed > 0:
            print(f"\nFailed tests:")
            for name, passed, msg in self.tests:
                if not passed:
                    print(f"  - {name}: {msg}")
        print(f"{'='*70}")
        return self.failed == 0


def run_command(cmd, check=True):
    """Run a shell command and return the result."""
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True
    )
    if check and result.returncode != 0:
        print(f"Command failed: {cmd}")
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        raise RuntimeError(f"Command failed with code {result.returncode}")
    return result


def cleanup_test_artifacts():
    """Clean up test artifacts from previous runs."""
    artifacts = [
        "SMA_TEST001.db",
        "Input",
        "Output",
        "test_data",
    ]
    for artifact in artifacts:
        if os.path.exists(artifact):
            if os.path.isdir(artifact):
                shutil.rmtree(artifact)
            else:
                os.remove(artifact)


def test_mkdb(results):
    """Test mkdb.py script."""
    print("\n--- Testing mkdb.py ---")
    
    # Test 1: Database creation
    try:
        result = run_command("python mkdb.py TEST001")
        db_exists = os.path.exists("SMA_TEST001.db")
        results.record(
            "mkdb: Database file created",
            db_exists,
            "Database file not found" if not db_exists else ""
        )
    except Exception as e:
        results.record("mkdb: Database file created", False, str(e))
        return
    
    # Test 2: Check database schema
    try:
        conn = sqlite3.connect("SMA_TEST001.db")
        cursor = conn.cursor()
        
        # Check tables exist
        cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )
        tables = {row[0] for row in cursor.fetchall()}
        expected_tables = {"Reads", "Mods", "Exp", "Refseq"}
        has_all_tables = expected_tables.issubset(tables)
        
        results.record(
            "mkdb: All required tables created",
            has_all_tables,
            f"Missing tables: {expected_tables - tables}" if not has_all_tables else ""
        )
        
        # Check Mods table is populated
        cursor.execute("SELECT COUNT(*) FROM Mods")
        mod_count = cursor.fetchone()[0]
        results.record(
            "mkdb: Mods table populated",
            mod_count > 0,
            f"Mods table has {mod_count} entries" if mod_count > 0 else "Mods table is empty"
        )
        
        conn.close()
    except Exception as e:
        results.record("mkdb: Database schema check", False, str(e))


def test_inputInit(results):
    """Test inputInit.py script."""
    print("\n--- Testing inputInit.py ---")
    
    # First create test data
    try:
        run_command("python test_setup.py")
    except Exception as e:
        results.record("inputInit: Test data creation", False, str(e))
        return
    
    # Test 1: Input standardization
    try:
        result = run_command(
            "python inputInit.py --bam test_data/TEST001_h_v5.2.0_1_6mA.bam "
            "--pod5 test_data/pod5 --ref test_data/reference.fa"
        )
        
        # Check Input directory created
        input_exists = os.path.exists("Input")
        results.record(
            "inputInit: Input directory created",
            input_exists,
            "Input directory not found" if not input_exists else ""
        )
        
        # Check symlinks created
        expected_links = [
            "Input/TEST001_h_v5.2.0_1_6mA.bam",
            "Input/TEST001_pod5",
            "Input/TEST001.fa"
        ]
        
        for link in expected_links:
            exists = os.path.exists(link)
            results.record(
                f"inputInit: Created {Path(link).name}",
                exists,
                f"{link} not found" if not exists else ""
            )
    
    except Exception as e:
        results.record("inputInit: Input standardization", False, str(e))


def test_ingest(results):
    """Test ingest.py script."""
    print("\n--- Testing ingest.py ---")
    
    # Test 1: Read processing
    try:
        result = run_command("python ingest.py TEST001")
        
        # Check Output directory created
        output_exists = os.path.exists("Output")
        results.record(
            "ingest: Output directory created",
            output_exists,
            "Output directory not found" if not output_exists else ""
        )
        
        # Check output BAM created
        output_bam = "Output/TEST001.bam"
        bam_exists = os.path.exists(output_bam)
        results.record(
            "ingest: Output BAM created",
            bam_exists,
            f"{output_bam} not found" if not bam_exists else ""
        )
        
    except Exception as e:
        results.record("ingest: Read processing", False, str(e))
        return
    
    # Test 2: Database population
    try:
        conn = sqlite3.connect("SMA_TEST001.db")
        cursor = conn.cursor()
        
        # Check reads were inserted
        cursor.execute("SELECT COUNT(*) FROM Reads")
        read_count = cursor.fetchone()[0]
        results.record(
            "ingest: Reads inserted into database",
            read_count > 0,
            f"Found {read_count} reads" if read_count > 0 else "No reads in database"
        )
        
        # Check reference sequences populated
        cursor.execute("SELECT COUNT(*) FROM Refseq")
        refseq_count = cursor.fetchone()[0]
        results.record(
            "ingest: Reference sequences populated",
            refseq_count == 2,
            f"Expected 2 refseqs, found {refseq_count}"
        )
        
        # Check reference matching
        cursor.execute(
            "SELECT refseq_id, COUNT(*) as cnt FROM Reads "
            "GROUP BY refseq_id"
        )
        ref_distribution = dict(cursor.fetchall())
        has_matches = any(k for k in ref_distribution.keys() if k is not None and k != '')
        results.record(
            "ingest: Reference matching performed",
            has_matches,
            f"Reference distribution: {ref_distribution}"
        )
        
        # Check quality metrics calculated
        cursor.execute(
            "SELECT COUNT(*) FROM Reads WHERE q_bc IS NOT NULL AND q_bc > 0"
        )
        qbc_count = cursor.fetchone()[0]
        results.record(
            "ingest: Quality metrics (q_bc) calculated",
            qbc_count > 0,
            f"{qbc_count} reads have q_bc scores"
        )
        
        # Check Levenshtein metrics for matched reads
        cursor.execute(
            "SELECT COUNT(*) FROM Reads WHERE refseq_id IS NOT NULL "
            "AND refseq_id != '' AND ed IS NOT NULL"
        )
        ed_count = cursor.fetchone()[0]
        results.record(
            "ingest: Levenshtein distance calculated for matched reads",
            ed_count > 0,
            f"{ed_count} matched reads have edit distance"
        )
        
        conn.close()
        
    except Exception as e:
        results.record("ingest: Database validation", False, str(e))


def test_integration(results):
    """Test complete pipeline integration."""
    print("\n--- Testing Pipeline Integration ---")
    
    try:
        conn = sqlite3.connect("SMA_TEST001.db")
        cursor = conn.cursor()
        
        # Verify data integrity
        cursor.execute(
            "SELECT r.uniq_id, r.readlen, r.q_bc, r.ed, rs.reflen "
            "FROM Reads r "
            "LEFT JOIN Refseq rs ON r.refseq_id = rs.refseq_id "
            "LIMIT 1"
        )
        sample_read = cursor.fetchone()
        
        results.record(
            "integration: Database contains complete read data",
            sample_read is not None,
            "No reads with complete data found"
        )
        
        # Check unique ID format
        if sample_read:
            uniq_id = sample_read[0]
            has_valid_format = (
                "TEST001" in uniq_id and
                "h" in uniq_id and
                "520" in uniq_id
            )
            results.record(
                "integration: Unique ID format correct",
                has_valid_format,
                f"Sample ID: {uniq_id}"
            )
        
        # Verify foreign key relationships
        cursor.execute(
            "SELECT COUNT(*) FROM Reads r "
            "LEFT JOIN Mods m ON r.mod_bitflag = m.mod_bitflag "
            "WHERE m.mod_bitflag IS NULL"
        )
        orphan_mods = cursor.fetchone()[0]
        results.record(
            "integration: Modification foreign keys valid",
            orphan_mods == 0,
            f"Found {orphan_mods} reads with invalid mod_bitflag"
        )
        
        conn.close()
        
    except Exception as e:
        results.record("integration: Data integrity checks", False, str(e))


def main():
    """Run all tests."""
    print("=" * 70)
    print("ONT-SMA-seq Pipeline Test Suite")
    print("=" * 70)
    
    # Clean up from previous runs
    cleanup_test_artifacts()
    
    # Run tests
    results = TestResult()
    
    test_mkdb(results)
    test_inputInit(results)
    test_ingest(results)
    test_integration(results)
    
    # Print summary
    success = results.summary()
    
    # Clean up after tests
    print("\nCleaning up test artifacts...")
    cleanup_test_artifacts()
    
    # Exit with appropriate code
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
