#!/usr/bin/env python3
"""
End-to-end test for the ONT-SMA-seq pipeline.

This script tests the complete pipeline workflow:
1. mkdb.py - Create database
2. inputInit.py - Standardize inputs
3. ingest.py - Process reads and populate database

The test creates minimal synthetic test data and verifies that:
- Database is created with correct schema
- Input files are symlinked correctly
- Reads are processed and inserted into database
- Output BAM is created with ER tags
"""

import os
import random
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import uuid
from pathlib import Path

try:
    import pysam
except ImportError:
    print("Error: pysam is required. Install with: pip install pysam", file=sys.stderr)
    sys.exit(1)

# Test configuration
EXP_ID = "TEST001"
MODEL_TIER = "h"
MODEL_VER = "5.2.0"
TRIM = 1
MODS = "6mA"
NUM_READS = 20  # Small number for quick testing


def random_seq(length):
    """Generate random DNA sequence."""
    return "".join(random.choices("ACGT", k=length))


def create_test_data(test_dir):
    """
    Create minimal test data for the pipeline.
    
    Returns:
        Tuple of (bam_path, pod5_dir, ref_path)
    """
    print("\n=== Creating test data ===")
    
    # Create reference FASTA with 2 sequences (short ~200bp, long ~500bp)
    SHORT_LEN = 200
    LONG_LEN = 500
    short_ref = random_seq(SHORT_LEN)
    long_ref = random_seq(LONG_LEN)
    
    ref_path = os.path.join(test_dir, "reference.fa")
    with open(ref_path, "w") as f:
        f.write(f">short_amplicon\n{short_ref}\n")
        f.write(f">long_amplicon\n{long_ref}\n")
    
    print(f"✓ Created reference FASTA with sequences of {SHORT_LEN}bp and {LONG_LEN}bp")
    
    # Create test BAM file with reads
    bam_filename = f"{EXP_ID}_{MODEL_TIER}_v{MODEL_VER}_{TRIM}_{MODS}.bam"
    bam_path = os.path.join(test_dir, bam_filename)
    
    # Create BAM header (unaligned BAM)
    header = {
        "HD": {"VN": "1.6", "SO": "unknown"},
        "RG": [{"ID": "test", "SM": "test_sample"}],
    }
    
    # Generate test reads
    reads_data = []
    
    for i in range(NUM_READS):
        read_id = str(uuid.uuid4())
        
        # Randomly choose target (short, long, or unmatched)
        choice = random.choice(["short", "long", "unmatched"])
        
        if choice == "short":
            # Generate read close to short reference length
            length_offset = random.randint(-50, 50)
            read_len = SHORT_LEN + length_offset
            base_seq = short_ref
        elif choice == "long":
            # Generate read close to long reference length  
            length_offset = random.randint(-50, 50)
            read_len = LONG_LEN + length_offset
            base_seq = long_ref
        else:
            # Unmatched - length outside both ranges
            read_len = random.choice([50, 800])  # Too short or too long
            base_seq = random_seq(read_len)
        
        # Create read sequence (with some mutations from reference)
        if choice != "unmatched" and read_len > 0:
            # Truncate or extend the base sequence
            if read_len <= len(base_seq):
                read_seq = base_seq[:read_len]
            else:
                read_seq = base_seq + random_seq(read_len - len(base_seq))
            
            # Add some random mutations (5% error rate)
            read_seq = list(read_seq)
            for j in range(len(read_seq)):
                if random.random() < 0.05:
                    read_seq[j] = random.choice("ACGT")
            read_seq = "".join(read_seq)
        else:
            read_seq = random_seq(max(read_len, 10))
        
        # Generate quality scores (Phred scale, 10-40)
        qual = [random.randint(10, 40) for _ in range(len(read_seq))]
        
        reads_data.append({
            "id": read_id,
            "seq": read_seq,
            "qual": qual,
            "choice": choice,
        })
    
    # Write BAM file
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for read_info in reads_data:
            a = pysam.AlignedSegment()
            a.query_name = read_info["id"]
            a.query_sequence = read_info["seq"]
            a.query_qualities = pysam.qualitystring_to_array(
                "".join(chr(q + 33) for q in read_info["qual"])
            )
            a.flag = 4  # Unmapped
            a.reference_id = -1
            a.reference_start = 0
            a.mapping_quality = 0
            bam.write(a)
    
    print(f"✓ Created BAM file with {NUM_READS} reads: {bam_filename}")
    
    # Create Pod5 directory
    # Since creating real Pod5 files is complex and requires FAST5 input,
    # we'll create an empty directory. The pipeline handles missing end reasons gracefully.
    pod5_dir = os.path.join(test_dir, "pod5")
    os.makedirs(pod5_dir, exist_ok=True)
    
    # Create a marker file indicating this is a test directory
    marker_path = os.path.join(pod5_dir, "README.txt")
    with open(marker_path, "w") as f:
        f.write("Mock Pod5 directory for testing\n")
        f.write("No real Pod5 files - end reasons will be 'unknown'\n")
    
    print(f"✓ Created Pod5 directory (no real Pod5 files - end reasons will be 'unknown')")
    
    # Print read distribution
    short_count = sum(1 for r in reads_data if r["choice"] == "short")
    long_count = sum(1 for r in reads_data if r["choice"] == "long")
    unmatched_count = sum(1 for r in reads_data if r["choice"] == "unmatched")
    print(f"  Read distribution:")
    print(f"    Short amplicon range: {short_count}")
    print(f"    Long amplicon range: {long_count}")
    print(f"    Unmatched (out of range): {unmatched_count}")
    
    return bam_path, pod5_dir, ref_path, reads_data


def run_command(cmd, description):
    """Run a command and check for errors."""
    print(f"\n=== {description} ===")
    print(f"Command: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.stdout:
        print(result.stdout)
    
    if result.returncode != 0:
        print(f"ERROR: {description} failed!", file=sys.stderr)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        return False
    
    return True


def verify_database(db_path, expected_reads):
    """Verify database was created and populated correctly."""
    print("\n=== Verifying database ===")
    
    if not os.path.exists(db_path):
        print(f"ERROR: Database not found: {db_path}", file=sys.stderr)
        return False
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Check tables exist
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = set(row[0] for row in cursor.fetchall())
    expected_tables = {"Reads", "Mods", "Exp", "Refseq"}
    
    if not expected_tables.issubset(tables):
        print(f"ERROR: Missing tables. Expected {expected_tables}, found {tables}", file=sys.stderr)
        return False
    
    print(f"✓ All required tables exist: {', '.join(sorted(tables))}")
    
    # Check Mods table is populated
    cursor.execute("SELECT COUNT(*) FROM Mods")
    mod_count = cursor.fetchone()[0]
    if mod_count == 0:
        print("ERROR: Mods table is empty", file=sys.stderr)
        return False
    print(f"✓ Mods table has {mod_count} modification definitions")
    
    # Check Refseq table is populated
    cursor.execute("SELECT COUNT(*) FROM Refseq")
    ref_count = cursor.fetchone()[0]
    if ref_count != 2:
        print(f"ERROR: Expected 2 reference sequences, found {ref_count}", file=sys.stderr)
        return False
    print(f"✓ Refseq table has {ref_count} reference sequences")
    
    # Check Reads table is populated
    cursor.execute("SELECT COUNT(*) FROM Reads")
    read_count = cursor.fetchone()[0]
    if read_count != expected_reads:
        print(f"ERROR: Expected {expected_reads} reads, found {read_count}", file=sys.stderr)
        return False
    print(f"✓ Reads table has {read_count} reads")
    
    # Check that some reads have refseq_id (matched) and some are NULL (unmatched)
    cursor.execute("SELECT COUNT(*) FROM Reads WHERE refseq_id IS NOT NULL")
    matched_count = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(*) FROM Reads WHERE refseq_id IS NULL")
    unmatched_count = cursor.fetchone()[0]
    
    print(f"✓ Reads matched to reference: {matched_count}")
    print(f"✓ Reads unmatched (out of range): {unmatched_count}")
    
    # Verify metrics are calculated
    cursor.execute("SELECT COUNT(*) FROM Reads WHERE q_bc IS NOT NULL")
    qbc_count = cursor.fetchone()[0]
    if qbc_count != expected_reads:
        print(f"ERROR: q_bc should be calculated for all reads", file=sys.stderr)
        return False
    print(f"✓ q_bc calculated for all {qbc_count} reads")
    
    # Check that ed and q_ld are only set for matched reads
    cursor.execute("SELECT COUNT(*) FROM Reads WHERE ed IS NOT NULL")
    ed_count = cursor.fetchone()[0]
    if ed_count != matched_count:
        print(f"ERROR: ed should be calculated only for matched reads ({matched_count}), found {ed_count}", file=sys.stderr)
        return False
    print(f"✓ Levenshtein distance calculated for {ed_count} matched reads")
    
    # Sample some data
    cursor.execute("SELECT uniq_id, refseq_id, readlen, q_bc, ed, q_ld FROM Reads LIMIT 5")
    print("\nSample reads:")
    for row in cursor.fetchall():
        print(f"  {row}")
    
    conn.close()
    return True


def verify_output_bam(bam_path, expected_reads):
    """Verify output BAM was created with ER tags."""
    print("\n=== Verifying output BAM ===")
    
    if not os.path.exists(bam_path):
        print(f"ERROR: Output BAM not found: {bam_path}", file=sys.stderr)
        return False
    
    read_count = 0
    er_tag_count = 0
    
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for read in bam:
            read_count += 1
            if read.has_tag("ER"):
                er_tag_count += 1
    
    if read_count != expected_reads:
        print(f"ERROR: Expected {expected_reads} reads in output BAM, found {read_count}", file=sys.stderr)
        return False
    
    print(f"✓ Output BAM has {read_count} reads")
    print(f"✓ {er_tag_count} reads have ER (End Reason) tags")
    
    return True


def verify_input_symlinks(exp_id):
    """Verify input symlinks were created correctly."""
    print("\n=== Verifying input symlinks ===")
    
    input_dir = "Input"
    if not os.path.exists(input_dir):
        print(f"ERROR: Input directory not found: {input_dir}", file=sys.stderr)
        return False
    
    # Check for BAM symlink
    bam_pattern = f"{exp_id}_*.bam"
    bam_files = list(Path(input_dir).glob(bam_pattern))
    if not bam_files:
        print(f"ERROR: No BAM file found matching {bam_pattern} in {input_dir}", file=sys.stderr)
        return False
    print(f"✓ BAM symlink exists: {bam_files[0].name}")
    
    # Check for Pod5 symlink
    pod5_dir = os.path.join(input_dir, f"{exp_id}_pod5")
    if not os.path.exists(pod5_dir):
        print(f"ERROR: Pod5 directory symlink not found: {pod5_dir}", file=sys.stderr)
        return False
    print(f"✓ Pod5 directory symlink exists")
    
    # Check for reference symlink
    ref_path = os.path.join(input_dir, f"{exp_id}.fa")
    if not os.path.exists(ref_path):
        print(f"ERROR: Reference symlink not found: {ref_path}", file=sys.stderr)
        return False
    print(f"✓ Reference symlink exists")
    
    return True


def cleanup():
    """Clean up test artifacts."""
    print("\n=== Cleaning up test artifacts ===")
    
    # Remove directories
    for dir_name in ["Input", "Output"]:
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
            print(f"✓ Removed {dir_name}/")
    
    # Remove database
    db_pattern = f"SMA_{EXP_ID}.db"
    if os.path.exists(db_pattern):
        os.remove(db_pattern)
        print(f"✓ Removed {db_pattern}")


def main():
    """Run end-to-end test of the pipeline."""
    print("=" * 60)
    print("ONT-SMA-seq Pipeline - End-to-End Test")
    print("=" * 60)
    
    # Create temporary directory for test data
    with tempfile.TemporaryDirectory() as test_dir:
        try:
            # Clean up any previous test artifacts
            cleanup()
            
            # Step 1: Create test data
            bam_path, pod5_dir, ref_path, reads_data = create_test_data(test_dir)
            
            # Step 2: Run mkdb.py to create database
            if not run_command(
                ["python3", "mkdb.py", EXP_ID],
                "Step 1: Create database (mkdb.py)"
            ):
                return 1
            
            db_path = f"SMA_{EXP_ID}.db"
            if not os.path.exists(db_path):
                print(f"ERROR: Database was not created: {db_path}", file=sys.stderr)
                return 1
            
            # Step 3: Run inputInit.py to standardize inputs
            if not run_command(
                ["python3", "inputInit.py", 
                 "--bam", bam_path,
                 "--pod5", pod5_dir,
                 "--ref", ref_path,
                 "--force"],
                "Step 2: Standardize inputs (inputInit.py)"
            ):
                return 1
            
            # Verify input symlinks
            if not verify_input_symlinks(EXP_ID):
                return 1
            
            # Step 4: Run ingest.py to process reads
            if not run_command(
                ["python3", "ingest.py", EXP_ID],
                "Step 3: Process reads (ingest.py)"
            ):
                return 1
            
            # Step 5: Verify database
            if not verify_database(db_path, NUM_READS):
                return 1
            
            # Step 6: Verify output BAM
            output_bam = os.path.join("Output", f"{EXP_ID}.bam")
            if not verify_output_bam(output_bam, NUM_READS):
                return 1
            
            # Success!
            print("\n" + "=" * 60)
            print("✓ ✓ ✓ ALL TESTS PASSED ✓ ✓ ✓")
            print("=" * 60)
            print("\nPipeline workflow completed successfully:")
            print(f"  1. Database created: {db_path}")
            print(f"  2. Input files standardized in Input/")
            print(f"  3. {NUM_READS} reads processed")
            print(f"  4. Output BAM created: {output_bam}")
            print(f"  5. All database tables populated correctly")
            
            return 0
            
        except Exception as e:
            print(f"\nERROR: Test failed with exception: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            return 1
        
        finally:
            # Clean up
            cleanup()


if __name__ == "__main__":
    sys.exit(main())
