#!/usr/bin/env python3
"""
Create test data for the ONT-SMA-seq pipeline.
"""

import os
import random
import string
import uuid

import pysam

# Test experiment ID
EXP_ID = "TEST001"
MODEL_TIER = "h"
MODEL_VER = "5.2.0"
TRIM = 1
MODS = "6mA"

# Create test directories
os.makedirs("test_data/pod5", exist_ok=True)

# Generate random DNA sequence
def random_seq(length):
    return "".join(random.choices("ACGT", k=length))

# Create reference FASTA with 2 sequences (short ~200bp, long ~500bp)
SHORT_LEN = 200
LONG_LEN = 500
short_ref = random_seq(SHORT_LEN)
long_ref = random_seq(LONG_LEN)

with open("test_data/reference.fa", "w") as f:
    f.write(f">short_amplicon\n{short_ref}\n")
    f.write(f">long_amplicon\n{long_ref}\n")

print(f"Created reference FASTA with sequences of {SHORT_LEN}bp and {LONG_LEN}bp")

# Create test BAM file with reads
bam_filename = f"{EXP_ID}_{MODEL_TIER}_v{MODEL_VER}_{TRIM}_{MODS}.bam"
bam_path = f"test_data/{bam_filename}"

# Create BAM header (unaligned BAM)
header = {
    "HD": {"VN": "1.6", "SO": "unknown"},
    "RG": [{"ID": "test", "SM": "test_sample"}],
}

# Generate test reads
NUM_READS = 50
reads_data = []

for i in range(NUM_READS):
    read_id = str(uuid.uuid4())

    # Randomly choose target (short or long) or make it unmatched
    choice = random.choice(["short", "long", "unmatched"])

    if choice == "short":
        # Generate read close to short reference length
        length_offset = random.randint(-100, 100)
        read_len = SHORT_LEN + length_offset
        base_seq = short_ref
    elif choice == "long":
        # Generate read close to long reference length
        length_offset = random.randint(-100, 100)
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

print(f"Created BAM file with {NUM_READS} reads: {bam_path}")

# Create a mock Pod5 "directory" with a placeholder
# Since creating real Pod5 files is complex, we'll create a simple marker
# and modify ingest.py to handle empty Pod5 gracefully
pod5_marker = "test_data/pod5/.pod5_test_marker"
with open(pod5_marker, "w") as f:
    f.write("# Mock Pod5 directory for testing\n")
    f.write("# Read IDs:\n")
    for read_info in reads_data:
        f.write(f"{read_info['id']}\n")

print(f"Created Pod5 directory marker (no real Pod5 files - end reasons will be 'unknown')")

# Create a simple Pod5 file using the pod5 library
try:
    import pod5
    from pod5.tools import pod5_convert  # Check if available
    print("Pod5 library available - but creating mock Pod5 files is complex")
    print("The pipeline will handle missing end reasons gracefully")
except Exception as e:
    print(f"Note: Could not create real Pod5 files: {e}")

print(f"\nTest data created in test_data/")
print(f"BAM filename follows convention: {bam_filename}")
print("\nRead distribution:")
short_count = sum(1 for r in reads_data if r["choice"] == "short")
long_count = sum(1 for r in reads_data if r["choice"] == "long")
unmatched_count = sum(1 for r in reads_data if r["choice"] == "unmatched")
print(f"  Short amplicon range: {short_count}")
print(f"  Long amplicon range: {long_count}")
print(f"  Unmatched (out of range): {unmatched_count}")
