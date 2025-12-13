#!/usr/bin/env python3
"""
ingest.py - Main Processing Script

Parses inputs, calculates metrics, tags BAMs, and populates the database.
This is the core data ingestion pipeline for ONT-SMA-seq.

Usage:
    python ingest.py <exp_id> [--db <database_path>]

Prerequisites:
    - Database must be created with mkdb.py
    - Input files must be standardized with inputInit.py

Input Structure:
    Input/
    ├── {exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam
    ├── {exp_id}_pod5/
    └── {exp_id}.fa

Output:
    Output/
    └── {exp_id}.bam  (with ER tags added)
"""

import argparse
import glob
import hashlib
import math
import os
import re
import sqlite3
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import pysam
except ImportError:
    print("Error: pysam is required. Install with: pip install pysam", file=sys.stderr)
    sys.exit(1)

try:
    import pod5
except ImportError:
    print("Error: pod5 is required. Install with: pip install pod5", file=sys.stderr)
    sys.exit(1)

try:
    import edlib
except ImportError:
    print("Error: edlib is required. Install with: pip install edlib", file=sys.stderr)
    sys.exit(1)


# Modification bitflag mapping
MOD_BITFLAGS = {
    "non": 0,
    "6mA": 1,
    "5mCG_5hmCG": 2,
    "5mC_5hmC": 4,
    "4mC_5mC": 8,
    "5mC": 16,
}

# Length tolerance for reference matching
LENGTH_TOLERANCE = 150


def parse_bam_filename(bam_path: str) -> Optional[Dict]:
    """
    Parse metadata from BAM filename.

    Expected format: {exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam
    """
    filename = os.path.basename(bam_path)

    if not filename.endswith(".bam"):
        return None

    base = filename[:-4]
    pattern = r"^(.+)_([shf])_v([\d.]+)_([01])_(.+)$"
    match = re.match(pattern, base)

    if not match:
        return None

    exp_id, model_tier, model_ver, trim, mods_str = match.groups()

    # Parse modification string to bitflag
    mod_bitflag = 0
    if mods_str != "non":
        for mod in mods_str.split("+"):
            mod = mod.strip()
            if mod in MOD_BITFLAGS:
                mod_bitflag += MOD_BITFLAGS[mod]

    return {
        "exp_id": exp_id,
        "model_tier": model_tier,
        "model_ver": model_ver,
        "trim": int(trim),
        "mods_str": mods_str,
        "mod_bitflag": mod_bitflag,
    }


def parse_reference_fasta(fasta_path: str) -> Dict[str, Dict]:
    """
    Parse reference FASTA file (expects exactly 2 sequences: long and short).

    Returns:
        Dictionary mapping defline to {seq, length, range_min, range_max}
    """
    references = {}

    with open(fasta_path, "r") as f:
        current_defline = None
        current_seq = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Save previous sequence
                if current_defline is not None:
                    seq = "".join(current_seq)
                    seq_len = len(seq)
                    references[current_defline] = {
                        "seq": seq,
                        "length": seq_len,
                        "range_min": seq_len - LENGTH_TOLERANCE,
                        "range_max": seq_len + LENGTH_TOLERANCE,
                    }

                # Start new sequence
                current_defline = line[1:].split()[0]  # Take first word after >
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget last sequence
        if current_defline is not None:
            seq = "".join(current_seq)
            seq_len = len(seq)
            references[current_defline] = {
                "seq": seq,
                "length": seq_len,
                "range_min": seq_len - LENGTH_TOLERANCE,
                "range_max": seq_len + LENGTH_TOLERANCE,
            }

    return references


def extract_pod5_metadata(pod5_dir: str) -> Dict[str, str]:
    """
    Extract end reason metadata from Pod5 files.

    Returns:
        Dictionary mapping read_id to end_reason string
    """
    end_reasons = {}
    pod5_files = glob.glob(os.path.join(pod5_dir, "*.pod5"))

    if not pod5_files:
        print(f"Warning: No .pod5 files found in {pod5_dir}", file=sys.stderr)
        return end_reasons

    for pod5_path in pod5_files:
        try:
            with pod5.Reader(pod5_path) as reader:
                for read in reader.reads():
                    read_id = str(read.read_id)
                    end_reason = read.end_reason.name if read.end_reason else "unknown"
                    end_reasons[read_id] = end_reason
        except Exception as e:
            print(f"Warning: Error reading {pod5_path}: {e}", file=sys.stderr)

    return end_reasons


def calculate_q_bc(quality_scores: List[int]) -> float:
    """
    Calculate probability-averaged basecall quality.

    q_bc = -10 * log10(mean(10^(-Q_base/10)))
    """
    if not quality_scores:
        return 0.0

    # Convert to error probabilities, average, then convert back
    error_probs = [10 ** (-q / 10) for q in quality_scores]
    mean_error = sum(error_probs) / len(error_probs)

    # Avoid log(0)
    if mean_error <= 0:
        return 60.0  # Cap at Q60

    q_bc = -10 * math.log10(mean_error)
    return round(q_bc, 4)


def calculate_levenshtein(read_seq: str, ref_seq: str) -> int:
    """
    Calculate Levenshtein distance between read and reference sequences.
    """
    result = edlib.align(read_seq, ref_seq, task="distance")
    return result["editDistance"]


def calculate_q_ld(edit_distance: int, ref_length: int) -> float:
    """
    Calculate Levenshtein quality score.

    q_ld = -10 * log10(min(max(1/L^2, ed/L), 1))
    """
    if ref_length <= 0:
        return 0.0

    L = ref_length
    ed_ratio = edit_distance / L
    min_ratio = 1 / (L * L)

    # max(1/L^2, ed/L)
    inner = max(min_ratio, ed_ratio)
    # min(..., 1)
    clamped = min(inner, 1.0)

    # Avoid log(0)
    if clamped <= 0:
        return 60.0

    q_ld = -10 * math.log10(clamped)
    return round(q_ld, 4)


def generate_unique_id(exp_id: str, model_tier: str, model_ver: str,
                       trim: int, mod_bitflag: int, read_id: str) -> str:
    """
    Generate unique identifier for a read.

    Format: {exp_id}{tier}{ver}t{trim}m{mod_flag}_{read_hash}
    """
    # Create short hash of read_id
    read_hash = hashlib.md5(read_id.encode()).hexdigest()[:8]

    # Clean version string (remove dots)
    ver_clean = model_ver.replace(".", "")

    return f"{exp_id}{model_tier}{ver_clean}t{trim}m{mod_bitflag}_{read_hash}"


def match_reference(read_length: int, references: Dict[str, Dict]) -> Optional[str]:
    """
    Match read to reference based on length.

    Returns:
        Reference ID (defline) if match found, None otherwise
    """
    for ref_id, ref_data in references.items():
        if ref_data["range_min"] <= read_length <= ref_data["range_max"]:
            return ref_id
    return None


def find_bam_file(input_dir: str, exp_id: str) -> Optional[str]:
    """
    Find the BAM file for the given experiment ID in the Input directory.
    """
    bam_pattern = os.path.join(input_dir, f"{exp_id}_*.bam")
    bam_files = glob.glob(bam_pattern)

    if not bam_files:
        return None

    # Return first match (should typically be only one)
    return bam_files[0]


def process_reads(exp_id: str, db_path: str, input_dir: str = "Input", output_dir: str = "Output"):
    """
    Main processing function that streams reads and populates the database.
    """
    # Validate paths
    input_path = Path(input_dir)
    if not input_path.exists():
        print(f"Error: Input directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)

    # Find BAM file
    bam_path = find_bam_file(input_dir, exp_id)
    if not bam_path:
        print(f"Error: No BAM file found for experiment {exp_id} in {input_dir}", file=sys.stderr)
        sys.exit(1)

    # Parse BAM filename for metadata
    metadata = parse_bam_filename(bam_path)
    if not metadata:
        print(f"Error: Could not parse BAM filename: {bam_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Processing experiment: {exp_id}")
    print(f"  BAM: {bam_path}")
    print(f"  Model: {metadata['model_tier']}_v{metadata['model_ver']}")
    print(f"  Trim: {metadata['trim']}, Mods: {metadata['mods_str']}")

    # Paths
    ref_path = os.path.join(input_dir, f"{exp_id}.fa")
    pod5_path = os.path.join(input_dir, f"{exp_id}_pod5")

    if not os.path.exists(ref_path):
        print(f"Error: Reference FASTA not found: {ref_path}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(pod5_path):
        print(f"Error: Pod5 directory not found: {pod5_path}", file=sys.stderr)
        sys.exit(1)

    # Parse reference sequences
    print("\nParsing reference sequences...")
    references = parse_reference_fasta(ref_path)
    print(f"  Found {len(references)} reference sequence(s)")
    for ref_id, ref_data in references.items():
        print(f"    {ref_id}: {ref_data['length']} bp (range: {ref_data['range_min']}-{ref_data['range_max']})")

    # Extract Pod5 metadata
    print("\nExtracting Pod5 metadata...")
    end_reasons = extract_pod5_metadata(pod5_path)
    print(f"  Found end reasons for {len(end_reasons)} reads")

    # Connect to database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Insert reference sequences into Refseq table
    print("\nPopulating Refseq table...")
    for ref_id, ref_data in references.items():
        cursor.execute(
            "INSERT OR REPLACE INTO Refseq (refseq_id, refseq, reflen) VALUES (?, ?, ?)",
            (ref_id, ref_data["seq"], ref_data["length"])
        )
    conn.commit()

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    output_bam_path = output_path / f"{exp_id}.bam"

    # Process reads
    print("\nProcessing reads...")
    processed = 0
    matched = 0
    unmatched = 0

    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as in_bam:
        with pysam.AlignmentFile(str(output_bam_path), "wb", header=in_bam.header) as out_bam:

            for read in in_bam:
                read_id = read.query_name
                read_seq = read.query_sequence or ""
                read_len = len(read_seq)
                quality_scores = list(read.query_qualities) if read.query_qualities else []

                # Get end reason from Pod5
                end_reason = end_reasons.get(read_id, "unknown")

                # Add ER tag to read
                read.set_tag("ER", end_reason, value_type="Z")

                # Write tagged read to output BAM
                out_bam.write(read)

                # Match to reference
                ref_id = match_reference(read_len, references)

                # Calculate metrics
                q_bc = calculate_q_bc(quality_scores)

                if ref_id:
                    # Calculate Levenshtein distance and quality
                    ref_seq = references[ref_id]["seq"]
                    ref_len = references[ref_id]["length"]
                    ed = calculate_levenshtein(read_seq, ref_seq)
                    q_ld = calculate_q_ld(ed, ref_len)
                    matched += 1
                else:
                    ed = None
                    q_ld = None
                    unmatched += 1

                # Generate unique ID
                uniq_id = generate_unique_id(
                    exp_id, metadata["model_tier"], metadata["model_ver"],
                    metadata["trim"], metadata["mod_bitflag"], read_id
                )

                # Extract additional read attributes (if available)
                # These are typically set in the BAM tags or can be derived
                channel = None
                well = None
                pore_type = None
                num_samples = None
                start_sample = None
                median_before = None
                scale = None
                offset = None
                forced = None

                # Try to get channel from tags
                try:
                    channel = read.get_tag("ch") if read.has_tag("ch") else None
                except:
                    pass

                # Insert into database
                cursor.execute("""
                    INSERT INTO Reads (
                        uniq_id, exp_id, refseq_id, read_id, readseq, readlen,
                        model_tier, model_ver, trim, mod_bitflag,
                        ed, q_bc, q_ld, ER, forced, channel, well, pore_type,
                        num_samples, start_sample, median_before, scale, offset
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    uniq_id, exp_id, ref_id, read_id, read_seq, read_len,
                    metadata["model_tier"], metadata["model_ver"],
                    metadata["trim"], metadata["mod_bitflag"],
                    ed, q_bc, q_ld, end_reason, forced, channel, well, pore_type,
                    num_samples, start_sample, median_before, scale, offset
                ))

                processed += 1
                if processed % 10000 == 0:
                    print(f"  Processed {processed} reads...")
                    conn.commit()

    # Final commit
    conn.commit()
    conn.close()

    print(f"\nProcessing complete!")
    print(f"  Total reads: {processed}")
    print(f"  Matched to reference: {matched}")
    print(f"  Unmatched (length out of range): {unmatched}")
    print(f"  Output BAM: {output_bam_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Process ONT-SMA-seq data and populate database",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python ingest.py EXP001
    python ingest.py EXP001 --db custom_path.db
    python ingest.py EXP001 --input-dir /data/Input --output-dir /data/Output
        """
    )
    parser.add_argument(
        "exp_id",
        type=str,
        help="Experiment identifier"
    )
    parser.add_argument(
        "--db",
        type=str,
        default=None,
        help="Path to SQLite database (default: SMA_{exp_id}.db)"
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default="Input",
        help="Input directory containing standardized files (default: Input)"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="Output",
        help="Output directory for processed BAM (default: Output)"
    )

    args = parser.parse_args()

    # Determine database path
    db_path = args.db or f"SMA_{args.exp_id}.db"

    if not os.path.exists(db_path):
        print(f"Error: Database not found: {db_path}", file=sys.stderr)
        print(f"Please run: python mkdb.py {args.exp_id}", file=sys.stderr)
        sys.exit(1)

    # Run processing
    process_reads(args.exp_id, db_path, args.input_dir, args.output_dir)


if __name__ == "__main__":
    main()
