#!/usr/bin/env python3
# ingest.py
# Main processing script: BAM tagging, Metric Calculation, DB Ingestion

import argparse
import sys
import math
import sqlite3
import pysam
import edlib
import pandas as pd
import numpy as np
from pathlib import Path

############
# argparse #
############

parser = argparse.ArgumentParser(description="Ingest BAM, calculate metrics, and populate DB")
parser.add_argument("-e", "--expid", required=True,
	help="Experiment ID (used for tagging and logic)")
parser.add_argument("-b", "--bam", required=True,
	help="Path to input uBAM file")
parser.add_argument("-s", "--summary", required=True,
	help="Path to input Pod5 summary TSV")
parser.add_argument("-r", "--ref", required=True,
	help="Path to Reference FASTA")
parser.add_argument("-d", "--database", required=True,
	help="Path to target SQLite database file")
parser.add_argument("-o", "--output_bam", required=True,
	help="Path to output tagged BAM file")
parser.add_argument("-k", "--tolerance", type=int, default=150,
	help="Length tolerance +/- bp for RefSeq matching [%(default)s]")

args = parser.parse_args()


##########
# config #
##########

EXP_ID = args.expid
INPUT_BAM = Path(args.bam)
SUMMARY_TSV = Path(args.summary)
REF_FASTA = Path(args.ref)
DB_PATH = Path(args.database)
OUTPUT_BAM = Path(args.output_bam)
TOLERANCE = args.tolerance

# Validate Inputs
if not INPUT_BAM.exists():
	sys.exit(f"[ingest] Error: BAM file {INPUT_BAM} not found.")
if not SUMMARY_TSV.exists():
	sys.exit(f"[ingest] Error: Summary TSV {SUMMARY_TSV} not found.")
if not REF_FASTA.exists():
	sys.exit(f"[ingest] Error: Reference FASTA {REF_FASTA} not found.")
if not DB_PATH.exists():
	sys.exit(f"[ingest] Error: Database {DB_PATH} not found. Run mkdb.py first.")

# Ensure Output directory exists
if not OUTPUT_BAM.parent.exists():
	print(f"[ingest] Creating output directory: {OUTPUT_BAM.parent}")
	OUTPUT_BAM.parent.mkdir(parents=True, exist_ok=True)

print(f"[ingest] Starting ingestion for Experiment: {EXP_ID}")
print(f"[ingest] Length Tolerance (-k): +/- {TOLERANCE} bp")


#################
# DB Connection #
#################

conn = sqlite3.connect(DB_PATH)
c = conn.cursor()
# Performance tuning for bulk inserts
c.execute("PRAGMA synchronous = OFF")
c.execute("PRAGMA journal_mode = MEMORY")


#####################
# Load Metadata TSV #
#####################

print(f"[ingest] Loading metadata from {SUMMARY_TSV}...")
# Use pandas for fast C-engine parsing
meta_df = pd.read_csv(SUMMARY_TSV, sep='\t')

# QC: Check uniqueness
if not meta_df['read_id'].is_unique:
	print(f"[ingest] Warning: Found {meta_df['read_id'].duplicated().sum()} duplicates in TSV. Dropping duplicates.")
	meta_df = meta_df.drop_duplicates(subset='read_id')

# Convert to Dictionary for O(1) Lookup: {read_id: end_reason}
meta_lookup = pd.Series(
	meta_df.end_reason.values,
	index=meta_df.read_id
).to_dict()

print(f"[ingest] Loaded metadata for {len(meta_lookup)} reads.")


#########################
# Parse Reference FASTA #
#########################

print(f"[ingest] Parsing Reference FASTA: {REF_FASTA}")

# Store RefSeq info: {id: {'seq': str, 'len': int, 'range': (min, max)}}
ref_data = {}

with pysam.FastaFile(str(REF_FASTA)) as fa:
	for ref_name in fa.references:
		seq = fa.fetch(ref_name).upper()
		length = len(seq)

		# Range logic: Length +/- TOLERANCE
		min_len = length - TOLERANCE
		max_len = length + TOLERANCE

		ref_data[ref_name] = {
			'seq': seq,
			'len': length,
			'range': (min_len, max_len)
		}

		# DB Action: Insert into Refseq table
		c.execute("INSERT OR REPLACE INTO Refseq (refseq_id, refseq, reflen) VALUES (?, ?, ?)",
				  (ref_name, seq, length))

conn.commit()
print(f"[ingest] Loaded {len(ref_data)} references.")


###########################
# Parse BAM Filename Meta #
###########################

print(f"[ingest] Processing BAM: {INPUT_BAM.name}")

filename_root = INPUT_BAM.stem # removes .bam extension
parts = filename_root.split("_")

# Convention: ... {tier}_{ver}_{trim}_{mods}
try:
	mods = int(parts[-1])
	trim = int(parts[-2][-1])
	ver_str = parts[-3].lstrip("v") # "5.2.0"
	tier = parts[-4]

	# Clean version string for ID generation (remove dots)
	ver_clean = ver_str.replace(".", "")

except (IndexError, ValueError) as e:
	sys.exit(f"[ingest] Error parsing metadata from filename {INPUT_BAM.name}: {e}")

print(f"[ingest] BAM metadata:\ntier - [{tier}]\nver  - [{ver_str}]\ntrim - [{trim}]\nmods - [{mods}]")

print(f"[ingest] Processed BAM filename metadata")

########################
# Main Processing Loop #
########################

# Prepared Statements for Speed
insert_sql = '''
	INSERT INTO Reads (
		uniq_id, exp_id, refseq_id, read_id, readseq, readlen,
		model_tier, model_ver, trim, mod_bitflag,
		ed, q_bc, q_ld, ER
	) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
'''

# Counters
count_total = 0
count_tagged = 0
count_mapped = 0

print(f"[ingest] Ingesting Reads into {DB_PATH} and writing tagged BAM to {OUTPUT_BAM}...")

with pysam.AlignmentFile(INPUT_BAM, "rb", check_sq=False) as in_bam, \
	 pysam.AlignmentFile(OUTPUT_BAM, "wb", template=in_bam) as out_bam:

	for read in in_bam:
		count_total += 1

		# A. Get basic read info
		r_id = read.query_name
		seq = read.query_sequence

		if not seq:
			continue

		r_len = len(seq)

		# B. End Reason Lookup & Tagging
		er_val = meta_lookup.get(r_id, "unknown")
		read.set_tag("ER", er_val, value_type="Z")
		count_tagged += 1

		# C. RefSeq Logic (Length-based matching)
		assigned_ref_id = None

		# Check against loaded references
		for ref_id, data in ref_data.items():
			min_l, max_l = data['range']
			if min_l <= r_len <= max_l:
				assigned_ref_id = ref_id
				break

		# D. Metrics Calculation

		# 1. q_bc (Basecall Quality)
		q_scores = read.query_qualities
		if q_scores and len(q_scores) > 0:
			probs = np.power(10, -np.array(q_scores) / 10.0)
			avg_prob = np.mean(probs)
			if avg_prob > 0:
				q_bc = -10 * math.log10(avg_prob)
			else:
				q_bc = 0.0
		else:
			q_bc = 0.0

		# 2. Alignment Metrics (ed, q_ld)
		ed = None
		q_ld = None

		if assigned_ref_id:
			count_mapped += 1
			ref_seq = ref_data[assigned_ref_id]['seq']
			ref_len = ref_data[assigned_ref_id]['len']

			result = edlib.align(seq, ref_seq, mode="NW", task="distance")
			ed = result["editDistance"]

			# q_ld calculation
			L = float(ref_len)
			term1 = 1.0 / (L * L)
			term2 = float(ed) / L
			inner_max = max(term1, term2)
			bounded_val = min(inner_max, 1.0)
			q_ld = -10 * math.log10(bounded_val)

		# E. Unique ID Construction
		# Format: {exp_id}{tier}{ver}t{trim}m{mod_flag}_{read_id}
		uniq_id = f"{EXP_ID}{tier}{ver_clean}t{trim}m{mods}_{r_id}"

		# F. DB Insertion
		row_data = (
			uniq_id, EXP_ID, assigned_ref_id, r_id, seq, r_len,
			tier, ver_str, trim, mods,
			ed, q_bc, q_ld, er_val
		)

		c.execute(insert_sql, row_data)
		out_bam.write(read)

# Cleanup
conn.commit()
conn.close()

print(f"[ingest] Ingested Reads into {DB_PATH} and wrote tagged BAM to {OUTPUT_BAM}")

print(f"[ingest] Complete.")
print(f"  - Total Reads:  {count_total}")
print(f"  - Tagged:       {count_tagged}")
print(f"  - Ref Matches:  {count_mapped}")
print(f"  - Output DB:    {DB_PATH}")
print(f"  - Output BAM:   {OUTPUT_BAM}")
