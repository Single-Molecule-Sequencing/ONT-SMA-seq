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
from typing import Tuple, Any


#########
# funcs #
#########

def parse_bam_filename(bam_path: Path) -> Tuple[str, str, int, int]:
	"""
	Parses metadata from BAM filename (Right-to-Left parsing).
	Expected Format: {exp_id}_{tier}_v{ver}_{trim}_{mods}.bam
	"""
	filename = bam_path.resolve().stem
	parts = filename.split('_')

	if len(parts) < 5:
		sys.exit(f"[ingest] Error: Metadata parse failed for '{filename}'.")

	try:
		# Parse from right to left to handle underscores in exp_id
		mods = int(parts[-1])
		trim = int(parts[-2][-1])  # 'trim0' -> 0
		ver_str = parts[-3].lstrip("v")
		tier = parts[-4]
		return tier, ver_str, trim, mods
	except Exception as e:
		sys.exit(f"[ingest] Error: Metadata parse failed for '{filename}'. {e}")


def calculate_q_bc(quality_scores: Any) -> float:
	"""Calculates probability-averaged Phred quality score."""
	if len(quality_scores) == 0:
		return 0.0

	q_arr = np.array(quality_scores, dtype=np.float64)
	probs = np.power(10, -q_arr / 10.0)
	avg_prob = np.mean(probs)

	if avg_prob > 0:
		return -10 * math.log10(avg_prob)
	return 0.0


def calculate_q_ld(ed: int, ref_len: int) -> float:
	"""
	Calculates Levenshtein Quality.
	Formula: -10 * log10( min( max(1/L^2, ed/L), 1 ) )
	"""
	if ref_len == 0:
		return 0.0

	L = float(ref_len)
	term1 = 1.0 / (L * L)
	term2 = float(ed) / L

	# Bounded error probability
	prob = min(max(term1, term2), 1.0)

	return -10 * math.log10(prob)


def insert_target(cursor: sqlite3.Cursor, tgt_id: str, tgt_seq: str, tgt_len: int) -> None:
	"""Inserts target sequence metadata. strictly DB IO."""
	cursor.execute('''
		INSERT OR REPLACE INTO Target (tgt_id, tgt_refseq, tgt_reflen)
		VALUES (?, ?, ?)
	''', (tgt_id, tgt_seq, tgt_len))


def insert_read(cursor: sqlite3.Cursor, data: Tuple) -> None:
	"""Inserts a single processed read record."""
	cursor.execute('''
		INSERT OR IGNORE INTO Reads (
			uniq_id, exp_id, tgt_id, read_id, readseq, readlen,
			model_tier, model_ver, trim, mod_bitflag,
			ed, q_bc, q_ld, ER
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
	''', data)


############
# argparse #
############

parser = argparse.ArgumentParser(description="Ingest BAM, calculate metrics, and populate DB")
parser.add_argument("-e", "--expid", required=True,
	help="Experiment ID")
parser.add_argument("-b", "--bam", required=True,
	help="Input uBAM file")
parser.add_argument("-s", "--summary", required=True,
	help="Pod5 summary TSV")
parser.add_argument("-r", "--ref", required=True,
	help="Target FASTA")
parser.add_argument("-d", "--database", required=True,
	help="Target SQLite DB")
parser.add_argument("-o", "--output_bam", required=True,
	help="Output tagged BAM")
args = parser.parse_args()


#########
# param #
#########

EXP_ID = args.expid
INPUT_BAM = Path(args.bam)
SUMMARY_TSV = Path(args.summary)
REF_FASTA = Path(args.ref)
DB_PATH = Path(args.database)
OUTPUT_BAM = Path(args.output_bam)

MODEL_TIER, MODEL_VER, TRIM, MOD_BITFLAG = parse_bam_filename(INPUT_BAM)

print(f"[ingest] Run Metadata:\n  Exp: {EXP_ID}\n  Tier: {MODEL_TIER}\n  Ver: {MODEL_VER}\n  Trim: {TRIM}\n  Mods: {MOD_BITFLAG}")


########
# main #
########

tgt_id = None
tgt_seq_parts = []

try:
	with open(REF_FASTA, 'r') as f:
		for line in f:
			line = line.strip()
			if not line: continue
			if line.startswith('>'):
				tgt_id = line[1:]
			else:
				tgt_seq_parts.append(line)

	tgt_seq = "".join(tgt_seq_parts)
	if not tgt_id or not tgt_seq:
		raise ValueError("Empty ID or Sequence")

	tgt_len = len(tgt_seq)
	print(f"[ingest] Target Loaded: {tgt_id} ({tgt_len} bp)")

except Exception as e:
	sys.exit(f"[ingest] Error loading FASTA: {e}")

print(f"[ingest] Loading Pod5 Summary from {SUMMARY_TSV}...")

try:
	summary_df = pd.read_csv(SUMMARY_TSV, sep='\t', usecols=['read_id', 'end_reason'])
	end_reason_map = pd.Series(summary_df.end_reason.values, index=summary_df.read_id).to_dict()
	print(f"[ingest] Loaded {len(end_reason_map)} records.")
except Exception as e:
	sys.exit(f"[ingest] Error loading summary TSV: {e}")

with sqlite3.connect(DB_PATH) as conn:
	c = conn.cursor()

	insert_target(c, tgt_id, tgt_seq, tgt_len)
	conn.commit()
	save_verbosity = pysam.set_verbosity(0)

	with pysam.AlignmentFile(INPUT_BAM, "rb", check_sq=False) as bam_in, \
			pysam.AlignmentFile(OUTPUT_BAM, "wb", template=bam_in) as bam_out:

		pysam.set_verbosity(save_verbosity)
		print("[ingest] Processing reads...")

		count_total = 0
		count_processed = 0
		count_skipped = 0
		count_unknown_er = 0

		for read in bam_in:
			count_total += 1
			seq = read.query_sequence
			quals = read.query_qualities

			if not seq or not quals:
				count_skipped += 1
				continue

			read_id = read.query_name
			read_len = len(seq)

			# Metrics
			q_bc = calculate_q_bc(quals)

			# Alignment
			align_res = edlib.align(seq, tgt_seq, mode="NW", task="distance")
			ed = align_res["editDistance"]
			q_ld = calculate_q_ld(ed, tgt_len)

			# Metadata
			er_val = end_reason_map.get(read_id, "unknown")
			if er_val == "unknown":
				count_unknown_er += 1
			uniq_id = f"{EXP_ID}_{MODEL_TIER}{MODEL_VER}t{TRIM}m{MOD_BITFLAG}_{read_id}"

			# DB Insert
			read_data = (
				uniq_id, EXP_ID, tgt_id, read_id, seq, read_len,
				MODEL_TIER, MODEL_VER, TRIM, MOD_BITFLAG,
				ed, q_bc, q_ld, er_val
			)
			insert_read(c, read_data)

			# BAM Write
			read.set_tag("ER", er_val, value_type="Z")
			bam_out.write(read)

			count_processed += 1
			if count_processed % 1000 == 0:
				print(f"  ...processed {count_processed} reads", end='\r')

	conn.commit()

print(f"[ingest] Complete.")
print(f"  - Total Reads:  {count_total}")
print(f"  - Processed:    {count_processed}")
print(f"  - Skipped:      {count_skipped}")
print(f"  - Unknown ER:   {count_unknown_er}")
print(f"  - Output DB:    {DB_PATH}")
print(f"  - Output BAM:   {OUTPUT_BAM}")
