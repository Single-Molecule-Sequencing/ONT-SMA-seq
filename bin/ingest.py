#!/usr/bin/env python3
# ingest.py
# Main processing script: BAM tagging, Metric Calculation, DB Ingestion
# Supports two modes:
#   1. Classification mode  (-ss / -rd) : barcode demultiplexing + per-target alignment
#   2. Single-target mode   (-r)        : original behaviour, backward compatible

import argparse
import sys
import math
import sqlite3
import pysam
import edlib
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple, Any, Dict, Optional

from barcodes import BARCODES, BARCODE_LENGTH, classify_barcode, reverse_complement
from construct import parse_construct_toml
from report_analysis import find_flank_position
from sample_sheet import parse_sample_sheet, detect_barcode_ambiguity


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SEGMENT_LEN = 100  # bp window for barcode search at read start/end


#########
# funcs #
#########


def classify_truncation(
	read_seq: str,
	bc_start_conf: float,
	bc_end_conf: float,
	bc_start_ed: int,
	bc_end_ed: int,
	flank_front: str,
	flank_rear: str,
	flank_rev_front: str,
	segment_len: int,
	confidence,
	min_target_length: int,
) -> str:
	"""Classify the truncation level of a read using construct geometry.

	Classification levels (evaluated in order):
	  - adapter_only    : start barcode confidence below threshold
	  - bc1_only        : read too short to contain target after front flank
	  - full_length     : end barcode confidence high AND rear flank found
	  - bc1_target_bc2  : end barcode has moderate confidence (0.3 <= conf < threshold)
	  - bc1_target      : has bc1 and target but no detectable bc2

	Parameters
	----------
	read_seq : str
		Full read sequence.
	bc_start_conf : float
		Start barcode confidence score.
	bc_end_conf : float
		End barcode confidence score.
	bc_start_ed : int
		Start barcode edit distance.
	bc_end_ed : int
		End barcode edit distance.
	flank_front : str
		Front flanking sequence (mask1_rear) to search after bc1.
	flank_rear : str
		Rear flanking sequence (mask2_front).
	flank_rev_front : str
		Reverse complement of mask2_front, to search near end of read.
	segment_len : int
		Search window length for barcodes at read start/end.
	confidence : ConfidenceConfig
		Confidence thresholds (full_length_threshold, start_barcode_min).
	min_target_length : int
		Minimum target length to consider a read as containing target.

	Returns
	-------
	str
		One of: 'adapter_only', 'bc1_only', 'full_length', 'bc1_target_bc2',
		'bc1_target'.
	"""
	# Step 1: Check start barcode confidence
	if bc_start_conf < confidence.start_barcode_min:
		return "adapter_only"

	# Step 2: Find front flank (mask1_rear) in first segment_len + flank_len region
	flank_len = len(flank_front)
	read_len = len(read_seq)
	search_end = min(segment_len + flank_len, read_len)
	front_hit = find_flank_position(read_seq, flank_front, 0, search_end)

	# Step 3: If read too short for target after flank -> bc1_only
	if front_hit is not None:
		target_start = front_hit["end"]
	else:
		# No flank found; estimate target starts after segment_len
		target_start = segment_len

	remaining = read_len - target_start
	if remaining < min_target_length:
		return "bc1_only"

	# Step 4: Find rear flank (RC of mask2_front) in last segment_len region
	rear_search_start = max(0, read_len - segment_len)
	rear_hit = find_flank_position(
		read_seq, flank_rev_front, rear_search_start, read_len
	)

	# Step 5: full_length check
	if bc_end_conf >= confidence.full_length_threshold and rear_hit is not None:
		return "full_length"

	# Step 6: bc1_target_bc2 (moderate end barcode confidence)
	if confidence.start_barcode_min <= bc_end_conf < confidence.full_length_threshold:
		return "bc1_target_bc2"

	# Step 7: default
	return "bc1_target"

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
	"""Inserts a single processed read record (21 columns)."""
	cursor.execute('''
		INSERT OR IGNORE INTO Reads (
			uniq_id, exp_id, tgt_id, read_id, readseq, readlen,
			model_tier, model_ver, trim, mod_bitflag,
			ed, q_bc, q_ld, ER,
			bc_start_id, bc_start_ed, bc_start_conf,
			bc_end_id, bc_end_ed, bc_end_conf,
			trunc_level
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
	''', data)


def load_references(ref_dir: Path) -> Dict[str, Tuple[str, int]]:
	"""Load all FASTA files from a directory.

	Globs for *.fasta and *.fa files. Each file should contain a single
	sequence with a header line (>name).

	Parameters
	----------
	ref_dir : Path
		Directory containing reference FASTA files.

	Returns
	-------
	dict
		Mapping of alias (header name) to (sequence, length).
	"""
	refs: Dict[str, Tuple[str, int]] = {}

	for pattern in ("*.fasta", "*.fa"):
		for fasta_path in ref_dir.glob(pattern):
			tgt_id = None
			seq_parts = []

			with open(fasta_path, 'r') as f:
				for line in f:
					line = line.strip()
					if not line:
						continue
					if line.startswith('>'):
						tgt_id = line[1:]
					else:
						seq_parts.append(line)

			if tgt_id and seq_parts:
				seq = "".join(seq_parts)
				refs[tgt_id] = (seq, len(seq))

	return refs


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
parser.add_argument("-r", "--ref", required=False, default=None,
	help="Target FASTA (single-target mode)")
parser.add_argument("-d", "--database", required=True,
	help="Target SQLite DB")
parser.add_argument("-o", "--output_bam", required=True,
	help="Output tagged BAM")
parser.add_argument("-ss", "--sample-sheet", required=False, default=None,
	help="MinKNOW sample sheet CSV (enables classification mode)")
parser.add_argument("-rd", "--ref-dir", required=False, default=None,
	help="Directory of per-target reference FASTAs")
parser.add_argument("--full-construct", action="store_true", default=False,
	help="Force full construct alignment for all reads")
parser.add_argument("--split-bams", required=False, default=None,
	help="Output directory for optional per-barcode BAMs")
parser.add_argument("-c", "--construct",
	help="Construct TOML file (enables truncation detection)")
parser.add_argument("--tag", action="store_true", default=False,
	help="Write BC tags (BS, BE, BA) to output BAM")
args = parser.parse_args()


#########
# param #
#########

EXP_ID = args.expid
INPUT_BAM = Path(args.bam)
SUMMARY_TSV = Path(args.summary)
DB_PATH = Path(args.database)
OUTPUT_BAM = Path(args.output_bam)

# Determine mode
CLASSIFICATION_MODE = args.sample_sheet is not None and args.ref_dir is not None

if not CLASSIFICATION_MODE and args.ref is None:
	sys.exit("[ingest] Error: Either -r (single-target) or -ss/-rd (classification) must be provided.")

MODEL_TIER, MODEL_VER, TRIM, MOD_BITFLAG = parse_bam_filename(INPUT_BAM)

print(f"[ingest] Run Metadata:\n  Exp: {EXP_ID}\n  Tier: {MODEL_TIER}\n  Ver: {MODEL_VER}\n  Trim: {TRIM}\n  Mods: {MOD_BITFLAG}")


########
# main #
########

if CLASSIFICATION_MODE:
	# -----------------------------------------------------------------------
	# Classification mode
	# -----------------------------------------------------------------------
	SS_PATH = Path(args.sample_sheet)
	REF_DIR = Path(args.ref_dir)

	print(f"[ingest] Classification mode enabled")
	print(f"  Sample sheet: {SS_PATH}")
	print(f"  Ref dir:      {REF_DIR}")

	# Parse sample sheet -> barcode pair to alias mapping
	barcode_pair_to_alias = parse_sample_sheet(SS_PATH)
	print(f"[ingest] Sample sheet: {len(barcode_pair_to_alias)} barcode pairs loaded")

	# Check ambiguity
	is_ambiguous = detect_barcode_ambiguity(barcode_pair_to_alias)
	if is_ambiguous or args.full_construct:
		print("[ingest] Barcode ambiguity detected or --full-construct forced: using full construct alignment")

	# Build expected_barcodes: only barcodes appearing in sample sheet
	used_barcodes = set()
	for (bc_up, bc_down) in barcode_pair_to_alias.keys():
		used_barcodes.add(bc_up)
		used_barcodes.add(bc_down)

	expected_barcodes = {bc_id: BARCODES[bc_id] for bc_id in used_barcodes if bc_id in BARCODES}
	expected_barcodes_rc = {bc_id: reverse_complement(seq) for bc_id, seq in expected_barcodes.items()}

	print(f"[ingest] Expected barcodes: {sorted(expected_barcodes.keys())}")

	# Load reference FASTAs
	references = load_references(REF_DIR)
	print(f"[ingest] Loaded {len(references)} references from {REF_DIR}")

	# Validate: every alias in sample sheet must have a matching reference
	aliases = set(barcode_pair_to_alias.values())
	missing_refs = aliases - set(references.keys())
	if missing_refs:
		sys.exit(f"[ingest] Error: Missing reference FASTAs for aliases: {missing_refs}")

	# --- Parse construct TOML if provided ---
	construct_cfg = None
	flank_front = None
	flank_rear = None
	flank_rev_front = None

	if args.construct:
		construct_cfg = parse_construct_toml(Path(args.construct))
		print(f"[ingest] Construct TOML loaded: {construct_cfg.arrangement.name}")

		# Use TOML targets for barcode_pair_to_alias if not already set from -ss
		# (-ss takes precedence when both are provided)
		if not barcode_pair_to_alias:
			barcode_pair_to_alias = construct_cfg.barcode_pair_to_alias()
			print(f"[ingest] Using TOML targets for barcode pairing: {len(barcode_pair_to_alias)} pairs")

		# Use TOML flanks
		flank_front = construct_cfg.flank_front    # mask1_rear
		flank_rear = construct_cfg.flank_rear      # mask2_front
		flank_rev_front = reverse_complement(flank_rear) if flank_rear else ""
		print(f"[ingest] Truncation detection enabled (flank_front={flank_front}, flank_rear={flank_rear})")

	# Load summary
	print(f"[ingest] Loading Pod5 Summary from {SUMMARY_TSV}...")
	try:
		summary_df = pd.read_csv(SUMMARY_TSV, sep='\t', usecols=['read_id', 'end_reason'])
		end_reason_map = pd.Series(summary_df.end_reason.values, index=summary_df.read_id).to_dict()
		print(f"[ingest] Loaded {len(end_reason_map)} records.")
	except Exception as e:
		sys.exit(f"[ingest] Error loading summary TSV: {e}")

	# Insert ALL targets into DB
	with sqlite3.connect(DB_PATH) as conn:
		c = conn.cursor()

		for alias, (ref_seq, ref_len) in references.items():
			insert_target(c, alias, ref_seq, ref_len)
		conn.commit()

		# Optional: split BAMs
		split_writers: Dict[str, pysam.AlignmentFile] = {}
		split_dir: Optional[Path] = None
		if args.split_bams:
			split_dir = Path(args.split_bams)
			split_dir.mkdir(parents=True, exist_ok=True)

		save_verbosity = pysam.set_verbosity(0)

		with pysam.AlignmentFile(INPUT_BAM, "rb", check_sq=False) as bam_in, \
				pysam.AlignmentFile(OUTPUT_BAM, "wb", template=bam_in) as bam_out:

			pysam.set_verbosity(save_verbosity)
			print("[ingest] Processing reads (classification mode)...")

			count_total = 0
			count_processed = 0
			count_skipped = 0
			count_unknown_er = 0
			count_matched = 0
			count_unmatched = 0

			# Per-alias counts for classification summary
			alias_counts: Dict[str, int] = {}

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

				# --- Barcode classification ---
				# Start: first SEGMENT_LEN bp
				start_segment = seq[:SEGMENT_LEN]
				start_result = classify_barcode(start_segment, expected_barcodes)
				bc_start_id = start_result["barcode_id"]
				bc_start_ed = start_result["edit_distance"]
				bc_start_conf = start_result["confidence"]

				# End: last SEGMENT_LEN bp, classify against RC barcodes
				end_segment = seq[-SEGMENT_LEN:]
				end_result = classify_barcode(end_segment, expected_barcodes_rc)
				bc_end_id = end_result["barcode_id"]
				bc_end_ed = end_result["edit_distance"]
				bc_end_conf = end_result["confidence"]

				# --- Target assignment from barcode pair ---
				pair_key = (bc_start_id, bc_end_id)
				alias = barcode_pair_to_alias.get(pair_key)

				if alias is not None:
					# Matched pair -> align against assigned target
					ref_seq, ref_len = references[alias]
					align_res = edlib.align(seq, ref_seq, mode="NW", task="distance")
					ed = align_res["editDistance"]
					q_ld = calculate_q_ld(ed, ref_len)
					tgt_id = alias
					count_matched += 1
					alias_counts[alias] = alias_counts.get(alias, 0) + 1
				else:
					# Unmatched pair
					tgt_id = f"unmatched_{bc_start_id}_{bc_end_id}"
					ed = None
					q_ld = None
					count_unmatched += 1

				# Metadata
				er_val = end_reason_map.get(read_id, "unknown")
				if er_val == "unknown":
					count_unknown_er += 1
				uniq_id = f"{EXP_ID}_{MODEL_TIER}{MODEL_VER}t{TRIM}m{MOD_BITFLAG}_{read_id}"

				# Truncation classification
				trunc_level = None
				if construct_cfg:
					trunc_level = classify_truncation(
						seq, bc_start_conf, bc_end_conf,
						bc_start_ed, bc_end_ed,
						flank_front, flank_rear, flank_rev_front,
						SEGMENT_LEN, construct_cfg.confidence,
						construct_cfg.truncation.min_target_length,
					)

				# DB Insert (21 columns)
				read_data = (
					uniq_id, EXP_ID, tgt_id, read_id, seq, read_len,
					MODEL_TIER, MODEL_VER, TRIM, MOD_BITFLAG,
					ed, q_bc, q_ld, er_val,
					bc_start_id, bc_start_ed, bc_start_conf,
					bc_end_id, bc_end_ed, bc_end_conf,
					trunc_level,
				)
				insert_read(c, read_data)

				# BAM tags
				read.set_tag("ER", er_val, value_type="Z")
				if args.tag:
					read.set_tag("BS", bc_start_id, value_type="Z")
					read.set_tag("BE", bc_end_id, value_type="Z")
					if alias is not None:
						read.set_tag("BA", alias, value_type="Z")
					else:
						read.set_tag("BA", tgt_id, value_type="Z")

				bam_out.write(read)

				# Optional split BAMs
				if split_dir is not None and alias is not None:
					if alias not in split_writers:
						split_bam_path = split_dir / f"{alias}.bam"
						split_writers[alias] = pysam.AlignmentFile(
							str(split_bam_path), "wb", template=bam_in
						)
					split_writers[alias].write(read)

				count_processed += 1
				if count_processed % 1000 == 0:
					print(f"  ...processed {count_processed} reads", end='\r')

		# Close split BAM writers
		for writer in split_writers.values():
			writer.close()

		conn.commit()

	# Classification summary
	print(f"\n[ingest] Complete (classification mode).")
	print(f"  - Total Reads:   {count_total}")
	print(f"  - Processed:     {count_processed}")
	print(f"  - Skipped:       {count_skipped}")
	print(f"  - Unknown ER:    {count_unknown_er}")
	print(f"  - Matched:       {count_matched}")
	print(f"  - Unmatched:     {count_unmatched}")
	print(f"  - Per-target breakdown:")
	for alias_name in sorted(alias_counts.keys()):
		print(f"      {alias_name}: {alias_counts[alias_name]}")
	print(f"  - Output DB:     {DB_PATH}")
	print(f"  - Output BAM:    {OUTPUT_BAM}")

else:
	# -----------------------------------------------------------------------
	# Single-target mode (backward compatible)
	# -----------------------------------------------------------------------
	REF_FASTA = Path(args.ref)

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

				# DB Insert (21 columns, barcode fields are NULL)
				read_data = (
					uniq_id, EXP_ID, tgt_id, read_id, seq, read_len,
					MODEL_TIER, MODEL_VER, TRIM, MOD_BITFLAG,
					ed, q_bc, q_ld, er_val,
					None, None, None,
					None, None, None,
					None
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
