# ingest.py
# Main logic for ingesting BAM reads into the SQLite database.

import edlib
import math
import pandas as pd
import pysam
import sqlite3
from pathlib import Path

# Maps full BAM filename tier token → single-char abbreviation used in DB and uniq_id
TIER_MAP = {"fast": "f", "hac": "h", "sup": "s"}

# Tables that must exist before ingestion can proceed
EXPECTED_TABLES = {"Reads", "Mods", "Exp", "Target"}

# Number of reads to process between intermediate commits
BATCH_SIZE = 10_000


def parse_bam_filename(bam_path: Path):
	"""
	Parses metadata from BAM filename (Right-to-Left parsing).
	Expected Format: {exp_id}_{tier}_v{ver}_{trim}_{mods}.bam
	"""
	filename = bam_path.resolve().stem
	parts = filename.split('_')

	if len(parts) < 5:
		raise ValueError(
			f"Metadata parse failed for '{filename}': Filename format incorrect.")

	try:
		mods = int(parts[-1])
		trim = int(parts[-2][-1])  # 'trim0' -> 0
		ver_str = parts[-3].lstrip("v")
		tier = parts[-4]
		exp_id = "_".join(parts[:-4])

		return exp_id, tier, ver_str, trim, mods
	except Exception as e:
		raise ValueError(f"Metadata parse failed for '{filename}': {e}")


def calculate_q_bc(quality_scores) -> float:
	"""Calculates probability-averaged Phred quality score."""
	if not quality_scores:
		return 0.0

	avg_prob = sum(10 ** (-q / 10.0) for q in quality_scores) / len(quality_scores)
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

	prob = min(max(term1, term2), 1.0)

	return -10 * math.log10(prob)


def load_end_reasons(tsv_path):
	"""Loads Read ID -> End Reason mapping from TSV."""
	print(f"[ingest] Loading Pod5 Summary from {tsv_path}...")
	df = pd.read_csv(tsv_path, sep='\t', usecols=['read_id', 'end_reason'],
					 dtype=str, engine='c')
	return pd.Series(df.end_reason.values, index=df.read_id).to_dict()


def run(bam_path, db_path, meta_path, batch_size=BATCH_SIZE):
	"""
	Main ingestion logic.

	Raises:
		ValueError: If the BAM filename cannot be parsed or the tier token is unknown.
		RuntimeError: If the DB is missing expected tables, Target is empty, or
		              mod_bitflag is not present in the Mods table.
		Any sqlite3 or pysam exception propagates after rollback.
	"""
	bam_file = Path(bam_path)
	db_file = Path(db_path)

	# Parse BAM filename metadata
	exp_id, tier_raw, ver, trim, mods = parse_bam_filename(bam_file)

	if tier_raw not in TIER_MAP:
		raise ValueError(
			f"[ingest] Unknown model tier '{tier_raw}'. Expected one of {list(TIER_MAP)}.")
	tier = TIER_MAP[tier_raw]

	print(f"[ingest] Run Metadata:\n  Exp: {exp_id}\n  Tier: {tier}\n  Ver: {ver}\n  Trim: {trim}\n  Mods: {mods}")

	conn = sqlite3.connect(db_file)
	# FK enforcement left OFF (SQLite default); referential integrity is
	# guaranteed by mkdb/init pre-flight validation and mod_bitflag check above.
	c = conn.cursor()

	try:
		# Pre-flight: verify DB was fully initialized by mkdb
		c.execute("SELECT name FROM sqlite_master WHERE type='table'")
		existing_tables = {row[0] for row in c.fetchall()}
		missing = EXPECTED_TABLES - existing_tables
		if missing:
			raise RuntimeError(
				f"[ingest] DB is missing expected tables: {missing}. Run 'mkdb' first.")

		# Load Target reference
		c.execute("SELECT tgt_id, tgt_refseq, tgt_reflen FROM Target LIMIT 1")
		row = c.fetchone()
		if not row:
			raise RuntimeError("[ingest] No Target reference in DB. Run 'init' first.")
		tgt_id, tgt_seq, tgt_len = row
		print(f"[ingest] Target Loaded: {tgt_id} ({tgt_len} bp)")

		# Validate mod_bitflag against Mods table
		c.execute("SELECT 1 FROM Mods WHERE mod_bitflag = ?", (mods,))
		if not c.fetchone():
			raise RuntimeError(
				f"[ingest] mod_bitflag={mods} not found in Mods table. "
				f"Check filename or update Mods table.")

		# Load end-reason metadata
		er_map = load_end_reasons(meta_path)

		# Stream BAM
		samfile = pysam.AlignmentFile(bam_file, "rb", check_sq=False)

		print("[ingest] Processing reads...")
		count_total = 0
		count_processed = 0
		count_skipped = 0
		count_non_primary = 0
		count_unknown_er = 0

		sql = '''INSERT INTO Reads (
					uniq_id, exp_id, tgt_id, read_id, readseq, readlen,
					model_tier, model_ver, trim, mod_bitflag,
					ed, q_bc, q_ld, ER
				 ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)'''

		for read in samfile:
			count_total += 1

			# Supplementary/secondary records are alignment artifacts of the same
			# physical molecule — skip to avoid duplicate uniq_id on the same read_id.
			if read.is_supplementary or read.is_secondary:
				count_non_primary += 1
				continue

			seq = read.query_sequence
			quals = read.query_qualities

			if not seq or not quals:
				count_skipped += 1
				continue

			# read_id is the raw sequencer UUID; retained in full for traceability
			read_id = read.query_name
			read_len = len(seq)

			q_bc = calculate_q_bc(quals)

			align_res = edlib.align(seq, tgt_seq, mode="NW", task="distance")
			ed = align_res["editDistance"]
			q_ld = calculate_q_ld(ed, tgt_len)

			er_val = er_map.get(read_id, "unknown")
			if er_val == "unknown":
				count_unknown_er += 1

			uniq_id = f"{exp_id}_{tier}{ver}t{trim}m{mods}_{read_id}"

			c.execute(sql, (
				uniq_id, exp_id, tgt_id, read_id, seq, read_len,
				tier, ver, trim, mods,
				ed, q_bc, q_ld, er_val
			))

			count_processed += 1
			if count_processed % batch_size == 0:
				conn.commit()
				print(f"  ...committed {count_processed} reads", end='\r')

		# Final commit for the last partial batch
		conn.commit()
		print(f"\n[ingest] Complete.")
		print(f"  - Total Records: {count_total}")
		print(f"  - Non-primary:   {count_non_primary}")
		print(f"  - Processed:     {count_processed}")
		print(f"  - Skipped:       {count_skipped}")
		print(f"  - Unknown ER:    {count_unknown_er}")

	except Exception:
		conn.rollback()
		raise
	finally:
		if 'samfile' in locals():
			samfile.close()
		conn.close()
