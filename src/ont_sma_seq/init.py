# init.py
# Lock a reference FASTA into the Target table.

import sqlite3
from pathlib import Path


def parse_reference(fasta_path):
	"""
	Validates that the FASTA contains exactly one sequence.
	Returns tuple: (header, sequence)
	"""
	path = Path(fasta_path)
	if not path.exists():
		raise FileNotFoundError(f"Reference file not found: {fasta_path}")

	header = None
	sequence = []

	with open(path, 'r') as f:
		for line in f:
			line = line.strip()
			if not line:
				continue

			if line.startswith(">"):
				if header is not None:
					raise ValueError("Input FASTA must contain exactly ONE sequence.")
				header = line[1:]
			else:
				if header is None:
					raise ValueError("FASTA file seems malformed (sequence found before header).")
				sequence.append(line)

	if header is None:
		raise ValueError("No sequence found in FASTA file.")

	return header, "".join(sequence)


def run(db_path, ref_fasta):
	"""
	Validates the reference FASTA and locks it into the Target table.
	Uses INSERT OR IGNORE so re-running against an existing DB is a safe no-op.

	Raises:
		FileNotFoundError: If the FASTA path does not exist.
		ValueError: If the FASTA is malformed or contains more than one sequence.
		Any sqlite3 exception propagates after rollback.
	"""
	print(f"[init] Validating reference: {ref_fasta}")
	header, sequence = parse_reference(ref_fasta)

	ref_len = len(sequence)
	print(f"[init] Reference '{header}' (Length: {ref_len})")

	conn = sqlite3.connect(db_path)
	try:
		c = conn.cursor()
		c.execute(
			"INSERT OR IGNORE INTO Target (tgt_id, tgt_refseq, tgt_reflen) VALUES (?, ?, ?)",
			(header, sequence, ref_len)
		)
		conn.commit()
	except Exception:
		conn.rollback()
		conn.close()
		raise
	conn.close()
	print(f"[init] Successfully locked reference into {db_path}")
