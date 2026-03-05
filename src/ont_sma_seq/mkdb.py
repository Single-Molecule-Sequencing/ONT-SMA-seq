# mkdb.py
# Initializes the SQLite database and schema for a new experiment.

import os
import sqlite3
import time


def run(exp_id, out_dir):
	"""
	Initializes the SQLite DB, parses Experiment ID, and populates static tables.
	"""

	# Parse experiment ID components
	exp_parts = exp_id.split('_')

	if len(exp_parts) >= 3:
		flow_cell_id = exp_parts[0]
		sample_id = exp_parts[1]
		alias = "_".join(exp_parts[2:])
	else:
		raise ValueError(
			f"[mkdb] Experiment ID '{exp_id}' does not follow 'FlowCell_Sample_Alias' convention.")

	print(f"[mkdb] Parsed Experiment Metadata:")
	print(f"  - Exp ID:    {exp_id}")
	print(f"  - Flow Cell: {flow_cell_id}")
	print(f"  - Sample:    {sample_id}")
	print(f"  - Alias:     {alias}")

	base_db_name = f"SMA_{exp_id}.db"
	db_path = os.path.join(out_dir, base_db_name)

	os.makedirs(out_dir, exist_ok=True)

	# Avoid overwriting an existing DB — append timestamp instead
	if os.path.exists(db_path):
		timestamp = time.strftime("%Y%m%d-%H%M%S")
		base_db_name = f"SMA_{exp_id}_{timestamp}.db"
		db_path = os.path.join(out_dir, base_db_name)
		print(f"[mkdb] Warning: Database already exists. Creating new database as {base_db_name} instead.")

	conn = sqlite3.connect(db_path)
	c = conn.cursor()

	# DDL in dependency order: Mods → Exp → Target → Reads

	c.execute('''
		CREATE TABLE Mods (
			mod_bitflag INTEGER PRIMARY KEY,
			mods TEXT
		)
	''')

	c.execute('''
		CREATE TABLE Exp (
			exp_id TEXT PRIMARY KEY,
			flow_cell_id TEXT,
			sample_id TEXT,
			alias TEXT,
			exp_desc TEXT
		)
	''')

	c.execute('''
		CREATE TABLE Target (
			tgt_id TEXT PRIMARY KEY,
			tgt_refseq TEXT,
			tgt_reflen INTEGER
		)
	''')

	c.execute('''
		CREATE TABLE Reads (
			uniq_id TEXT PRIMARY KEY,
			exp_id TEXT,
			tgt_id TEXT,
			read_id TEXT,
			readseq TEXT,
			readlen INTEGER,
			model_tier TEXT,
			model_ver TEXT,
			trim INTEGER,
			mod_bitflag INTEGER,
			ed INTEGER,
			q_bc REAL,
			q_ld REAL,
			ER TEXT,
			FOREIGN KEY(exp_id) REFERENCES Exp(exp_id),
			FOREIGN KEY(tgt_id) REFERENCES Target(tgt_id),
			FOREIGN KEY(mod_bitflag) REFERENCES Mods(mod_bitflag)
		)
	''')

	# --- Data Population ---

	# Static modification bitflag definitions
	# 4-bit encoding: bit 3 = 6mA (independent), bits 2-0 = C-mod enum (0=none, 1-4 mutually exclusive)
	mods_data = [
		(0,  "non"),
		(8,  "6mA"),
		(1,  "5mCG_5hmCG"),
		(2,  "5mC_5hmC"),
		(3,  "4mC_5mC"),
		(4,  "5mC"),
		(9,  "6mA,5mCG_5hmCG"),
		(10, "6mA,5mC_5hmC"),
		(11, "6mA,4mC_5mC"),
		(12, "6mA,5mC")
	]

	c.executemany(
		"INSERT OR IGNORE INTO Mods (mod_bitflag, mods) VALUES (?, ?)", mods_data
	)

	c.execute(
		"INSERT INTO Exp (exp_id, flow_cell_id, sample_id, alias, exp_desc) VALUES (?, ?, ?, ?, ?)",
		(exp_id, flow_cell_id, sample_id, alias, "Initialized via mkdb")
	)

	conn.commit()
	conn.close()

	print(f"[mkdb] Database {db_path} created successfully.")
	return db_path
