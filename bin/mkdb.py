#!/usr/bin/env python3
# mkdb.py
# initialize SQLite DB and populates the static Mods table

import argparse
import os
import sqlite3
import time


############
# argparse #
############

parser = argparse.ArgumentParser(description="Initialize SMA Database")
parser.add_argument("-e", "--expid", required=True,
	help="Experiment ID")
parser.add_argument("-o", "--outdir", default=".",
	help="Output directory [%(default)s]")
args = parser.parse_args()


##########
# SQLite #
##########

# Construct the initial filename
base_db_name = f"SMA_{args.expid}.db"
db_path = os.path.join(args.outdir, base_db_name)

# Check if DB exists; if so, timestamp the new one to avoid overwriting
if os.path.exists(db_path):
	timestamp = time.strftime("%Y%m%d-%H%M%S")
	base_db_name = f"SMA_{args.expid}_{timestamp}.db"
	db_path = os.path.join(args.outdir, base_db_name)
	print(f"[mkdb] Warning: Database {db_path} already exists. Creating new database as {base_db_name} instead.")

# Ensure output directory exists
os.makedirs(args.outdir, exist_ok=True)

conn = sqlite3.connect(db_path)
c = conn.cursor()

# Reads Table
c.execute('''
	CREATE TABLE Reads (
		uniq_id TEXT PRIMARY KEY,
		exp_id TEXT,
		refseq_id TEXT,
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
		FOREIGN KEY(refseq_id) REFERENCES Refseq(refseq_id),
		FOREIGN KEY(mod_bitflag) REFERENCES Mods(mod_bitflag)
	)
''')

# Mods Table
c.execute('''
	CREATE TABLE Mods (
		mod_bitflag INTEGER PRIMARY KEY,
		mods TEXT
	)
''')

# Exp Table
c.execute('''
	CREATE TABLE Exp (
		exp_id TEXT PRIMARY KEY,
		exp_desc TEXT
	)
''')

# Refseq Table
c.execute('''
	CREATE TABLE Refseq (
		refseq_id TEXT PRIMARY KEY,
		refseq TEXT,
		reflen INTEGER
	)
''')


# Mods Table Population
mods_data = [
	(0, "non"),
	(1, "6mA"),
	(2, "5mCG_5hmCG"),
	(4, "5mC_5hmC"),
	(8, "4mC_5mC"),
	(16, "5mC"),
	(3, "6mA,5mCG_5hmCG"),
	(5, "6mA,5mC_5hmC"),
	(9, "6mA,4mC_5mC"),
	(17, "6mA,5mC")
]

c.executemany(
	"INSERT OR IGNORE INTO Mods (mod_bitflag, mods) VALUES (?, ?)", mods_data)

# --- Commit & Close ---
conn.commit()
conn.close()

print(f"[mkdb] Database {db_path} created successfully.")
