# merge.py
# Merges multiple SMA SQLite databases into a single master database.

import sqlite3
import os
import re
import time

# Canonical table order: dependencies before dependents
TABLES = ["Mods", "Exp", "Target", "Reads"]


def run(output_db, input_dbs):
	"""
	Merges one or more source databases into a master database.

	If the output DB already exists, data is merged into it (pre-existing rows
	are preserved via INSERT OR IGNORE). If it does not exist, it is created.

	Args:
		output_db (str): Path to the master output database.
		input_dbs (list[str]): Paths to source .db files to merge.

	Raises:
		ValueError: If input_dbs is empty.
		RuntimeError: If schema extraction fails or a source DB cannot be read.
		Any sqlite3.OperationalError that is not a transient busy/lock error
		propagates after all retry attempts are exhausted.
	"""
	if not input_dbs:
		raise ValueError("[merge] No input databases provided.")

	output_path = os.path.abspath(output_db)

	if os.path.exists(output_path):
		print(f"[merge] Warning: Output DB '{output_db}' exists. Merging into it.")

	print(f"[merge] Output Master DB: {output_db}")

	# isolation_level=None enables explicit transaction control via BEGIN/COMMIT/ROLLBACK
	conn = sqlite3.connect(output_path, timeout=60.0, isolation_level=None)
	cur = conn.cursor()

	cur.execute("PRAGMA busy_timeout = 60000")  # milliseconds
	cur.execute("PRAGMA main.synchronous = OFF")
	cur.execute("PRAGMA main.foreign_keys = OFF")

	# --- PHASE 1: EXTRACT AND APPLY SCHEMA ---
	first_db = os.path.abspath(input_dbs[0])
	print(f"[merge] Extracting schema from {input_dbs[0]}...")

	cur.execute(f"ATTACH DATABASE 'file:{first_db}?mode=ro' AS schema_src")
	try:
		cur.execute(
			"SELECT sql FROM schema_src.sqlite_master "
			"WHERE type='table' AND name NOT LIKE 'sqlite_%' AND sql IS NOT NULL"
		)
		schema_statements = [row[0] for row in cur.fetchall()]
	finally:
		cur.execute("DETACH DATABASE schema_src")

	if not schema_statements:
		raise RuntimeError(
			f"[merge] No table definitions found in schema source: {input_dbs[0]}")

	# Inject IF NOT EXISTS so this is safe when merging into an existing master DB
	safe_statements = [
		re.sub(r"(?i)\bCREATE TABLE\b", "CREATE TABLE IF NOT EXISTS", stmt)
		for stmt in schema_statements
	]

	cur.execute("BEGIN")
	try:
		for stmt in safe_statements:
			cur.execute(stmt)
		cur.execute("COMMIT")
	except Exception:
		cur.execute("ROLLBACK")
		raise

	print("[merge] Schema initialized.")
	print(f"[merge] Merging {len(input_dbs)} source file(s)...")

	# --- PHASE 2: MERGE DATA ---
	for i, db_file in enumerate(input_dbs):
		abs_path = os.path.abspath(db_file)

		if abs_path == output_path:
			print(f"  -> Skipping {db_file} (same as output)")
			continue

		print(f"  -> Processing {i + 1}/{len(input_dbs)}: {db_file}")

		alias = f"src_{i}"  # unique alias per source to avoid collisions
		attached = False

		for attempt in range(1, 6):
			try:
				cur.execute(f"ATTACH DATABASE 'file:{abs_path}?mode=ro' AS {alias}")
				attached = True

				# Per-source transaction so each DB is committed independently
				cur.execute("BEGIN")
				try:
					for table in TABLES:
						cur.execute(
							f"INSERT OR IGNORE INTO main.{table} SELECT * FROM {alias}.{table}"
						)
					cur.execute("COMMIT")
				except Exception:
					cur.execute("ROLLBACK")
					raise

				break  # success — move to next source DB

			except sqlite3.OperationalError as e:
				msg = str(e).lower()
				if "locked" in msg or "busy" in msg:
					wait_s = min(2 * attempt, 10)
					print(f"[merge] Busy/locked. Retry {attempt}/5 in {wait_s}s: {e}")
					time.sleep(wait_s)
					continue
				raise

			finally:
				if attached:
					try:
						cur.execute(f"DETACH DATABASE {alias}")
					except Exception:
						pass
					attached = False

	conn.close()
	print(f"[merge] Success! All databases merged into {output_db}")
