#!/usr/bin/env python3
"""Initialize the central unified nanopore database."""

import argparse
import sys
from pathlib import Path

# Add lib to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from lib.db_schema import create_central_db
from lib.db_lookups import populate_lookups


def main():
    parser = argparse.ArgumentParser(
        description="Initialize the central unified nanopore database"
    )
    parser.add_argument(
        "-o", "--output",
        default="nanopore_unified.db",
        help="Output database path [%(default)s]"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing database"
    )
    args = parser.parse_args()

    db_path = Path(args.output)

    if db_path.exists() and not args.force:
        print(f"[init_central_db] Error: {db_path} already exists. Use --force to overwrite.")
        sys.exit(1)

    if db_path.exists() and args.force:
        db_path.unlink()
        print(f"[init_central_db] Removed existing database: {db_path}")

    print(f"[init_central_db] Creating database: {db_path}")
    conn = create_central_db(str(db_path))

    print("[init_central_db] Populating lookup tables...")
    populate_lookups(conn)

    # Report table counts
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = cursor.fetchall()
    print(f"[init_central_db] Created {len(tables)} tables")

    # Count lookup entries
    lk_tables = [t[0] for t in tables if t[0].startswith('lk_')]
    total_lookups = 0
    for table in lk_tables:
        cursor.execute(f"SELECT COUNT(*) FROM {table}")
        count = cursor.fetchone()[0]
        total_lookups += count
    print(f"[init_central_db] Populated {total_lookups} lookup entries across {len(lk_tables)} lookup tables")

    conn.close()
    print(f"[init_central_db] Done. Database ready at: {db_path}")


if __name__ == "__main__":
    main()
