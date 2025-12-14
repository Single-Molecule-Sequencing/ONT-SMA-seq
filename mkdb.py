#!/usr/bin/env python3
"""
mkdb.py - Database Initialization Script

Creates the SQLite database schema and populates static lookup tables
for the ONT-SMA-seq pipeline.

Usage:
    python mkdb.py <exp_id>

Output:
    SMA_{exp_id}.db
"""

import argparse
import sqlite3
import os
import sys


# Modification bitflag definitions
MODIFICATIONS = [
    (0, "non"),           # No modifications
    (1, "6mA"),           # Independent
    (2, "5mCG_5hmCG"),    # Mutually exclusive C-Mod
    (4, "5mC_5hmC"),      # Mutually exclusive C-Mod
    (8, "4mC_5mC"),       # Mutually exclusive C-Mod
    (16, "5mC"),          # Mutually exclusive C-Mod
]


def create_database(exp_id: str) -> str:
    """
    Create the SQLite database with the required schema.

    Args:
        exp_id: Experiment identifier

    Returns:
        Path to the created database file
    """
    db_path = f"SMA_{exp_id}.db"

    # Check if database already exists
    if os.path.exists(db_path):
        print(f"Warning: Database {db_path} already exists. It will be overwritten.")
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Create Mods table (lookup table for modification bitflags)
    cursor.execute("""
        CREATE TABLE Mods (
            mod_bitflag INTEGER PRIMARY KEY,
            mods TEXT NOT NULL
        )
    """)

    # Create Exp table (experiment metadata)
    cursor.execute("""
        CREATE TABLE Exp (
            exp_id TEXT PRIMARY KEY,
            exp_desc TEXT
        )
    """)

    # Create Refseq table (reference sequences)
    cursor.execute("""
        CREATE TABLE Refseq (
            refseq_id TEXT PRIMARY KEY,
            refseq TEXT NOT NULL,
            reflen INTEGER NOT NULL
        )
    """)

    # Create Reads table (primary data table)
    cursor.execute("""
        CREATE TABLE Reads (
            uniq_id TEXT PRIMARY KEY,
            exp_id TEXT,
            refseq_id TEXT,
            read_id TEXT NOT NULL,
            readseq TEXT NOT NULL,
            readlen INTEGER NOT NULL,
            model_tier TEXT NOT NULL,
            model_ver TEXT NOT NULL,
            trim INTEGER NOT NULL,
            mod_bitflag INTEGER NOT NULL,
            ed INTEGER,
            q_bc REAL,
            q_ld REAL,
            ER TEXT,
            forced INTEGER,
            channel INTEGER,
            well INTEGER,
            pore_type TEXT,
            num_samples INTEGER,
            start_sample INTEGER,
            median_before REAL,
            scale REAL,
            offset REAL,
            FOREIGN KEY (exp_id) REFERENCES Exp(exp_id),
            FOREIGN KEY (refseq_id) REFERENCES Refseq(refseq_id),
            FOREIGN KEY (mod_bitflag) REFERENCES Mods(mod_bitflag)
        )
    """)

    # Create indices for common queries
    cursor.execute("CREATE INDEX idx_reads_exp_id ON Reads(exp_id)")
    cursor.execute("CREATE INDEX idx_reads_refseq_id ON Reads(refseq_id)")
    cursor.execute("CREATE INDEX idx_reads_read_id ON Reads(read_id)")

    # Populate Mods table with standard bitflag definitions
    cursor.executemany(
        "INSERT INTO Mods (mod_bitflag, mods) VALUES (?, ?)",
        MODIFICATIONS
    )

    # Also insert common combinations
    # 6mA can coexist with any C-mod
    combinations = [
        (1 + 2, "6mA+5mCG_5hmCG"),
        (1 + 4, "6mA+5mC_5hmC"),
        (1 + 8, "6mA+4mC_5mC"),
        (1 + 16, "6mA+5mC"),
    ]
    cursor.executemany(
        "INSERT INTO Mods (mod_bitflag, mods) VALUES (?, ?)",
        combinations
    )

    conn.commit()
    conn.close()

    return db_path


def main():
    parser = argparse.ArgumentParser(
        description="Initialize SQLite database for ONT-SMA-seq pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python mkdb.py EXP001
    python mkdb.py my_experiment_2024
        """
    )
    parser.add_argument(
        "exp_id",
        type=str,
        help="Experiment identifier (used to name the database file)"
    )

    args = parser.parse_args()

    # Validate exp_id (basic sanity check)
    if not args.exp_id or args.exp_id.isspace():
        print("Error: exp_id cannot be empty or whitespace", file=sys.stderr)
        sys.exit(1)

    # Create the database
    db_path = create_database(args.exp_id)
    print(f"Database created: {db_path}")

    # Print summary
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM Mods")
    mod_count = cursor.fetchone()[0]
    conn.close()

    print(f"  - Tables created: Reads, Mods, Exp, Refseq")
    print(f"  - Mods table populated with {mod_count} modification definitions")


if __name__ == "__main__":
    main()
