#!/usr/bin/env python3
"""
Database Initialization Script for ONT-SMA-seq Pipeline

Creates SQLite database schema and populates static lookup tables.
"""

import argparse
import sqlite3
import sys
from pathlib import Path


def create_tables(conn):
    """Create database schema tables."""
    cursor = conn.cursor()
    
    # Create Exp table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS Exp (
            exp_id TEXT PRIMARY KEY,
            exp_desc TEXT
        )
    ''')
    
    # Create Refseq table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS Refseq (
            refseq_id TEXT PRIMARY KEY,
            refseq TEXT NOT NULL,
            reflen INTEGER NOT NULL
        )
    ''')
    
    # Create Mods table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS Mods (
            mod_bitflag INTEGER PRIMARY KEY,
            mods TEXT NOT NULL
        )
    ''')
    
    # Create Reads table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS Reads (
            uniq_id TEXT PRIMARY KEY,
            exp_id TEXT,
            refseq_id TEXT,
            read_id TEXT NOT NULL,
            readseq TEXT NOT NULL,
            readlen INTEGER NOT NULL,
            model_tier TEXT NOT NULL,
            model_ver TEXT NOT NULL,
            trim INTEGER NOT NULL,
            mod_bitflag INTEGER,
            ed INTEGER,
            q_bc REAL NOT NULL,
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
    ''')
    
    conn.commit()


def populate_mods_table(conn):
    """Populate the Mods table with standard modification bitflags."""
    cursor = conn.cursor()
    
    # Modification bitflags as specified in the documentation
    # Note: 6mA (1) can coexist with C-mods, but C-mods (2,4,8,16) are mutually exclusive
    # Only including valid/common combinations
    modifications = [
        (0, 'non'),
        (1, '6mA'),
        (2, '5mCG_5hmCG'),
        (3, '6mA+5mCG_5hmCG'),
        (4, '5mC_5hmC'),
        (5, '6mA+5mC_5hmC'),
        (8, '4mC_5mC'),
        (9, '6mA+4mC_5mC'),
        (16, '5mC'),
        (17, '6mA+5mC'),
    ]
    
    cursor.executemany(
        'INSERT OR IGNORE INTO Mods (mod_bitflag, mods) VALUES (?, ?)',
        modifications
    )
    
    conn.commit()


def main():
    """Main entry point for database initialization."""
    parser = argparse.ArgumentParser(
        description='Initialize SQLite database for ONT-SMA-seq pipeline'
    )
    parser.add_argument(
        'exp_id',
        type=str,
        help='Experiment ID for database naming'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='.',
        help='Output directory for database file (default: current directory)'
    )
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Database file path
    db_path = output_dir / f'SMA_{args.exp_id}.db'
    
    # Check if database already exists
    if db_path.exists():
        print(f'Warning: Database {db_path} already exists. Overwriting...', file=sys.stderr)
        db_path.unlink()
    
    # Create database connection
    try:
        conn = sqlite3.connect(str(db_path))
        
        # Create tables
        create_tables(conn)
        print(f'Created database schema in {db_path}')
        
        # Populate static lookup tables
        populate_mods_table(conn)
        print('Populated Mods table with modification bitflags')
        
        conn.close()
        print(f'Successfully initialized database: {db_path}')
        return 0
        
    except sqlite3.Error as e:
        print(f'Database error: {e}', file=sys.stderr)
        return 1
    except Exception as e:
        print(f'Error: {e}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
