# Unified Nanopore Database Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a federated SQLite database system that unifies experiment discovery, basecalling tracking, and per-read SMA-seq metrics with experimental design specifications.

**Architecture:** Central metadata database (`nanopore_unified.db`) stores experiments, basecall runs, library specs, and lookup tables. Per-experiment databases (`SMA_{exp_id}.db`) store per-read metrics. Python library module (`lib/`) provides database operations. Orchestrator script coordinates the pipeline.

**Tech Stack:** Python 3.x, SQLite3, pysam, edlib, pandas (existing dependencies)

---

## Phase 1: Central Database Schema

### Task 1: Create lib/ directory and database initialization module

**Files:**
- Create: `lib/__init__.py`
- Create: `lib/db_schema.py`

**Step 1: Create lib directory and __init__.py**

```bash
mkdir -p lib
```

```python
# lib/__init__.py
"""ONT-SMA-seq database library."""

from .db_schema import create_central_db, create_sma_db

__all__ = ['create_central_db', 'create_sma_db']
```

**Step 2: Create db_schema.py with central database DDL**

```python
# lib/db_schema.py
"""Database schema definitions and creation functions."""

import sqlite3
from pathlib import Path

CENTRAL_SCHEMA = """
-- ============================================================
-- CORE TABLES
-- ============================================================

CREATE TABLE IF NOT EXISTS experiments (
    exp_id          TEXT PRIMARY KEY,
    experiment_path TEXT UNIQUE,
    instrument      TEXT,
    flow_cell_id    TEXT,
    sample_id       TEXT,
    protocol        TEXT,
    started         TEXT,
    pod5_count      INTEGER,
    pod5_dir        TEXT,
    created_at      TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS basecall_runs (
    run_id          TEXT PRIMARY KEY,
    exp_id          TEXT REFERENCES experiments(exp_id),
    model_tier      TEXT,
    model_version   TEXT,
    trim            INTEGER,
    mod_bitflag     INTEGER,
    dorado_version  TEXT,
    dorado_args     TEXT,
    batch_size      INTEGER,
    emit_moves      INTEGER,
    gpu_model       TEXT,
    slurm_job_id    TEXT,
    runtime_seconds REAL,
    bam_path        TEXT,
    sma_db_path     TEXT,
    status          TEXT DEFAULT 'pending',
    created_at      TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS references (
    ref_id    TEXT PRIMARY KEY,
    ref_path  TEXT,
    ref_seq   TEXT,
    ref_len   INTEGER
);

-- ============================================================
-- LIBRARY SPECIFICATION TABLES
-- ============================================================

CREATE TABLE IF NOT EXISTS library_specs (
    lib_id              TEXT PRIMARY KEY,
    lib_name            TEXT NOT NULL,
    lib_version         TEXT DEFAULT '1.0',
    description         TEXT,
    created_by          TEXT,
    project             TEXT,
    pore_type           TEXT,
    flow_cell_type      TEXT,
    sequencing_kit      TEXT,
    status              TEXT DEFAULT 'draft',
    created_at          TEXT DEFAULT CURRENT_TIMESTAMP,
    updated_at          TEXT
);

CREATE TABLE IF NOT EXISTS lib_constructs (
    construct_id        TEXT PRIMARY KEY,
    lib_id              TEXT REFERENCES library_specs(lib_id),
    name                TEXT NOT NULL,
    aliases             TEXT,
    construct_type      TEXT,
    sequence            TEXT,
    length              INTEGER,
    gc_content          REAL,
    source_file         TEXT,
    source_format       TEXT,
    genbank_accession   TEXT,
    annotations_json    TEXT,
    organism            TEXT,
    created_at          TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS lib_fragments (
    fragment_id         TEXT PRIMARY KEY,
    construct_id        TEXT REFERENCES lib_constructs(construct_id),
    name                TEXT NOT NULL,
    fragment_type       TEXT,
    sequence            TEXT,
    expected_len        INTEGER,
    len_min             INTEGER,
    len_max             INTEGER,
    len_tolerance       INTEGER DEFAULT 150,
    prep_method         TEXT,
    restriction_sites   TEXT,
    primer_fwd          TEXT,
    primer_rev          TEXT,
    expected_strand     TEXT,
    strand_ratio        REAL DEFAULT 0.5,
    adapter_5p          TEXT,
    adapter_3p          TEXT,
    ligation_method     TEXT,
    is_sequenceable     INTEGER DEFAULT 1,
    sequenceable_notes  TEXT,
    expected_fraction   REAL,
    created_at          TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS lib_modifications (
    mod_id              INTEGER PRIMARY KEY,
    fragment_id         TEXT REFERENCES lib_fragments(fragment_id),
    mod_type            TEXT NOT NULL,
    mod_code            INTEGER,
    expected_status     TEXT NOT NULL,
    expected_fraction   REAL,
    position_type       TEXT,
    positions_json      TEXT,
    motif               TEXT,
    modification_source TEXT,
    enzyme_used         TEXT,
    notes               TEXT
);

CREATE TABLE IF NOT EXISTS lib_barcodes (
    barcode_id          TEXT PRIMARY KEY,
    lib_id              TEXT REFERENCES library_specs(lib_id),
    barcode_kit         TEXT,
    barcode_name        TEXT,
    barcode_sequence    TEXT,
    sample_name         TEXT NOT NULL,
    sample_description  TEXT,
    expected_fragments  TEXT,
    expected_constructs TEXT,
    min_reads           INTEGER,
    min_purity          REAL,
    replicate_group     TEXT,
    condition           TEXT,
    notes               TEXT
);

CREATE TABLE IF NOT EXISTS lib_sample_prep (
    prep_id             TEXT PRIMARY KEY,
    lib_id              TEXT REFERENCES library_specs(lib_id),
    input_type          TEXT,
    extraction_method   TEXT,
    extraction_kit      TEXT,
    input_mass_ng       REAL,
    qubit_conc_ng_ul    REAL,
    nanodrop_260_280    REAL,
    nanodrop_260_230    REAL,
    size_selection      TEXT,
    size_cutoff_bp      INTEGER,
    fragmentation       TEXT,
    target_length_bp    INTEGER,
    library_kit         TEXT,
    library_protocol    TEXT,
    end_repair          INTEGER DEFAULT 1,
    a_tailing           INTEGER DEFAULT 1,
    adapter_ligation    INTEGER DEFAULT 1,
    final_mass_ng       REAL,
    final_conc_ng_ul    REAL,
    bioanalyzer_file    TEXT,
    prep_date           TEXT,
    prepped_by          TEXT,
    notes               TEXT
);

CREATE TABLE IF NOT EXISTS lib_qc_thresholds (
    lib_id              TEXT PRIMARY KEY REFERENCES library_specs(lib_id),
    min_mean_qscore     REAL DEFAULT 10.0,
    min_median_qscore   REAL DEFAULT 12.0,
    target_qscore       REAL DEFAULT 20.0,
    min_purity          REAL DEFAULT 0.8,
    target_purity       REAL DEFAULT 0.95,
    max_error_rate      REAL DEFAULT 0.10,
    target_error_rate   REAL DEFAULT 0.02,
    min_n50             INTEGER,
    target_n50          INTEGER,
    min_signal_positive REAL,
    max_unblock         REAL,
    min_reads           INTEGER,
    target_reads        INTEGER,
    min_bases           INTEGER,
    target_bases        INTEGER,
    custom_thresholds   TEXT,
    notes               TEXT
);

CREATE TABLE IF NOT EXISTS lib_custom_fields (
    id                  INTEGER PRIMARY KEY,
    lib_id              TEXT REFERENCES library_specs(lib_id),
    entity_type         TEXT,
    entity_id           TEXT,
    field_name          TEXT NOT NULL,
    field_value         TEXT,
    field_type          TEXT DEFAULT 'text',
    created_at          TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS lib_experiment_link (
    lib_id              TEXT REFERENCES library_specs(lib_id),
    exp_id              TEXT REFERENCES experiments(exp_id),
    relationship        TEXT DEFAULT 'primary',
    run_purpose         TEXT,
    notes               TEXT,
    PRIMARY KEY (lib_id, exp_id)
);

-- ============================================================
-- INDEXES
-- ============================================================

CREATE INDEX IF NOT EXISTS idx_basecall_runs_exp_id ON basecall_runs(exp_id);
CREATE INDEX IF NOT EXISTS idx_basecall_runs_status ON basecall_runs(status);
CREATE INDEX IF NOT EXISTS idx_lib_constructs_lib_id ON lib_constructs(lib_id);
CREATE INDEX IF NOT EXISTS idx_lib_fragments_construct_id ON lib_fragments(construct_id);
"""


def create_central_db(db_path: str) -> sqlite3.Connection:
    """Create central database with full schema.

    Args:
        db_path: Path to database file

    Returns:
        sqlite3.Connection to the created database
    """
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.executescript(CENTRAL_SCHEMA)
    conn.commit()
    return conn


SMA_SCHEMA = """
CREATE TABLE IF NOT EXISTS reads (
    read_id         TEXT,
    run_id          TEXT,
    readseq         TEXT,
    readlen         INTEGER,
    q_bc            REAL,
    ed              INTEGER,
    q_ld            REAL,
    ref_id          TEXT,
    end_reason      TEXT,
    PRIMARY KEY (read_id, run_id)
);

CREATE INDEX IF NOT EXISTS idx_reads_run_id ON reads(run_id);
CREATE INDEX IF NOT EXISTS idx_reads_ref_id ON reads(ref_id);
CREATE INDEX IF NOT EXISTS idx_reads_end_reason ON reads(end_reason);
CREATE INDEX IF NOT EXISTS idx_reads_ed ON reads(ed) WHERE ed IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_reads_qbc ON reads(q_bc);

CREATE TABLE IF NOT EXISTS run_summaries (
    run_id          TEXT PRIMARY KEY,
    total_reads     INTEGER,
    matched_reads   INTEGER,
    mean_qbc        REAL,
    median_qbc      REAL,
    mean_ed         REAL,
    median_ed       REAL,
    mean_qld        REAL,
    er_signal_pos   INTEGER,
    er_signal_neg   INTEGER,
    er_unblock      INTEGER,
    er_other        INTEGER
);
"""


def create_sma_db(db_path: str) -> sqlite3.Connection:
    """Create per-experiment SMA database.

    Args:
        db_path: Path to database file

    Returns:
        sqlite3.Connection to the created database
    """
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.executescript(SMA_SCHEMA)
    conn.commit()
    return conn
```

**Step 3: Verify module imports work**

```bash
cd /nfs/turbo/umms-atheylab/gregfar/SMS/ONT-SMA-seq/.worktrees/unified-db
python3 -c "from lib.db_schema import create_central_db, create_sma_db; print('OK')"
```

Expected: `OK`

**Step 4: Commit**

```bash
git add lib/
git commit -m "feat: add database schema module with central and SMA schemas"
```

---

### Task 2: Create lookup tables module

**Files:**
- Create: `lib/db_lookups.py`

**Step 1: Create lookup tables population module**

```python
# lib/db_lookups.py
"""Lookup table definitions and population functions."""

import sqlite3

LOOKUP_TABLES = """
-- ============================================================
-- PORE TYPES
-- ============================================================

CREATE TABLE IF NOT EXISTS lk_pore_types (
    pore_type       TEXT PRIMARY KEY,
    description     TEXT,
    release_year    INTEGER,
    is_current      INTEGER DEFAULT 1
);

INSERT OR IGNORE INTO lk_pore_types VALUES
    ('R9.4.1', 'Legacy pore, good accuracy', 2017, 0),
    ('R10.4.1', 'Current high-accuracy pore', 2022, 1),
    ('R10.4.2', 'Latest iteration', 2024, 1);

-- ============================================================
-- FLOW CELLS
-- ============================================================

CREATE TABLE IF NOT EXISTS lk_flow_cells (
    flow_cell_type  TEXT PRIMARY KEY,
    pore_type       TEXT,
    device          TEXT,
    pore_count      INTEGER,
    description     TEXT
);

INSERT OR IGNORE INTO lk_flow_cells VALUES
    ('FLO-FLG114', 'R10.4.1', 'Flongle', 126, 'Flongle R10.4.1'),
    ('FLO-MIN114', 'R10.4.1', 'MinION', 2048, 'MinION R10.4.1'),
    ('FLO-PRO114M', 'R10.4.1', 'PromethION', 2675, 'PromethION R10.4.1'),
    ('FLO-MIN106', 'R9.4.1', 'MinION', 2048, 'MinION R9.4.1 (legacy)');

-- ============================================================
-- SEQUENCING KITS
-- ============================================================

CREATE TABLE IF NOT EXISTS lk_sequencing_kits (
    kit_id          TEXT PRIMARY KEY,
    kit_name        TEXT,
    kit_type        TEXT,
    barcoding       INTEGER DEFAULT 0,
    barcode_count   INTEGER,
    pore_compat     TEXT,
    description     TEXT
);

INSERT OR IGNORE INTO lk_sequencing_kits VALUES
    ('SQK-LSK114', 'Ligation Sequencing Kit V14', 'ligation', 0, NULL, '["R10.4.1"]', 'Standard ligation prep'),
    ('SQK-NBD114.24', 'Native Barcoding Kit 24 V14', 'ligation', 1, 24, '["R10.4.1"]', '24-plex native barcoding'),
    ('SQK-NBD114.96', 'Native Barcoding Kit 96 V14', 'ligation', 1, 96, '["R10.4.1"]', '96-plex native barcoding'),
    ('SQK-RBK114.24', 'Rapid Barcoding Kit 24 V14', 'rapid', 1, 24, '["R10.4.1"]', '24-plex rapid barcoding'),
    ('SQK-RBK114.96', 'Rapid Barcoding Kit 96 V14', 'rapid', 1, 96, '["R10.4.1"]', '96-plex rapid barcoding'),
    ('SQK-RAD114', 'Rapid Sequencing Kit V14', 'rapid', 0, NULL, '["R10.4.1"]', 'Rapid transposase prep'),
    ('SQK-LSK109', 'Ligation Sequencing Kit', 'ligation', 0, NULL, '["R9.4.1"]', 'Legacy R9 ligation');

-- ============================================================
-- BASECALLING MODELS
-- ============================================================

CREATE TABLE IF NOT EXISTS lk_basecall_models (
    model_id        TEXT PRIMARY KEY,
    pore_type       TEXT,
    speed           TEXT,
    tier            TEXT,
    version         TEXT,
    dorado_min_ver  TEXT,
    description     TEXT
);

INSERT OR IGNORE INTO lk_basecall_models VALUES
    ('dna_r10.4.1_e8.2_400bps_fast@v5.0.0', 'R10.4.1', '400bps', 'fast', 'v5.0.0', '0.5.0', 'Fast model'),
    ('dna_r10.4.1_e8.2_400bps_hac@v5.0.0', 'R10.4.1', '400bps', 'hac', 'v5.0.0', '0.5.0', 'High accuracy'),
    ('dna_r10.4.1_e8.2_400bps_sup@v5.0.0', 'R10.4.1', '400bps', 'sup', 'v5.0.0', '0.5.0', 'Super accuracy'),
    ('dna_r10.4.1_e8.2_400bps_sup@v5.2.0', 'R10.4.1', '400bps', 'sup', 'v5.2.0', '0.8.0', 'Latest super accuracy');

-- ============================================================
-- MODIFICATIONS
-- ============================================================

CREATE TABLE IF NOT EXISTS lk_modifications (
    mod_code        INTEGER PRIMARY KEY,
    mod_type        TEXT,
    mod_name        TEXT,
    dorado_flag     TEXT,
    description     TEXT
);

INSERT OR IGNORE INTO lk_modifications VALUES
    (0, 'none', 'No modifications', NULL, 'Canonical bases only'),
    (1, '6mA', '6-methyladenine', '6mA', 'Bacterial/mitochondrial'),
    (2, '5mCG_5hmCG', '5mC and 5hmC in CpG', '5mCG_5hmCG', 'Mammalian CpG methylation'),
    (4, '5mC_5hmC', '5mC and 5hmC all contexts', '5mC_5hmC', 'All-context cytosine mods'),
    (8, '4mC_5mC', '4mC and 5mC', '4mC_5mC', 'Bacterial cytosine mods'),
    (16, '5mC', '5-methylcytosine only', '5mC', '5mC without hydroxymethyl'),
    (3, '6mA+5mCG_5hmCG', 'Combined', '6mA,5mCG_5hmCG', 'Adenine + CpG cytosine'),
    (5, '6mA+5mC_5hmC', 'Combined', '6mA,5mC_5hmC', 'Adenine + all-context cytosine');

-- ============================================================
-- CONSTRUCT TYPES
-- ============================================================

CREATE TABLE IF NOT EXISTS lk_construct_types (
    construct_type  TEXT PRIMARY KEY,
    description     TEXT
);

INSERT OR IGNORE INTO lk_construct_types VALUES
    ('plasmid', 'Circular plasmid DNA'),
    ('genome', 'Genomic DNA (bacterial, viral, etc.)'),
    ('amplicon', 'PCR amplification product'),
    ('synthetic', 'Synthesized DNA construct'),
    ('cdna', 'Complementary DNA from RNA'),
    ('chromosome', 'Eukaryotic chromosome/fragment');

-- ============================================================
-- FRAGMENT TYPES
-- ============================================================

CREATE TABLE IF NOT EXISTS lk_fragment_types (
    fragment_type   TEXT PRIMARY KEY,
    description     TEXT
);

INSERT OR IGNORE INTO lk_fragment_types VALUES
    ('whole_plasmid', 'Uncut circular plasmid (rolling circle)'),
    ('linearized', 'Single-cut linearized plasmid'),
    ('restriction_digest', 'Multi-enzyme restriction digest'),
    ('pcr_amplicon', 'PCR amplification product'),
    ('sheared', 'Mechanically sheared DNA'),
    ('enzymatic_fragment', 'Enzymatic fragmentation'),
    ('native', 'Native/unprocessed');

-- ============================================================
-- END REASONS
-- ============================================================

CREATE TABLE IF NOT EXISTS lk_end_reasons (
    end_reason      TEXT PRIMARY KEY,
    category        TEXT,
    is_good         INTEGER,
    description     TEXT
);

INSERT OR IGNORE INTO lk_end_reasons VALUES
    ('signal_positive', 'complete', 1, 'Read completed normally (good)'),
    ('signal_negative', 'complete', 1, 'Read completed, no signal'),
    ('unblock_mux_change', 'rejected', 0, 'Adaptive sampling rejection'),
    ('mux_change', 'technical', 0, 'Mux scan during read'),
    ('data_service_unblock_mux_change', 'rejected', 0, 'Software-triggered unblock'),
    ('unknown', 'unknown', 0, 'Unknown termination');

-- ============================================================
-- EXTRACTION METHODS
-- ============================================================

CREATE TABLE IF NOT EXISTS lk_extraction_methods (
    method          TEXT PRIMARY KEY,
    description     TEXT
);

INSERT OR IGNORE INTO lk_extraction_methods VALUES
    ('phenol_chloroform', 'Traditional organic extraction'),
    ('column_silica', 'Silica column purification'),
    ('magnetic_beads', 'Magnetic bead-based extraction'),
    ('alkaline_lysis', 'Alkaline lysis (plasmid miniprep)'),
    ('kit_qiagen', 'Qiagen kit extraction'),
    ('kit_zymo', 'Zymo Research kit extraction'),
    ('kit_monarch', 'NEB Monarch kit extraction'),
    ('direct', 'Direct input (no extraction)');

-- ============================================================
-- SIZE SELECTION METHODS
-- ============================================================

CREATE TABLE IF NOT EXISTS lk_size_selection (
    method          TEXT PRIMARY KEY,
    description     TEXT
);

INSERT OR IGNORE INTO lk_size_selection VALUES
    ('none', 'No size selection'),
    ('ampure_0.4x', 'AMPure XP 0.4x ratio (>1kb)'),
    ('ampure_0.5x', 'AMPure XP 0.5x ratio (>500bp)'),
    ('ampure_0.6x', 'AMPure XP 0.6x ratio (>300bp)'),
    ('ampure_0.8x', 'AMPure XP 0.8x ratio (cleanup)'),
    ('bluepippin', 'BluePippin gel extraction'),
    ('gel_extraction', 'Manual gel extraction'),
    ('sspei', 'Short fragment eliminator');

-- ============================================================
-- LIBRARY STATUSES
-- ============================================================

CREATE TABLE IF NOT EXISTS lk_library_status (
    status          TEXT PRIMARY KEY,
    description     TEXT
);

INSERT OR IGNORE INTO lk_library_status VALUES
    ('draft', 'In design, not finalized'),
    ('active', 'Currently in use'),
    ('archived', 'No longer active'),
    ('deprecated', 'Replaced by newer version');

-- ============================================================
-- RUN PURPOSES
-- ============================================================

CREATE TABLE IF NOT EXISTS lk_run_purposes (
    purpose         TEXT PRIMARY KEY,
    description     TEXT
);

INSERT OR IGNORE INTO lk_run_purposes VALUES
    ('production', 'Production data collection'),
    ('test', 'Test/validation run'),
    ('optimization', 'Protocol optimization'),
    ('troubleshooting', 'Debugging issues'),
    ('training', 'Training/demonstration');
"""


def populate_lookups(conn: sqlite3.Connection) -> None:
    """Populate lookup tables with predefined values.

    Args:
        conn: Database connection
    """
    conn.executescript(LOOKUP_TABLES)
    conn.commit()


def get_valid_options(conn: sqlite3.Connection, table: str) -> list:
    """Get valid options from a lookup table.

    Args:
        conn: Database connection
        table: Lookup table name (e.g., 'lk_pore_types')

    Returns:
        List of primary key values
    """
    cursor = conn.execute(f"SELECT * FROM {table}")
    return cursor.fetchall()
```

**Step 2: Update lib/__init__.py**

```python
# lib/__init__.py
"""ONT-SMA-seq database library."""

from .db_schema import create_central_db, create_sma_db
from .db_lookups import populate_lookups, get_valid_options

__all__ = [
    'create_central_db',
    'create_sma_db',
    'populate_lookups',
    'get_valid_options'
]
```

**Step 3: Verify lookups work**

```bash
python3 -c "
from lib.db_schema import create_central_db
from lib.db_lookups import populate_lookups, get_valid_options
import tempfile
import os

with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
    db_path = f.name

conn = create_central_db(db_path)
populate_lookups(conn)
kits = get_valid_options(conn, 'lk_sequencing_kits')
print(f'Loaded {len(kits)} sequencing kits')
conn.close()
os.unlink(db_path)
"
```

Expected: `Loaded 7 sequencing kits`

**Step 4: Commit**

```bash
git add lib/
git commit -m "feat: add lookup tables with predefined ONT options"
```

---

### Task 3: Create database initialization CLI

**Files:**
- Create: `bin/init_central_db.py`

**Step 1: Create initialization script**

```python
#!/usr/bin/env python3
# bin/init_central_db.py
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
```

**Step 2: Test the script**

```bash
python3 bin/init_central_db.py -o /tmp/test_unified.db
sqlite3 /tmp/test_unified.db "SELECT COUNT(*) FROM lk_sequencing_kits"
rm /tmp/test_unified.db
```

Expected:
- Script outputs table/lookup counts
- Query returns `7`

**Step 3: Commit**

```bash
git add bin/init_central_db.py
git commit -m "feat: add central database initialization CLI"
```

---

## Phase 2: Database Operations Module

### Task 4: Create database operations module

**Files:**
- Create: `lib/db_ops.py`

**Step 1: Create operations module**

```python
# lib/db_ops.py
"""Database CRUD operations for the unified nanopore database."""

import sqlite3
from typing import Optional, Dict, List, Any
from pathlib import Path


class CentralDB:
    """Operations for the central nanopore_unified.db database."""

    def __init__(self, db_path: str):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row

    def close(self):
        self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    # --- Experiments ---

    def insert_experiment(self, exp_id: str, **kwargs) -> None:
        """Insert or update an experiment record."""
        cols = ['exp_id'] + list(kwargs.keys())
        placeholders = ','.join(['?'] * len(cols))
        values = [exp_id] + list(kwargs.values())

        self.conn.execute(
            f"INSERT OR REPLACE INTO experiments ({','.join(cols)}) VALUES ({placeholders})",
            values
        )
        self.conn.commit()

    def get_experiment(self, exp_id: str) -> Optional[Dict]:
        """Get experiment by ID."""
        cursor = self.conn.execute(
            "SELECT * FROM experiments WHERE exp_id = ?", (exp_id,)
        )
        row = cursor.fetchone()
        return dict(row) if row else None

    def list_experiments(self) -> List[Dict]:
        """List all experiments."""
        cursor = self.conn.execute("SELECT * FROM experiments ORDER BY created_at DESC")
        return [dict(row) for row in cursor.fetchall()]

    # --- Basecall Runs ---

    def insert_basecall_run(self, run_id: str, exp_id: str, **kwargs) -> None:
        """Insert a basecall run record."""
        cols = ['run_id', 'exp_id'] + list(kwargs.keys())
        placeholders = ','.join(['?'] * len(cols))
        values = [run_id, exp_id] + list(kwargs.values())

        self.conn.execute(
            f"INSERT OR REPLACE INTO basecall_runs ({','.join(cols)}) VALUES ({placeholders})",
            values
        )
        self.conn.commit()

    def get_basecall_run(self, run_id: str) -> Optional[Dict]:
        """Get basecall run by ID."""
        cursor = self.conn.execute(
            "SELECT * FROM basecall_runs WHERE run_id = ?", (run_id,)
        )
        row = cursor.fetchone()
        return dict(row) if row else None

    def list_basecall_runs(self, exp_id: Optional[str] = None,
                          status: Optional[str] = None) -> List[Dict]:
        """List basecall runs with optional filters."""
        query = "SELECT * FROM basecall_runs WHERE 1=1"
        params = []

        if exp_id:
            query += " AND exp_id = ?"
            params.append(exp_id)
        if status:
            query += " AND status = ?"
            params.append(status)

        query += " ORDER BY created_at DESC"
        cursor = self.conn.execute(query, params)
        return [dict(row) for row in cursor.fetchall()]

    def update_basecall_run_status(self, run_id: str, status: str,
                                   sma_db_path: Optional[str] = None) -> None:
        """Update basecall run status and optionally SMA DB path."""
        if sma_db_path:
            self.conn.execute(
                "UPDATE basecall_runs SET status = ?, sma_db_path = ? WHERE run_id = ?",
                (status, sma_db_path, run_id)
            )
        else:
            self.conn.execute(
                "UPDATE basecall_runs SET status = ? WHERE run_id = ?",
                (status, run_id)
            )
        self.conn.commit()

    # --- Library Specs ---

    def insert_library_spec(self, lib_id: str, lib_name: str, **kwargs) -> None:
        """Insert a library specification."""
        cols = ['lib_id', 'lib_name'] + list(kwargs.keys())
        placeholders = ','.join(['?'] * len(cols))
        values = [lib_id, lib_name] + list(kwargs.values())

        self.conn.execute(
            f"INSERT OR REPLACE INTO library_specs ({','.join(cols)}) VALUES ({placeholders})",
            values
        )
        self.conn.commit()

    def get_library_spec(self, lib_id: str) -> Optional[Dict]:
        """Get library spec by ID."""
        cursor = self.conn.execute(
            "SELECT * FROM library_specs WHERE lib_id = ?", (lib_id,)
        )
        row = cursor.fetchone()
        return dict(row) if row else None

    def link_library_experiment(self, lib_id: str, exp_id: str,
                                relationship: str = 'primary',
                                run_purpose: Optional[str] = None) -> None:
        """Link a library spec to an experiment."""
        self.conn.execute(
            """INSERT OR REPLACE INTO lib_experiment_link
               (lib_id, exp_id, relationship, run_purpose) VALUES (?, ?, ?, ?)""",
            (lib_id, exp_id, relationship, run_purpose)
        )
        self.conn.commit()

    # --- Pending Work Queries ---

    def get_pending_sma_analysis(self) -> List[Dict]:
        """Get basecall runs that are complete but not yet analyzed."""
        cursor = self.conn.execute("""
            SELECT br.*, e.pod5_dir, e.sample_id
            FROM basecall_runs br
            JOIN experiments e ON br.exp_id = e.exp_id
            WHERE br.status = 'complete' AND br.sma_db_path IS NULL
            ORDER BY br.created_at
        """)
        return [dict(row) for row in cursor.fetchall()]


class SMADB:
    """Operations for per-experiment SMA_{exp_id}.db databases."""

    def __init__(self, db_path: str):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        # Performance tuning for bulk inserts
        self.conn.execute("PRAGMA synchronous = OFF")
        self.conn.execute("PRAGMA journal_mode = MEMORY")

    def close(self):
        self.conn.commit()
        self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def insert_read(self, read_id: str, run_id: str, **kwargs) -> None:
        """Insert a single read record."""
        cols = ['read_id', 'run_id'] + list(kwargs.keys())
        placeholders = ','.join(['?'] * len(cols))
        values = [read_id, run_id] + list(kwargs.values())

        self.conn.execute(
            f"INSERT OR REPLACE INTO reads ({','.join(cols)}) VALUES ({placeholders})",
            values
        )

    def insert_reads_batch(self, reads: List[tuple]) -> None:
        """Insert multiple reads in a batch.

        Args:
            reads: List of tuples matching reads table columns
        """
        self.conn.executemany(
            """INSERT OR REPLACE INTO reads
               (read_id, run_id, readseq, readlen, q_bc, ed, q_ld, ref_id, end_reason)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            reads
        )

    def insert_run_summary(self, run_id: str, **kwargs) -> None:
        """Insert or update run summary."""
        cols = ['run_id'] + list(kwargs.keys())
        placeholders = ','.join(['?'] * len(cols))
        values = [run_id] + list(kwargs.values())

        self.conn.execute(
            f"INSERT OR REPLACE INTO run_summaries ({','.join(cols)}) VALUES ({placeholders})",
            values
        )
        self.conn.commit()

    def get_run_summary(self, run_id: str) -> Optional[Dict]:
        """Get summary for a specific run."""
        cursor = self.conn.execute(
            "SELECT * FROM run_summaries WHERE run_id = ?", (run_id,)
        )
        row = cursor.fetchone()
        return dict(row) if row else None

    def compute_run_summary(self, run_id: str) -> Dict:
        """Compute summary statistics for a run from reads table."""
        cursor = self.conn.execute("""
            SELECT
                COUNT(*) as total_reads,
                SUM(CASE WHEN ref_id IS NOT NULL THEN 1 ELSE 0 END) as matched_reads,
                AVG(q_bc) as mean_qbc,
                AVG(ed) as mean_ed,
                AVG(q_ld) as mean_qld,
                SUM(CASE WHEN end_reason = 'signal_positive' THEN 1 ELSE 0 END) as er_signal_pos,
                SUM(CASE WHEN end_reason = 'signal_negative' THEN 1 ELSE 0 END) as er_signal_neg,
                SUM(CASE WHEN end_reason = 'unblock_mux_change' THEN 1 ELSE 0 END) as er_unblock
            FROM reads
            WHERE run_id = ?
        """, (run_id,))
        row = cursor.fetchone()
        return dict(row) if row else {}
```

**Step 2: Update lib/__init__.py**

```python
# lib/__init__.py
"""ONT-SMA-seq database library."""

from .db_schema import create_central_db, create_sma_db
from .db_lookups import populate_lookups, get_valid_options
from .db_ops import CentralDB, SMADB

__all__ = [
    'create_central_db',
    'create_sma_db',
    'populate_lookups',
    'get_valid_options',
    'CentralDB',
    'SMADB'
]
```

**Step 3: Test operations**

```bash
python3 -c "
from lib import create_central_db, populate_lookups, CentralDB
import tempfile
import os

with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
    db_path = f.name

conn = create_central_db(db_path)
populate_lookups(conn)
conn.close()

with CentralDB(db_path) as db:
    db.insert_experiment('test_exp', sample_id='sample1', instrument='P2S')
    exp = db.get_experiment('test_exp')
    print(f'Experiment: {exp[\"exp_id\"]} - {exp[\"sample_id\"]}')

    db.insert_basecall_run('test_run', 'test_exp', model_tier='sup', status='complete')
    runs = db.list_basecall_runs(status='complete')
    print(f'Complete runs: {len(runs)}')

os.unlink(db_path)
"
```

Expected:
```
Experiment: test_exp - sample1
Complete runs: 1
```

**Step 4: Commit**

```bash
git add lib/
git commit -m "feat: add database operations module with CentralDB and SMADB classes"
```

---

## Phase 3: Migration from Existing Database

### Task 5: Create migration script from nanopore_experiments.db

**Files:**
- Create: `bin/migrate_experiments.py`

**Step 1: Create migration script**

```python
#!/usr/bin/env python3
# bin/migrate_experiments.py
"""Migrate data from existing nanopore_experiments.db to unified database."""

import argparse
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from lib import CentralDB


def migrate_experiments(source_db: str, target_db: str, dry_run: bool = False):
    """Migrate experiments from source to target database.

    Args:
        source_db: Path to existing nanopore_experiments.db
        target_db: Path to new nanopore_unified.db
        dry_run: If True, don't actually insert data
    """
    source_conn = sqlite3.connect(source_db)
    source_conn.row_factory = sqlite3.Row

    # Get experiments from source
    cursor = source_conn.execute("""
        SELECT
            e.id,
            e.experiment_path,
            e.unique_id,
            e.instrument,
            e.flow_cell_id,
            e.sample_id,
            e.protocol,
            e.started,
            e.pod5_files as pod5_count
        FROM experiments e
    """)

    experiments = cursor.fetchall()
    print(f"[migrate] Found {len(experiments)} experiments in source database")

    if dry_run:
        print("[migrate] Dry run - not inserting data")
        for exp in experiments[:5]:
            print(f"  - {exp['sample_id'] or exp['unique_id']}: {exp['experiment_path']}")
        if len(experiments) > 5:
            print(f"  ... and {len(experiments) - 5} more")
        source_conn.close()
        return

    with CentralDB(target_db) as db:
        migrated = 0
        for exp in experiments:
            # Generate exp_id from unique_id or path
            exp_id = exp['unique_id'] or Path(exp['experiment_path']).name

            # Extract pod5_dir from experiment_path
            pod5_dir = str(Path(exp['experiment_path']))

            db.insert_experiment(
                exp_id=exp_id,
                experiment_path=exp['experiment_path'],
                instrument=exp['instrument'],
                flow_cell_id=exp['flow_cell_id'],
                sample_id=exp['sample_id'],
                protocol=exp['protocol'],
                started=exp['started'],
                pod5_count=exp['pod5_count'],
                pod5_dir=pod5_dir
            )
            migrated += 1

        print(f"[migrate] Migrated {migrated} experiments")

    source_conn.close()


def main():
    parser = argparse.ArgumentParser(
        description="Migrate experiments from nanopore_experiments.db to unified database"
    )
    parser.add_argument(
        "-s", "--source",
        default="/nfs/turbo/umms-atheylab/nanopore_experiments.db",
        help="Source database path [%(default)s]"
    )
    parser.add_argument(
        "-t", "--target",
        default="nanopore_unified.db",
        help="Target database path [%(default)s]"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be migrated without actually doing it"
    )
    args = parser.parse_args()

    if not Path(args.source).exists():
        print(f"[migrate] Error: Source database not found: {args.source}")
        sys.exit(1)

    if not Path(args.target).exists():
        print(f"[migrate] Error: Target database not found: {args.target}")
        print("[migrate] Run bin/init_central_db.py first")
        sys.exit(1)

    migrate_experiments(args.source, args.target, args.dry_run)


if __name__ == "__main__":
    main()
```

**Step 2: Test with dry run**

```bash
python3 bin/migrate_experiments.py --dry-run
```

Expected: Lists first 5 experiments without inserting

**Step 3: Commit**

```bash
git add bin/migrate_experiments.py
git commit -m "feat: add migration script for existing nanopore_experiments.db"
```

---

## Phase 4: Integration with Existing Scripts

### Task 6: Update ingest.py to use unified database

**Files:**
- Modify: `bin/ingest.py`

**Step 1: Add central database integration to ingest.py**

Add after the existing imports at the top of `bin/ingest.py`:

```python
# Add after existing imports
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Try to import unified DB (optional integration)
try:
    from lib import CentralDB, SMADB, create_sma_db
    UNIFIED_DB_AVAILABLE = True
except ImportError:
    UNIFIED_DB_AVAILABLE = False
```

Add new arguments to argparse section:

```python
parser.add_argument("--central-db",
    help="Path to central unified database (enables integration)")
parser.add_argument("--run-id",
    help="Basecall run ID for tracking (required if --central-db is set)")
```

Add after processing loop (before final print statements):

```python
# Update central database if provided
if args.central_db and UNIFIED_DB_AVAILABLE:
    if not args.run_id:
        print("[ingest] Warning: --central-db requires --run-id. Skipping central DB update.")
    else:
        print(f"[ingest] Updating central database: {args.central_db}")
        with CentralDB(args.central_db) as central:
            central.update_basecall_run_status(
                args.run_id,
                status='analyzed',
                sma_db_path=str(DB_PATH)
            )
        print(f"[ingest] Updated run {args.run_id} status to 'analyzed'")
```

**Step 2: Test that existing functionality still works**

```bash
python3 bin/ingest.py --help
```

Expected: Help shows new --central-db and --run-id options

**Step 3: Commit**

```bash
git add bin/ingest.py
git commit -m "feat: add unified database integration to ingest.py"
```

---

### Task 7: Create orchestrator script

**Files:**
- Create: `bin/orchestrate.py`

**Step 1: Create orchestrator**

```python
#!/usr/bin/env python3
# bin/orchestrate.py
"""Orchestrate the unified nanopore pipeline.

Watches for:
1. New experiments → ready for basecalling
2. Completed basecall runs → ready for SMA-seq analysis
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from lib import CentralDB


def show_status(db_path: str):
    """Show current pipeline status."""
    with CentralDB(db_path) as db:
        experiments = db.list_experiments()
        print(f"\n{'='*60}")
        print("PIPELINE STATUS")
        print(f"{'='*60}")
        print(f"\nExperiments: {len(experiments)}")

        # Count by status
        all_runs = db.list_basecall_runs()
        pending = len([r for r in all_runs if r['status'] == 'pending'])
        running = len([r for r in all_runs if r['status'] == 'running'])
        complete = len([r for r in all_runs if r['status'] == 'complete'])
        analyzed = len([r for r in all_runs if r['status'] == 'analyzed'])

        print(f"\nBasecall Runs:")
        print(f"  Pending:  {pending}")
        print(f"  Running:  {running}")
        print(f"  Complete: {complete}")
        print(f"  Analyzed: {analyzed}")

        # Show pending SMA analysis
        pending_sma = db.get_pending_sma_analysis()
        if pending_sma:
            print(f"\nReady for SMA Analysis ({len(pending_sma)}):")
            for run in pending_sma[:5]:
                print(f"  - {run['run_id']}: {run['sample_id']} ({run['model_tier']})")
            if len(pending_sma) > 5:
                print(f"  ... and {len(pending_sma) - 5} more")


def list_pending_sma(db_path: str):
    """List basecall runs ready for SMA-seq analysis."""
    with CentralDB(db_path) as db:
        pending = db.get_pending_sma_analysis()

        if not pending:
            print("No basecall runs pending SMA analysis.")
            return

        print(f"{'Run ID':<40} {'Exp ID':<20} {'Tier':<6} {'BAM Path'}")
        print("-" * 100)
        for run in pending:
            print(f"{run['run_id']:<40} {run['exp_id']:<20} {run['model_tier']:<6} {run['bam_path'] or 'N/A'}")


def main():
    parser = argparse.ArgumentParser(
        description="Orchestrate the unified nanopore pipeline"
    )
    parser.add_argument(
        "-d", "--database",
        default="nanopore_unified.db",
        help="Path to central database [%(default)s]"
    )

    subparsers = parser.add_subparsers(dest='command', help='Commands')

    # Status command
    subparsers.add_parser('status', help='Show pipeline status')

    # List pending SMA
    subparsers.add_parser('pending-sma', help='List runs ready for SMA analysis')

    args = parser.parse_args()

    if not Path(args.database).exists():
        print(f"[orchestrate] Error: Database not found: {args.database}")
        sys.exit(1)

    if args.command == 'status':
        show_status(args.database)
    elif args.command == 'pending-sma':
        list_pending_sma(args.database)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
```

**Step 2: Test orchestrator**

```bash
# First create a test database
python3 bin/init_central_db.py -o /tmp/test_orch.db

# Then test status command
python3 bin/orchestrate.py -d /tmp/test_orch.db status

# Cleanup
rm /tmp/test_orch.db
```

Expected: Shows status with 0 experiments, 0 runs

**Step 3: Commit**

```bash
git add bin/orchestrate.py
git commit -m "feat: add pipeline orchestrator script"
```

---

## Phase 5: Documentation and Final Integration

### Task 8: Update CLAUDE.md with unified database info

**Files:**
- Modify: `CLAUDE.md`

**Step 1: Add unified database section to CLAUDE.md**

Add the following section after the "Dependencies" section:

```markdown
## Unified Database System

The unified database integrates experiment discovery, basecalling tracking, and SMA-seq metrics.

### Database Files

- **`nanopore_unified.db`** - Central metadata (experiments, basecall runs, library specs)
- **`SMA_{exp_id}.db`** - Per-experiment read data (created by ingest.py)

### Commands

```bash
# Initialize central database
python bin/init_central_db.py -o nanopore_unified.db

# Migrate existing experiments
python bin/migrate_experiments.py -s /path/to/nanopore_experiments.db -t nanopore_unified.db

# Check pipeline status
python bin/orchestrate.py -d nanopore_unified.db status

# List runs ready for SMA analysis
python bin/orchestrate.py -d nanopore_unified.db pending-sma

# Run ingest with central DB integration
python bin/ingest.py \
  -e <exp_id> \
  -b input.bam \
  -s summary.tsv \
  -r reference.fa \
  -d Output/SMA_<exp_id>.db \
  -o Output/<exp_id>_tagged.bam \
  --central-db nanopore_unified.db \
  --run-id <run_id>
```

### Library Module

```python
from lib import CentralDB, SMADB, create_central_db, create_sma_db

# Central database operations
with CentralDB('nanopore_unified.db') as db:
    db.insert_experiment('exp1', sample_id='sample1')
    db.insert_basecall_run('run1', 'exp1', model_tier='sup')
    pending = db.get_pending_sma_analysis()

# Per-experiment SMA database
with SMADB('SMA_exp1.db') as sma:
    sma.insert_read('read1', 'run1', q_bc=25.0, ed=5)
    summary = sma.compute_run_summary('run1')
```
```

**Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: add unified database documentation to CLAUDE.md"
```

---

### Task 9: Final verification and merge preparation

**Step 1: Verify all modules work together**

```bash
# Full integration test
python3 -c "
from lib import (
    create_central_db, create_sma_db,
    populate_lookups, get_valid_options,
    CentralDB, SMADB
)
import tempfile
import os

# Create central DB
central_path = tempfile.mktemp(suffix='.db')
conn = create_central_db(central_path)
populate_lookups(conn)
conn.close()

# Test operations
with CentralDB(central_path) as db:
    db.insert_experiment('test_exp', sample_id='test_sample')
    db.insert_basecall_run('test_run', 'test_exp', model_tier='sup', status='complete')
    pending = db.get_pending_sma_analysis()
    assert len(pending) == 1, 'Should have 1 pending run'

# Create SMA DB
sma_path = tempfile.mktemp(suffix='.db')
conn = create_sma_db(sma_path)
conn.close()

with SMADB(sma_path) as sma:
    sma.insert_read('read1', 'test_run', readseq='ACGT', readlen=4, q_bc=20.0)
    sma.conn.commit()
    summary = sma.compute_run_summary('test_run')
    assert summary['total_reads'] == 1, 'Should have 1 read'

# Update central DB
with CentralDB(central_path) as db:
    db.update_basecall_run_status('test_run', 'analyzed', sma_path)
    pending = db.get_pending_sma_analysis()
    assert len(pending) == 0, 'Should have 0 pending runs'

os.unlink(central_path)
os.unlink(sma_path)

print('All integration tests passed!')
"
```

Expected: `All integration tests passed!`

**Step 2: Review changes**

```bash
git log --oneline feature/unified-database ^main
```

**Step 3: Commit any final changes**

```bash
git status
# If any uncommitted changes:
git add -A
git commit -m "chore: final cleanup for unified database feature"
```

---

## Summary

After completing all tasks, the unified database system provides:

1. **Central database** (`nanopore_unified.db`) with:
   - Experiments table (from existing nanopore_experiments.db)
   - Basecall runs tracking (dorado-bench integration)
   - Library specifications (experimental design)
   - 11 lookup tables with predefined ONT options

2. **Per-experiment databases** (`SMA_{exp_id}.db`) with:
   - Per-read metrics
   - Run summaries

3. **CLI tools**:
   - `bin/init_central_db.py` - Initialize central database
   - `bin/migrate_experiments.py` - Migrate from existing DB
   - `bin/orchestrate.py` - Pipeline status and coordination
   - `bin/ingest.py` - Updated with central DB integration

4. **Python library** (`lib/`):
   - `db_schema.py` - Schema definitions
   - `db_lookups.py` - Lookup table management
   - `db_ops.py` - CRUD operations

To merge: `git checkout main && git merge feature/unified-database`
