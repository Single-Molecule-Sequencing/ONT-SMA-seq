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

CREATE TABLE IF NOT EXISTS "references" (
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
