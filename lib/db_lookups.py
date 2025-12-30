"""Lookup tables with predefined ONT options."""

import sqlite3
from typing import List, Tuple, Any

LOOKUP_TABLES = """
-- ============================================================
-- LOOKUP TABLES - Predefined ONT Options
-- ============================================================

-- Pore types
CREATE TABLE IF NOT EXISTS lk_pore_types (
    pore_type   TEXT PRIMARY KEY,
    description TEXT
);

INSERT OR IGNORE INTO lk_pore_types (pore_type, description) VALUES
    ('R9.4.1', 'Legacy pore, widely used through 2023'),
    ('R10.4.1', 'Current standard pore with improved accuracy'),
    ('R10.4.2', 'Latest pore revision with enhanced performance');

-- Flow cell types
CREATE TABLE IF NOT EXISTS lk_flow_cells (
    flow_cell_type  TEXT PRIMARY KEY,
    pore_type       TEXT,
    device_type     TEXT,
    description     TEXT
);

INSERT OR IGNORE INTO lk_flow_cells (flow_cell_type, pore_type, device_type, description) VALUES
    ('FLO-FLG114', 'R10.4.1', 'Flongle', 'Flongle flow cell, R10.4.1 pores'),
    ('FLO-MIN114', 'R10.4.1', 'MinION', 'MinION flow cell, R10.4.1 pores'),
    ('FLO-PRO114M', 'R10.4.1', 'PromethION', 'PromethION flow cell, R10.4.1 pores'),
    ('FLO-MIN106', 'R9.4.1', 'MinION', 'MinION flow cell, R9.4.1 pores (legacy)');

-- Sequencing kits
CREATE TABLE IF NOT EXISTS lk_sequencing_kits (
    kit_id      TEXT PRIMARY KEY,
    kit_type    TEXT,
    barcoded    INTEGER DEFAULT 0,
    barcode_count INTEGER,
    description TEXT
);

INSERT OR IGNORE INTO lk_sequencing_kits (kit_id, kit_type, barcoded, barcode_count, description) VALUES
    ('SQK-LSK114', 'ligation', 0, NULL, 'Ligation Sequencing Kit V14'),
    ('SQK-NBD114.24', 'native_barcoding', 1, 24, 'Native Barcoding Kit 24 V14'),
    ('SQK-NBD114.96', 'native_barcoding', 1, 96, 'Native Barcoding Kit 96 V14'),
    ('SQK-RBK114.24', 'rapid_barcoding', 1, 24, 'Rapid Barcoding Kit 24 V14'),
    ('SQK-RBK114.96', 'rapid_barcoding', 1, 96, 'Rapid Barcoding Kit 96 V14'),
    ('SQK-RAD114', 'rapid', 0, NULL, 'Rapid Sequencing Kit V14'),
    ('SQK-LSK109', 'ligation', 0, NULL, 'Ligation Sequencing Kit (legacy R9.4.1)');

-- Basecall models
CREATE TABLE IF NOT EXISTS lk_basecall_models (
    model_name      TEXT PRIMARY KEY,
    pore_type       TEXT,
    chemistry       TEXT,
    speed           TEXT,
    tier            TEXT,
    modifications   TEXT,
    description     TEXT
);

INSERT OR IGNORE INTO lk_basecall_models (model_name, pore_type, chemistry, speed, tier, modifications, description) VALUES
    ('dna_r10.4.1_e8.2_400bps_fast@v5.0.0', 'R10.4.1', 'e8.2', '400bps', 'fast', NULL, 'Fast model for R10.4.1'),
    ('dna_r10.4.1_e8.2_400bps_hac@v5.0.0', 'R10.4.1', 'e8.2', '400bps', 'hac', NULL, 'High accuracy model for R10.4.1'),
    ('dna_r10.4.1_e8.2_400bps_sup@v5.0.0', 'R10.4.1', 'e8.2', '400bps', 'sup', NULL, 'Super accuracy model for R10.4.1'),
    ('dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG', 'R10.4.1', 'e8.2', '400bps', 'hac', '5mCG_5hmCG', 'HAC with 5mC and 5hmC calling'),
    ('dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG', 'R10.4.1', 'e8.2', '400bps', 'sup', '5mCG_5hmCG', 'SUP with 5mC and 5hmC calling'),
    ('dna_r10.4.1_e8.2_400bps_hac@v5.0.0_4mC_5mC', 'R10.4.1', 'e8.2', '400bps', 'hac', '4mC_5mC', 'HAC with 4mC and 5mC calling'),
    ('dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC', 'R10.4.1', 'e8.2', '400bps', 'sup', '4mC_5mC', 'SUP with 4mC and 5mC calling'),
    ('dna_r10.4.1_e8.2_400bps_hac@v5.0.0_6mA', 'R10.4.1', 'e8.2', '400bps', 'hac', '6mA', 'HAC with 6mA calling'),
    ('dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA', 'R10.4.1', 'e8.2', '400bps', 'sup', '6mA', 'SUP with 6mA calling');

-- Modification bitflags
CREATE TABLE IF NOT EXISTS lk_modifications (
    mod_code    INTEGER PRIMARY KEY,
    mod_name    TEXT NOT NULL,
    description TEXT
);

INSERT OR IGNORE INTO lk_modifications (mod_code, mod_name, description) VALUES
    (0, 'none', 'No modifications called'),
    (1, '5mC', '5-methylcytosine'),
    (2, '5hmC', '5-hydroxymethylcytosine'),
    (4, '6mA', 'N6-methyladenine'),
    (8, '4mC', 'N4-methylcytosine'),
    (16, '5fC', '5-formylcytosine'),
    (3, '5mCG_5hmCG', '5mC and 5hmC in CpG context (1+2)'),
    (5, '5mC_6mA', '5mC and 6mA combined (1+4)');

-- Construct types
CREATE TABLE IF NOT EXISTS lk_construct_types (
    construct_type  TEXT PRIMARY KEY,
    description     TEXT
);

INSERT OR IGNORE INTO lk_construct_types (construct_type, description) VALUES
    ('plasmid', 'Circular plasmid DNA'),
    ('genome', 'Genomic DNA'),
    ('amplicon', 'PCR amplification product'),
    ('synthetic', 'Synthetic DNA construct'),
    ('cdna', 'Complementary DNA from RNA'),
    ('chromosome', 'Full or partial chromosome');

-- Fragment types
CREATE TABLE IF NOT EXISTS lk_fragment_types (
    fragment_type   TEXT PRIMARY KEY,
    description     TEXT
);

INSERT OR IGNORE INTO lk_fragment_types (fragment_type, description) VALUES
    ('whole_plasmid', 'Complete circular plasmid'),
    ('linearized', 'Linearized plasmid by single cut'),
    ('restriction_digest', 'Fragment from restriction enzyme digestion'),
    ('pcr_amplicon', 'PCR amplification product'),
    ('sheared', 'Mechanically sheared fragment'),
    ('enzymatic_fragment', 'Fragment from enzymatic processing'),
    ('native', 'Native unprocessed DNA');

-- End reasons (from POD5 metadata)
CREATE TABLE IF NOT EXISTS lk_end_reasons (
    end_reason  TEXT PRIMARY KEY,
    category    TEXT,
    description TEXT
);

INSERT OR IGNORE INTO lk_end_reasons (end_reason, category, description) VALUES
    ('signal_positive', 'complete', 'Read completed normally on positive strand'),
    ('signal_negative', 'complete', 'Read completed normally on negative strand'),
    ('unblock_mux_change', 'rejected', 'Read rejected and unblocked'),
    ('mux_change', 'interrupted', 'Multiplexer channel change'),
    ('data_service_unblock_mux_change', 'rejected', 'Data service triggered unblock'),
    ('unknown', 'unknown', 'Unknown or unrecognized end reason');

-- Extraction methods
CREATE TABLE IF NOT EXISTS lk_extraction_methods (
    extraction_method   TEXT PRIMARY KEY,
    category            TEXT,
    description         TEXT
);

INSERT OR IGNORE INTO lk_extraction_methods (extraction_method, category, description) VALUES
    ('phenol_chloroform', 'manual', 'Phenol-chloroform extraction'),
    ('column_silica', 'column', 'Silica column-based extraction'),
    ('magnetic_beads', 'beads', 'Magnetic bead-based extraction'),
    ('alkaline_lysis', 'manual', 'Alkaline lysis for plasmids'),
    ('kit_qiagen', 'kit', 'Qiagen extraction kit'),
    ('kit_zymo', 'kit', 'Zymo Research extraction kit'),
    ('kit_monarch', 'kit', 'NEB Monarch extraction kit'),
    ('direct', 'none', 'Direct use without extraction');

-- Size selection methods
CREATE TABLE IF NOT EXISTS lk_size_selection (
    size_selection  TEXT PRIMARY KEY,
    method_type     TEXT,
    description     TEXT
);

INSERT OR IGNORE INTO lk_size_selection (size_selection, method_type, description) VALUES
    ('none', 'none', 'No size selection performed'),
    ('ampure_0.4x', 'beads', 'AMPure XP beads at 0.4x ratio (>1kb)'),
    ('ampure_0.5x', 'beads', 'AMPure XP beads at 0.5x ratio (>800bp)'),
    ('ampure_0.6x', 'beads', 'AMPure XP beads at 0.6x ratio (>500bp)'),
    ('ampure_0.8x', 'beads', 'AMPure XP beads at 0.8x ratio (>300bp)'),
    ('bluepippin', 'electrophoresis', 'BluePippin automated size selection'),
    ('gel_extraction', 'electrophoresis', 'Gel extraction and purification'),
    ('sspei', 'electrophoresis', 'SageELF or similar electrophoresis');

-- Library status
CREATE TABLE IF NOT EXISTS lk_library_status (
    status      TEXT PRIMARY KEY,
    description TEXT
);

INSERT OR IGNORE INTO lk_library_status (status, description) VALUES
    ('draft', 'Library specification in progress'),
    ('active', 'Library specification finalized and in use'),
    ('archived', 'Library specification archived, no longer active'),
    ('deprecated', 'Library specification deprecated, do not use');

-- Run purposes
CREATE TABLE IF NOT EXISTS lk_run_purposes (
    run_purpose TEXT PRIMARY KEY,
    description TEXT
);

INSERT OR IGNORE INTO lk_run_purposes (run_purpose, description) VALUES
    ('production', 'Production sequencing run'),
    ('test', 'Test or validation run'),
    ('optimization', 'Protocol optimization run'),
    ('troubleshooting', 'Troubleshooting run for issue diagnosis'),
    ('training', 'Training or demonstration run');
"""


def populate_lookups(conn: sqlite3.Connection) -> None:
    """Populate all lookup tables with predefined values.

    Args:
        conn: Active database connection

    Note:
        Uses INSERT OR IGNORE to avoid duplicates on re-population.
    """
    conn.executescript(LOOKUP_TABLES)
    conn.commit()


def get_valid_options(conn: sqlite3.Connection, table: str) -> List[Tuple[Any, ...]]:
    """Get all rows from a lookup table.

    Args:
        conn: Active database connection
        table: Name of lookup table (must start with 'lk_')

    Returns:
        List of tuples containing all rows from the table

    Raises:
        ValueError: If table name doesn't start with 'lk_'
    """
    if not table.startswith('lk_'):
        raise ValueError(f"Invalid lookup table name: {table}. Must start with 'lk_'")

    # Sanitize table name to prevent SQL injection
    # Only allow alphanumeric and underscore characters
    if not all(c.isalnum() or c == '_' for c in table):
        raise ValueError(f"Invalid characters in table name: {table}")

    cursor = conn.execute(f"SELECT * FROM {table}")
    return cursor.fetchall()
