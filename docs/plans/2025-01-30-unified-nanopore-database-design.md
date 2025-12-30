# Unified Nanopore Database Design

**Date:** 2025-01-30
**Status:** Draft
**Author:** Athey Lab

## Overview

This design unifies three existing systems into a single coherent database architecture:

1. **nanopore_experiments.db** - Experiment discovery and tracking
2. **dorado-bench** - Basecalling benchmarking
3. **ONT-SMA-seq** - Per-read metrics and analysis

The unified system enables end-to-end tracking from experimental design through basecalling to per-read quality metrics, with flexible querying at all levels.

## Architecture

### Federated SQLite Approach

```
nanopore_unified.db              # Central: metadata, configs, references
├── experiments                  # Raw sequencing runs
├── basecall_runs               # Dorado configurations and outputs
├── library_specs               # Experimental designs
├── references                  # Shared reference sequences
└── lookup tables               # Standardized options
    │
    └── links to per-experiment databases:
        SMA_exp1.db             # Per-read data for experiment 1
        SMA_exp2.db             # Per-read data for experiment 2
        ...
```

**Rationale:**
- Parallelizes perfectly (pipeline jobs write to separate DBs)
- No lock contention between experiments
- Can archive/delete individual experiments
- Cross-experiment queries via SQLite ATTACH

## Integration Flow

```
┌─────────────────────────────────────────────────────────────────────────┐
│                     EXPERIMENTAL DESIGN                                 │
│  library_specs → lib_constructs → lib_fragments → lib_modifications    │
│  Define what you expect BEFORE sequencing                              │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                        DISCOVERY                                        │
│  Scanner populates experiments table from sequencing runs               │
│  Links experiments to library_specs via lib_experiment_link             │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                        BASECALLING                                      │
│  dorado-bench queries experiments for POD5 paths                        │
│  Generates commands, runs basecalling                                   │
│  INSERT INTO basecall_runs (status='complete', bam_path=...)            │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                        SMA-SEQ ANALYSIS                                 │
│  Query basecall_runs WHERE status='complete' AND sma_db_path IS NULL    │
│  Run ingest.py → produce SMA_{exp_id}.db                                │
│  UPDATE basecall_runs SET sma_db_path=..., status='analyzed'            │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                        QUERIES & ANALYSIS                               │
│  Cross-experiment via ATTACH or Python multi-DB loading                 │
│  Compare configs, evaluate against expectations, generate reports       │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Central Database Schema

### Core Tables

```sql
-- ============================================================
-- EXPERIMENTS (raw sequencing runs from MinKNOW)
-- ============================================================

CREATE TABLE experiments (
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

-- ============================================================
-- BASECALL RUNS (dorado executions)
-- ============================================================

CREATE TABLE basecall_runs (
    run_id          TEXT PRIMARY KEY,
    exp_id          TEXT REFERENCES experiments(exp_id),
    -- Model config
    model_tier      TEXT,              -- 'fast', 'hac', 'sup'
    model_version   TEXT,              -- '5.2.0'
    trim            INTEGER,           -- 0 or 1
    mod_bitflag     INTEGER,           -- 0, 1, 2, 4, etc.
    -- Dorado details
    dorado_version  TEXT,
    dorado_args     TEXT,              -- Full CLI args as JSON
    batch_size      INTEGER,
    emit_moves      INTEGER,
    -- Execution details
    gpu_model       TEXT,
    slurm_job_id    TEXT,
    runtime_seconds REAL,
    -- Outputs
    bam_path        TEXT,
    sma_db_path     TEXT,
    status          TEXT DEFAULT 'pending',
    created_at      TEXT DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================
-- REFERENCES (shared across experiments)
-- ============================================================

CREATE TABLE references (
    ref_id    TEXT PRIMARY KEY,
    ref_path  TEXT,
    ref_seq   TEXT,
    ref_len   INTEGER
);
```

### Experimental Design Tables

```sql
-- ============================================================
-- LIBRARY SPECIFICATIONS
-- ============================================================

CREATE TABLE library_specs (
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

-- ============================================================
-- CONSTRUCTS (source sequences)
-- ============================================================

CREATE TABLE lib_constructs (
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

-- ============================================================
-- FRAGMENTS (sequenceable portions)
-- ============================================================

CREATE TABLE lib_fragments (
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

-- ============================================================
-- MODIFICATIONS (expected methylation)
-- ============================================================

CREATE TABLE lib_modifications (
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

-- ============================================================
-- BARCODING
-- ============================================================

CREATE TABLE lib_barcodes (
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

-- ============================================================
-- SAMPLE PREP METADATA
-- ============================================================

CREATE TABLE lib_sample_prep (
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

-- ============================================================
-- QC THRESHOLDS
-- ============================================================

CREATE TABLE lib_qc_thresholds (
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

-- ============================================================
-- CUSTOM FIELDS (extensibility)
-- ============================================================

CREATE TABLE lib_custom_fields (
    id                  INTEGER PRIMARY KEY,
    lib_id              TEXT REFERENCES library_specs(lib_id),
    entity_type         TEXT,
    entity_id           TEXT,
    field_name          TEXT NOT NULL,
    field_value         TEXT,
    field_type          TEXT DEFAULT 'text',
    created_at          TEXT DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================
-- LIBRARY-EXPERIMENT LINK
-- ============================================================

CREATE TABLE lib_experiment_link (
    lib_id              TEXT REFERENCES library_specs(lib_id),
    exp_id              TEXT REFERENCES experiments(exp_id),
    relationship        TEXT DEFAULT 'primary',
    run_purpose         TEXT,
    notes               TEXT,
    PRIMARY KEY (lib_id, exp_id)
);
```

---

## Lookup Tables

### Pore Types

```sql
CREATE TABLE lk_pore_types (
    pore_type       TEXT PRIMARY KEY,
    description     TEXT,
    release_year    INTEGER,
    is_current      INTEGER DEFAULT 1
);

INSERT INTO lk_pore_types VALUES
    ('R9.4.1', 'Legacy pore, good accuracy', 2017, 0),
    ('R10.4.1', 'Current high-accuracy pore', 2022, 1),
    ('R10.4.2', 'Latest iteration', 2024, 1);
```

### Flow Cells

```sql
CREATE TABLE lk_flow_cells (
    flow_cell_type  TEXT PRIMARY KEY,
    pore_type       TEXT,
    device          TEXT,
    pore_count      INTEGER,
    description     TEXT
);

INSERT INTO lk_flow_cells VALUES
    ('FLO-FLG114', 'R10.4.1', 'Flongle', 126, 'Flongle R10.4.1'),
    ('FLO-MIN114', 'R10.4.1', 'MinION', 2048, 'MinION R10.4.1'),
    ('FLO-PRO114M', 'R10.4.1', 'PromethION', 2675, 'PromethION R10.4.1'),
    ('FLO-MIN106', 'R9.4.1', 'MinION', 2048, 'MinION R9.4.1 (legacy)');
```

### Sequencing Kits

```sql
CREATE TABLE lk_sequencing_kits (
    kit_id          TEXT PRIMARY KEY,
    kit_name        TEXT,
    kit_type        TEXT,
    barcoding       INTEGER DEFAULT 0,
    barcode_count   INTEGER,
    pore_compat     TEXT,
    description     TEXT
);

INSERT INTO lk_sequencing_kits VALUES
    ('SQK-LSK114', 'Ligation Sequencing Kit V14', 'ligation', 0, NULL, '["R10.4.1"]', 'Standard ligation prep'),
    ('SQK-NBD114.24', 'Native Barcoding Kit 24 V14', 'ligation', 1, 24, '["R10.4.1"]', '24-plex native barcoding'),
    ('SQK-NBD114.96', 'Native Barcoding Kit 96 V14', 'ligation', 1, 96, '["R10.4.1"]', '96-plex native barcoding'),
    ('SQK-RBK114.24', 'Rapid Barcoding Kit 24 V14', 'rapid', 1, 24, '["R10.4.1"]', '24-plex rapid barcoding'),
    ('SQK-RBK114.96', 'Rapid Barcoding Kit 96 V14', 'rapid', 1, 96, '["R10.4.1"]', '96-plex rapid barcoding'),
    ('SQK-RAD114', 'Rapid Sequencing Kit V14', 'rapid', 0, NULL, '["R10.4.1"]', 'Rapid transposase prep'),
    ('SQK-LSK109', 'Ligation Sequencing Kit', 'ligation', 0, NULL, '["R9.4.1"]', 'Legacy R9 ligation');
```

### Basecalling Models

```sql
CREATE TABLE lk_basecall_models (
    model_id        TEXT PRIMARY KEY,
    pore_type       TEXT,
    speed           TEXT,
    tier            TEXT,
    version         TEXT,
    dorado_min_ver  TEXT,
    description     TEXT
);

INSERT INTO lk_basecall_models VALUES
    ('dna_r10.4.1_e8.2_400bps_fast@v5.0.0', 'R10.4.1', '400bps', 'fast', 'v5.0.0', '0.5.0', 'Fast model'),
    ('dna_r10.4.1_e8.2_400bps_hac@v5.0.0', 'R10.4.1', '400bps', 'hac', 'v5.0.0', '0.5.0', 'High accuracy'),
    ('dna_r10.4.1_e8.2_400bps_sup@v5.0.0', 'R10.4.1', '400bps', 'sup', 'v5.0.0', '0.5.0', 'Super accuracy'),
    ('dna_r10.4.1_e8.2_400bps_sup@v5.2.0', 'R10.4.1', '400bps', 'sup', 'v5.2.0', '0.8.0', 'Latest super accuracy');
```

### Modifications

```sql
CREATE TABLE lk_modifications (
    mod_code        INTEGER PRIMARY KEY,
    mod_type        TEXT,
    mod_name        TEXT,
    dorado_flag     TEXT,
    description     TEXT
);

INSERT INTO lk_modifications VALUES
    (0, 'none', 'No modifications', NULL, 'Canonical bases only'),
    (1, '6mA', '6-methyladenine', '6mA', 'Bacterial/mitochondrial'),
    (2, '5mCG_5hmCG', '5mC and 5hmC in CpG', '5mCG_5hmCG', 'Mammalian CpG methylation'),
    (4, '5mC_5hmC', '5mC and 5hmC all contexts', '5mC_5hmC', 'All-context cytosine mods'),
    (8, '4mC_5mC', '4mC and 5mC', '4mC_5mC', 'Bacterial cytosine mods'),
    (16, '5mC', '5-methylcytosine only', '5mC', '5mC without hydroxymethyl'),
    (3, '6mA+5mCG_5hmCG', 'Combined', '6mA,5mCG_5hmCG', 'Adenine + CpG cytosine'),
    (5, '6mA+5mC_5hmC', 'Combined', '6mA,5mC_5hmC', 'Adenine + all-context cytosine');
```

### Construct Types

```sql
CREATE TABLE lk_construct_types (
    construct_type  TEXT PRIMARY KEY,
    description     TEXT
);

INSERT INTO lk_construct_types VALUES
    ('plasmid', 'Circular plasmid DNA'),
    ('genome', 'Genomic DNA (bacterial, viral, etc.)'),
    ('amplicon', 'PCR amplification product'),
    ('synthetic', 'Synthesized DNA construct'),
    ('cdna', 'Complementary DNA from RNA'),
    ('chromosome', 'Eukaryotic chromosome/fragment');
```

### Fragment Types

```sql
CREATE TABLE lk_fragment_types (
    fragment_type   TEXT PRIMARY KEY,
    description     TEXT
);

INSERT INTO lk_fragment_types VALUES
    ('whole_plasmid', 'Uncut circular plasmid (rolling circle)'),
    ('linearized', 'Single-cut linearized plasmid'),
    ('restriction_digest', 'Multi-enzyme restriction digest'),
    ('pcr_amplicon', 'PCR amplification product'),
    ('sheared', 'Mechanically sheared DNA'),
    ('enzymatic_fragment', 'Enzymatic fragmentation'),
    ('native', 'Native/unprocessed');
```

### End Reasons

```sql
CREATE TABLE lk_end_reasons (
    end_reason      TEXT PRIMARY KEY,
    category        TEXT,
    is_good         INTEGER,
    description     TEXT
);

INSERT INTO lk_end_reasons VALUES
    ('signal_positive', 'complete', 1, 'Read completed normally (good)'),
    ('signal_negative', 'complete', 1, 'Read completed, no signal'),
    ('unblock_mux_change', 'rejected', 0, 'Adaptive sampling rejection'),
    ('mux_change', 'technical', 0, 'Mux scan during read'),
    ('data_service_unblock_mux_change', 'rejected', 0, 'Software-triggered unblock'),
    ('unknown', 'unknown', 0, 'Unknown termination');
```

### Extraction Methods

```sql
CREATE TABLE lk_extraction_methods (
    method          TEXT PRIMARY KEY,
    description     TEXT
);

INSERT INTO lk_extraction_methods VALUES
    ('phenol_chloroform', 'Traditional organic extraction'),
    ('column_silica', 'Silica column purification'),
    ('magnetic_beads', 'Magnetic bead-based extraction'),
    ('alkaline_lysis', 'Alkaline lysis (plasmid miniprep)'),
    ('kit_qiagen', 'Qiagen kit extraction'),
    ('kit_zymo', 'Zymo Research kit extraction'),
    ('kit_monarch', 'NEB Monarch kit extraction'),
    ('direct', 'Direct input (no extraction)');
```

### Size Selection Methods

```sql
CREATE TABLE lk_size_selection (
    method          TEXT PRIMARY KEY,
    description     TEXT
);

INSERT INTO lk_size_selection VALUES
    ('none', 'No size selection'),
    ('ampure_0.4x', 'AMPure XP 0.4x ratio (>1kb)'),
    ('ampure_0.5x', 'AMPure XP 0.5x ratio (>500bp)'),
    ('ampure_0.6x', 'AMPure XP 0.6x ratio (>300bp)'),
    ('ampure_0.8x', 'AMPure XP 0.8x ratio (cleanup)'),
    ('bluepippin', 'BluePippin gel extraction'),
    ('gel_extraction', 'Manual gel extraction'),
    ('sspei', 'Short fragment eliminator');
```

### Library Statuses

```sql
CREATE TABLE lk_library_status (
    status          TEXT PRIMARY KEY,
    description     TEXT
);

INSERT INTO lk_library_status VALUES
    ('draft', 'In design, not finalized'),
    ('active', 'Currently in use'),
    ('archived', 'No longer active'),
    ('deprecated', 'Replaced by newer version');
```

### Run Purposes

```sql
CREATE TABLE lk_run_purposes (
    purpose         TEXT PRIMARY KEY,
    description     TEXT
);

INSERT INTO lk_run_purposes VALUES
    ('production', 'Production data collection'),
    ('test', 'Test/validation run'),
    ('optimization', 'Protocol optimization'),
    ('troubleshooting', 'Debugging issues'),
    ('training', 'Training/demonstration');
```

---

## Per-Experiment Database Schema

Each experiment gets its own `SMA_{exp_id}.db` file:

```sql
-- ============================================================
-- READS (per-read metrics)
-- ============================================================

CREATE TABLE reads (
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

CREATE INDEX idx_reads_run_id ON reads(run_id);
CREATE INDEX idx_reads_ref_id ON reads(ref_id);
CREATE INDEX idx_reads_end_reason ON reads(end_reason);
CREATE INDEX idx_reads_ed ON reads(ed) WHERE ed IS NOT NULL;
CREATE INDEX idx_reads_qbc ON reads(q_bc);

-- ============================================================
-- RUN SUMMARIES (pre-computed aggregates)
-- ============================================================

CREATE TABLE run_summaries (
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
```

---

## Example Queries

### 1. Find experiments ready for SMA-seq analysis

```sql
SELECT
    e.exp_id,
    e.sample_id,
    br.model_tier,
    br.model_version,
    br.bam_path,
    ls.lib_name
FROM basecall_runs br
JOIN experiments e ON br.exp_id = e.exp_id
LEFT JOIN lib_experiment_link lel ON e.exp_id = lel.exp_id
LEFT JOIN library_specs ls ON lel.lib_id = ls.lib_id
WHERE br.status = 'complete'
  AND br.sma_db_path IS NULL;
```

### 2. Compare basecaller performance

```sql
SELECT
    br.model_tier,
    br.model_version,
    AVG(r.ed) as mean_ed,
    AVG(r.q_bc) as mean_qbc,
    COUNT(*) as read_count
FROM reads r
JOIN basecall_runs br ON r.run_id = br.run_id
JOIN experiments e ON br.exp_id = e.exp_id
JOIN lib_experiment_link lel ON e.exp_id = lel.exp_id
WHERE lel.lib_id = 'PVU2_CYP_panel_v1'
  AND r.ref_id IS NOT NULL
GROUP BY br.model_tier, br.model_version
ORDER BY mean_ed ASC;
```

### 3. QC check against expectations

```sql
SELECT
    e.exp_id,
    ls.lib_name,
    (SELECT COUNT(*) FROM reads r
     JOIN basecall_runs br ON r.run_id = br.run_id
     WHERE br.exp_id = e.exp_id AND r.ref_id IS NOT NULL) * 100.0 /
    (SELECT COUNT(*) FROM reads r
     JOIN basecall_runs br ON r.run_id = br.run_id
     WHERE br.exp_id = e.exp_id) as actual_purity,
    qt.min_purity * 100 as expected_min_purity,
    qt.target_purity * 100 as target_purity,
    CASE
        WHEN actual_purity >= qt.target_purity * 100 THEN 'PASS'
        WHEN actual_purity >= qt.min_purity * 100 THEN 'MARGINAL'
        ELSE 'FAIL'
    END as purity_status
FROM experiments e
JOIN lib_experiment_link lel ON e.exp_id = lel.exp_id
JOIN library_specs ls ON lel.lib_id = ls.lib_id
JOIN lib_qc_thresholds qt ON ls.lib_id = qt.lib_id;
```

### 4. Fragment analysis

```sql
SELECT
    lf.name as fragment_name,
    lf.expected_len,
    lf.expected_fraction * 100 as expected_pct,
    COUNT(r.read_id) as actual_reads,
    COUNT(r.read_id) * 100.0 / SUM(COUNT(r.read_id)) OVER() as actual_pct,
    AVG(r.ed) as mean_ed,
    AVG(r.q_bc) as mean_qbc
FROM lib_fragments lf
JOIN lib_constructs lc ON lf.construct_id = lc.construct_id
JOIN library_specs ls ON lc.lib_id = ls.lib_id
JOIN lib_experiment_link lel ON ls.lib_id = lel.lib_id
JOIN basecall_runs br ON lel.exp_id = br.exp_id
LEFT JOIN reads r ON br.run_id = r.run_id
    AND r.ref_id = lf.fragment_id
WHERE lel.exp_id = '30_40_PVU2'
  AND br.model_tier = 'sup'
GROUP BY lf.fragment_id, lf.name, lf.expected_len, lf.expected_fraction
ORDER BY actual_reads DESC;
```

### 5. End reason distribution

```sql
SELECT
    e.exp_id,
    r.end_reason,
    ler.category,
    ler.is_good,
    COUNT(*) as count,
    COUNT(*) * 100.0 / SUM(COUNT(*)) OVER(PARTITION BY e.exp_id) as pct
FROM reads r
JOIN basecall_runs br ON r.run_id = br.run_id
JOIN experiments e ON br.exp_id = e.exp_id
LEFT JOIN lk_end_reasons ler ON r.end_reason = ler.end_reason
WHERE br.model_tier = 'sup'
GROUP BY e.exp_id, r.end_reason, ler.category, ler.is_good;
```

### 6. Best config per library type

```sql
SELECT
    lc.construct_type,
    br.model_tier,
    br.model_version,
    COUNT(DISTINCT e.exp_id) as experiments,
    AVG(rs.mean_qbc) as avg_qbc,
    AVG(rs.mean_ed) as avg_ed,
    AVG(rs.matched_reads * 100.0 / rs.total_reads) as avg_purity
FROM run_summaries rs
JOIN basecall_runs br ON rs.run_id = br.run_id
JOIN experiments e ON br.exp_id = e.exp_id
JOIN lib_experiment_link lel ON e.exp_id = lel.exp_id
JOIN library_specs ls ON lel.lib_id = ls.lib_id
JOIN lib_constructs lc ON ls.lib_id = lc.lib_id
GROUP BY lc.construct_type, br.model_tier, br.model_version
ORDER BY lc.construct_type, avg_ed ASC;
```

### 7. Orphan reads analysis

```sql
SELECT
    e.exp_id,
    r.end_reason,
    AVG(r.readlen) as mean_len,
    AVG(r.q_bc) as mean_qbc,
    COUNT(*) as orphan_count
FROM reads r
JOIN basecall_runs br ON r.run_id = br.run_id
JOIN experiments e ON br.exp_id = e.exp_id
WHERE r.ref_id IS NULL
  AND br.model_tier = 'sup'
GROUP BY e.exp_id, r.end_reason
ORDER BY orphan_count DESC;
```

---

## Implementation Plan

### Phase 1: Central Database Setup
1. Create `nanopore_unified.db` with all tables
2. Populate lookup tables with predefined options
3. Migrate existing data from `nanopore_experiments.db`

### Phase 2: Library Spec Tooling
1. CLI/GUI for defining library specifications
2. Import from SnapGene/GenBank files
3. Barcode scheme configuration

### Phase 3: Integration Scripts
1. `sma_orchestrator.py` - Watches for new experiments, triggers workflows
2. Modify `dorado-bench` to read/write central DB
3. Modify `ONT-SMA-seq/ingest.py` to update central DB status

### Phase 4: Reporting & Visualization
1. QC dashboard comparing actual vs expected
2. Basecaller benchmarking reports
3. Cross-experiment analysis tools

---

## File Locations

```
/nfs/turbo/umms-atheylab/
├── nanopore_unified.db           # Central metadata database
├── sma_databases/                # Per-experiment SMA databases
│   ├── SMA_30_40_PVU2.db
│   ├── SMA_30_40_SRF1.db
│   └── ...
├── ONT-SMA-seq/                  # Pipeline code
│   ├── bin/
│   ├── lib/                      # New: shared library code
│   │   ├── db_central.py         # Central DB operations
│   │   ├── db_sma.py             # SMA DB operations
│   │   └── orchestrator.py       # Workflow orchestration
│   └── ...
└── dorado-bench/                 # Basecalling benchmarking
    └── ...
```
