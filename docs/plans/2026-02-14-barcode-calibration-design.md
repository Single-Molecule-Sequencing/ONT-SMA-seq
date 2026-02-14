# SMA-seq Barcode Calibration & Threshold Visualization Tool

**Date:** 2026-02-14
**Status:** Approved

## Overview

Data-driven barcode calibration system for SMA-seq experiments. Replaces hardcoded classification thresholds with empirically derived values based on observed barcode edit distance and confidence distributions from ground-truth experiments. Provides interactive visualization for threshold selection and cross-experiment comparison.

## Architecture

**Approach:** CLI pipeline + separate browser-based visualization app (Approach B)

Two distinct tools sharing a common SQLite database schema:

- **`sma_calibrate.py` CLI**: Auto-discovers MinKNOW output directories, validates basecalling consistency, merges BAMs, ingests reads with SMA-seq demultiplexing, extracts signal duration from sequencing summaries.
- **`sma_calibrate_viz` browser app**: FastAPI + HTMX + D3.js interactive tool for exploring barcode distributions, dragging thresholds, viewing confusion matrices, and comparing experiments. Exports self-contained static HTML and publication SVG figures.

```
sma_calibrate.py (CLI)                    sma_calibrate_viz (Browser)
+------------------------------+          +---------------------------------+
| discover:                    |          | Distributions:                  |
|   Walk MinKNOW output dirs   |          |   KDE: signal_dur, read_len     |
|   Parse final_summary_*.txt  |          |   Per-experiment / per-barcode   |
|   Parse sequencing_summary   |          |   Per-target / per-end_reason    |
|   Parse sample_sheet_*.csv   |          |   Overlay multiple experiments   |
|                              |          |                                 |
| validate:                    |          | Thresholds:                     |
|   Basecalling model match    |    DB    |   Draggable on distributions     |
|   No trim / no demux         |-------->|   Live confusion matrix update   |
|   Same flow cell + config    |          |   Classification counts table    |
|                              |          |   Affected reads list            |
| merge:                       |          |                                 |
|   Concatenate per-barcode    |          | Comparison:                     |
|   and per-run BAMs           |          |   Multi-DB overlay mode          |
|                              |          |   Side-by-side panel mode        |
| ingest:                      |          |                                 |
|   SMA-seq demultiplexing     |          | Export:                         |
|   Signal duration from       |          |   Static HTML (self-contained)   |
|   sequencing_summary         |          |   Publication KDE/SVG figures    |
|   End reason from summary    |          |                                 |
+------------------------------+          +---------------------------------+
```

## Ground Truth Experiments

Four experiments with increasing construct complexity, all with Sanger-verified target DNA, TapeStation-verified product sizes, and unique barcode-per-strand-per-target assignments:

| Experiment | Adapters | Products | End Barcode? | Path |
|---|---|---|---|---|
| `no_trim` | All missing 5'P | Truncated only (4 molecules) | No | `12282025_IF_DoubleBC_SMA_seq_no_trim` |
| `one_nick` | One per target missing 5'P | Rev: truncated; Fwd: up to full-length | Rev: no; Fwd: yes | `12302025_IF_DoubleBC_SMA_seq_one_nick` |
| `Ggallphos_extended` | All phosphorylated | Full-length (multiple short targets) | Yes | `20260129_IF_GG_Part5_Odd/.../extended` |
| `Ggallphos_extended_1` | All phosphorylated | Full-length (single long target) | Yes | `20260129_IF_GG_Part5_Odd/.../extended_1` |

## MinKNOW Output Structure (ONT Specification)

Directory pattern: `{protocol_group_id}/{sample_id}/{start_time}_{device_id}_{flow_cell_id}_{short_protocol_run_id}/`

Key files per run:
- `final_summary_{flow_cell_id}_{short_protocol_run_id}_{short_run_id}.txt` — run metadata
- `sequencing_summary_{flow_cell_id}_{short_protocol_run_id}_{short_run_id}.txt` — per-read TSV with `duration`, `end_reason`, `barcode_arrangement`, `barcode_score`, `barcode_front_score`, `barcode_rear_score`, `sequence_length_template`, `mean_qscore_template`
- `sample_sheet_{flow_cell_id}_{daq_start_time}_{short_protocol_run_id}.csv` — barcode-to-alias mappings
- `bam_pass/` — BAM files, optionally in per-alias subdirectories
- `pod5_pass/` — raw signal data (used only if sequencing summary unavailable)

Signal duration comes from the sequencing summary `duration` column (seconds), not POD5 extraction.

## CLI Pipeline

### Discovery

Walks root directory recursively, identifies MinKNOW output directories by presence of `final_summary_*.txt`. For each run:
1. Parse `final_summary_*.txt` for flow cell ID, device, protocol run ID, sample ID
2. Parse `sequencing_summary_*.txt` for per-read metadata
3. Parse `sample_sheet_*.csv` for barcode-to-alias mappings
4. Locate BAM files in `bam_pass/` (flat or per-alias subdirectories)

### Merge Logic

Group runs by `(flow_cell_id, sample_id, kit)`. Within a group, runs sorted by start time. Runs with same flow cell and sample within a configurable window (default 24h) are merged. Rules:
- **Merge**: Multiple sequencing runs of same sample on same flow cell (restarts with same config)
- **Don't merge**: Post-wash runs, different samples, different flow cells
- **Validate**: All BAMs in a merge group must share basecalling model, no trimming, no prior demultiplexing

### Validation

Before merging, validate all runs in a group share:
- Same basecalling model (from BAM `@PG` header)
- No trimming applied (check `@PG` for `--trim` flags)
- No demultiplexing by Dorado (warn if `barcode_arrangement` populated)

Failures produce clear error messages identifying which run disagrees.

### Merge

1. Concatenate all BAM files across per-barcode subdirectories and per-run directories into `merged.bam`
2. Record provenance: source BAM paths, checksums, run IDs
3. Store in `RunMetadata` table

### Ingest

Wraps existing `ingest.py` logic, enhanced:
- Reads signal duration from sequencing summary (joined by `read_id`)
- Stores `signal_duration_s` in Reads table
- Performs SMA-seq barcode classification (ignoring MinKNOW assignments)
- Stores all barcode metrics

### Config Generation

For experiments without SMA-seq configs:
- Generate `construct.toml` from user-provided barcode-target mappings
- Generate sample sheet CSV
- Generate reference FASTA files from provided target sequences

### CLI Interface

```
sma_calibrate discover <root_dir>
sma_calibrate validate <root_dir>
sma_calibrate ingest <root_dir> -c <construct.toml>
sma_calibrate viz <db_path> [<db_path>...]
sma_calibrate export <db_path> [--output <dir>]
sma_calibrate export <db1> <db2> [--compare] [--output <dir>]
```

## Database Schema

### Extended Reads Table

Two new columns added to existing schema:

```sql
signal_duration_s REAL    -- from sequencing_summary duration column (seconds)
mean_qscore REAL          -- from sequencing_summary mean_qscore_template
```

### New RunMetadata Table

```sql
CREATE TABLE RunMetadata (
    run_id TEXT PRIMARY KEY,
    flow_cell_id TEXT,
    device_id TEXT,
    sample_id TEXT,
    experiment_id TEXT,
    kit TEXT,
    protocol_run_id TEXT,
    start_time TEXT,
    basecall_model TEXT,
    source_bam_count INTEGER,
    source_bam_paths TEXT,       -- JSON array
    merge_timestamp TEXT
);
```

### New ReadRun Table

```sql
CREATE TABLE ReadRun (
    read_id TEXT PRIMARY KEY,
    run_id TEXT,
    FOREIGN KEY (run_id) REFERENCES RunMetadata(run_id)
);
```

## Visualization App

### Page 1: Experiment Overview (`/`)

Summary cards per loaded database: experiment ID, flow cell, kit, basecalling model, total reads, pass/fail counts, mean Q-score, barcode pair count, target count, run count, merge provenance.

### Page 2: Length Distributions (`/distributions`)

Two paired KDE plots (publication-quality, Savitzky-Golay peak detection):
- **Signal duration** (seconds) from `signal_duration_s`
- **Read length** (bases) from `readlen`

Grouping modes (tabs):
- All reads (single distribution)
- Per-barcode (one curve per `bc_start_id` or `bc_end_id`)
- Per-target (one curve per `tgt_id`)
- Per-end_reason (one curve per `ER` value)

Each plot includes:
- Draggable vertical threshold lines (D3 drag)
- Legend with read counts per group
- Expected product sizes as reference markers
- End reason subplots below each main KDE

Cross-experiment: overlay distributions from multiple databases on same axes, or side-by-side panel mode.

### Page 3: Barcode Confidence (`/confidence`)

Distribution plots:
- KDE of `bc_start_conf` across all reads, colored by true/false assignment
- KDE of `bc_end_conf`, same treatment
- Scatter: `bc_start_conf` vs `bc_end_conf`, colored by `trunc_level`

Interactive thresholds:
- Draggable lines for `start_barcode_min` and `full_length_threshold`
- Live-updating panels:
  1. Classification counts table (reads per trunc_level)
  2. Confusion matrix (true barcode vs. classified barcode)
  3. Affected reads list (scrollable table of reads that change classification)

Per-barcode detail: click barcode ID to see its edit distance distribution against all other barcodes, showing separation between true-match peak and next-best-match peak.

### Page 4: Barcode Separation (`/separation`)

Pairwise edit distance heatmap (D3): all experiment barcodes, color-coded by minimum edit distance, clickable cells for underlying distributions.

Separation metrics table: per-barcode mean true-match ED, mean next-best ED, separation gap, estimated error rate at current threshold. Sortable by separation gap.

Barcode set analysis: expected false-positive rate per barcode, cross-assignment rate between pairs, overall demultiplexing accuracy estimate.

### Page 5: Threshold Optimization (`/thresholds`)

ROC-style curves:
- `bc_start_conf` threshold: TPR vs FPR across values
- `bc_end_conf` threshold: same
- `full_length_threshold`: correct full-length classification vs truncated misclassification

Threshold recommendation: given user-selectable target false-positive rate, compute optimal thresholds. Display recommended values with table showing threshold, sensitivity, specificity, F1, reads affected.

Export: write selected thresholds back to `construct.toml` confidence section.

### Page 6: Cross-Experiment Comparison (`/compare`)

Overlay mode: select 2-4 databases, overlay distributions on shared axes (length, confidence, end reason).

Panel mode: same experiments in grid layout for side-by-side inspection.

Comparison table: per-experiment summary metrics in columns (total reads, mean length, mean signal duration, demux accuracy, mean barcode confidence, truncation proportions, end reason proportions).

## Static HTML Export

### Output Structure

```
exports/{experiment_id}_{timestamp}/
  index.html              # Overview + navigation
  distributions.html      # Signal + read length KDEs with end reason subplots
  confidence.html         # Barcode confidence distributions + confusion matrix
  separation.html         # Pairwise edit distance heatmap + separation table
  thresholds.html         # ROC curves + threshold recommendations
  figures/
    signal_length_kde.svg
    read_length_kde.svg
    confidence_scatter.svg
    confusion_matrix.svg
    separation_heatmap.svg
```

All HTML files self-contained with inline CSS, JS (D3 bundled), and data as JSON script blocks. No external dependencies. Print-friendly CSS media queries.

### Triggers

- CLI: `sma_calibrate export <db_path> [--output <dir>]`
- Browser: `/api/export` endpoint
- Cross-experiment: `sma_calibrate export <db1> <db2> --compare`

## File Structure

```
bin/
  calibrate/
    __init__.py
    cli.py              # CLI entry point (argparse subcommands)
    discover.py         # MinKNOW output directory auto-detection
    merge.py            # BAM merge + basecalling validation
    signal.py           # Signal duration from sequencing_summary or POD5
    config_gen.py       # Generate construct TOML + sample sheet
  calibrate_viz/
    __init__.py
    app.py              # FastAPI app, startup, launch browser
    api.py              # REST endpoints for distributions, thresholds, comparison
    distributions.py    # KDE computation, peak detection, threshold analysis
    confusion.py        # Barcode confusion matrix computation
    comparison.py       # Cross-experiment overlay logic
    export.py           # Static HTML + SVG export generator
    templates/
      base.html
      overview.html
      distributions.html
      confidence.html
      separation.html
      thresholds.html
      compare.html
    static/
      distributions.js  # D3 KDE plots with draggable thresholds
      confidence.js     # D3 scatter + confidence distributions
      separation.js     # D3 heatmap
      thresholds.js     # D3 ROC curves
      style.css
```
