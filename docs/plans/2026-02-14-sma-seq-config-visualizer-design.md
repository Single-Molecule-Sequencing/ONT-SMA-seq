# SMA-seq Experiment Configuration Visualizer

**Date:** 2026-02-14
**Status:** Approved

## Overview

Interactive browser-based tool for inspecting, creating, and managing SMA-seq experiment configurations. Provides full CRUD for sample sheets, barcode designs, adapter structure, construct layout, target sequences, and pipeline assumptions. Generates self-contained static HTML exports for sharing and archival.

## Architecture

**Approach:** FastAPI + HTMX + D3.js (Hybrid server-rendered with targeted interactivity)

- **Backend:** FastAPI (Python 3), Jinja2 templates, uvicorn
- **Frontend:** HTMX (CDN) for CRUD interactivity, D3.js (CDN) for construct diagrams, minimal CSS framework (CDN)
- **Data:** ONT-native config files as source of truth + minimal SMA-seq overlay TOML
- **Launch:** `python -m bin.viz.app` opens browser at `localhost:8050`

### File Structure

```
bin/
  viz/
    __init__.py
    app.py              # FastAPI app, startup, launch browser
    api.py              # REST endpoints for CRUD operations
    models.py           # Pydantic models for config validation
    config_store.py     # Read/write ONT config files from disk
    export.py           # Static HTML export generator
    templates/
      base.html         # Layout shell (nav, sidebar)
      dashboard.html    # Overview landing page
      sample_sheet.html # Sample sheet editor
      barcodes.html     # Barcode design viewer/editor
      construct.html    # Read construct diagram + truncation
      targets.html      # Target sequence manager
      assumptions.html  # Pipeline assumptions explainer
    static/
      construct.js      # D3 construct diagram rendering
      style.css         # Minimal custom styles
```

## Config Files (ONT-Native Formats)

The app reads and writes official ONT file formats wherever possible:

| File | ONT Format | Purpose |
|------|------------|---------|
| `sample_sheet.csv` | MinKNOW CSV | Barcode-to-alias mapping with experiment_id, kit, flow_cell_id, type |
| `arrangement.toml` | Dorado custom barcodes | Barcode flanks, scoring params, pattern definitions |
| `barcodes.fasta` | Dorado barcode sequences | Actual barcode sequences (FASTA, `BC%02i` headers) |
| `primers.fasta` | Dorado primer format | Custom primer sequences for trimming |
| `references/*.fasta` | Standard FASTA | Per-target reference sequences |

### SMA-seq Overlay (`sma_experiment.toml`)

Minimal metadata file for SMA-seq-specific context not covered by ONT formats:

```toml
[experiment]
description = "CYP2D6 variant library with duplexed native barcoding"

[construct]
adapter_5prime = "Ligation adapter"
adapter_3prime = "Ligation adapter"
insert_type = "amplicon"  # amplicon | genomic | synthetic
notes = "Upstream barcode is in forward orientation, downstream is RC"

[demultiplexing]
mode = "dual"  # "start_only" | "end_only" | "dual"

[demultiplexing.start_barcode]
assignments = { "nb05" = "library_A", "nb10" = "library_B" }

[demultiplexing.end_barcode]
assignments = { "nb10" = "forward", "nb05" = "reverse" }

[demultiplexing.pairs]
entries = [
    { start = "nb05", end = "nb10", alias = "CYP2D6_v04_fwd" },
    { start = "nb10", end = "nb05", alias = "CYP2D6_v04_rev" },
]

[classification]
barcode_search_window = 100
confidence_formula = "1.0 - (ed / barcode_length)"
ambiguity_triggers_full_construct = true

[quality]
q_bc = "Probability-averaged Phred: -10*log10(mean(10^(-Qi/10)))"
q_ld = "-10*log10(min(max(1/L^2, ed/L), 1.0))"

[truncation]
min_barcode_confidence = 0.5
min_target_fraction = 0.1
adapter_search_window = 50

[truncation.rules]
full            = "assign_to_target"
trunc_3prime    = "assign_to_target"
trunc_target    = "assign_by_start_barcode"
trunc_barcode   = "assign_by_start_barcode"
adapter_only    = "unclassified"
chimeric        = "flag_for_review"

[truncation.auto_references]
enabled = true

[assumptions]
entries = [
    { key = "all_reads_ingested", text = "All reads captured; filtering is post-hoc via SQL", why = "Avoids premature filtering bias" },
    { key = "untrimmed_bams", text = "BAMs must be untrimmed for barcode classification", why = "Barcodes are in the read sequence" },
    { key = "native_barcodes_only", text = "Only ONT native 24bp barcodes supported", why = "Hardcoded in barcodes.py" },
    { key = "duplexed_required", text = "Duplexed barcode pairs required", why = "Single barcodes can't uniquely identify fwd/rev targets" },
]
```

## Construct Model & Truncation Classification

### Full Library Product

```
5'-[seq_adapter]-[barcode1]-[target]-[RC_barcode2]-[RC_seq_adapter]-3'
```

### Truncation Classes

| Class | Structure | BC1 | BC2 | Assignment Rule |
|-------|-----------|:---:|:---:|-----------------|
| **Full** | adapter + bc1 + target + RC_bc2 + RC_adapter | yes | yes | assign_to_target |
| **Trunc-3prime** | adapter + bc1 + target + RC_bc2 | yes | yes | assign_to_target |
| **Trunc-target** | adapter + bc1 + target (partial/full) | yes | no | assign_by_start_barcode |
| **Trunc-barcode** | adapter + bc1 | yes | no | assign_by_start_barcode |
| **Adapter-only** | adapter | no | no | unclassified |
| **Chimeric** | unexpected barcode combinations | mismatch | mismatch | flag_for_review |

### Barcode Assignment Modes

Three independent strategies, selectable per experiment:

- **`start_only`**: Classify by upstream barcode only
- **`end_only`**: Classify by downstream barcode only
- **`dual`**: Require matching barcode pair (start, end) -> alias

Each barcode position supports different barcodes. When a barcode is undetected, truncation rules determine assignment.

### Auto-Generated Truncated References

When enabled, for each target alias the app generates:
- `{alias}_full.fasta`: adapter + bc1 + target + RC_bc2 + RC_adapter
- `{alias}_no_adapter.fasta`: bc1 + target + RC_bc2
- `{alias}_no_end_bc.fasta`: adapter + bc1 + target
- `{alias}_bc_only.fasta`: adapter + bc1

Output to: `references/auto/{alias}_*.fasta`

## UI Pages

### 1. Dashboard (`/`)
- Experiment summary card: exp_id, flow_cell, kit, description
- File status indicators: which ONT files are present/missing
- Truncation class pie chart (if DB exists)
- Quick actions: Create New, Import from MinKNOW, Export HTML

### 2. Sample Sheet (`/sample-sheet`)
- CRUD table with HTMX inline editing
- MinKNOW CSV format validation
- Import/export MinKNOW-compatible CSV
- Barcode ambiguity warnings
- Used vs available barcode indicators

### 3. Barcodes (`/barcodes`)
- Kit Barcodes tab (arrangement.toml / built-in ONT set)
- Custom Barcodes tab (user-defined)
- Sequence table: ID, sequence, RC, GC%, Tm
- Arrangement editor: flanks, scoring params
- FASTA import/export
- Pairwise edit distance heatmap (D3)

### 4. Construct (`/construct`)
- D3 annotated molecule diagram with clickable regions
- Truncation class ladder (6 classes, color-coded)
- Barcode assignment mode selector + tables
- Auto-reference generation toggle + preview
- Rule editor for truncation assignment

### 5. Targets (`/targets`)
- CRUD for reference FASTA files
- Per-target stats: length, GC%, linked barcode pair
- Upload or paste sequences
- Auto-generated truncated refs shown as children
- Header-to-alias validation

### 6. Assumptions (`/assumptions`)
- Editable key/description/rationale list
- Pre-populated defaults for new experiments
- Add/edit/remove entries

## API Endpoints

```
GET  /api/experiment              -> experiment metadata
PUT  /api/experiment              -> update metadata

GET  /api/sample-sheet            -> parsed sample sheet
PUT  /api/sample-sheet            -> update (rewrites CSV)
POST /api/sample-sheet/import     -> upload MinKNOW CSV
GET  /api/sample-sheet/export     -> download MinKNOW CSV

GET  /api/barcodes                -> barcode sequences + arrangement
PUT  /api/barcodes/arrangement    -> update arrangement.toml
PUT  /api/barcodes/sequences      -> update barcodes.fasta
GET  /api/barcodes/distances      -> pairwise edit distance matrix

GET  /api/construct               -> construct definition + truncation rules
PUT  /api/construct               -> update construct config
POST /api/construct/auto-refs     -> generate truncated reference FASTAs

GET  /api/targets                 -> all target references
POST /api/targets                 -> add new target FASTA
PUT  /api/targets/{alias}         -> update target
DELETE /api/targets/{alias}       -> remove target

GET  /api/assumptions             -> all assumptions
PUT  /api/assumptions             -> update assumptions

POST /api/export                  -> generate static HTML bundle
GET  /api/validate                -> validate all configs for consistency
```

## Static HTML Export

### Trigger Methods
1. UI button on dashboard
2. CLI: `python -m bin.viz.export <experiment_dir>/ [--output <dir>]`
3. Post-pipeline hook after ingest.py

### Output Structure

```
exports/{exp_id}_{timestamp}/
  index.html              # Overview with links to all sub-pages
  sample_sheet.html       # Barcode-to-alias table, validation results
  barcodes.html           # Sequences, arrangement, distance heatmap
  construct.html          # Molecule diagram, truncation classes, assignment rules
  targets.html            # Reference sequences, auto-generated truncated refs
  assumptions.html        # Pipeline logic documentation
  validation_report.html  # Cross-config consistency check results
```

### Export Characteristics
- Fully self-contained: inline CSS, JS (D3 bundled), data as JSON script blocks
- No external dependencies; viewable by double-clicking in any browser
- Print-friendly CSS media queries
- index.html links to all sub-pages

## Cross-Cutting Features

- **Live validation:** Every edit triggers `/api/validate` via HTMX, inline warnings
- **Undo:** Server keeps last 10 versions of each config file (timestamped copies)
- **MinKNOW import:** Auto-discover configs from MinKNOW output directory structure
- **Consistent nav:** Sidebar with page links, validation badge showing error count
