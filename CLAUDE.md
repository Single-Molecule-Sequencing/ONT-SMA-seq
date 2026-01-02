# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ONT-SMA-seq is a pure Python pipeline for processing Oxford Nanopore single-molecule sequencing data with SQLite database storage. It transforms uBAM files and POD5 metadata into a structured database with per-read quality metrics, end reason tagging, and length-based reference matching.

## Architecture

Four-stage pipeline with separate scripts for each concern:

```
mkdb.py → inputInit.py → extractMeta.py → ingest.py
   ↓           ↓              ↓              ↓
 Create     Symlink &      POD5 view      Tag BAM,
 SQLite     validate       CLI for        calculate
 schema     inputs         end reasons    metrics, DB insert
```

**Data flow:**
- Input: uBAM file (strict naming convention), POD5 directory, reference FASTA (2 sequences)
- Output: Tagged BAM with ER tags, SQLite database with Reads/Mods/Refseq/Exp tables

## Commands

### Environment Setup
```bash
conda env create -f env/env.yml
conda activate ont-sma-seq
```

### Full Pipeline Execution
```bash
# 1. Initialize database
python bin/mkdb.py -e <exp_id> -o Output

# 2. Standardize inputs (creates Input/ directory with symlinks)
python bin/inputInit.py \
  -b /path/to/<exp_id>_<tier>_v<ver>_trim<0|1>_<mods>.bam \
  -p /path/to/pod5_directory \
  -r /path/to/reference.fa

# 3. Extract end reason metadata from POD5
python bin/extractMeta.py \
  -i Input/<exp_id>_pod5 \
  -o Input/<exp_id>_summary.tsv

# 4. Process BAM, calculate metrics, populate database
python bin/ingest.py \
  -e <exp_id> \
  -b Input/<exp_id>_<tier>_v<ver>_trim<0|1>_<mods>.bam \
  -s Input/<exp_id>_summary.tsv \
  -r Input/<exp_id>.fa \
  -d Output/SMA_<exp_id>.db \
  -o Output/<exp_id>_tagged.bam \
  -k 150  # length tolerance in bp
```

### Script Help
```bash
python bin/<script>.py --help
```

## BAM Naming Convention

**Required format:** `{exp_id}_{tier}_v{version}_{trim}_{mods}.bam`

| Field | Values |
|-------|--------|
| tier | `fast`, `hac`, `sup` |
| version | e.g., `v5.2.0` |
| trim | `trim0` or `trim1` |
| mods | Bitflag integer (see below) |

**Modification Bitflags:**
- 0: none, 1: 6mA, 2: 5mCG_5hmCG, 4: 5mC_5hmC, 8: 4mC_5mC, 16: 5mC
- Combine via addition: 6mA + 5mC_5hmC = 1 + 4 = 5

**Example:** `sample1_hac_v5.2.0_trim1_5.bam`

## Database Schema

**Reads table** - primary output, one row per read:
- `uniq_id` (PK): `{exp_id}{tier}{ver}t{trim}m{mods}_{read_id}`
- `refseq_id`: Matched reference (NULL if length out of range)
- `ed`: Levenshtein edit distance (NULL if no match)
- `q_bc`: Probability-averaged basecall quality
- `q_ld`: Levenshtein-derived quality (NULL if no match)
- `ER`: End reason from POD5

**Mods table** - pre-populated bitflag lookups

**Refseq table** - reference sequences with length ranges

## Key Implementation Details

**Reference matching:** Length-based with configurable tolerance (`-k`). Read must fall within `ref_length +/- tolerance` to be assigned.

**Quality formulas:**
- `q_bc = -10 * log10(mean(10^(-Q/10)))` for all base Q-scores
- `q_ld = -10 * log10(min(max(1/L^2, ed/L), 1))` where L = reference length

**Performance optimizations:**
- `PRAGMA synchronous=OFF` and `journal_mode=MEMORY` for bulk inserts
- Pandas Series->dict for O(1) metadata lookup
- pysam streaming (never loads full BAM)
- `check_sq=False` for unmapped BAMs

## Dependencies

Core: `pysam`, `edlib`, `pandas`, `numpy`, `pod5` (CLI and library)

Visualization (future): `matplotlib`, `seaborn`

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
