# ONT-SMA-seq

[![DOI](https://zenodo.org/badge/1110728589.svg)](https://doi.org/10.5281/zenodo.18872468)

The Single-Molecule-Accuracy-seq protocol for Oxford Nanopore Technology experiments, implemented as a Python package with a SQLite database backend. This workflow processes an unaligned BAM file and its parent Pod5 files for a single target sequence, storing per-read metrics and metadata for downstream analysis.

## Installation

### 1. Create and activate the conda environment

```bash
conda env create -f env/env.yml
conda activate ont-sma-seq
```

**Note:** `edlib` may fail to import despite proper conda installation. If encountered, reinstall via pip:

```bash
pip install edlib --force-reinstall --no-cache-dir
```

### 2. Install the package

```bash
pip install -e .
```

This exposes the `ont-sma` command. To also install downstream analysis dependencies (`matplotlib`, `seaborn`):

```bash
pip install -e ".[analysis]"
```

---

## Manifest

All modules are located in `src/ont_sma_seq/`.

* `mkdb.py`: Initializes the SQLite database schema, parses the Experiment ID, and populates static lookup tables (`Mods`, `Exp`).
* `init.py`: Validates the reference FASTA (exactly one sequence) and inserts it into the `Target` table.
* `meta.py`: Wraps `pod5 view` to extract `read_id` and `end_reason` into a TSV.
* `ingest.py`: Streams BAM reads, computes per-read quality and Levenshtein metrics, and populates the `Reads` table in batches.
* `merge.py`: Merges one or more per-run databases into a master database using `INSERT OR IGNORE`.
* `cli.py`: Argparse entry point exposing all subcommands and the `run` pipeline orchestrator.

## Usage

### Push-button run

Edit `config.yml` with your experiment paths, then run:

```bash
ont-sma run [-c config.yml]
```

This chains all steps sequentially: `mkdb` → `init` → `meta` → `ingest`. If any step fails, the pipeline stops and reports which step failed.

### Individual subcommands

Each step can also be run independently:

```bash
ont-sma mkdb    -e <exp_id> [-o <outdir>]
ont-sma init    -d <db> -r <ref>
ont-sma meta    -i <pod5_dir> [-o <summary.tsv>]
ont-sma ingest  -b <bam> -d <db> -m <summary.tsv>
ont-sma merge   -o <master.db> <input_db> [<input_db> ...]
```

#### `mkdb`

Initializes the SQLite database, parses the Experiment ID, and populates static tables (`Mods`, `Exp`).

| Flag             | Description                                                        |
| ---------------- | ------------------------------------------------------------------ |
| `-e`, `--expid`  | Experiment ID in `FlowCellID_SampleID_Alias` format (**required**) |
| `-o`, `--outdir` | Output directory (default: `Output`)                               |

```bash
ont-sma mkdb -e FAL12345_20260129-IF_SMA -o Output
```

#### `init`

Validates the reference FASTA (must contain exactly one sequence) and inserts it into the `Target` table.

| Flag          | Description                                |
| ------------- | ------------------------------------------ |
| `-d`, `--db`  | Path to the SQLite database (**required**) |
| `-r`, `--ref` | Path to reference FASTA (**required**)     |

```bash
ont-sma init -d Output/SMA_FAL12345_20260129-IF_SMA.db -r target.fa
```

#### `meta`

Wraps `pod5 view` to extract `read_id` and `end_reason` into a TSV.

| Flag             | Description                                     |
| ---------------- | ----------------------------------------------- |
| `-i`, `--input`  | Pod5 input directory (**required**)             |
| `-o`, `--output` | Output TSV path (default: `Output/summary.tsv`) |

```bash
ont-sma meta -i raw_pod5/ -o Output/summary.tsv
```

#### `ingest`

Streams BAM reads, computes per-read metrics (basecall quality, Levenshtein distance), and inserts rows into the `Reads` table. Commits in batches of 10 000 reads.

| Flag             | Description                                                           |
| ---------------- | --------------------------------------------------------------------- |
| `-b`, `--bam`    | Input BAM file (**required**)                                         |
| `-d`, `--db`     | Target SQLite database (**required**)                                 |
| `-m`, `--meta`   | Metadata TSV from `meta` step (**required**)                          |
| `--len-min MULT` | Minimum read length as a multiplier of target length (default: `0.5`) |
| `--len-max MULT` | Maximum read length as a multiplier of target length (default: `2.0`) |

```bash
ont-sma ingest \
  -b FAL12345_20260129-IF_SMA_sup_v5.2.0_trim1_0.bam \
  -d Output/SMA_FAL12345_20260129-IF_SMA.db \
  -m Output/summary.tsv
```

The BAM filename encodes experiment metadata and must follow the format:

```txt
{exp_id}_{tier}_v{model_ver}_trim{0|1}_{mod_bitflag}.bam
```

#### `merge`

Merges one or more per-run databases into a master database using `INSERT OR IGNORE`.

| Flag             | Description                                   |
| ---------------- | --------------------------------------------- |
| `-o`, `--output` | Master output database (**required**)         |
| `inputs`         | One or more source `.db` files (**required**) |

```bash
ont-sma merge -o master.db run1.db run2.db run3.db
```

---

## `config.yml`

Used by `ont-sma run`. Copy and edit `config.yml`:

```yaml
exp_id:       "FAL12345_20260129-IF_SMA"   # FlowCellID_SampleID_Alias
bam:          "path/to/reads.bam"          # Raw uBAM
model_tier:   "sup"                        # Basecalling model tier (e.g. "sup", "hac", "fast")
model_ver:    "5.0.0"                      # Basecalling model version
trim:         "yes"                        # Set trimming status on/off ("yes"/"no")
mods_bitflag: "10"
pod5_dir:     "path/to/pod5/"              # Raw Pod5 directory
ref:          "path/to/target.fa"          # Single-sequence reference FASTA

outdir:       "Output"                     # Directory for DB output
summary_tsv:  "Output/summary.tsv"         # Intermediate: pod5 end-reason metadata

len_min_mult: 0.5
len_max_mult: 2.0
```

---

## Workflow

```txt
Input:
  Unaligned BAM   e.g. FAL12345_20260129-IF_SMA_sup_v5.2.0_trim1_0.bam
  Pod5 Directory  e.g. raw_pod5/
  Reference FASTA e.g. target.fa
        │
        ▼
┌────────────────────────┐
│  1. mkdb               │  Creates Output/SMA_<EXP_ID>.db
└────────────────────────┘
        │
        ▼
┌────────────────────────┐
│  2. init               │  Inserts target sequence into DB
└────────────────────────┘
        │
        ▼
┌────────────────────────┐
│  3. meta               │  Extracts end_reason → Output/summary.tsv
└────────────────────────┘
        │
        ▼
┌────────────────────────┐
│  4. ingest             │  Streams reads, computes metrics, populates DB
└────────────────────────┘
        │
        ▼
  Output/SMA_<EXP_ID>.db
```

---

## Database Schema

### `Reads` Table

Contains metrics for every read processed.

| Column        | Type          | Description                                   |
| ------------- | ------------- | --------------------------------------------- |
| `uniq_id`     | **TEXT (PK)** | Composite unique identifier.                  |
| `exp_id`      | TEXT (FK)     | Experiment ID.                                |
| `tgt_id`      | TEXT (FK)     | Target Sequence ID.                           |
| `read_id`     | TEXT          | Original ONT Read UUID.                       |
| `readlen`     | INT           | Length of the Read.                           |
| `readseq`     | TEXT          | The Basecalled Read Sequence.                 |
| `readscr`     | TEXT          | Phred+33 ASCII quality string.                |
| `model_tier`  | TEXT          | Basecaller tier (`f`=fast, `h`=hac, `s`=sup). |
| `model_ver`   | TEXT          | Basecaller version (e.g., 5.2.0).             |
| `trim`        | INT           | Barcode trimming status (0 or 1).             |
| `mod_bitflag` | INT (FK)      | Integer sum of modification flags.            |
| `ed`          | INT           | Levenshtein Distance vs Target.               |
| `q_bc`        | REAL          | Probability-averaged basecall quality.        |
| `q_ld`        | REAL          | Levenshtein quality vs Target.                |
| `ER`          | TEXT          | End Reason (from `pod5 view`).                |

### `Target` Table

Stores the single target sequence for this database.

| Column       | Type          | Description                    |
| ------------ | ------------- | ------------------------------ |
| `tgt_id`     | **TEXT (PK)** | Target Sequence ID (defline).  |
| `tgt_refseq` | TEXT          | The Target Sequence content.   |
| `tgt_reflen` | INT           | Length of the Target Sequence. |

### `Mods` Table

Static lookup table.

| Column        | Type         | Description                                                     |
| ------------- | ------------ | --------------------------------------------------------------- |
| `mod_bitflag` | **INT (PK)** | 4-bit flag: bit 3 = `6mA`, bits 2–0 = C-mod enum (0=none, 1–4). |
| `mods`        | TEXT         | Modifications present.                                          |

### `Exp` Table

Stores experiment metadata parsed from the `exp_id`.

| Column         | Type          | Description                                         |
| -------------- | ------------- | --------------------------------------------------- |
| `exp_id`       | **TEXT (PK)** | Experiment ID. {flow_cell_id}\_{sample_id}\_{alias} |
| `flow_cell_id` | TEXT          | Flow Cell ID.                                       |
| `sample_id`    | TEXT          | Library Prep Info (YYYYMMDD-INITIALS).              |
| `alias`        | TEXT          | Experiment Alias.                                   |
| `exp_desc`     | TEXT          | Experiment description.                             |
