# ONT-SMA-seq

The Single-Molecule-Accuracy-seq protocol for Oxford Nanopore Technology experiments, implemented in pure Python with a SQLite database backend. This workflow processes an unaligned BAM file and its parent Pod5 files for a single target sequence within an overarching experiment, storing read metrics and metadata into a structured database for analysis.

## Setup

Set up and activate the conda environment using the provided `env/env.yml` file.

```bash
conda env create -f env/env.yml
conda activate ont-sma-seq
```

Note: `edlib` might fail to be imported despite proper conda installation. If this issue is encountered, reinstall it with `pip`.

Within activated conda env:

```bash
pip install edlib --force-reinstall --no-cache-dir
```

## Manifest

All scripts are located in the `bin/` directory.

* `mkdb.py`: Initializes the SQLite database schema, parses the Experiment ID (FlowCell/Sample), and populates static tables.
* `inputInit.py`: Standardizes inputs. Sanitizes the Reference FASTA (ensuring only **one** sequence exists) and creates symlinks for BAMs and Pod5s.
* `extractMeta.py`: Extracts lightweight metadata (Read ID and End Reason) from Pod5 files using the `pod5` CLI.
* `ingest.py`: The core processing script. It parses inputs, calculates quality and alignment metrics (Levenshtein) for all reads against the target, tags BAMs with End Reasons, and populates the database.

## Usage

### `mkdb.py`

Creates the initial database file and parses the Experiment ID.
**Note:** `exp_id` should follow the format `FlowCellID_SampleID` (e.g., `FAL12345_20260129_IF`).

* `-e`, `--expid`: Experiment ID (Required).
* `-o`, `--outdir`: Output directory (Default: `Output`).

```bash
python3 bin/mkdb.py -e <EXP_ID> -o <PATH_TO_OUTPUT_DIR>
```

### `inputInit.py`

Standardizes file paths and sanitizes the reference.

* `-b`, `--bam`: Path to raw uBAM file.
* `-p`, `--pod5_dir`: Path to raw Pod5 directory.
* `-r`, `--ref`: Path to Reference FASTA.

```bash
python3 bin/inputInit.py -b <PATH_TO_BAM> -p <PATH_TO_POD5_DIR> -r <PATH_TO_REF_FASTA>
```

### `extractMeta.py`

Wrapper for `pod5 view`. Extracts End Reason metadata required for the DB.

* `-i`, `--input`: Path to `Input/pod5` (Default).
* `-o`, `--output`: Path to `Input/summary.tsv` (Default).

```bash
python3 bin/extractMeta.py -i <PATH_TO_INPUT_DIR> -o <PATH_TO_OUTPUT_TSV>
```

### `prepare.py`

Discovers MinKNOW runs, merges POD5/BAM data, initializes the database, aligns all reads against target references, and generates diagnostic QC plots.

* `-d`, `--expdir`: Path to experiment directory (Required).
* `-e`, `--expid`: Experiment ID (Required).
* `-r`, `--ref`: Multi-sequence FASTA with target sequences (Required).
* `-o`, `--outdir`: Output directory (Default: `Output`).
* `--force-rebasecall`: Always re-basecall from POD5.
* `--dry-run`: Print merge plan without executing.

```bash
python3 bin/prepare.py -d /path/to/experiment -e EXP_ID -r targets.fa -o Output/
```

### `ingest.py`

Calculates metrics and ingests data. Requires the database created by `mkdb.py`.

* `-e`, `--expid`: Experiment ID (must match `mkdb.py`).
* `-b`, `--bam`: Path to `Input/reads.bam`.
* `-s`, `--summary`: Path to `Input/summary.tsv`.
* `-r`, `--ref`: Path to `Input/target.fa`.
* `-d`, `--database`: Path to the SQLite DB file.
* `-o`, `--output_bam`: Path to output tagged BAM file.

```bash
python3 bin/ingest.py -e <EXP_ID> -b <PATH_TO_BAM> -s <PATH_TO_SUMMARY_TSV> -r <PATH_TO_REF_FASTA> -d <PATH_TO_DB> -o <PATH_TO_OUTPUT_BAM>
```

## Workflow

The workflow ingests raw ONT data from a single experiment into a SQLite database.

```txt
Input:
- Unaligned BAM (e.g., 'FAL12345_20260129_IF_sup_v5.2.0_trim1_0.bam')
- Pod5 Directory (e.g., 'raw_pod5/')
- Single-Seq FASTA (e.g., 'target.fa')
  │
  ▼
┌──────────────────────────┐
│  1. Initialize Database  │
└──────────────────────────┘
  │
  │ Input: Experiment ID (FlowCell_Sample)
  │ Output: 'Output/SMA_<EXP_ID>.db'
  │
  └─► python3 bin/mkdb.py -e <EXP_ID> -o Output
  │
  ▼
┌─────────────────────────┐
│  2. Standardize Inputs  │
└─────────────────────────┘
  │
  │ Input: Raw BAM, Pod5 Dir, Target FASTA
  │ Output: 'Input/' with symlinks & sanitized 'target.fa'
  │
  └─► python3 bin/inputInit.py -b <BAM_FILE> -p <POD5_DIR> -r <TARGET_FASTA>
  │
  ▼
┌───────────────────────┐
│  3. Extract Metadata  │
└───────────────────────┘
  │
  │ Input: 'Input/pod5/' (Symlinked)
  │ Output: 'Input/summary.tsv'
  │
  └─► python3 bin/extractMeta.py -i Input/pod5/ -o Input/summary.tsv
  │
  ▼
┌──────────────────┐
│  4. Ingest Data  │
└──────────────────┘
  │
  │ Input: 'Input/reads.bam', 'Input/summary.tsv', 'Input/target.fa'
  │ Output: Populated DB, Tagged BAM
  │
  └─► python3 bin/ingest.py -e <EXP_ID> -b Input/reads.bam -s Input/summary.tsv -r Input/target.fa -d Output/SMA_<EXP_ID>.db -o Output/tagged.bam
  │
  ▼
Final Output:
  * Populated SQLite DB ('Output/SMA_<EXP_ID>.db')
  * End-Reason Tagged BAM ('Output/tagged.bam')
```

## Database Schema

### `Reads` Table

Contains metrics for every read processed.

| Column        | Type          | Description                            |
| ------------- | ------------- | -------------------------------------- |
| `uniq_id`     | **TEXT (PK)** | Composite unique identifier.           |
| `exp_id`      | TEXT (FK)     | Experiment ID.                         |
| `tgt_id`      | TEXT (FK)     | Target Sequence ID.                    |
| `read_id`     | TEXT          | Original ONT Read UUID.                |
| `readseq`     | TEXT          | The Basecalled Read Sequence.          |
| `readlen`     | INT           | Length of the Read.                    |
| `model_tier`  | TEXT          | Basecaller model tier (fast/hac/sup).  |
| `model_ver`   | TEXT          | Basecaller version (e.g., 5.2.0).      |
| `trim`        | INT           | Barcode trimming status (0 or 1).      |
| `mod_bitflag` | INT (FK)      | Integer sum of modification flags.     |
| `ed`          | INT           | Levenshtein Distance vs Target.        |
| `q_bc`        | REAL          | Probability-averaged basecall quality. |
| `q_ld`        | REAL          | Levenshtein quality vs Target.         |
| `ER`          | TEXT          | End Reason (from `pod5 view`).         |

### `Target` Table

Stores the single target sequence for this database.

| Column       | Type          | Description                    |
| ------------ | ------------- | ------------------------------ |
| `tgt_id`     | **TEXT (PK)** | Target Sequence ID (defline).  |
| `tgt_refseq` | TEXT          | The Target Sequence content.   |
| `tgt_reflen` | INT           | Length of the Target Sequence. |

### `Mods` Table

Static lookup table.

| Column        | Type         | Description                   |
| ------------- | ------------ | ----------------------------- |
| `mod_bitflag` | **INT (PK)** | Sum of modification flags (). |
| `mods`        | TEXT         | Modifications present.        |

### `Exp` Table

Stores experiment metadata parsed from the `exp_id`.

| Column         | Type          | Description                                         |
| -------------- | ------------- | --------------------------------------------------- |
| `exp_id`       | **TEXT (PK)** | Experiment ID. {flow_cell_id}\_{sample_id}\_{alias} |
| `flow_cell_id` | TEXT          | Flow Cell ID.                                       |
| `sample_id`    | TEXT          | Library Prep Info (YYYYMMDD-INITIALS).              |
| `alias`        | TEXT          | Experiemnt Alias.                                   |
| `exp_desc`     | TEXT          | Experiment description.                             |

## Future Plans

* **Script Renaming**: Removing `.py` extensions for cleaner CLI usage.
* **Database Merging**: Scripts to merge multiple Per-Target databases into a Per-Experiment database.
* **Multithreading**: Parallel creation of DBs based on a config file.
* **Querying & Plotting**: Post-hoc SQL analysis scripts to filter reads and generate reports.
