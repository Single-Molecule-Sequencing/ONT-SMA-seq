# ONT-SMA-seq

The Single-Molecule-Accuracy-seq protocol for Oxford Nanopore Technology experiments, implemented in pure Python with a SQLite database backend. This workflow processes an unaligned BAM file and its parent Pod5 files from a single ONT experiment, storing read metrics and metadata into a structured database for analysis.

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

* `mkdb.py`: Initializes the SQLite database schema and populates the static modification bitflag lookup table.
* `inputInit.py`: Standardizes input file paths and directory structures by creating symlinks for BAMs, Pod5s, and reference sequences.
* `extractMeta.py`: Extracts lightweight metadata (Read ID and End Reason) from Pod5 files using the `pod5` CLI.
* `ingest.py`: The core processing script. It parses inputs, calculates quality and alignment metrics, tags BAMs with End Reasons, and populates the database.

## Usage

### `mkdb.py`

Creates the initial database file.

* `-e`, `--expid`: Experiment ID (Required).
* `-o`, `--outdir`: Output directory (Default: `Output`).

```bash
python3 bin/mkdb.py -e <EXP_ID> -o <PATH_TO_OUTPUT_DIR>
```

### `inputInit.py`

Sets up the `Input/` directory with standardized symlinks.

* `-b`, `--bam`: Path to raw uBAM file. Must follow naming convention: `{exp_id}_{model}_{ver}_{trim}_{mods}.bam` (Required).
* `-p`, `--pod5_dir`: Path to raw Pod5 directory (Required).
* `-r`, `--ref`: Path to Reference FASTA containing two sequences of distinct lengths (Required).

```bash
python3 bin/inputInit.py -b <PATH_TO_BAM> -p <PATH_TO_POD5_DIR> -r <PATH_TO_REF_FASTA>
```

### `extractMeta.py`

Generates a summary TSV of read metadata.

* `-i`, `--input`: Path to standardized input directory containing `.pod5` files (Required).
* `-o`, `--output`: Path to output summary TSV file (Default: `./pod5_ER_summary.tsv`).

```bash
python3 bin/extractMeta.py -i <PATH_TO_INPUT_DIR> -o <PATH_TO_OUTPUT_TSV>
```

### `ingest.py`

Populates the database with read metrics and generates a tagged BAM.

* `-e`, `--expid`: Experiment ID (Required).
* `-b`, `--bam`: Path to standardized input uBAM file (Required).
* `-s`, `--summary`: Path to input Pod5 summary TSV (Required).
* `-r`, `--ref`: Path to standardized Reference FASTA (Required).
* `-d`, `--database`: Path to target SQLite database file (Required).
* `-o`, `--output_bam`: Path to output tagged BAM file (Required).
* `-k`, `--tolerance`: Length tolerance (+/- bp) for RefSeq matching (Default: 150).

```bash
python3 bin/ingest.py -e <EXP_ID> -b <PATH_TO_BAM> -s <PATH_TO_SUMMARY_TSV> -r <PATH_TO_REF_FASTA> -d <PATH_TO_DB> -o <PATH_TO_OUTPUT_BAM> -k <TOLERANCE>
```

## Workflow

The workflow ingests raw ONT data from a single experiment into a SQLite database.

```txt
Input:
- Unaligned BAM (e.g., '<EXP_ID>_sup_v5.2.0_trim1_0.bam')
- Pod5 Directory (e.g., 'raw_pod5/')
- Reference FASTA (e.g., 'reference.fa')
  │
  ▼
┌──────────────────────────┐
│  1. Initialize Database  │
└──────────────────────────┘
  │
  │ Input: Experiment ID
  │ Output: 'Output/SMA_<EXP_ID>.db'
  │
  └─► python3 bin/mkdb.py -e <EXP_ID> -o Output
  │
  ▼
┌─────────────────────────┐
│  2. Standardize Inputs  │
└─────────────────────────┘
  │
  │ Input: Raw BAM, Pod5 Dir, Ref.fa
  │ Output: 'Input/' directory with symlinks
  │
  └─► python3 bin/inputInit.py -b <BAM_FILE> -p <POD5_DIR> -r <REF_FASTA>
  │
  ▼
┌───────────────────────┐
│  3. Extract Metadata  │
└───────────────────────┘
  │
  │ Input: 'Input/<EXP_ID>_pod5/'
  │ Output: 'Output/<EXP_ID>_summary.tsv'
  │
  └─► python3 bin/extractMeta.py -i Input/<EXP_ID>_pod5/ -o Output/<EXP_ID>_summary.tsv
  │
  ▼
┌──────────────────┐
│  4. Ingest Data  │
└──────────────────┘
  │
  │ Input: Symlinked BAM, Summary TSV, Symlinked FASTA, Database
  │ Output: Populated DB, Tagged BAM
  │
  └─► python3 bin/ingest.py -e <EXP_ID> -b Input/<BAM_FILE> -s Output/<SUMMARY_TSV> -r Input/<REF_FASTA> -d Output/<DATABASE> -o Output/<TAGGED_BAM>
  │
  ▼
Final Output:
- Populated SQLite DB ('Output/SMA_<EXP_ID>.db')
- End-Reason Tagged BAM ('Output/<EXP_ID>_ER.bam')
```

## Database Schema

### `Reads` Table

| Column        | Type          | Description                                |
| :------------ | :------------ | :----------------------------------------- |
| `uniq_id`     | **TEXT (PK)** | Composite unique identifier.               |
| `exp_id`      | TEXT (FK)     | Experiment ID.                             |
| `refseq_id`   | TEXT (FK)     | RefSeq ID. `NULL` if length out of range.  |
| `read_id`     | TEXT          | Original ONT Read UUID.                    |
| `readseq`     | TEXT          | The Basecalled Read Sequence.              |
| `readlen`     | INT           | Length of the Read.                        |
| `model_tier`  | TEXT          | 's', 'h', or 'f'.                          |
| `model_ver`   | TEXT          | e.g., '5.2.0'.                             |
| `trim`        | INT           | 1 (True) or 0 (False).                     |
| `mod_bitflag` | INT (FK)      | Integer sum of modification flags.         |
| `ed`          | INT           | Levenshtein Distance. `NULL` if no RefSeq. |
| `q_bc`        | REAL          | Probability-averaged basecall quality.     |
| `q_ld`        | REAL          | Levenshtein quality. `NULL` if no RefSeq.  |
| `ER`          | TEXT          | End Reason (from `pod5 view`).             |

### `Mods` Table

| Column        | Type         | Description                        |
| :------------ | :----------- | :--------------------------------- |
| `mod_bitflag` | **INT (PK)** | Sum of modification flags ($2^n$). |
| `mods`        | TEXT         | Modifications present.             |

### `Exp` Table

| Column     | Type          | Description             |
| :--------- | :------------ | :---------------------- |
| `exp_id`   | **TEXT (PK)** | Experiment ID.          |
| `exp_desc` | TEXT          | Experiment description. |

### `Refseq` Table

| Column      | Type          | Description                      |
| :---------- | :------------ | :------------------------------- |
| `refseq_id` | **TEXT (PK)** | Reference Sequence ID (defline). |
| `refseq`    | TEXT          | Reference Sequence.              |
| `reflen`    | INT           | Reference Sequence Length.       |

## Future Plans

* **Script Renaming**: Removing `.py` extensions for cleaner CLI usage.
* **Query Scripts**: Developing scripts to query the database for specific metrics and data extraction.
* **Database Merging**: Tools to merge multiple experiment databases into a single master database.
* **Plotting**: Visualization scripts to generate reports directly from database queries.
