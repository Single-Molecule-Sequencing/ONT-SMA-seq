# ONT-SMA-seq

The SMA-seq protocol for Oxford Nanopore Technology, in pure Python and SQLite database.

## Overview

ONT-SMA-seq is a Python-based pipeline for processing Oxford Nanopore sequencing data, calculating quality metrics, and storing results in an SQLite database. The pipeline consists of three main scripts:

1. **`mkdb.py`** - Database initialization
2. **`inputInit.py`** - Input file standardization
3. **`ingest.py`** - Main processing and data ingestion

## Installation

Install the required dependencies:

```bash
pip install -r requirements.txt
```

For development and testing:

```bash
pip install -r requirements-dev.txt
```

## Pipeline Usage

### Step 1: Create Database

Initialize the SQLite database for your experiment:

```bash
python mkdb.py <exp_id> [--output-dir <dir>]
```

This creates `SMA_{exp_id}.db` with the required schema and populates static lookup tables.

### Step 2: Standardize Input Files

Create symlinks to standardize input file locations:

```bash
python inputInit.py <bam_file> <pod5_dir> <reference_fasta> [--output-dir Input]
```

**BAM file naming convention:**
`{exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam`

Example: `exp001_s_v5.2.0_1_6mA.bam`

### Step 3: Process Reads

Process reads, calculate metrics, and populate the database:

```bash
python ingest.py <exp_id> [--input-dir Input] [--output-dir Output] [--db-dir .]
```

This script:
- Parses reference sequences
- Extracts Pod5 metadata
- Calculates quality metrics (q_bc, q_ld)
- Tags BAM files with end reasons
- Populates the database with all read information

## Testing

Run the test suite:

```bash
# Run all tests
pytest tests/

# Run specific test files
pytest tests/test_mkdb.py -v
pytest tests/test_inputInit.py -v
pytest tests/test_ingest.py -v
pytest tests/test_integration.py -v

# Run with coverage
pytest tests/ --cov=. --cov-report=html
```

## Database Schema

The pipeline creates four main tables:

- **Exp** - Experiment metadata
- **Refseq** - Reference sequences
- **Mods** - Modification bitflags
- **Reads** - Read-level data and metrics

See `overhaul.md` for detailed schema documentation.

## Modification Bitflags

Modifications are stored as integer bitflags:

| Modification | Bit Value |
|--------------|-----------|
| non          | 0         |
| 6mA          | 1         |
| 5mCG_5hmCG   | 2         |
| 5mC_5hmC     | 4         |
| 4mC_5mC      | 8         |
| 5mC          | 16        |

Combinations are sums (e.g., 6mA + 5mC_5hmC = 5).

## License

See LICENSE file for details.
