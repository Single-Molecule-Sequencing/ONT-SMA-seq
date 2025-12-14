# ONT-SMA-seq

The SMA-seq protocol for Oxford Nanopore Technology, in pure Python and SQLite database.

## Overview

This pipeline processes Oxford Nanopore Technology (ONT) sequencing data through a three-step workflow:

1. **mkdb.py** - Initialize SQLite database schema
2. **inputInit.py** - Standardize input file paths and extract metadata
3. **ingest.py** - Process reads, calculate metrics, and populate database

## Installation

Install required dependencies:

```bash
pip install -r requirements.txt
```

Required packages:
- pysam >= 0.20.0
- pod5 >= 0.3.0
- edlib >= 1.3.9

## Usage

### Step 1: Create Database

```bash
python mkdb.py <exp_id>
```

Creates `SMA_{exp_id}.db` with the required schema and lookup tables.

### Step 2: Standardize Inputs

```bash
python inputInit.py --bam <path_to_bam> --pod5 <path_to_pod5_dir> --ref <path_to_ref_fasta>
```

BAM filename must follow the naming convention:
`{exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam`

Creates standardized symlinks in the `Input/` directory.

### Step 3: Process Reads

```bash
python ingest.py <exp_id>
```

Processes reads, calculates metrics (quality scores, Levenshtein distance), and populates the database. Creates tagged output BAM in `Output/` directory.

## Testing

Run the end-to-end test to verify the pipeline:

```bash
./run_test.sh
# or
python test_e2e.py
```

This test:
- Creates synthetic test data (reference sequences and BAM reads)
- Runs all three pipeline steps
- Verifies database population and output files
- Cleans up test artifacts

Expected output: All tests should pass with `✓ ✓ ✓ ALL TESTS PASSED ✓ ✓ ✓`

See [TESTING.md](TESTING.md) for detailed testing documentation.

## Documentation

See [overhaul.md](overhaul.md) for detailed documentation on:
- Database schema
- Modification bitflags
- Metric calculations
- Input file formats
