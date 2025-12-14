# ONT-SMA-seq

The SMA-seq protocol for Oxford Nanopore Technology, in pure Python and SQLite database.

## Documentation

- **[SUMMARY.md](SUMMARY.md)** - Quick reference with ASCII diagrams and text-based flowcharts
- **[WORKFLOW.md](WORKFLOW.md)** - Comprehensive workflow flowcharts, database schema diagrams, and data flow documentation
- **[DIAGRAMS.md](DIAGRAMS.md)** - Quick reference diagrams for pipeline architecture and data structures
- **[TESTING.md](TESTING.md)** - Testing guide and test suite documentation
- **[overhaul.md](overhaul.md)** - Detailed technical specification and implementation notes

## Quick Overview

This pipeline processes Oxford Nanopore sequencing data through three main scripts:

1. **mkdb.py** - Initialize SQLite database with schema
2. **inputInit.py** - Standardize input files and extract metadata
3. **ingest.py** - Process reads, calculate metrics, and populate database

See [SUMMARY.md](SUMMARY.md) for a text-based overview, [WORKFLOW.md](WORKFLOW.md) for detailed flowcharts, and [DIAGRAMS.md](DIAGRAMS.md) for visual architecture.

## Quick Start

```bash
# 1. Create database
python mkdb.py EXP001

# 2. Standardize inputs
python inputInit.py --bam data/EXP001_h_v5.2.0_1_6mA.bam --pod5 data/pod5/ --ref data/reference.fa

# 3. Process data
python ingest.py EXP001
```

## Dependencies

```bash
pip install pysam pod5 edlib
```
