# ONT-SMA-seq

The SMA-seq protocol for Oxford Nanopore Technology, in pure Python and SQLite database.

## Documentation

- **[WORKFLOW.md](WORKFLOW.md)** - Comprehensive workflow flowcharts, database schema diagrams, and data flow documentation
- **[DIAGRAMS.md](DIAGRAMS.md)** - Quick reference diagrams for pipeline architecture and data structures
- **[overhaul.md](overhaul.md)** - Detailed technical specification and implementation notes

## Quick Overview

This pipeline processes Oxford Nanopore sequencing data through three main scripts:

1. **mkdb.py** - Initialize SQLite database with schema
2. **inputInit.py** - Standardize input files and extract metadata
3. **ingest.py** - Process reads, calculate metrics, and populate database

See [WORKFLOW.md](WORKFLOW.md) for detailed flowcharts and [DIAGRAMS.md](DIAGRAMS.md) for visual architecture.
