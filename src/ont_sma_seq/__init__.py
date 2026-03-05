"""
ont_sma_seq — ONT Single-Molecule-Accuracy Sequencing pipeline.

Provides CLI entry point `ont-sma` with subcommands:
  mkdb    — initialise SQLite database and schema
  init    — lock a reference FASTA into the Target table
  meta    — extract end_reason metadata from Pod5 files
  ingest  — stream BAM reads and populate the database
  merge   — combine multiple per-run databases into a master DB
  run     — chain all steps from a config.yml
"""

__version__ = "0.2.0"
