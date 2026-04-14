# CLAUDE.md — ONT-SMA-seq

Reference implementation of the Single-Molecule-Accuracy-seq (SMA-seq) protocol
for Oxford Nanopore data. A Python package + SQLite backend that processes an
unaligned BAM + parent Pod5 files for a single target sequence and stores
per-read metrics for downstream analysis.

**Status:** Public, citable (Zenodo DOI 10.5281/zenodo.18872468).

## Repo layout

```
ONT-SMA-seq/
├── src/            # Python package (ont-sma CLI lives here)
├── config.yml      # Default protocol configuration
├── env/env.yml     # Conda environment spec
├── pyproject.toml  # Package metadata + entry points
└── README.md       # Installation + usage (authoritative)
```

## Setup

```bash
conda env create -f env/env.yml
conda activate ont-sma-seq
pip install -e .          # exposes `ont-sma` CLI
pip install edlib --force-reinstall --no-cache-dir   # if edlib import fails
```

## How this relates to the rest of the lab

- **`/sma-pipeline`** (ont-ecosystem skill) is the lab-wide orchestrator for
  SMA-seq runs; it dispatches per-experiment work through `ont-sma` installed
  from this repo. When the orchestrator needs bug fixes in the core SMA logic,
  they land here and propagate via `pip install -e .`.
- **`smaseq-qc`** (separate repo) is the QC-metrics companion and consumes the
  SQLite outputs from ont-sma.
- **`sma-seq-workspace`** is the analysis workspace (notebooks + scripts) that
  runs on top of ont-sma outputs.

## Conventions

- Python ≥3.9, pathlib (not os.path), ruff linter.
- Unit tests: pytest under `tests/` (run via `pytest -q`).
- Database schema changes must be backward compatible or include a migration.

## When NOT to edit here

- If the bug is in an analysis notebook or a figure script, fix it in
  `sma-seq-workspace` or `smaseq-qc` instead; this repo is the library layer.
- If the bug is in lab-wide orchestration (experiment registry, report
  generation), fix it in `ont-ecosystem/skills/sma-pipeline/` instead.
