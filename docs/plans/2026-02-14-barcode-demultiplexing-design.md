# Barcode Demultiplexing for SMA-seq

**Date:** 2026-02-14
**Status:** Approved
**Repos:** ONT-SMA-seq, Reference_Fasta_Generator

---

## Problem

`ingest.py` processes all reads against a single target reference. SMA-seq experiments use duplexed native barcoding where each read carries an upstream barcode near its start and a reverse-complemented downstream barcode near its end. Reads must be classified by their barcode pair to assign each read to the correct target sequence. Currently there is no classification step, so all reads get compared to one target regardless of which construct they belong to.

## Decision

Integrate barcode classification directly into `ingest.py` using a two-stage hybrid approach. No separate `demux.py` script. All reads are ingested using untrimmed basecalled sequences.

## Pipeline (Unchanged Flow)

```
mkdb.py -> inputInit.py -> extractMeta.py -> ingest.py (enhanced)
```

Only `ingest.py` changes. The other three scripts remain as-is.

## Classification Algorithm

### Stage 1: Dual-End Barcode Alignment (Always)

For each read:

1. **Start barcode:** Align first ~100bp against all expected barcode sequences using edlib semi-global alignment (`HW` mode). Record best match: `bc_start_id`, `bc_start_ed`, `bc_start_conf`.

2. **End barcode:** Align last ~100bp against reverse complement of all expected barcode sequences using edlib `HW` mode. Record best match: `bc_end_id`, `bc_end_ed`, `bc_end_conf`.

3. **Target assignment:** Look up `(bc_start_id, bc_end_id)` pair in the sample sheet to get the target alias. Use alias to select the correct reference FASTA for `ed`/`q_ld` calculation.

### Stage 2: Full Construct Alignment (Conditional)

Only triggered when:
- User passes `--full-construct` flag, OR
- The barcode design has a single barcode mapped to 2+ different target sequences (detected automatically from the sample sheet)

Aligns the entire read against all full construct reference FASTAs. Best match overrides the Stage 1 assignment.

### Confidence Metric

```
confidence = 1.0 - (barcode_ed / barcode_length)
```

Where `barcode_length` = 24bp (ONT native barcodes). Perfect match = 1.0. Interpretable for wet lab users as a percentage.

## Inputs

### New Arguments for `ingest.py`

| Argument | Short | Description | Required |
|----------|-------|-------------|----------|
| `--sample-sheet` | `-ss` | MinKNOW sample sheet CSV | Yes (for classification) |
| `--ref-dir` | `-rd` | Directory of per-target reference FASTAs | Yes (for classification) |
| `--full-construct` | | Force full construct alignment for all reads | No |
| `--split-bams` | | Output directory for per-barcode BAMs | No |
| `--tag` | | Write BC:Z: tags to output BAM | No |

Existing `-r / --ref` (single FASTA) still works for backward-compatible single-target mode.

### Sample Sheet Format

MinKNOW sample sheet CSV with duplexed barcode column:

```csv
flow_cell_id,kit,sample_id,experiment_id,barcode,alias
FAL12345,SQK-NBD114-96,20260129_IF,exp001,barcode05--barcode10,CYP2D6_v04_fwd
FAL12345,SQK-NBD114-96,20260129_IF,exp001,barcode10--barcode05,CYP2D6_v04_rev
```

The `barcode` column uses `barcodeNN--barcodeNN` format for duplexed pairs. The `alias` column matches FASTA filenames in the reference directory.

### Reference Directory

```
references/
  CYP2D6_v04_fwd.fasta    # Target-only sequence, header matches alias
  CYP2D6_v04_rev.fasta
  ...
```

Full construct FASTAs (from Reference_Fasta_Generator) in a separate directory, only needed when `--full-construct` is used or when barcode ambiguity is detected.

## DB Schema Changes

New columns on the `Reads` table:

| Column | Type | Description |
|--------|------|-------------|
| `bc_start_id` | TEXT | Upstream barcode assignment (e.g., nb05) |
| `bc_start_ed` | INT | Edit distance to assigned start barcode |
| `bc_start_conf` | REAL | Classification confidence (start) |
| `bc_end_id` | TEXT | Downstream barcode assignment (e.g., nb10) |
| `bc_end_ed` | INT | Edit distance to assigned end barcode |
| `bc_end_conf` | REAL | Classification confidence (end) |

The `Target` table now holds multiple targets (one per alias). No schema change needed -- it already supports multiple entries.

The `tgt_id` in `Reads` is set from the sample sheet lookup `(bc_start, bc_end) -> alias -> tgt_id`, not from a single reference.

## All Reads Are Classified

Every read is classified and ingested. No reads are discarded. Poorly matched reads (low confidence) remain in the DB and serve as diagnostics for run quality issues.

Query poorly classified reads:
```sql
SELECT * FROM Reads WHERE bc_start_conf < 0.7 OR bc_end_conf < 0.7;
```

## Optional Post-Classification Outputs

- `--split-bams DIR` -- Write per-barcode BAMs (one per target alias)
- `--tag` -- Add `BC:Z:nb05` and `BA:Z:alias` tags to the output BAM
- Classification summary printed to stdout: barcode counts, mean confidence, unmatched pair counts

## Error Handling

Clear, actionable error messages for wet lab users:
- "Alias 'CYP2D6_v04' in sample sheet has no matching FASTA in references/"
- "Barcode pair nb05--nb10 not found in sample sheet"
- "Sample sheet has no 'barcode' column -- expected MinKNOW format"

## Backward Compatibility

When `-ss` and `-rd` are not provided, `ingest.py` behaves exactly as before: single target from `-r`, no barcode classification. The new barcode columns are NULL in this case.

## Dependencies

No new dependencies. Uses existing:
- `edlib` -- barcode alignment (HW mode for semi-global)
- `pysam` -- BAM I/O
- `pandas` -- sample sheet CSV parsing
- `sqlite3` -- DB operations

Barcode sequences sourced from Reference_Fasta_Generator's `Make_Fasta.py` (96 ONT native barcodes, hardcoded).
