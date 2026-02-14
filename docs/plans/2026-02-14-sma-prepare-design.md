# Design: SMA-seq Experiment Preparation & Alignment Classification

**Date:** 2026-02-14
**Status:** Approved
**Repo:** `/tmp/ont-sma-seq/`

---

## Problem

The existing SMA-seq pipeline (`mkdb.py` -> `inputInit.py` -> `extractMeta.py` -> `ingest.py`) assumes a single pre-made BAM file as input. In practice, a single SMA-seq experiment produces raw MinKNOW output spread across multiple sequencing runs (stop/restart), multiple barcode subdirectories, and hundreds of POD5/BAM files. There is no tool to consolidate this raw output into a single, validated, consistently-basecalled BAM — nor to perform all-vs-all alignment classification against target sequences with diagnostic visualizations.

## Solution

A new `bin/prepare.py` CLI with a staged pipeline architecture, plus `bin/align.py` (alignment engine) and `bin/qc.py` (visualization module). The tool discovers MinKNOW runs within an experiment directory, builds a merge plan with user confirmation, consolidates POD5/BAM data, initializes the SMA-seq database with provenance, aligns every read against every target reference, and generates publication-quality diagnostic plots for threshold-based classification.

---

## Architecture: Staged Pipeline (Approach A)

Single entry point with 5 sequential stages. Each stage writes progress to a manifest JSON for resume capability.

```
prepare.py
  Stage 1: DISCOVER  -> Scan experiment dir, parse MinKNOW metadata
  Stage 2: PLAN      -> Build merge plan, detect basecalling, present to user
  Stage 3: MERGE     -> Merge POD5s, merge/re-basecall BAMs, sort+index
  Stage 4: INIT      -> Create SMA-seq SQLite DB with provenance
  Stage 5: ALIGN+QC  -> All-vs-all alignment metrics + diagnostic plots
```

### CLI Interface

```bash
python3 bin/prepare.py \
  -d /mnt/d/12282025_IF_DoubleBC_SMA_seq_no_trim \
  -e FAL12345_20251228_IF \
  -r targets.fa \
  -o Output/
```

**Arguments:**

| Flag | Required | Description |
|------|----------|-------------|
| `-d / --expdir` | Yes | Path to experiment directory on disk |
| `-e / --expid` | Yes | Experiment ID (`FlowCell_Sample_Alias` format) |
| `-r / --ref` | Yes | Multi-sequence FASTA with all target sequences |
| `-o / --outdir` | No | Output directory (default: `Output/`) |
| `--force-rebasecall` | No | Skip smart detection, always re-basecall from POD5 |
| `--dry-run` | No | Print merge plan without executing |

### New Files

```
bin/
  prepare.py     (~400 lines — Stages 1-4: discover, plan, merge, init)
  align.py       (~300 lines — Stage 5a: alignment engine + metrics)
  qc.py          (~400 lines — Stage 5b: all visualizations)
```

All existing scripts (`mkdb.py`, `inputInit.py`, `extractMeta.py`, `ingest.py`) remain unchanged and functional for single-BAM workflows.

---

## Stage 1: DISCOVER

Recursively scan `--expdir` for MinKNOW run folders. A valid run folder is identified by the presence of `pod5_pass/` (or `pod5/`) alongside a `sequencing_summary_*.txt`.

**Per-run metadata extracted:**

| Field | Source |
|-------|--------|
| `run_id` | Directory name (e.g., `20251228_2219_MD-100098_FBD69411_34fa833d`) |
| `flow_cell_id` | Parsed from directory name (e.g., `FBD69411`) |
| `device_id` | Parsed from directory name (e.g., `MD-100098`) |
| `start_time` | From sequencing summary header |
| `sample_sheet` | Parsed if present (kit, barcode aliases) |
| `basecall_model` | From BAM `@PG` header lines, if BAMs exist |
| `pod5_count` | Number of POD5 files |
| `pod5_bytes` | Total POD5 size |
| `bam_count` | Number of BAM files |
| `barcode_dirs` | Barcode subdirectories found in `bam_pass/` |

Output: List of `RunInfo` dataclass objects in memory.

---

## Stage 2: PLAN

Build a merge plan from discovered runs using three rules:

### Rule 1: Group by Flow Cell

Runs on the same `flow_cell_id` are candidates for merging (same physical experiment, stop/restart). Runs on different flow cells are flagged as parallel runs requiring separate handling.

### Rule 2: Detect Wash Boundaries

If multiple experiments used the same flow cell with a wash between them, the `protocol_run_id` or `sample_id` in the sample sheet will differ. These are NOT merged. User is prompted to confirm.

### Rule 3: Smart Basecalling Check

For each run's BAMs, parse `@PG` header to extract model name, version, and parameters:
- All BAMs `sup` tier? Same model version? Untrimmed (`--no-trim`)?
- If yes -> action = `REUSE` (concatenate existing BAMs)
- If no -> action = `REBASECALL` (re-basecall from POD5)

### Interactive Confirmation

```
[prepare] Merge Plan for FAL12345_20251228_IF
---------------------------------------------
Run 1: 20251228_2219_MD-100098_FBD69411_34fa833d
  Flow Cell: FBD69411 | POD5: 12 files (54 GB) | BAMs: sup v5.2.0 untrimmed
  Action: REUSE existing BAMs

Run 2: 20251229_1055_MD-100098_FBD69411_5b5c57a9
  Flow Cell: FBD69411 | POD5: 32 files (225 GB) | BAMs: sup v5.2.0 untrimmed
  Action: REUSE existing BAMs

Merge group: [Run 1 + Run 2] -> same flow cell, sequential
Output: Output/merged.bam (sorted, indexed)

Proceed? [Y/n]
```

User can confirm, exclude runs, or force re-basecall.

---

## Stage 3: MERGE

### 3a: POD5 Consolidation

Symlink (not copy) all POD5 files into an organized directory:

```
Output/pod5/
  run1_FBD69411_34fa833d/  -> symlinks to original POD5s
  run2_FBD69411_5b5c57a9/  -> symlinks to original POD5s
```

Then extract end-reason metadata:

```bash
pod5 view Output/pod5/**/*.pod5 --include "read_id,end_reason" --output Output/summary.tsv
```

### 3b: BAM Merge

**Path REUSE** (all BAMs consistent):
1. Collect all BAMs from `bam_pass/` across merge group (all barcodes + unclassified)
2. Concatenate with `pysam.merge()`
3. Sort by read name: `pysam.sort("-n", ...)`
4. Index: `pysam.index(...)`
5. Validate: merged read count == sum of source BAM read counts

**Path REBASECALL** (inconsistent basecalling):
1. Print detected inconsistency
2. Prompt for target model (default: latest sup, untrimmed, `--emit-moves`)
3. Run `dorado basecaller` on all POD5s
4. Sort and index output BAM

---

## Stage 4: INIT

Create `Output/SMA_{exp_id}.db` using the same schema as `mkdb.py`, plus populate `RunMetadata` with rich provenance from Stage 1.

**Tables created:**
- `Reads` — populated later by `ingest.py`
- `Mods` — static bitflag lookup (pre-populated, 10 entries)
- `Exp` — experiment metadata from `--expid`
- `Target` — populated later by `ingest.py`
- `RunMetadata` — populated now with per-run provenance
- `ReadRun` — populated later by `ingest.py`

**RunMetadata population:**

```sql
INSERT INTO RunMetadata (
  run_id, flow_cell_id, device_id, sample_id,
  experiment_id, kit, protocol_run_id, start_time,
  basecall_model, source_bam_count, source_bam_paths,
  merge_timestamp
) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
```

Each run in the merge group gets its own row for full traceability.

---

## Stage 5a: ALIGN — All-vs-All Alignment Classification

### Alignment Engine

Use `edlib` in `NW` (Needleman-Wunsch / global alignment) mode with `task="path"` to obtain the full CIGAR string.

```python
result = edlib.align(read_seq, ref_seq, mode="NW", task="path")
```

**Critical constraint:** Forward-strand only. No reverse complement. Each reference sequence in the FASTA already represents a specific strand orientation. Reads either match in the 5'->3' direction or they don't.

For every read in `merged.bam`, align against every target in `-r` and compute the full metric set.

### Metric Set (per read x reference pair)

#### Whole-Alignment Metrics

| Metric | Type | Description |
|--------|------|-------------|
| `ed` | INT | Raw edit distance |
| `ned` | REAL | Normalized edit distance: `ed / ref_len` |
| `identity` | REAL | `matches / alignment_length` from CIGAR |
| `ref_coverage` | REAL | Fraction of reference covered by non-gap aligned bases |
| `read_to_ref_ratio` | REAL | `read_len / ref_len` — detects gross size mismatch |

#### Segmented Metrics (reference divided into equal thirds)

The reference is split at `ref_len/3` and `2*ref_len/3`. The CIGAR is walked to map alignment columns to reference coordinates, then metrics are computed within each segment.

| Metric | Segments | Description |
|--------|----------|-------------|
| `seg_identity` | 5', mid, 3' | Matches / segment length |
| `seg_contiguity` | 5', mid, 3' | Longest run of consecutive matches (no indels) |
| `seg_gap_count` | 5', mid, 3' | Number of gaps within segment |

#### Indel Detection Metrics

| Metric | Type | Description |
|--------|------|-------------|
| `max_ins` | INT | Length of largest insertion |
| `max_del` | INT | Length of largest deletion |
| `n_sig_indels` | INT | Count of indels >= 5bp |
| `indel_positions` | LIST | `(type, size, ref_position)` for indels >= 3bp |

#### 5' Alignment Quality

| Metric | Type | Description |
|--------|------|-------------|
| `five_prime_offset` | INT | Reference bases at 5' end that are gapped (read doesn't cover) |
| `five_prime_identity_20` | REAL | Identity over first 20bp of reference |
| `five_prime_first_match` | INT | Position of first match — how far before alignment engages |

#### 3' Alignment Quality

| Metric | Type | Description |
|--------|------|-------------|
| `three_prime_offset` | INT | Reference bases at 3' end that are gapped |
| `three_prime_identity_20` | REAL | Identity over last 20bp of reference |

#### Classification Columns

| Metric | Type | Description |
|--------|------|-------------|
| `rank` | INT | Rank of this reference for this read (1 = best by `ned`) |
| `margin` | REAL | `ned(rank=2) - ned(rank=1)` — separation between best and runner-up |

### Output Files

**All metrics (one row per read x reference pair):**

```
Output/alignments.tsv
```

Columns: `read_id, ref_id, read_len, ref_len, ed, ned, identity, ref_coverage, read_to_ref_ratio, seg5_identity, seg5_contiguity, segM_identity, segM_contiguity, seg3_identity, seg3_contiguity, max_ins, max_del, n_sig_indels, five_prime_offset, five_prime_identity_20, three_prime_offset, three_prime_identity_20, rank, margin`

**Classification summary (one row per read):**

```
Output/classification.tsv
```

Columns: `read_id, read_len, assigned_ref, best_ned, second_ned, margin, confidence_flag`

Where `confidence_flag`:
- `HIGH` — margin > 0.1
- `LOW` — margin <= 0.1
- `POOR` — best ned > 0.5

These thresholds are initial defaults. The plots below are designed to refine them.

---

## Stage 5b: QC — Diagnostic Visualizations

All plots: matplotlib, `seaborn-v0_8-whitegrid` style, 8x5 inches, 300 DPI, saved as PNG + PDF. KDE uses Savitzky-Golay smoothing (window=51, polyorder=3).

### Plot 1: Classification Landscape Heatmap

`n_reads x n_refs` heatmap (subsampled to ~5000 reads if larger). Color = `ned`. Reads sorted by assigned target, then by `ned` within assignment. A clean block-diagonal pattern indicates good separability.

```
Output/qc/classification_heatmap.{png,pdf}
```

### Plot 2: Score Margin Distribution

KDE of `margin` values for all reads. Annotated threshold candidates. Bimodal distribution (peak near 0 = ambiguous, peak >> 0 = confident) confirms targets are distinguishable.

```
Output/qc/margin_kde.{png,pdf}
```

### Plot 3: NED Distribution by Assigned Target

**Top subplot:** One KDE curve per target, overlaid. X = `ned` against assigned reference.
**Bottom subplot:** Same but showing `ned` against the second-best reference.

The gap between these subplots is the classification margin.

```
Output/qc/ned_by_target.{png,pdf}
```

### Plot 4: Segmented Identity Violin Plots

For each assigned target, three violins: `seg5_identity`, `segM_identity`, `seg3_identity`. Reveals if a specific segment (e.g., 5' flank) is systematically worse, indicating adapter/barcode issues or truncation.

```
Output/qc/segmented_identity_violins.{png,pdf}
```

### Plot 5: Read Length vs NED Scatter

Scatter colored by assigned target. X = read length, Y = `ned`. Horizontal dashed lines at NED thresholds. Vertical dashed lines at expected target lengths. Outliers are misclassifications or truncated reads.

```
Output/qc/length_vs_ned_scatter.{png,pdf}
```

### Plot 6: 5' Alignment Quality

**Top panel:** KDE of `five_prime_offset` per target
**Bottom panel:** KDE of `five_prime_identity_20` per target

Poor 5' alignment suggests adapter ligation issues or systematic truncation.

```
Output/qc/five_prime_quality.{png,pdf}
```

### Plot 7: Significant Indel Profile

**Top panel:** KDE of `max_ins` and `max_del` per target (overlaid)
**Bottom panel:** Histogram of indel positions along the reference (binned by segment) — reveals hotspots

```
Output/qc/indel_profile.{png,pdf}
```

### Plot 8: Misclassification Threshold Sweep

For NED thresholds from 0.1 to 0.5 (step 0.02):
- **Left Y-axis:** Fraction of reads classified (above threshold = unassigned)
- **Right Y-axis:** Estimated misclassification rate (reads where margin < cutoff)

Key plot for picking thresholds — look for the elbow.

```
Output/qc/threshold_sweep.{png,pdf}
```

### Plot 9: Per-Read Multi-Target Profiles (sampled)

Random sample of ~50 reads. Small-multiples grid: each read's `ned` against all targets as a bar chart. Reveals whether reads have a single clear best match or are ambiguous.

```
Output/qc/per_read_profiles_sample.{png,pdf}
```

### Plot 10: End Reason Facets

Plots 2-5 regenerated with faceting by end reason (`signal_positive` vs `data_service_unblock_mux_change`). Adaptive sampling rejects may have systematically different alignment quality.

```
Output/qc/ned_by_target_by_endreason.{png,pdf}
Output/qc/length_vs_ned_by_endreason.{png,pdf}
Output/qc/segmented_identity_by_endreason.{png,pdf}
Output/qc/margin_by_endreason.{png,pdf}
```

---

## Manifest (Checkpoint / Provenance)

`Output/prepare_manifest.json` records all state for resume and traceability:

```json
{
  "exp_id": "FAL12345_20251228_IF",
  "expdir": "/mnt/d/12282025_IF_DoubleBC_SMA_seq_no_trim",
  "timestamp": "2026-02-14T10:30:00",
  "runs": [
    {
      "run_id": "20251228_2219_MD-100098_FBD69411_34fa833d",
      "flow_cell_id": "FBD69411",
      "device_id": "MD-100098",
      "pod5_count": 12,
      "pod5_bytes": 54000000000,
      "bam_count": 1500,
      "basecall_model": "dna_r10.4.1_e8.2_400bps_sup@v5.2.0",
      "action": "REUSE"
    }
  ],
  "merge_groups": [["run1", "run2"]],
  "output_bam": "Output/merged.bam",
  "output_bam_reads": 45832,
  "database_path": "Output/SMA_FAL12345_20251228_IF.db",
  "ref_path": "targets.fa",
  "ref_targets": ["V0-4.2_BC02_fwd", "V0-4.4_BC09_fwd", "V0-4.2_BC07_rev", "V0-4.4_BC04_rev"],
  "alignment_reads_processed": 45832,
  "classification_high": 42100,
  "classification_low": 2800,
  "classification_poor": 932,
  "stages_completed": ["discover", "plan", "merge", "init", "align", "qc"],
  "qc_plots": ["classification_heatmap.png", "margin_kde.png", "..."]
}
```

If `prepare.py` is re-run and finds a valid manifest, it offers to resume from the last completed stage.

---

## Output Layout

```
Output/
  merged.bam                    Sorted, indexed merged BAM
  merged.bam.bai                BAM index
  summary.tsv                   read_id -> end_reason from POD5
  SMA_{exp_id}.db               Initialized SMA-seq database
  prepare_manifest.json         Checkpoint + provenance
  alignments.tsv                All read x ref alignment metrics
  classification.tsv            Best assignment per read
  pod5/                         Organized POD5 symlinks
    run1_.../
    run2_.../
  qc/                           Publication-quality diagnostic plots
    classification_heatmap.{png,pdf}
    margin_kde.{png,pdf}
    ned_by_target.{png,pdf}
    segmented_identity_violins.{png,pdf}
    length_vs_ned_scatter.{png,pdf}
    five_prime_quality.{png,pdf}
    indel_profile.{png,pdf}
    threshold_sweep.{png,pdf}
    per_read_profiles_sample.{png,pdf}
    *_by_endreason.{png,pdf}
```

---

## Updated Workflow

```
# Step 0: Prepare (interactive — new)
python3 bin/prepare.py \
  -d /mnt/d/EXPERIMENT_DIR \
  -e EXP_ID \
  -r targets.fa \
  -o Output/

# Step 1: Ingest (batch — existing)
python3 bin/ingest.py \
  -e EXP_ID \
  -b Output/merged.bam \
  -s Output/summary.tsv \
  -r Input/target.fa \
  -d Output/SMA_EXP_ID.db \
  -o Output/tagged.bam
```

Steps 1-3 of the old workflow (`mkdb`, `inputInit`, `extractMeta`) are subsumed by `prepare.py`.

---

## Dependencies

All already available in `env/env.yml`:

| Package | Purpose |
|---------|---------|
| `pysam` | BAM I/O, merge, sort, index |
| `edlib` | Global alignment (NW mode with CIGAR) |
| `matplotlib` | All plots |
| `scipy` | Savitzky-Golay peak detection |
| `numpy` | Array operations for KDE and metrics |
| `pod5` | CLI for end-reason extraction |
| `sqlite3` | Database (stdlib) |

No new dependencies required.

---

## Design Decisions

1. **Approach A (Staged Pipeline)** chosen over Orchestrator+Workers (B) and Library+CLI (C) for simplicity. Single entry point, linear flow, resume via manifest.

2. **Forward-strand-only alignment** — each reference sequence represents a specific strand. Reads match 5'->3' or not at all. No reverse complement alignment.

3. **Global alignment (NW)** rather than semi-global — forces end-to-end comparison. A read that is grossly the wrong size still gets aligned, but `ned` and `ref_coverage` will reflect the mismatch.

4. **Segmented metrics in thirds** — simple and interpretable. The 5'/mid/3' split catches adapter issues (5'), target core quality (mid), and 3' truncation (3').

5. **Threshold sweep plot** as the primary classification tool — rather than hardcoding thresholds, the user examines the elbow plot and decides.

6. **Symlink POD5s** — avoids copying hundreds of GB. Original data stays read-only.

7. **Separate `align.py` and `qc.py`** — alignment is CPU-intensive and produces a reusable TSV. Plots are fast and may be re-run with adjusted parameters without re-computing alignments.
