# SMA-seq Experiment Preparation Tool — Design

## Problem

SMA-seq experiments produce raw data across multiple MinKNOW sequencing runs that may need merging. The existing MinKNOW outputs include demuxed BAMs and quality-filtered reads, but SMA-seq analysis requires a single untrimmed, undemuxed BAM file with all reads from an experiment. A tool is needed to:

1. Detect and merge runs from the same experiment (stop/restart on same flow cell)
2. Ensure consistent basecalling (same model, version, parameters)
3. Produce a single sorted/indexed BAM with all reads and metadata tags
4. Generate initial QC to validate library preparation before classification

## Architecture: Pipeline Stages

Three independent tools coordinated by CLI invocation:

```
sma_scan.py    →  manifest.json  →  sma_basecall.py  →  merged.bam  →  sma_init.py  →  DB + QC report
```

Each stage is independently testable and resumable. The manifest JSON between stages is human-inspectable and editable.

---

## Stage 1: `sma_scan.py` — Experiment Scanner

### Purpose

Scan an experiment directory, discover MinKNOW runs, detect merge candidates, validate consistency, output structured manifest.

### CLI

```
sma_scan.py <experiment_dir> -o manifest.json [--interactive] [--max-gap 24h]
```

### Discovery Logic

1. Recursively find MinKNOW run folders (identified by `final_summary*.txt`, `pod5*/`, or `pod5_pass/`)
2. For each run, extract metadata from multiple sources:

| Source | Key Fields |
|--------|-----------|
| final_summary.txt | flow_cell_id, instrument, protocol, protocol_run_id, started/stopped, file counts |
| sample_sheet.csv | experiment_id, kit, flow_cell_product_code, barcode→alias mapping |
| report.json | Software versions (basecaller, MinKNOW, bream), device type, channel count |
| sequencing_summary.txt | path recorded (not parsed at scan time — too large) |
| BAM headers | @PG basecaller model/version, @RG sample/run info |

3. Group runs by `flow_cell_id`
4. Check merge eligibility within each group:
   - Same protocol, same kit, same experiment_id, time gap < `--max-gap`
   - Flag potential wash events (gap > threshold, different sample_id, different experiment_id)

### Merge Decision Logic

- **Auto-merge** (default): Same flow cell + protocol + kit + experiment_id, gap < max_gap
- **Flag for review**: Gap > max_gap or mismatched sample/experiment
- **`--interactive`**: Always prompt before finalizing merge decisions

### Basecalling Consistency Validation

For `--from-bam` workflows: parse BAM headers from all files, verify same @PG basecaller line. Flag any mismatches in model, version, or trimming parameters.

### Output: `manifest.json`

```json
{
  "experiment_id": "12282025_IF_DoubleBC_SMA_seq_no_trim",
  "experiment_dir": "/mnt/d/12282025_IF_DoubleBC_SMA_seq_no_trim",
  "scan_timestamp": "2026-02-14T12:00:00",
  "kit": "SQK-NBD114-24",
  "flow_cell_product": "FLO-MIN114",
  "software": {
    "minknow": "25.09.16",
    "basecaller": "7.11.0+5d1db4a52",
    "bream": "8.8.3"
  },
  "barcode_aliases": {
    "barcode02": "V04_2",
    "barcode04": "V04_4",
    "barcode07": "V04_2_2",
    "barcode09": "V04_4_2"
  },
  "flow_cells": [
    {
      "flow_cell_id": "FBD69411",
      "instrument": "MD-100098",
      "merge_decision": "auto",
      "runs": [
        {
          "run_id": "34fa833d",
          "path": "no_sample_id/20251228_2219_.../",
          "started": "2025-12-28T22:19:00",
          "stopped": null,
          "pod5_dir": "pod5/",
          "pod5_count": 12,
          "bam_dir": "bam_pass/",
          "bam_demuxed": true,
          "bam_count": 5,
          "sequencing_summary": null,
          "final_summary": null
        },
        {
          "run_id": "5b5c57a9",
          "path": "no_sample_id/20251229_1055_.../",
          "started": "2025-12-29T10:56:12",
          "stopped": "2026-01-01T10:56:12",
          "pod5_dir": "pod5/",
          "pod5_count": 72,
          "bam_dir": "bam_pass/",
          "bam_demuxed": true,
          "bam_count": 2165,
          "sequencing_summary": "sequencing_summary_FBD69411_5b5c57a9_5452ebd1.txt",
          "final_summary": "final_summary_FBD69411_5b5c57a9_5452ebd1.txt"
        }
      ],
      "time_gap_hours": 12.6,
      "warnings": []
    }
  ],
  "action": {
    "mode": "from_pod5",
    "needs_rebasecall": true,
    "pod5_sources": [
      "no_sample_id/20251228_.../pod5/",
      "no_sample_id/20251229_.../pod5/"
    ],
    "sequencing_summary_paths": [
      "no_sample_id/20251229_.../sequencing_summary_FBD69411_5b5c57a9_5452ebd1.txt"
    ],
    "total_pod5_files": 84,
    "total_reads_estimated": null
  }
}
```

The `needs_rebasecall` flag is true when:
- BAMs are demuxed (subdirectories in bam_pass/)
- BAM headers show trimming was applied
- No BAMs exist
- User requested `--from-pod5`

---

## Stage 2: `sma_basecall.py` — Basecalling & BAM Merge

### Purpose

Take a manifest and produce a single sorted/indexed BAM with all reads — untrimmed, undemuxed, with all tags.

### CLI

```
sma_basecall.py manifest.json --from-pod5 -o output/ [--model sup] [--subsample N] [--hpc]
sma_basecall.py manifest.json --from-bam -o output/
```

### Mode 1: `--from-pod5` (re-basecall)

1. Collect all pod5 paths from manifest
2. Run Dorado: `dorado basecaller <model> <pod5_dir> --no-trim --emit-moves --emit-sam`
3. If multiple pod5 directories: basecall each, then merge
4. Sort and index with samtools
5. Validate: read count matches pod5 read count

### Mode 2: `--from-bam` (merge existing)

1. Collect all BAM files from all runs (including demuxed subdirectories)
2. Validate consistency (same @PG basecaller line)
3. Merge: `samtools merge`
4. Sort and index
5. Record original MinKNOW bin assignments as metadata

### Subsampling: `--subsample N`

For quick QC before full basecalling:

1. Read pod5 file headers to get total read count and read_id list
2. Randomly sample N read_ids (default: 10,000)
3. Extract sampled reads using `pod5 filter` into temporary pod5
4. Basecall subsample with same parameters
5. Tag reads with metadata from sequencing_summary.txt:
   - end_reason (ER tag)
   - signal_duration_s
   - mean_qscore
6. Output: `{experiment_id}_subsample_{N}.bam`

### HPC Dispatch: `--hpc`

1. Generate SLURM script for Dorado on Great Lakes
2. Submit via `sbatch`
3. Provide command to check status / poll for completion
4. Transfer results back when done
5. Continue with merge/sort/index locally

### Sequencing Summary Merge

Merge sequencing_summary.txt files from all runs into one combined TSV, preserving the `run_id` column. This provides per-read end_reason, duration, and mean_qscore without pod5 re-extraction.

### Output

```
output/
  {experiment_id}_merged.bam           # Single sorted, indexed BAM
  {experiment_id}_merged.bam.bai       # BAM index
  {experiment_id}_summary.tsv          # Merged sequencing summary
  {experiment_id}_basecall.json        # Provenance:
    {
      "mode": "from_pod5",
      "model": "dna_r10.4.1_e8.2_400bps_sup@v5.0.0",
      "dorado_version": "0.8.4",
      "parameters": {"no_trim": true, "emit_moves": true},
      "pod5_sources": [...],
      "total_reads": 171660,
      "total_bases": 48923456,
      "mean_qscore": 13.2,
      "subsample": null,
      "timestamp": "2026-02-14T..."
    }
```

### Validation Checks

- Read count in output BAM matches pod5 read count (or subsample size)
- No duplicate read_ids
- All reads have expected tags
- BAM header records basecaller model, version, parameters

---

## Stage 3: `sma_init.py` — Database Creation & QC Report

### Purpose

Create SMA-seq database from merged BAM, compute pre-classification QC metrics, generate interactive HTML report.

### CLI

```
sma_init.py \
  --bam merged.bam \
  --sequencing-summary merged_summary.tsv \
  --manifest manifest.json \
  -o output/ \
  [--construct construct.toml]
```

### Database Creation

Uses existing `mkdb.py` schema. For each read in BAM:
- read_id, readseq, readlen
- q_bc (from BAM qual scores, probability-averaged)
- model_tier, model_ver (from BAM header @PG)
- ER (end_reason from sequencing summary)
- signal_duration_s (from sequencing summary `duration` column)
- mean_qscore (from sequencing summary)

**No barcode classification at this stage.** Raw reads only. Classification uses the existing `ingest.py` with construct TOML.

### QC Metrics Computed

| Metric | Description | Expected Pattern |
|--------|-------------|-----------------|
| Read length distribution | KDE of all read lengths | Multimodal: peaks at each target size |
| Signal duration distribution | Raw signal duration histogram + scatter vs read length | Linear (translocation speed ~400bp/s) |
| End reason breakdown | Counts by end_reason category | Depends on adaptive sampling config |
| Q-score distribution | Basecaller Q-score histogram | Unimodal around Q12-Q15 for sup |
| Per-channel activity | Reads per channel (512 channels) | Uniform, identifies dead channels |
| Throughput over time | Reads/minute over run duration | Stable or graceful degradation |

### HTML QC Report

Interactive report (same framework as classification_report.html):

| Tab | Content |
|-----|---------|
| Overview | Experiment metadata from manifest, total reads/bases, mean Q, run duration, merge info |
| Read Lengths | KDE with vertical lines at expected sizes (if construct.toml provided), per-run overlay if merged |
| Signal & Duration | Duration vs length scatter, duration histogram, translocation speed |
| End Reasons | Stacked bar by end_reason, pie chart, per-run comparison if merged |
| Quality | Q-score KDE, Q vs read length scatter |
| Channel Activity | 512-channel heatmap, throughput over time curve |

### Validation Section (when construct.toml provided)

- Overlay expected target sizes on read length distribution
- Flag missing or unexpected peaks
- Compare observed size ratios against expected

### Output

```
output/
  {experiment_id}.db                    # SMA-seq SQLite database
  {experiment_id}_qc_report.html        # Interactive QC report
  {experiment_id}_init.json             # Provenance record
```

---

## Full Workflow Example

```bash
# 1. Scan experiment directory
sma_scan.py /mnt/d/12282025_IF_DoubleBC_SMA_seq_no_trim/ \
  -o /mnt/d/12282025_IF_DoubleBC_SMA_seq_no_trim/manifest.json

# 2. Quick QC with 10K subsample
sma_basecall.py manifest.json --from-pod5 --subsample 10000 \
  -o /mnt/d/12282025_IF_DoubleBC_SMA_seq_no_trim/Output/

# 3. Generate QC report from subsample
sma_init.py \
  --bam Output/no_trim_subsample_10000.bam \
  --sequencing-summary Output/no_trim_summary.tsv \
  --manifest manifest.json \
  -o Output/

# 4. Review QC report, then full basecall
sma_basecall.py manifest.json --from-pod5 --model sup \
  -o /mnt/d/12282025_IF_DoubleBC_SMA_seq_no_trim/Output/

# 5. Create full database + QC
sma_init.py \
  --bam Output/no_trim_merged.bam \
  --sequencing-summary Output/no_trim_summary.tsv \
  --manifest manifest.json \
  --construct construct.toml \
  -o Output/

# 6. Classification (existing tool)
python bin/ingest.py -d Output/no_trim.db -b Output/no_trim_merged.bam \
  --construct construct.toml --sample-sheet sample_sheet.csv
```

---

## Experiments This Tool Will Process

| Experiment | Dir | Runs | Pod5 | Status |
|-----------|-----|------|------|--------|
| no_trim | D:\12282025_IF_DoubleBC_SMA_seq_no_trim | 2 (merge) | 84 | BAMs demuxed, needs rebasecall |
| one_nick | D:\12302025_IF_DoubleBC_SMA_seq_one_nick | 1 | 437 | BAMs demuxed, needs rebasecall |
| extended | D:\20260129_IF_GG_Part5_Odd\..._extended | 1 | ~40 | BAMs demuxed, needs rebasecall |
| extended_1 | D:\20260129_IF_GG_Part5_Odd\..._extended_1 | 1 | ~36 | BAMs demuxed, needs rebasecall |
| level1 | D:\20260209_IF_GG_Part5_Level1 | 1 | 51 | No BAMs, needs basecalling |
| level0_odd | D:\20260209_IF_GG_Part5_Level0_Odd | 1 | 5 | BAMs demuxed, needs rebasecall |

### Special Considerations

- **no_trim + one_nick**: Use non-standard barcodes (22bp, 23bp, 24bp mixed). Custom demux config needed.
- **extended + level0_odd**: Multiple shorter targets, standard ONT 24bp barcodes.
- **extended_1 + level1**: Single long target, standard ONT 24bp barcodes.
- **no_trim**: Two runs on same flow cell (FBD69411) ~12h apart — auto-merge candidate.

---

## Dependencies

```
dorado >= 0.8.0     # Basecalling
samtools >= 1.19    # BAM operations
pod5 >= 0.3         # Pod5 read/filter
pysam >= 0.22       # BAM parsing in Python
```

## File Location

All three tools go in `/tmp/ont-sma-seq/bin/`:
- `bin/sma_scan.py`
- `bin/sma_basecall.py`
- `bin/sma_init.py`

Tests in `tests/test_sma_scan.py`, `tests/test_sma_basecall.py`, `tests/test_sma_init.py`.
