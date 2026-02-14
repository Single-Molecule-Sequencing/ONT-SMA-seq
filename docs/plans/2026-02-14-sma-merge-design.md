# SMA-seq Experiment Merge Tool Design

**Date**: 2026-02-14
**Status**: Approved

## Goal

Build `sma-merge`, a CLI tool that discovers MinKNOW run directories for an SMA-seq experiment, validates basecalling consistency, and produces clean merged or subsampled BAM files with full metadata tagging (including end_reason from POD5).

## Two Modes

1. **Full merge**: Discover -> Validate -> Merge POD5s -> Re-basecall with dorado (no trim, no demux) -> Tag end_reason -> Single clean BAM + merged POD5
2. **Subsample**: Discover -> Validate -> Random sample N reads from POD5 -> Basecall -> Tag -> Quick analysis BAM + subsampled POD5

## Architecture

CLI tool at `bin/sma_merge/` with modules:

```
bin/sma_merge/
    __init__.py
    cli.py          # argparse CLI entry point
    discover.py     # MinKNOW run discovery + metadata extraction
    validate.py     # Basecall consistency + run grouping
    merge.py        # POD5 merge + dorado basecalling
    subsample.py    # Random subsampling from POD5
    tag.py          # Post-basecall end_reason tagging
```

### CLI Interface

```bash
sma-merge discover <experiment_path>
sma-merge merge <experiment_path> -o output/
sma-merge subsample <experiment_path> -n 5000 -o output/ [--seed 42]
```

## Discovery & Validation

### Discovery

Scans a path for MinKNOW run directories by looking for:
```
<path>/
  <sample_id>/
    <YYYYMMDD_HHMM_device_flowcell_runid>/
      pod5_pass/
      final_summary_*.txt
```

Extracts from POD5 `run_info.context_tags` and `run_info.tracking_id`:
- `flow_cell_id` (e.g., FBD66244)
- `device_id` (e.g., MD-101527)
- `protocol_group_id` (e.g., 12302025_IF_DoubleBC_SMA_seq_one_nick)
- `basecall_model_simplex` (e.g., dna_r10.4.1_e8.2_400bps_hac@v5.2.0)
- `sample_id`, `run_id`

### Validation

Groups runs and checks:
1. Same basecall model across all runs in a group
2. Same flowcell ID = same physical run (merge OK)
3. Different flowcell/protocol_group = different experiment (don't merge)
4. Warns if sample_id is empty

Output: summary table + interactive confirmation.

## Full Merge

1. Merge all POD5 files with `pod5 merge`
2. Run dorado with auto-detected model, `--no-trim`, no `--kit-name`, `--emit-moves`
3. Post-process BAM: add `er` tag (end_reason) from POD5

Output: `{experiment}_merged.pod5` + `{experiment}_merged.bam`

## Subsample

1. Scan POD5 files, collect read IDs + metadata
2. Random sample N read IDs (seeded)
3. `pod5 filter` to extract subsampled POD5
4. Basecall with dorado (same flags)
5. Post-process: add `er` and `sl` tags from POD5

Output: `{experiment}_sub{N}.pod5` + `{experiment}_sub{N}.bam`

## End Reason Tagging

Post-basecalling step since dorado doesn't include end_reason:
1. Build lookup: `{read_id: (end_reason, num_samples)}` from POD5
2. Read dorado output BAM with pysam
3. Add `er` tag (end_reason string) and `sl` tag (raw signal length)
4. Write tagged BAM

Standard dorado tags already present: qs, du, ns, ts, ch, mx, rn, st, sm, sd, sv, dx, mv, MM, ML, RG.

## Key Technical Details

- **POD5 end_reason enum**: SIGNAL_POSITIVE, DATA_SERVICE_UNBLOCK_MUX_CHANGE, UNBLOCK_MUX_CHANGE, SIGNAL_NEGATIVE, MUX_CHANGE
- **Dorado model auto-detection**: Read from POD5 `context_tags['basecall_model_simplex']`
- **Raw signal length**: `num_samples` from POD5 read, also `ns` tag in dorado BAM
- **Duration**: `num_samples / sample_rate` (sample_rate from POD5 run_info)
- **Dependencies**: pod5, pysam, dorado (external binary), numpy

## Integration

The output BAM can feed into:
1. `ingest_bams.py` for barcode classification + calibration viz database
2. Calibration viz dashboard for interactive exploration
3. SMA-seq analysis pipeline
