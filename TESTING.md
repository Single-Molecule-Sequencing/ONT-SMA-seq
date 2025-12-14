# Testing the ONT-SMA-seq Pipeline

This document describes how to test the ONT-SMA-seq pipeline.

## Quick Test

Run the automated end-to-end test:

```bash
python test_e2e.py
```

This will:
1. Generate synthetic test data (20 reads, 2 reference sequences)
2. Run the complete pipeline workflow
3. Verify all outputs and database tables
4. Clean up test artifacts

**Expected output:** All tests pass with success message.

## What the Test Validates

### Database Creation (mkdb.py)
- ✓ Database file is created (`SMA_{exp_id}.db`)
- ✓ All required tables exist (Reads, Mods, Exp, Refseq)
- ✓ Mods table is pre-populated with modification bitflags
- ✓ Indices are created for common queries

### Input Standardization (inputInit.py)
- ✓ BAM filename is parsed correctly for metadata
- ✓ Symlinks are created in `Input/` directory
- ✓ Directory structure follows convention

### Read Processing (ingest.py)
- ✓ Reference sequences are parsed and stored
- ✓ Reads are matched to references by length
- ✓ Quality metrics are calculated (q_bc)
- ✓ Levenshtein distance is calculated for matched reads (ed, q_ld)
- ✓ Reads are tagged with End Reason (ER) in output BAM
- ✓ All reads are inserted into database
- ✓ NULL values for unmatched reads are handled correctly

## Test Data

The test creates:
- **Reference FASTA**: 2 sequences (short ~200bp, long ~500bp)
- **Input BAM**: 20 synthetic reads
  - ~8 reads matching short reference (±50bp)
  - ~7 reads matching long reference (±50bp)
  - ~5 reads outside both ranges (unmatched)
- **Pod5 directory**: Empty (end reasons will be 'unknown')

## Manual Testing

For manual testing with your own data:

### 1. Create Database
```bash
python mkdb.py MY_EXP
```

### 2. Prepare Input Files

Ensure your BAM file follows the naming convention:
```
{exp_id}_{model_tier}_v{model_ver}_{trim}_{mods}.bam
```

Example: `MY_EXP_h_v5.2.0_1_6mA.bam`

### 3. Standardize Inputs
```bash
python inputInit.py \
  --bam /path/to/MY_EXP_h_v5.2.0_1_6mA.bam \
  --pod5 /path/to/pod5_dir \
  --ref /path/to/reference.fa \
  --force
```

### 4. Process Reads
```bash
python ingest.py MY_EXP
```

### 5. Verify Outputs

Check database:
```bash
sqlite3 SMA_MY_EXP.db "SELECT COUNT(*) FROM Reads;"
```

Check output BAM:
```bash
samtools view Output/MY_EXP.bam | head
```

## Troubleshooting

### Missing Dependencies
If you see import errors, install dependencies:
```bash
pip install -r requirements.txt
```

### BAM Filename Format Error
Ensure your BAM filename follows the exact format:
```
{exp_id}_{s|h|f}_v{version}_{0|1}_{modifications}.bam
```

- Model tier: `s` (sup), `h` (hac), or `f` (fast)
- Trim: `0` (no trim) or `1` (trimmed)
- Modifications: `non`, `6mA`, `5mC`, or combinations like `6mA+5mC_5hmC`

### No Pod5 Files
The pipeline handles missing Pod5 files gracefully. End reasons will be set to "unknown" for all reads.

### Reference Matching
Reads are matched to references based on length ± 150bp tolerance. Reads outside this range will have NULL refseq_id, ed, and q_ld values.

## Test Coverage

Current test coverage:
- ✓ Full pipeline workflow (mkdb → inputInit → ingest)
- ✓ Database schema and population
- ✓ Metadata parsing from filenames
- ✓ Reference sequence matching
- ✓ Quality metric calculations
- ✓ BAM tagging with ER tags
- ✓ Handling of matched and unmatched reads

Not currently tested (requires real data):
- Pod5 end reason extraction (test uses empty Pod5 directory)
- Additional BAM tags extraction (channel, well, etc.)
- Large-scale performance (test uses only 20 reads)
