# Testing the ONT-SMA-seq Pipeline

This document describes how to test the ONT-SMA-seq pipeline.

## Quick Test

Run the comprehensive test suite:

```bash
python test_pipeline.py
```

This will:
1. Create test data using `test_setup.py`
2. Test database creation with `mkdb.py`
3. Test input standardization with `inputInit.py`
4. Test read processing with `ingest.py`
5. Verify database integrity and data quality
6. Clean up all test artifacts

## Test Coverage

The test suite (`test_pipeline.py`) includes **17 tests** covering:

### mkdb.py (3 tests)
- ✓ Database file creation
- ✓ Database schema (all tables present)
- ✓ Mods table population

### inputInit.py (4 tests)
- ✓ Input directory creation
- ✓ BAM symlink creation
- ✓ Pod5 directory symlink creation
- ✓ Reference FASTA symlink creation

### ingest.py (7 tests)
- ✓ Output directory creation
- ✓ Output BAM file creation
- ✓ Database read insertion
- ✓ Reference sequence population
- ✓ Reference matching logic
- ✓ Quality metric calculation (q_bc)
- ✓ Levenshtein distance calculation

### Integration (3 tests)
- ✓ Complete read data in database
- ✓ Unique ID format validation
- ✓ Foreign key relationship integrity

## Manual Testing

### 1. Create Test Data

```bash
python test_setup.py
```

This creates:
- `test_data/reference.fa` - Reference FASTA with 2 sequences (200bp and 500bp)
- `test_data/TEST001_h_v5.2.0_1_6mA.bam` - Test BAM with 50 reads
- `test_data/pod5/` - Empty Pod5 directory (end reasons will be 'unknown')

### 2. Run Pipeline

```bash
# Initialize database
python mkdb.py TEST001

# Standardize inputs
python inputInit.py --bam test_data/TEST001_h_v5.2.0_1_6mA.bam \
                    --pod5 test_data/pod5 \
                    --ref test_data/reference.fa

# Process reads
python ingest.py TEST001
```

### 3. Verify Results

Check database contents:

```bash
# Count total reads
sqlite3 SMA_TEST001.db "SELECT COUNT(*) FROM Reads;"

# View reference distribution
sqlite3 SMA_TEST001.db "SELECT refseq_id, COUNT(*) FROM Reads GROUP BY refseq_id;"

# Check quality metrics
sqlite3 SMA_TEST001.db "SELECT AVG(q_bc), AVG(q_ld) FROM Reads WHERE refseq_id IS NOT NULL;"
```

Check output BAM:

```bash
samtools view -c Output/TEST001.bam
```

### 4. Clean Up

```bash
rm -rf SMA_TEST001.db Input/ Output/ test_data/
```

## Test Results Example

```
======================================================================
ONT-SMA-seq Pipeline Test Suite
======================================================================

--- Testing mkdb.py ---
✓ mkdb: Database file created
✓ mkdb: All required tables created
✓ mkdb: Mods table populated

--- Testing inputInit.py ---
✓ inputInit: Input directory created
✓ inputInit: Created TEST001_h_v5.2.0_1_6mA.bam
✓ inputInit: Created TEST001_pod5
✓ inputInit: Created TEST001.fa

--- Testing ingest.py ---
✓ ingest: Output directory created
✓ ingest: Output BAM created
✓ ingest: Reads inserted into database
✓ ingest: Reference sequences populated
✓ ingest: Reference matching performed
✓ ingest: Quality metrics (q_bc) calculated
✓ ingest: Levenshtein distance calculated for matched reads

--- Testing Pipeline Integration ---
✓ integration: Database contains complete read data
✓ integration: Unique ID format correct
✓ integration: Modification foreign keys valid

======================================================================
Test Results: 17/17 passed
======================================================================
```

## Dependencies

Testing requires:
- Python 3.6+
- pysam
- edlib
- pod5 (optional - pipeline handles missing pod5 gracefully)

Install with:
```bash
pip install pysam edlib pod5
```
