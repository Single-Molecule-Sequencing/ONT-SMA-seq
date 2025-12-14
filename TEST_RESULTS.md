# ONT-SMA-seq Pipeline Test Results

## Test Execution Summary

**Date:** December 14, 2024  
**Pipeline Version:** Updated from claude/explain-code-functionality branch  
**Test Suite:** test_pipeline.py  
**Total Tests:** 17  
**Status:** ✅ ALL TESTS PASSING

---

## Test Results

### ✅ mkdb.py - Database Initialization (3/3 passed)

```
✓ mkdb: Database file created
✓ mkdb: All required tables created
✓ mkdb: Mods table populated
```

**Validates:**
- SQLite database creation
- Schema structure (Reads, Mods, Exp, Refseq tables)
- Modification bitflag lookup table population

---

### ✅ inputInit.py - Input Standardization (4/4 passed)

```
✓ inputInit: Input directory created
✓ inputInit: Created TEST001_h_v5.2.0_1_6mA.bam
✓ inputInit: Created TEST001_pod5
✓ inputInit: Created TEST001.fa
```

**Validates:**
- Input directory structure creation
- Symbolic link creation for BAM file
- Symbolic link creation for Pod5 directory
- Symbolic link creation for reference FASTA

---

### ✅ ingest.py - Read Processing (7/7 passed)

```
✓ ingest: Output directory created
✓ ingest: Output BAM created
✓ ingest: Reads inserted into database
✓ ingest: Reference sequences populated
✓ ingest: Reference matching performed
✓ ingest: Quality metrics (q_bc) calculated
✓ ingest: Levenshtein distance calculated for matched reads
```

**Validates:**
- Output directory creation
- Tagged BAM file generation
- Database population with read data
- Reference sequence insertion
- Length-based reference matching
- Basecall quality metric calculation
- Edit distance calculation for matched reads

---

### ✅ Integration Testing (3/3 passed)

```
✓ integration: Database contains complete read data
✓ integration: Unique ID format correct
✓ integration: Modification foreign keys valid
```

**Validates:**
- End-to-end pipeline execution
- Unique ID generation format (exp_id + model + version + trim + mods + hash)
- Database foreign key relationships
- Data integrity across all tables

---

## Sample Test Data

**Test Experiment:** TEST001  
**Model:** h_v5.2.0  
**Trim:** 1  
**Modifications:** 6mA  

**Reference Sequences:**
- short_amplicon: 200 bp (matching range: 50-350 bp)
- long_amplicon: 500 bp (matching range: 350-650 bp)

**Test Reads:** 50 total
- Matched to short_amplicon: ~14 reads
- Matched to long_amplicon: ~22 reads  
- Unmatched (out of range): ~14 reads

---

## Pipeline Verification

### Database Validation
- ✅ All reads inserted (50/50)
- ✅ Reference sequences populated (2/2)
- ✅ Quality metrics calculated for all reads
- ✅ Edit distance calculated for matched reads only
- ✅ Foreign key constraints maintained

### File Output Validation
- ✅ Output BAM created with ER tags
- ✅ Database schema matches specification
- ✅ Symlink structure correct

---

## How to Run Tests

```bash
# Quick test
python test_pipeline.py

# Manual test workflow
python test_setup.py
python mkdb.py TEST001
python inputInit.py --bam test_data/TEST001_h_v5.2.0_1_6mA.bam \
                    --pod5 test_data/pod5 \
                    --ref test_data/reference.fa
python ingest.py TEST001
```

---

## Dependencies

All required dependencies are installed and working:
- ✅ Python 3.x
- ✅ pysam (BAM I/O)
- ✅ edlib (Levenshtein distance)
- ✅ pod5 (Pod5 metadata - optional, gracefully handles absence)
- ✅ sqlite3 (built-in)

---

## Conclusion

**The ONT-SMA-seq pipeline is fully functional and tested.**

All components work correctly:
- Database initialization ✅
- Input standardization ✅
- Read processing and metric calculation ✅
- Data integrity and integration ✅

The pipeline successfully processes ONT sequencing data, calculates quality metrics, matches reads to references, and stores all data in a well-structured SQLite database.
