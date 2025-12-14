# ONT-SMA-seq Testing Documentation

## Test Coverage Summary

This pipeline has been fully tested with a comprehensive test suite consisting of 25 tests across 4 test modules.

### Test Results

```
✓ 25 tests passing
✓ 0 tests failing
✓ 100% success rate
```

## Test Modules

### 1. Unit Tests: `test_mkdb.py` (4 tests)

Tests for database initialization script:

- ✓ `test_create_tables` - Verifies all database tables are created correctly
- ✓ `test_populate_mods_table` - Verifies modification bitflags are populated
- ✓ `test_main_creates_database` - Tests end-to-end database creation
- ✓ `test_main_overwrites_existing_database` - Tests database overwrite functionality

### 2. Unit Tests: `test_inputInit.py` (10 tests)

Tests for input standardization script:

- ✓ `test_parse_bam_filename_valid` - Tests parsing of valid BAM filename
- ✓ `test_parse_bam_filename_complex` - Tests parsing with complex experiment IDs
- ✓ `test_parse_bam_filename_no_mods` - Tests parsing with no modifications
- ✓ `test_parse_bam_filename_invalid` - Tests error handling for invalid filenames
- ✓ `test_create_symlink` - Tests creating symbolic links to files
- ✓ `test_create_symlink_directory` - Tests creating symbolic links to directories
- ✓ `test_create_symlink_source_not_exists` - Tests error handling for missing sources
- ✓ `test_create_symlink_target_exists` - Tests error handling for existing targets
- ✓ `test_create_symlink_force_overwrite` - Tests force overwrite functionality
- ✓ `test_main_creates_symlinks` - Tests end-to-end symlink creation

### 3. Unit Tests: `test_ingest.py` (9 tests)

Tests for main processing script:

- ✓ `test_calculate_q_bc` - Tests basecall quality calculation
- ✓ `test_calculate_levenshtein` - Tests Levenshtein distance calculation
- ✓ `test_calculate_q_ld` - Tests Levenshtein quality calculation
- ✓ `test_match_reference` - Tests reference sequence matching
- ✓ `test_generate_uniq_id` - Tests unique ID generation
- ✓ `test_parse_bam_filename` - Tests BAM filename parsing
- ✓ `test_get_mod_bitflag` - Tests modification bitflag conversion
- ✓ `test_parse_reference_fasta` - Tests reference FASTA parsing
- ✓ `test_parse_reference_fasta_wrong_count` - Tests error handling for wrong sequence count

### 4. Integration Tests: `test_integration.py` (2 tests)

Tests for complete pipeline:

- ✓ `test_pipeline_integration` - Tests full pipeline workflow
- ✓ `test_database_schema` - Verifies database schema matches specification

## Running Tests

### Run all tests:
```bash
pytest tests/
```

### Run with verbose output:
```bash
pytest tests/ -v
```

### Run specific test module:
```bash
pytest tests/test_mkdb.py -v
pytest tests/test_inputInit.py -v
pytest tests/test_ingest.py -v
pytest tests/test_integration.py -v
```

### Run with coverage:
```bash
pytest tests/ --cov=. --cov-report=html
```

## Manual Verification

The scripts have been manually tested and verified to work correctly:

### mkdb.py
```bash
$ python mkdb.py demo_exp
Created database schema in SMA_demo_exp.db
Populated Mods table with modification bitflags
Successfully initialized database: SMA_demo_exp.db

$ sqlite3 SMA_demo_exp.db "SELECT * FROM Mods;"
0|non
1|6mA
2|5mCG_5hmCG
...
```

### inputInit.py
```bash
$ python inputInit.py demo_exp_s_v5.2.0_1_6mA.bam pod5/ ref.fa
Parsed metadata from BAM filename:
  Experiment ID: demo_exp
  Model Tier: s
  Model Version: 5.2.0
  Trim: 1
  Modifications: 6mA

Created directory: Input
Created symlink: Input/demo_exp_s_v5.2.0_1_6mA.bam -> demo_exp_s_v5.2.0_1_6mA.bam
Created symlink: Input/demo_exp_pod5 -> pod5/
Created symlink: Input/demo_exp.fa -> ref.fa

Successfully standardized inputs in Input
```

## Test Quality Metrics

- **Code Coverage**: Core functionality fully tested
- **Edge Cases**: Error handling and boundary conditions tested
- **Integration**: Full pipeline workflow verified
- **Documentation**: All tests have clear docstrings
- **Maintainability**: Tests are independent and reproducible

## Dependencies

Testing requires the following packages:
- pytest >= 7.0.0
- pysam >= 0.21.0
- edlib >= 1.3.9
- pod5 >= 0.2.0

Install with:
```bash
pip install -r requirements.txt
pip install -r requirements-dev.txt
```
