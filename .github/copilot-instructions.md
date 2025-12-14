# Copilot Instructions for ONT-SMA-seq

## Project Overview

ONT-SMA-seq is a Python-based bioinformatics pipeline for analyzing Oxford Nanopore Technology (ONT) sequencing data using the SMA-seq protocol. The pipeline processes nanopore sequencing reads, calculates quality metrics, and stores results in a SQLite database.

## Architecture

The project consists of three main Python scripts with distinct responsibilities:

1. **`mkdb.py`** - Database initialization
   - Creates SQLite database schema
   - Populates static lookup tables (Mods, Exp, Refseq)
   - Output: `SMA_{exp_id}.db`

2. **`inputInit.py`** - Input standardization
   - Parses BAM filename metadata
   - Creates symlinks for standardized input paths
   - Input structure: `{exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam`

3. **`ingest.py`** - Main processing pipeline
   - Parses reference FASTA and Pod5 metadata
   - Tags BAM files with end reason (ER)
   - Calculates quality metrics (q_bc, q_ld)
   - Computes Levenshtein distance for matched reads
   - Populates the Reads table

## Technology Stack

- **Language**: Python 3.x
- **Key Libraries**:
  - `pysam` - BAM file I/O
  - `pod5` - End Reason extraction from Pod5 files
  - `edlib` - Levenshtein distance calculation
  - `sqlite3` - Database operations
  - `argparse` - Command-line argument parsing

## Database Schema

The SQLite database contains four main tables:

### Reads Table
Primary storage for sequencing read data and metrics. Key columns:
- `uniq_id` (TEXT PK): Composite identifier in format `{exp_id}{tier}{ver}t{trim}m{mod_flag}_{read_hash}`
- `refseq_id` (TEXT FK): Reference sequence ID, NULL if length out of range
- `ed` (INT): Levenshtein distance, NULL if no RefSeq match
- `q_bc` (REAL): Probability-averaged basecall quality
- `q_ld` (REAL): Levenshtein quality, NULL if no RefSeq match

### Mods Table
Static lookup table for modification bitflags:
- `mod_bitflag` (INT PK): Sum of modification flags (2^n)
- `mods` (TEXT): Modification names

### Exp Table
Experiment metadata storage

### Refseq Table
Reference sequence information populated during ingestion

## File Naming Conventions

### BAM Files
Format: `{exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam`
- `exp_id`: Experiment identifier
- `bc_model_type`: Basecall model tier (s/h/f for standard/high/fast)
- `bc_model_version`: Version number (e.g., 5.2.0)
- `trim`: Binary flag (0/1)
- `modifications`: Modification bitflag value

### Modification Bitflags
Modifications are encoded as integer sums of powers of 2:
- `non` = 0 (no modifications)
- `6mA` = 1 (can coexist with others)
- `5mCG_5hmCG` = 2 (mutually exclusive C-mod)
- `5mC_5hmC` = 4 (mutually exclusive C-mod)
- `4mC_5mC` = 8 (mutually exclusive C-mod)
- `5mC` = 16 (mutually exclusive C-mod)

Example: `5mC_5hmC` + `6mA` = 4 + 1 = 5

## Quality Metrics

### Basecall Quality (q_bc)
Probability-averaged quality score:
```
q_bc = -10 * log10(sum(10^(-Q_base/10)) / n)
```

### Levenshtein Quality (q_ld)
Calculated only when read matches a reference sequence:
```
q_ld = -10 * log10(min(max(1/L^2, ed/L), 1))
```
where L is the length of the matched RefSeq and ed is the Levenshtein distance.

## Coding Conventions

1. **NULL Handling**: Use `None` (Python NULL) for:
   - `refseq_id` when read length is outside reference range
   - `ed` and `q_ld` when no RefSeq match exists

2. **RefSeq Matching**: Compare read length against `refseqrange` (Len_seq - 150, Len_seq + 150)

3. **BAM Tagging**: Add SAM-compliant tag `ER:Z:<end_reason>` to output BAM files

4. **File Organization**:
   - Input files: `Input/` directory with symlinks
   - Output files: `Output/{exp_id}/` directory
   - Database: Root directory as `SMA_{exp_id}.db`

5. **Database Operations**: Always commit transactions and close file handles properly

## Testing and Validation

When making changes:
1. Ensure database schema integrity is maintained
2. Validate filename parsing against the naming convention
3. Verify NULL handling for reads outside RefSeq ranges
4. Check quality metric calculations match formulas
5. Ensure BAM file tagging follows SAM specification

## Best Practices

- Use context managers (`with` statements) for file operations
- Validate input file formats before processing
- Handle edge cases: reads outside RefSeq range, missing Pod5 data
- Maintain transaction atomicity for database operations
- Use appropriate data types matching the schema (TEXT, INT, REAL)
- Follow Python PEP 8 style guidelines
