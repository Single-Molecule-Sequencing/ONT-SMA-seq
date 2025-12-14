# ONT-SMA-seq Copilot Instructions

## About This Project

ONT-SMA-seq is a bioinformatics pipeline for analyzing Oxford Nanopore Technology (ONT) sequencing data. This repository implements a three-stage Python pipeline that processes nanopore reads, calculates quality metrics, and stores results in a SQLite database.

## Development Guidelines

### Technology Choices
- **Language**: Python 3.x exclusively
- **Database**: SQLite3 for all data storage
- **Required Libraries**: 
  - `pysam` for BAM file operations
  - `pod5` for extracting End Reason metadata
  - `edlib` for Levenshtein distance calculations
  - `argparse` for CLI argument parsing

### Pipeline Architecture

The system consists of three independent Python scripts:

#### 1. mkdb.py - Database Initialization
Creates the SQLite schema and populates static reference tables.
- Input: Experiment ID (exp_id)
- Output: `SMA_{exp_id}.db` file
- Creates tables: Reads, Mods, Exp, Refseq
- Pre-populates Mods table with modification bitflag mappings

#### 2. inputInit.py - Input Standardization  
Parses metadata from filenames and creates standardized symlinks.
- Parses BAM filename: `{exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam`
- Creates `Input/` directory structure
- Symlinks: BAM file, Pod5 directory, Reference FASTA

#### 3. ingest.py - Read Processing & Database Population
Main pipeline that processes reads and populates the database.
- Parses reference FASTA (exactly 2 sequences: long and short)
- Builds Pod5 lookup dictionary for End Reason values
- Streams BAM reads with calculated metrics into database
- Tags output BAM with `ER:Z:<end_reason>` SAM tag

## Critical Implementation Details

### File Naming Convention
BAM files must follow this exact format:
```
{exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam
```
- `bc_model_type`: Single letter - 's' (standard), 'h' (high), or 'f' (fast)
- `bc_model_version`: Semantic version (e.g., "5.2.0")
- `trim`: Binary - 0 or 1
- `modifications`: Integer bitflag value

### Modification Bitflag Encoding
Modifications use power-of-2 encoding (can sum for combinations):
```
non          = 0   (no modifications)
6mA          = 1   (independent, can combine with C-mods)
5mCG_5hmCG   = 2   (mutually exclusive C-mod)
5mC_5hmC     = 4   (mutually exclusive C-mod)
4mC_5mC      = 8   (mutually exclusive C-mod)
5mC          = 16  (mutually exclusive C-mod)
```
Example: `5mC_5hmC` (4) + `6mA` (1) = 5

### Database Schema

#### Reads Table
Primary data storage. Key columns:
```
uniq_id       TEXT PRIMARY KEY   Format: {exp_id}{tier}{ver}t{trim}m{mod_flag}_{read_hash}
refseq_id     TEXT FOREIGN KEY   NULL if read length outside reference range
ed            INT                Levenshtein distance (NULL if refseq_id is NULL)
q_bc          REAL               Probability-averaged basecall quality
q_ld          REAL               Levenshtein quality (NULL if refseq_id is NULL)
ER            TEXT               End Reason from Pod5
```

#### Mods Table (Static Reference)
```
mod_bitflag   INT PRIMARY KEY    Sum of modification flags
mods          TEXT               Comma-separated modification names
```

#### Exp Table
```
exp_id        TEXT PRIMARY KEY   Experiment identifier
exp_desc      TEXT               Description
```

#### Refseq Table (Populated by ingest.py)
```
refseq_id     TEXT PRIMARY KEY   Sequence defline
refseq        TEXT               Full sequence
reflen        INT                Sequence length
```

### Quality Metric Formulas

**Basecall Quality (q_bc)** - Always calculated:
```python
q_bc = -10 * log10(sum(10**(-Q_base/10)) / n)
```

**Levenshtein Quality (q_ld)** - Only when refseq_id is matched:
```python
q_ld = -10 * log10(min(max(1/L**2, ed/L), 1))
```
where `L` = matched RefSeq length, `ed` = Levenshtein distance

### RefSeq Matching Logic

Reference sequences have a tolerance range:
```
refseqrange = (Len_seq - 150, Len_seq + 150)
```

When processing reads:
- If `read_length` within ANY refseqrange → assign that refseq_id, calculate ed and q_ld
- If `read_length` outside ALL ranges → set refseq_id = NULL, ed = NULL, q_ld = NULL

**Important**: Reads with NULL refseq_id are still inserted into the database.

## Coding Standards

### NULL Handling
Use Python `None` for SQL NULL values:
- `refseq_id = None` when read length doesn't match any reference
- `ed = None` and `q_ld = None` when refseq_id is None

### File Organization
```
Input/
  {exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam
  {exp_id}_pod5/
  {exp_id}.fa

Output/
  {exp_id}/
    (output BAM files)

SMA_{exp_id}.db
```

### BAM Processing
- Always use context managers (`with pysam.AlignmentFile(...)`)
- Add SAM tag `ER:Z:<end_reason>` to every output read
- Output BAM goes to `Output/{exp_id}.bam`

### Database Operations
- Use parameterized queries (no string interpolation)
- Commit transactions explicitly
- Close all database connections in finally blocks
- Match data types: TEXT for strings, INT for integers, REAL for floats

### Python Style
- Follow PEP 8
- Use type hints where beneficial
- Validate inputs before processing
- Handle missing Pod5 data gracefully
- Log progress for long-running operations

## Common Patterns

### Parsing BAM Filename
```python
# Extract from: exp123_h_v5.2.0_1_5.bam
parts = filename.replace('.bam', '').split('_')
exp_id = parts[0]
bc_model_type = parts[1]  # 's', 'h', or 'f'
bc_model_version = parts[2].replace('v', '')
trim = int(parts[3])
mod_bitflag = int(parts[4])
```

### Constructing uniq_id
```python
# Format: {exp_id}{tier}{ver}t{trim}m{mod_flag}_{read_hash}
# Example: exp123h5.2.0t1m5_abc123
import hashlib
read_hash = hashlib.md5(read_id.encode()).hexdigest()[:8]
uniq_id = f"{exp_id}{tier}{ver}t{trim}m{mod_flag}_{read_hash}"
```

### Pod5 End Reason Lookup
```python
import pod5
# Build lookup dictionary
end_reasons = {}
for pod5_file in pod5_directory:
    with pod5.Reader(pod5_file) as reader:
        for read in reader:
            end_reasons[read.read_id] = read.end_reason
```

## Testing Checklist

When implementing or modifying code:
- [ ] Filename parsing handles all valid format variations
- [ ] Database schema matches specification exactly
- [ ] NULL values inserted correctly for unmatched reads
- [ ] RefSeq range calculation: (length - 150, length + 150)
- [ ] Quality metrics use correct formulas
- [ ] BAM output includes ER tag
- [ ] All file handles and DB connections closed properly
- [ ] Edge cases handled: missing Pod5 data, reads outside all ranges
