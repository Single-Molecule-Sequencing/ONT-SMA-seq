# ONT-SMA-seq Pipeline Summary

## Pipeline Overview

```
┌─────────────┐      ┌──────────────┐      ┌─────────────┐
│  mkdb.py    │ ───> │ inputInit.py │ ───> │  ingest.py  │
│ (Database)  │      │  (Organize)  │      │  (Process)  │
└─────────────┘      └──────────────┘      └─────────────┘
       │                     │                     │
       ↓                     ↓                     ↓
 SMA_exp.db            Input/ dir           Output/ dir
  (Schema)             (Symlinks)          (Tagged BAM)
                                           + Full DB
```

## Inputs and Outputs

### Script 1: mkdb.py

```
INPUT:                          OUTPUT:
┌──────────────────┐           ┌──────────────────────────┐
│ exp_id (string)  │  ──────>  │ SMA_{exp_id}.db          │
└──────────────────┘           │   ├── Reads (empty)      │
                               │   ├── Mods (populated)   │
                               │   ├── Exp (empty)        │
                               │   └── Refseq (empty)     │
                               └──────────────────────────┘
```

### Script 2: inputInit.py

```
INPUTS:                                  OUTPUTS:
┌────────────────────────────┐          ┌─────────────────────────────┐
│ Raw BAM (with metadata in  │          │ Input/                      │
│ filename)                  │  ─────>  │   ├── {exp_id}*.bam ──> src │
│                            │          │   ├── {exp_id}_pod5/ ──> src│
│ Pod5 directory             │  ─────>  │   └── {exp_id}.fa ──> src   │
│                            │          └─────────────────────────────┘
│ Reference FASTA            │  ─────>
│ (2 sequences)              │
└────────────────────────────┘
```

### Script 3: ingest.py

```
INPUTS:                                  OUTPUTS:
┌────────────────────────────┐          ┌─────────────────────────────┐
│ SMA_{exp_id}.db (schema)   │  ─────>  │ SMA_{exp_id}.db (populated) │
│                            │          │   ├── Reads (all data)      │
│ Input/{exp_id}*.bam        │  ─────>  │   └── Refseq (2 refs)       │
│                            │          │                             │
│ Input/{exp_id}_pod5/       │  ─────>  │ Output/{exp_id}.bam         │
│                            │          │   (Tagged with ER field)    │
│ Input/{exp_id}.fa          │  ─────>  └─────────────────────────────┘
└────────────────────────────┘
```

## Database Schema (Simplified)

```
┌─────────────────────────────────────────────────────────────┐
│                        Reads Table                          │
├─────────────────────────────────────────────────────────────┤
│ Primary Key: uniq_id                                        │
├─────────────────────────────────────────────────────────────┤
│ Foreign Keys:                                               │
│   • exp_id ────────────────> Exp.exp_id                     │
│   • refseq_id ──────────────> Refseq.refseq_id (nullable)   │
│   • mod_bitflag ────────────> Mods.mod_bitflag              │
├─────────────────────────────────────────────────────────────┤
│ Core Data:                                                  │
│   • read_id, readseq, readlen                               │
│   • model_tier, model_ver, trim                             │
├─────────────────────────────────────────────────────────────┤
│ Calculated Metrics:                                         │
│   • ed (edit distance, NULL if no refseq match)             │
│   • q_bc (basecall quality)                                 │
│   • q_ld (Levenshtein quality, NULL if no refseq match)     │
├─────────────────────────────────────────────────────────────┤
│ Metadata:                                                   │
│   • ER (end reason from Pod5)                               │
│   • channel, well, pore_type, num_samples, etc.             │
└─────────────────────────────────────────────────────────────┘

┌──────────────────────┐  ┌──────────────────────┐
│    Mods Table        │  │    Refseq Table      │
├──────────────────────┤  ├──────────────────────┤
│ mod_bitflag (PK)     │  │ refseq_id (PK)       │
│ mods (name)          │  │ refseq (sequence)    │
└──────────────────────┘  │ reflen (length)      │
                          └──────────────────────┘

┌──────────────────────┐
│    Exp Table         │
├──────────────────────┤
│ exp_id (PK)          │
│ exp_desc             │
└──────────────────────┘
```

## Data Processing Flow

```
For each read in BAM:
  
  1. Extract Read Data
     ├── read_id (UUID)
     ├── readseq (ATCG...)
     ├── readlen (integer)
     └── quality_scores (array)
  
  2. Lookup Pod5 Metadata
     └── end_reason ── from {read_id → end_reason} map
  
  3. Add ER Tag
     └── Write read with ER:Z:{end_reason} to Output BAM
  
  4. Match to Reference
     ├── IF readlen in (reflen_long ± 150):
     │   └── refseq_id = "long"
     ├── ELSE IF readlen in (reflen_short ± 150):
     │   └── refseq_id = "short"
     └── ELSE:
         └── refseq_id = NULL
  
  5. Calculate Metrics
     ├── q_bc = -10*log10(mean(10^(-Q/10)))
     ├── IF refseq_id is NOT NULL:
     │   ├── ed = Levenshtein(readseq, refseq)
     │   └── q_ld = -10*log10(min(max(1/L², ed/L), 1))
     └── ELSE:
         ├── ed = NULL
         └── q_ld = NULL
  
  6. Generate Unique ID
     └── uniq_id = {exp_id}{tier}{ver}t{trim}m{mod}_{hash}
  
  7. Insert to Database
     └── INSERT INTO Reads (all columns) VALUES (...)
```

## Modification Bitflags

```
Base Values (2^n):
┌────────────────┬───────┬─────────────────────┐
│ Modification   │ Value │ Type                │
├────────────────┼───────┼─────────────────────┤
│ non            │   0   │ None                │
│ 6mA            │   1   │ Independent         │
│ 5mCG_5hmCG     │   2   │ C-Mod (exclusive)   │
│ 5mC_5hmC       │   4   │ C-Mod (exclusive)   │
│ 4mC_5mC        │   8   │ C-Mod (exclusive)   │
│ 5mC            │  16   │ C-Mod (exclusive)   │
└────────────────┴───────┴─────────────────────┘

Combinations (6mA + one C-Mod):
┌────────────────┬───────┐
│ 6mA+5mCG_5hmCG │   3   │  (1 + 2)
│ 6mA+5mC_5hmC   │   5   │  (1 + 4)
│ 6mA+4mC_5mC    │   9   │  (1 + 8)
│ 6mA+5mC        │  17   │  (1 + 16)
└────────────────┴───────┘
```

## Filename Convention

```
Format: {exp_id}_{model_tier}_v{model_ver}_{trim}_{modifications}.bam

Example: EXP001_h_v5.2.0_1_6mA.bam
         ──┬──  ┬  ──┬──   ┬  ─┬─
           │    │    │      │   └─── Modifications (6mA)
           │    │    │      └─────── Trim (1 = yes, 0 = no)
           │    │    └────────────── Model version (5.2.0)
           │    └─────────────────── Model tier (h = high accuracy)
           └──────────────────────── Experiment ID

Parsed to:
  • exp_id: "EXP001"
  • model_tier: "h"
  • model_ver: "5.2.0"
  • trim: 1
  • mods_str: "6mA"
  • mod_bitflag: 1
```

## Unique ID Format

```
Format: {exp_id}{tier}{ver}t{trim}m{mod}_{hash}

Components:
  • exp_id: EXP001
  • tier: h
  • ver: 520 (5.2.0 with dots removed)
  • trim: 1
  • mod: 1 (bitflag)
  • hash: a3f2d8e1 (MD5 of read_id, first 8 chars)

Result: EXP001h520t1m1_a3f2d8e1
```

## Metric Formulas

### Basecall Quality (q_bc)
```
q_bc = -10 × log₁₀(Σ(10^(-Qᵢ/10)) / n)

Where:
  Qᵢ = individual base quality scores
  n = number of bases
```

### Levenshtein Distance (ed)
```
ed = LevenshteinDistance(read_sequence, reference_sequence)

Only calculated when read length matches reference (± 150 bp)
```

### Levenshtein Quality (q_ld)
```
q_ld = -10 × log₁₀(min(max(1/L², ed/L), 1))

Where:
  ed = edit distance
  L = reference length
```

## Quick Start Commands

```bash
# 1. Create database
python mkdb.py EXP001

# 2. Standardize inputs
python inputInit.py \
  --bam data/EXP001_h_v5.2.0_1_6mA.bam \
  --pod5 data/pod5/ \
  --ref data/reference.fa

# 3. Process data
python ingest.py EXP001
```

## Output Files

```
project/
├── SMA_EXP001.db                    # Complete SQLite database
│   ├── Reads: All processed reads with metrics
│   ├── Mods: Modification lookup table
│   ├── Exp: Experiment metadata
│   └── Refseq: Reference sequences
│
├── Input/                           # Standardized inputs (symlinks)
│   ├── EXP001_h_v5.2.0_1_6mA.bam   # → original BAM
│   ├── EXP001_pod5/                 # → original pod5 dir
│   └── EXP001.fa                    # → original reference
│
└── Output/
    └── EXP001.bam                   # Tagged BAM with ER field
```

## Dependencies

```bash
pip install pysam pod5 edlib
```

- **pysam**: BAM file I/O
- **pod5**: End reason extraction from Pod5 files
- **edlib**: Levenshtein distance calculation
- **sqlite3**: Database (Python built-in)

## Key Features

✓ Single-pass database ingestion
✓ Efficient BAM streaming (no full load to memory)
✓ Reference matching by length (± 150 bp tolerance)
✓ Quality metrics calculation (q_bc, q_ld)
✓ Edit distance for matched reads
✓ NULL handling for unmatched reads
✓ Modification bitflag system
✓ Pod5 metadata integration
✓ Tagged output BAM with end reasons
