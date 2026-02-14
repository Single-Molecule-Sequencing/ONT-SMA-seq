# ONT-Native Configuration & Truncated Molecule Detection Design

## Goal

Replace ad-hoc CLI flags and custom CSV parsing with ONT-native configuration
formats (Dorado TOML arrangement + MinKNOW sample sheet). Add structural
detection of truncated library molecules with auto-generated truncated reference
sequences. Provide an interactive wizard for building construct definitions.

## Architecture

Use Dorado's TOML arrangement format as the base layer, extended with an `[sma]`
section for barcode pairing, confidence thresholds, and truncation rules. TOML
ignores unknown sections, so the file remains valid Dorado input. A CLI wizard
generates the TOML interactively. A reference generator tool creates truncated
reference FASTAs from the construct definition. The ingest pipeline uses
structural element detection combined with confidence gating to classify each
read into one of five truncation levels.

## Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Config format | Extended Dorado TOML | Single source of truth; ONT-compatible |
| Config architecture | Single file (not companion) | Atomic construct definition |
| Truncation detection | Structural + confidence gating | Most rigorous classification |
| Truncated references | Pre-generated tool (mkrefs.py) | Explicit, inspectable, reusable |
| Barcode modes | Dual independent + start-only | Covers SMA-seq use cases |
| Wizard | Interactive CLI (construct_wizard.py) | User-friendly TOML generation |

---

## 1. Construct TOML Format

### Standard Dorado Sections

```toml
[arrangement]
name = "SMA_CYP2D6_dual"
kit = "SQK-NBD114-96"

# Flanking sequences around barcode 1 (5' end)
mask1_front = "AAGGTTAA"       # leading flank before barcode1
mask1_rear  = "CAGCACCT"       # trailing flank after barcode1

# Flanking sequences around barcode 2 (3' end, forward orientation)
mask2_front = "AGGTGCTG"       # leading flank before RC(barcode2)
mask2_rear  = "TTAACCTTAGCAAT" # trailing flank after RC(barcode2)

# Barcode naming pattern (matches FASTA headers)
barcode1_pattern = "NB%02i"
barcode2_pattern = "NB%02i"
first_index = 1
last_index = 96

[scoring]
max_barcode_penalty = 11
min_barcode_penalty_dist = 3
front_barcode_window = 100
rear_barcode_window = 100
min_flank_score = 0.5
```

### SMA Extension Section

```toml
[sma]
mode = "dual_independent"      # "dual_independent" | "start_only"
barcode_fasta = "barcodes.fasta"  # optional: custom barcode sequences

[sma.confidence]
full_length_threshold = 0.75   # bc_end_conf >= this = full-length
start_barcode_min = 0.6        # bc_start_conf below this = adapter_only
flank_max_error_rate = 0.5     # flank ED / flank_len > this = not found

[sma.truncation]
auto_generate_refs = true
min_target_length = 20         # min bases of target to count as present

[[sma.targets]]
barcode1 = "NB05"
barcode2 = "NB10"
alias = "CYP2D6_v04_fwd"
reference = "references/CYP2D6_v04_fwd.fasta"

[[sma.targets]]
barcode1 = "NB10"
barcode2 = "NB05"
alias = "CYP2D6_v04_rev"
reference = "references/CYP2D6_v04_rev.fasta"
```

### Field Reference

**[arrangement]** (Dorado standard):

| Field | Required | Description |
|-------|----------|-------------|
| `name` | Yes | Arrangement identifier |
| `kit` | Yes | ONT barcode kit (e.g., SQK-NBD114-96) |
| `mask1_front` | Yes | Leading flank before barcode1 |
| `mask1_rear` | Yes | Trailing flank after barcode1 |
| `mask2_front` | Dual only | Leading flank before RC(barcode2) |
| `mask2_rear` | Dual only | Trailing flank after RC(barcode2) |
| `barcode1_pattern` | Yes | FASTA header pattern for front barcodes |
| `barcode2_pattern` | Dual only | FASTA header pattern for rear barcodes |
| `first_index` | Yes | First barcode index |
| `last_index` | Yes | Last barcode index |

**[scoring]** (Dorado standard):

| Field | Default | Description |
|-------|---------|-------------|
| `max_barcode_penalty` | 11 | Maximum edit distance for classification |
| `min_barcode_penalty_dist` | 3 | Min gap between top-2 candidates |
| `front_barcode_window` | 100 | Search window at 5' end (bp) |
| `rear_barcode_window` | 100 | Search window at 3' end (bp) |
| `min_flank_score` | 0.5 | Min alignment score for flanking sequences |

**[sma]** (SMA extension):

| Field | Default | Description |
|-------|---------|-------------|
| `mode` | `dual_independent` | `dual_independent` or `start_only` |
| `barcode_fasta` | built-in | Path to custom barcode FASTA |

**[sma.confidence]**:

| Field | Default | Description |
|-------|---------|-------------|
| `full_length_threshold` | 0.75 | bc_end_conf >= this = full-length |
| `start_barcode_min` | 0.6 | bc_start below this = adapter_only |
| `flank_max_error_rate` | 0.5 | Flank ED / len > this = not found |

**[sma.truncation]**:

| Field | Default | Description |
|-------|---------|-------------|
| `auto_generate_refs` | true | mkrefs.py generates truncated refs |
| `min_target_length` | 20 | Min target bp to count as present |

**[[sma.targets]]** (array):

| Field | Required | Description |
|-------|----------|-------------|
| `barcode1` | Yes | Upstream barcode ID (e.g., NB05) |
| `barcode2` | Dual only | Downstream barcode ID |
| `alias` | Yes | Target name |
| `reference` | Yes | Path to reference FASTA |

---

## 2. Truncated Molecule Detection

### Library Structure

Full-length SMA-seq library molecule (5' → 3'):

```
adapter ── mask1_front ── BC1 ── mask1_rear ── TARGET ── mask2_front ── RC(BC2) ── mask2_rear ── RC(adapter)
```

Truncated molecules are missing elements from the 3' end:

```
Full:        adapter─BC1─flank─TARGET─flank─RC(BC2)─adapter
bc1_target:  adapter─BC1─flank─TARGET
bc1_only:    adapter─BC1─flank
adapter_only: adapter
```

### Detection Algorithm

For each read, detect structural elements and apply confidence gates:

1. **Find BC1**: edlib HW alignment of all barcodes against first
   `front_barcode_window` bp. Best match = bc_start_id, bc_start_conf.

2. **Find front flank** (mask1_rear): edlib HW search in region after BC1
   position. Confirms BC1 boundary.

3. **Find rear flank** (mask2_front): edlib HW search in last
   `rear_barcode_window` bp region. If found, target end is located.

4. **Find BC2**: edlib HW alignment of RC(barcodes) against last
   `rear_barcode_window` bp. Best match = bc_end_id, bc_end_conf.

5. **Measure target region**: Distance between front flank end and rear flank
   start. Must be >= `min_target_length` to count as target present.

6. **Classify truncation level** (evaluated in order):

| Level | Structural test | Confidence gate |
|-------|----------------|-----------------|
| `full_length` | BC1 + front_flank + target + rear_flank + BC2 | bc_start >= 0.6 AND bc_end >= 0.75 |
| `bc1_target_bc2` | BC1 + front_flank + target + BC2 present but low | bc_start >= 0.6 AND 0.3 <= bc_end < 0.75 |
| `bc1_target` | BC1 + front_flank + target, no rear_flank/BC2 | bc_start >= 0.6, bc_end not detected |
| `bc1_only` | BC1 found, read too short for target | bc_start >= 0.6, len < bc1 + flank + min_target |
| `adapter_only` | BC1 below threshold | bc_start < 0.6 |

### Target Assignment

After truncation level is determined:

- **full_length / bc1_target_bc2**: Look up (bc_start, bc_end) pair in targets.
  If found → alias. If not → `unmatched_{bc_start}_{bc_end}`.
- **bc1_target**: Assign `{bc_start}_target_only`. If bc_start maps to exactly
  one target pair, infer the alias: `{alias}_truncated`.
- **bc1_only**: Assign `{bc_start}_bc_only`.
- **adapter_only**: Assign `unclassified`.

---

## 3. Reference Generation Tool (mkrefs.py)

### CLI

```
python bin/mkrefs.py -c construct.toml -o references/
```

### Generated Files

For each `[[sma.targets]]` entry with `auto_generate_refs = true`:

```
references/truncated/
  {alias}_full.fasta          # mask1_F + BC1 + mask1_R + target + mask2_F + RC(BC2) + mask2_R
  {alias}_bc1_target.fasta    # mask1_F + BC1 + mask1_R + target
  {bc1_id}_bc_only.fasta      # mask1_F + BC1 + mask1_R  (one per unique bc1)
```

### Validation

- Verifies all referenced barcode IDs exist in kit/FASTA
- Verifies all reference FASTAs exist and are valid
- Reports construct lengths for each truncation level
- Writes a manifest: `references/truncated/manifest.tsv`

---

## 4. Interactive Construct Wizard (construct_wizard.py)

### CLI

```
python bin/construct_wizard.py -o construct.toml [-rd references/]
python bin/construct_wizard.py -i existing.toml -o updated.toml  # edit mode
```

### Step-by-Step Flow

7 interactive steps with diagrams and explanations:

1. **Library Design Mode** — Select dual_independent or start_only. Shows
   construct diagrams for each mode.

2. **Barcode Kit** — Select from known ONT kits or provide custom FASTA.
   Shows barcode count and length.

3. **Flanking Sequences** — Enter mask1_front, mask1_rear (and mask2_front,
   mask2_rear for dual mode). Shows position diagram.

4. **Barcode Pair → Target Mapping** — Auto-discovers reference FASTAs in
   provided directory. User selects BC1 number, BC2 number, and reference
   file for each pair. Shows construct diagram per pair.

5. **Confidence Thresholds** — Defaults shown, user can override. Brief
   explanation of what each threshold controls.

6. **Truncation Settings** — Auto-generate refs toggle, min target length.

7. **Review & Save** — Full summary with construct diagram, pair table,
   truncation levels. Confirm and write TOML.

### Features

- Auto-discovers reference FASTAs in `-rd` directory
- Shows ASCII construct diagrams at each step
- Validates all inputs (barcode range, FASTA existence, sequence chars)
- Sensible defaults for all optional parameters
- Edit mode loads existing TOML for modification
- Preview before saving

---

## 5. Database Schema Changes

### Reads Table

Add column:

```sql
trunc_level TEXT  -- 'full_length', 'bc1_target_bc2', 'bc1_target', 'bc1_only', 'adapter_only'
```

NULL when not using construct TOML (backward compatible).

### No Other Schema Changes

Existing barcode columns (bc_start_id, bc_start_ed, bc_start_conf, bc_end_id,
bc_end_ed, bc_end_conf) and tgt_id are sufficient.

---

## 6. File Plan

| Action | File | Purpose | ~Lines |
|--------|------|---------|--------|
| Create | `bin/construct.py` | Parse/validate construct TOML | ~200 |
| Create | `bin/construct_wizard.py` | Interactive TOML builder | ~400 |
| Create | `bin/mkrefs.py` | Generate truncated references | ~150 |
| Modify | `bin/barcodes.py` | Load from FASTA or built-in | ~30 delta |
| Modify | `bin/sample_sheet.py` | Accept TOML targets as alt | ~20 delta |
| Modify | `bin/ingest.py` | Construct config + truncation | ~80 delta |
| Modify | `bin/mkdb.py` | Add trunc_level column | ~5 delta |
| Modify | `bin/report_analysis.py` | Truncation level stats | ~40 delta |
| Modify | `bin/report_template.py` | Truncation breakdown | ~60 delta |
| Create | `tests/test_construct.py` | TOML parsing/validation | ~150 |
| Create | `tests/test_mkrefs.py` | Ref generation tests | ~100 |
| Modify | `tests/test_ingest.py` | Truncation classification | ~60 delta |

---

## 7. Updated Workflow

```bash
# Step 0: Build construct definition (interactive wizard)
python bin/construct_wizard.py -o construct.toml -rd references/

# Step 1: Generate truncated references
python bin/mkrefs.py -c construct.toml -o references/

# Step 2: Create database
python bin/mkdb.py -e FAL99999_20260214_TEST -o Output/

# Step 3: Ingest reads
python bin/ingest.py -e FAL99999_20260214_TEST \
  -b reads.bam -s summary.tsv -d Output/SMA_*.db \
  -o tagged.bam -c construct.toml -rd references/

# Step 4: Generate report
python bin/report.py -d Output/SMA_*.db -c construct.toml -o report.html
```

---

## 8. Backward Compatibility

- `-c construct.toml` is optional everywhere. Without it, existing behavior
  is preserved.
- `-ss sample_sheet.csv` still works. Gets converted internally to the same
  pairing dict that the TOML produces.
- `--flank-front`/`--flank-rear` CLI flags still work as overrides.
- DB schema adds `trunc_level` column (NULL when not using construct TOML).
- Old databases without `trunc_level` still work with report.py (graceful
  fallback).
- `barcodes.py` still has built-in 96 barcodes. FASTA loading is additive.

---

## 9. ONT Format Compatibility

| ONT Format | How We Use It |
|------------|---------------|
| Dorado TOML arrangement | Base format for construct.toml |
| Dorado TOML scoring | Scoring parameters in [scoring] section |
| Dorado barcode FASTA | Optional custom barcode sequences |
| MinKNOW sample sheet CSV | Accepted via `-ss` (backward compat) |
| Sequencing summary TSV | End reason extraction (existing) |
| MinKNOW output structure | BAM/POD5 discovery patterns |

The `[sma]` section is transparent to Dorado (TOML ignores unknown sections).
The construct.toml can be used directly with `dorado demux --barcode-arrangement`
for the standard demultiplexing fields.
