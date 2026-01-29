# Overhaul Plan: Single-Target Database Ingestion

## Architecture Overview

The pipeline is now restructured to operate on a **Per-Target Sequence** basis. For every target sequence in an overarching experiment, a dedicated SQLite database is generated. This creates a modular "building block" architecture where individual Target DBs can be merged into a master database later.

**Language:** Python 3.x (with Shell/CLI for Pod5 extraction)

**Key Libraries:** `pysam` (BAM I/O), `pandas` (CSV/TSV parsing), `edlib` (Levenshtein), `sqlite3`, `argparse`.

**CLI Tools:** `pod5` (via `pod5 view`).

**Scope:** 1 Input BAM + 1 Input Pod5 Dir + 1 Reference FASTA (containing **single** sequence) -> 1 SQLite DB.

**Philosophy:** "Capture All, Filter Later." All reads are ingested and compared against the single target reference. Length filtering and quality thresholds are applied during post-hoc SQL analysis instead of during ingestion.

**Mergeability:** Primary Keys (`uniq_id`, `tgt_id`, `exp_id`) are designed to be distinct or consistent to allow safe merging of multiple Target DBs.

---

## **1. Script Breakdown & Workflow**

### **Script A: `mkdb.py` (Database Initialization)**

**Purpose:** Initializes the SQLite database schema and populates static lookup tables.

**Input:** `exp_id` (via argparse).

**Output:** `SMA_{exp_id}_{tgt_id}.db` (Naming convention suggested to avoid overwrites if running in parallel, though user defines output path).

**Logic:**

1. Create the database file.
2. Execute DDL to create tables: `Reads`, `Mods`, `Exp`, `Target`.
3. **Populate `Mods` Table:** Insert the standard bitflag definitions (see *Modification Bitflags* section).

### **Script B: `inputInit.py` (Input Standardization)**

**Purpose:** Standardizes input file paths and directory structures.

**Inputs (via argparse):**

1. Path to raw uBAM.
2. Path to raw Pod5 directory.
3. Path to Reference FASTA (**Must contain exactly 1 sequence**).

**Logic:**

1. **Parse BAM Filename:** Extract metadata components (`exp_id`, `model`, `ver`, `trim`, `mods`).
2. **Directory Setup:** Create `Input/` directory.
3. **Symlink Inputs:**
   1. BAM  `Input/reads.bam` (Simplified name since metadata is in DB).
   2. Pod5 Dir  `Input/pod5/`
   3. RefSeq  `Input/target.fa`

### **Script C: `extractMeta.py` (Metadata Extraction)**

**Purpose:** Generates the lightweight metadata summary required for ingestion using the optimized `pod5` CLI.

**Inputs (via argparse):** Path to Input Pod5 directory.

**Logic:**

1. **Execute Extraction:** Run `pod5 view` on `Input/pod5/`.
   1. **Command:** `pod5 view Input/pod5/*.pod5 --include "read_id,end_reason" --output Input/summary.tsv`

### **Script D: `ingest.py` (Main Processing)**

**Purpose:** Parses inputs, calculates metrics against the single target, tags BAMs, and populates the database.

**Logic Flow:**

   1. **Parse Target (FASTA):**
      * Read `Input/target.fa`.
      * **Validation:** Ensure only **1 sequence** exists.
      * Extract `tgt_id` (defline header) and `tgt_refseq` (sequence).
      * **DB Action:** Insert entry into the `Target` table.
   2. **Load Metadata Summary:**
      * Read `Input/summary.tsv` into memory (Map: `{read_id: end_reason}`).
   3. **Stream & Process Reads:**
      * Open Input BAM and Output BAM.
      * **Loop through every read:**
        1. **ER Lookup:** Retrieve End Reason.
        2. **Tagging:** Write `ER:Z:<end_reason>` tag to read.
        3. **ID Generation:** Construct `uniq_id`.
        4. **Metric Calculation (All Reads):**
           * Calculate  (Basecall Quality).
           * Calculate  (Levenshtein Distance) vs `tgt_refseq`.
           * Calculate  (Levenshtein Quality).
   4. **DB Action:** Insert record into `Reads` table.
      * `tgt_id` is always the single ID extracted from FASTA.
   5. **Cleanup:** Commit DB transactions.

---

## **2. Filename Parsing & Metadata**

**Naming Convention:** `{exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam`

* **Standard Fields:** Extract `exp_id`, `model_tier` (s/h/f), `model_ver` (e.g., 5.2.0), `trim` (0/1).
* **Modification Bitflags:**
    Combinations of modifications are stored as integers (Sum of $2^n$).

| Modification   | Bit Value ($2^n$) | Logic                                 |
| :------------- | :---------------- | :------------------------------------ |
| **non**        | **0**             | No Modifications                      |
| **6mA**        | **1**             | Independent (Can coexist with others) |
| **5mCG_5hmCG** | **2**             | Mutually Exclusive C-Mod              |
| **5mC_5hmC**   | **4**             | Mutually Exclusive C-Mod              |
| **4mC_5mC**    | **8**             | Mutually Exclusive C-Mod              |
| **5mC**        | **16**            | Mutually Exclusive C-Mod              |

* *Example:* `5mC_5hmC` + `6mA` = $4 + 1 = 5$.

---

## **3. ID Construction & Metrics**

**`uniq_id` Format:** `{exp_id}{tier}{ver}t{trim}m{mod_flag}_{read_hash}`

**Metric Formulas:**

1. **$q_{bc}$ (Prob-Avg):**
    $$q_{bc} = -10 \log_{10} \left( \frac{\sum 10^{-Q_{base}/10}}{n} \right)$$
2. **Levenshtein Distance (`ed`):**
    $$ed = \text{Levenshtein}(ReadSeq, RefSeq)$$
    *(Calculated for ALL reads regardless of length)*
3. **$q_{ld}$ (Levenshtein Quality):**
    $$q_{ld} = -10 \log_{10} \left( \min \left( \max \left( \frac{1}{L^2}, \frac{ed}{L} \right), 1 \right) \right)$$
    *(Where $L$ is the length of the Target Reference)*

---

## **4. Database Schema**

**Reads** table

| Column        | Type          | Description                            |
| ------------- | ------------- | -------------------------------------- |
| `uniq_id`     | **TEXT (PK)** | Composite unique identifier.           |
| `exp_id`      | TEXT (FK)     | Experiment ID.                         |
| `tgt_id`      | TEXT (FK)     | Target Sequence ID.                    |
| `read_id`     | TEXT          | Original ONT Read UUID.                |
| `readseq`     | TEXT          | The Basecalled Read Sequence.          |
| `readlen`     | INT           | Length of the Read.                    |
| `model_tier`  | TEXT          | 's', 'h', or 'f'.                      |
| `model_ver`   | TEXT          | e.g., '5.2.0'.                         |
| `trim`        | INT           | 1 (True) or 0 (False).                 |
| `mod_bitflag` | INT (FK)      | Integer sum of modification flags.     |
| `ed`          | INT           | Levenshtein Distance vs Target.        |
| `q_bc`        | REAL          | Probability-averaged basecall quality. |
| `q_ld`        | REAL          | Levenshtein quality.                   |
| `ER`          | TEXT          | End Reason (from `pod5 view`).         |

**Target** table (Replaces Refseq)

| Column       | Type          | Description                    |
| ------------ | ------------- | ------------------------------ |
| `tgt_id`     | **TEXT (PK)** | Target Sequence ID (defline).  |
| `tgt_refseq` | TEXT          | The Target Sequence content.   |
| `tgt_reflen` | INT           | Length of the Target Sequence. |

**Exp** table

| Column         | Type          | Description                               |
| -------------- | ------------- | ----------------------------------------- |
| `exp_id`       | **TEXT (PK)** | Experiment ID {flow_cell_id}_{sample_id}. |
| `flow_cell_id` | TEXT          | Flow Cell ID.                             |
| `sample_id`    | TEXT          | Library Prep Info YYYYMMDD_INITIALS.      |
| `exp_desc`     | TEXT          | Experiment Description.                   |

**Mods** table (Static Lookup)

| Column        | Type         | Description                  |
| ------------- | ------------ | ---------------------------- |
| `mod_bitflag` | **INT (PK)** | Sum of modification flags () |
| `mods`        | TEXT         | Modifications present        |

---

## **5. Future**

Since every run produces a self-contained per-target DB with a shared Schema:

**Merging:** Can be done by attaching databases in SQLite into a master DB.

**Multithreading** might be on the table, too.
