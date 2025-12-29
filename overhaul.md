# Overhaul Plan: Single-Pass Database Ingestion

## **Architecture Overview**

The overhaul turns the original NextFlow pipeline into an efficient database structure using distinct scripts to separate concerns: database creation, input standardization, metadata extraction, and core ingestion logic.

* **Language:** Python 3.x (with Shell/CLI for Pod5 extraction)
* **Key Libraries:** `pysam` (BAM I/O), `pandas` (CSV/TSV parsing), `edlib` (Levenshtein), `sqlite3`, `argparse`.
* **CLI Tools:** `pod5` (via `pod5 view`).

---

## **1. Script Breakdown & Workflow**

### **Script A: `mkdb.py` (Database Initialization)**

* **Purpose:** Initializes the SQLite database schema and populates static lookup tables.
* **Input:** `exp_id` (via argparse).
* **Output:** `SMA_{exp_id}.db`.
* **Logic:**
    1. Create the database file `SMA_{exp_id}.db`.
    2. Execute DDL to create tables: `Reads`, `Mods`, `Exp`, `Refseq`.
    3. **Populate `Mods` Table:** Insert the standard bitflag definitions (see *Modification Bitflags* section).

### **Script B: `inputInit.py` (Input Standardization)**

* **Purpose:** Standardizes input file paths and directory structures.
* **Inputs (via argparse):**
    1. Path to raw uBAM (Must follow naming convention: `{exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam`).
    2. Path to raw Pod5 directory.
    3. Path to Reference FASTA.
* **Logic:**
    1. **Parse BAM Filename:** Extract metadata components (`exp_id`, `model`, `ver`, `trim`, `mods`) directly from the input BAM filename.
    2. **Directory Setup:** Create `Input/` directory if it doesn't exist.
    3. **Symlink Inputs:**
        * BAM $\to$ `Input/{exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam`
        * Pod5 Dir $\to$ `Input/{exp_id}_pod5/`
        * RefSeq $\to$ `Input/{exp_id}.fa`

### **Script C: `extractMeta.py` (Metadata Extraction)**

* **Purpose:** Generates the lightweight metadata summary required for ingestion using the optimized pod5 CLI.
* **Inputs (via argparse):**
    1. `exp_id` (used to locate the standardized input directory).
* **Logic:**
    1. **Locate Input:** Target `Input/{exp_id}_pod5/`.
    2. **Execute Extraction:** Run the `pod5 view` command via `subprocess` or shell wrapper.
        * **Command:** `pod5 view Input/{exp_id}_pod5/*.pod5 --include "read_id,end_reason" --output Input/{exp_id}_summary.tsv`
        * *Note:* Only extracting `read_id` and `end_reason` to keep I/O minimal.

### **Script D: `ingest.py` (Main Processing)**

* **Purpose:** Parses inputs, calculates metrics, tags BAMs, and populates the database.
* **Logic Flow:**
    1. **Parse Reference (FASTA):**
        * Read `Input/{exp_id}.fa` (contains exactly 2 sequences: long and short).
        * Calculate `refseqrange` for each: $(Len_{seq} - 150, Len_{seq} + 150)$.
        * **DB Action:** Insert entries into the `Refseq` table.
    2. **Load Metadata Summary (Fast):**
        * Read `Input/{exp_id}_summary.tsv` (generated in Script C) into a Pandas DataFrame.
        * **Map:** `{read_id: end_reason}`.
    3. **Stream & Process Reads:**
        * Create `Output/{exp_id}` directory.
        * Open Input BAM (Read) and `Output/{exp_id}.bam` (Write).
        * **Loop through every read:**
            1. **ER Lookup:** Retrieve End Reason from the loaded summary dictionary.
            2. **Tagging:** Write read to Output BAM with new SAM-compliant tag `ER:Z:<end_reason>`.
            3. **RefSeq Matching:** Compare `readlen` against the calculated `refseqrange` of the two references.
                * *If inside range:* Assign specific `refseq_id`.
                * *If outside range:* Set `refseq_id` to `NULL` (None).
            4. **ID Generation:** Construct `uniq_id`.
            5. **Metric Calculation:**
                * Calculate $q_{bc}$ (Basecall Quality).
                * **Conditional Metric:**
                    * *If `refseq_id` is identified:* Calculate $ed$ (Levenshtein) and $q_{ld}$.
                    * *If `refseq_id` is NULL:* Set $ed$ and $q_{ld}$ to `NULL` (None).
            6. **DB Action:** Insert the record into `Reads` table (including those with NULL RefSeq).
    4. **Cleanup:** Commit DB transactions and close file handles.

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

* *Note: read_hash is a shortened hash of the read_id/sequence to ensure uniqueness.*

**Metric Formulas:**

1. **$q_{bc}$ (Prob-Avg):**
    $$q_{bc} = -10 \log_{10} \left( \frac{\sum 10^{-Q_{base}/10}}{n} \right)$$
2. **Levenshtein Distance (`ed`):**
    $$ed = \text{Levenshtein}(ReadSeq, RefSeq)$$
    *(Only calculated if read fits in a RefSeq length range)*
3. **$q_{ld}$ (Levenshtein Quality):**
    $$q_{ld} = -10 \log_{10} \left( \min \left( \max \left( \frac{1}{L^2}, \frac{ed}{L} \right), 1 \right) \right)$$
    *(Where $L$ is the length of the matched RefSeq)*

---

## **4. Database Schema**

**Reads** table

| Column          | Type          | Description                                  |
| :-------------- | :------------ | :------------------------------------------- |
| `uniq_id`       | **TEXT (PK)** | Composite unique identifier.                 |
| `exp_id`        | TEXT (FK)     | Experiment ID.                               |
| `refseq_id`     | TEXT (FK)     | RefSeq ID. **NULL** if length out of range.  |
| `read_id`       | TEXT          | Original ONT Read UUID.                      |
| `readseq`       | TEXT          | The Basecalled Read Sequence.                |
| `readlen`       | INT           | Length of the Read.                          |
| `model_tier`    | TEXT          | 's', 'h', or 'f'.                            |
| `model_ver`     | TEXT          | e.g., '5.2.0'.                               |
| `trim`          | INT           | 1 (True) or 0 (False).                       |
| `mod_bitflag`   | INT (FK)      | Integer sum of modification flags.           |
| `ed`            | INT           | Levenshtein Distance. **NULL** if no RefSeq. |
| `q_bc`          | REAL          | Probability-averaged basecall quality.       |
| `q_ld`          | REAL          | Levenshtein quality. **NULL** if no RefSeq.  |
| `ER`            | TEXT          | End Reason (from `pod5 view`).               |

**Mods** table (Populated by `mkdb.py`)

| Column        | Type         | Description                       |
| :------------ | :----------- | :-------------------------------- |
| `mod_bitflag` | **INT (PK)** | Sum of modification flags ($2^n$) |
| `mods`        | TEXT         | Modifications present             |

**Exp** table

| Column     | Type          | Description             |
| :--------- | :------------ | :---------------------- |
| `exp_id`   | **TEXT (PK)** | Experiment ID.          |
| `exp_desc` | TEXT          | Experiment description. |

**Refseq** table (Populated by `ingest.py`)

| Column      | Type          | Description                      |
| :---------- | :------------ | :------------------------------- |
| `refseq_id` | **TEXT (PK)** | Reference Sequence ID (defline). |
| `refseq`    | TEXT          | Reference Sequence.              |
| `reflen`    | INT           | Reference Sequence Length.       |
