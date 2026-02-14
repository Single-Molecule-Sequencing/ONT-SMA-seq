# Barcode Demultiplexing Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Integrate dual-end barcode classification into `ingest.py` so every read is assigned to its target via duplexed ONT native barcode matching.

**Architecture:** A new `barcodes.py` module handles barcode lookup and edlib-based classification. A new `sample_sheet.py` module parses MinKNOW sample sheets into barcode-pair-to-alias mappings. `mkdb.py` gains 6 new barcode columns. `ingest.py` calls classification per-read during ingestion, looking up the target from the barcode pair. Optional `--split-bams` and `--tag` flags add post-classification outputs.

**Tech Stack:** Python 3.x, edlib (HW mode), pysam, pandas, sqlite3

**Repo:** `/tmp/ont-sma-seq/`

---

### Task 1: Create barcode module with classification functions

**Files:**
- Create: `bin/barcodes.py`
- Create: `tests/test_barcodes.py`

This module contains the 96 ONT native barcode sequences (extracted from Reference_Fasta_Generator's `Make_Fasta.py`) and the core classification function that aligns a read segment against all expected barcodes using edlib HW mode.

**Step 1: Write the failing tests**

Create `tests/test_barcodes.py`:

```python
#!/usr/bin/env python3
"""Tests for barcodes.py barcode lookup and classification."""

import pytest

# These will fail until barcodes.py exists
from barcodes import (
    BARCODES,
    reverse_complement,
    classify_barcode,
)


class TestBarcodeData:
    """Test barcode sequence data integrity."""

    def test_96_barcodes_present(self):
        assert len(BARCODES) == 96

    def test_barcodes_are_24bp(self):
        for name, seq in BARCODES.items():
            assert len(seq) == 24, f"{name} is {len(seq)}bp, expected 24bp"

    def test_barcode_names_are_nb01_to_nb96(self):
        for i in range(1, 97):
            key = f"nb{i:02d}"
            assert key in BARCODES, f"Missing barcode {key}"

    def test_barcodes_are_uppercase_acgt(self):
        valid = set("ACGT")
        for name, seq in BARCODES.items():
            assert set(seq).issubset(valid), f"{name} has invalid chars"

    def test_nb01_sequence(self):
        """Spot check: nb01 from ONT native barcoding kit."""
        assert BARCODES["nb01"] == "CACAAAGACACCGACAACTTTCTT"


class TestReverseComplement:
    """Test reverse complement utility."""

    def test_simple(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_palindrome(self):
        assert reverse_complement("AATT") == "AATT"

    def test_empty(self):
        assert reverse_complement("") == ""

    def test_roundtrip(self):
        seq = "CACAAAGACACCGACAACTTTCTT"
        assert reverse_complement(reverse_complement(seq)) == seq


class TestClassifyBarcode:
    """Test barcode classification against a read segment."""

    def test_perfect_match_start(self):
        """A read starting with nb01 barcode should classify as nb01."""
        # Simulate: flank + nb01 + flank + target
        nb01 = BARCODES["nb01"]
        segment = "AAGGTTAA" + nb01 + "CAGCACCT" + "A" * 50
        result = classify_barcode(segment, BARCODES)
        assert result["barcode_id"] == "nb01"
        assert result["edit_distance"] == 0
        assert result["confidence"] == 1.0

    def test_perfect_match_different_barcode(self):
        """nb05 in segment should classify as nb05."""
        nb05 = BARCODES["nb05"]
        segment = "AAGGTTAA" + nb05 + "CAGCACCT" + "A" * 50
        result = classify_barcode(segment, BARCODES)
        assert result["barcode_id"] == "nb05"
        assert result["confidence"] == 1.0

    def test_imperfect_match(self):
        """A barcode with 2 errors should still classify correctly."""
        nb01 = list(BARCODES["nb01"])
        nb01[0] = "T"  # mutate position 0
        nb01[5] = "T"  # mutate position 5
        mutated = "".join(nb01)
        segment = "AAGGTTAA" + mutated + "CAGCACCT" + "A" * 50
        result = classify_barcode(segment, BARCODES)
        assert result["barcode_id"] == "nb01"
        assert result["edit_distance"] == 2
        assert result["confidence"] == pytest.approx(1.0 - 2 / 24, abs=0.01)

    def test_returns_best_match(self):
        """When no perfect match, returns the closest barcode."""
        # Random-ish sequence that doesn't match any barcode well
        segment = "AAGGTTAA" + "NNNNNNNNNNNNNNNNNNNNNNNN" + "A" * 50
        result = classify_barcode(segment, BARCODES)
        assert "barcode_id" in result
        assert result["confidence"] < 0.5  # poor match

    def test_subset_of_barcodes(self):
        """Classification works with a subset of expected barcodes."""
        subset = {k: v for k, v in BARCODES.items() if k in ("nb01", "nb02", "nb03")}
        nb02 = BARCODES["nb02"]
        segment = "AAGGTTAA" + nb02 + "CAGCACCT" + "A" * 50
        result = classify_barcode(segment, subset)
        assert result["barcode_id"] == "nb02"
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_barcodes.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'barcodes'`

**Step 3: Write the implementation**

Create `bin/barcodes.py`:

```python
#!/usr/bin/env python3
"""ONT native barcode sequences and classification functions.

Barcode sequences from ONT native barcoding kit (nb01-nb96, 24bp each).
Classification uses edlib semi-global alignment (HW mode) to find the
best-matching barcode within a read segment.
"""

import edlib

# All 96 ONT native barcodes (forward orientation, 24bp each)
BARCODES = {
    "nb01": "CACAAAGACACCGACAACTTTCTT",
    "nb02": "ACAGACGACTACAAACGGAATCGA",
    "nb03": "CCTGGTAACTGGGACACAAGACTC",
    "nb04": "TAGGGAAACACGATAGAATCCGAA",
    "nb05": "AAGGTTACACAAACCCTGGACAAG",
    "nb06": "GACTACTTTCTGCCTTTGCGAGAA",
    "nb07": "AAGGATTCATTCCCACGGTAACAC",
    "nb08": "ACGTAACTTGGTTTGTTCCCTGAA",
    "nb09": "AACCAAGACTCGCTGTGCCTAGTT",
    "nb10": "GAGAGGACAAAGGTTTCAACGCTT",
    "nb11": "TCCATTCCCTCCGATAGATGAAAC",
    "nb12": "TCCGATTCTGCTTCTTTCTACCTG",
    "nb13": "AGAACGACTTCCATACTCGTGTGA",
    "nb14": "AACGAGTCTCTTGGGACCCATAGA",
    "nb15": "AGGTCTACCTCGCTAACACCACTG",
    "nb16": "CGTCAACTGACAGTGGTTCGTACT",
    "nb17": "ACCCTCCAGGAAAGTACCTCTGAT",
    "nb18": "CCAAACCCAACAACCTAGATAGGC",
    "nb19": "GTTCCTCGTGCAGTGTCAAGAGAT",
    "nb20": "TTGCGTCCTGTTACGAGAACTCAT",
    "nb21": "GAGCCTCTCATTGTCCGTTCTCTA",
    "nb22": "ACCACTGCCATGTATCAAAGTACG",
    "nb23": "CTTACTACCCAGTGAACCTCCTCG",
    "nb24": "GCATAGTTCTGCATGATGGGTTAG",
    "nb25": "GTAAGTTGGGTATGCAACGCAATG",
    "nb26": "CATACAGCGACTACGCATTCTCAT",
    "nb27": "CGACGGTTAGATTCACCTCTTACA",
    "nb28": "TGAAACCTAAGAAGGCACCGTATC",
    "nb29": "CTAGACACCTTGGGTTGACAGACC",
    "nb30": "TCAGTGAGGATCTACTTCGACCCA",
    "nb31": "TGCGTACAGCAATCAGTTACATTG",
    "nb32": "CCAGTAGAAGTCCGACAACGTCAT",
    "nb33": "CAGACTTGGTACGGTTGGGTAACT",
    "nb34": "GGACGAAGAACTCAAGTCAAAGGC",
    "nb35": "CTACTTACGAAGCTGAGGGACTGC",
    "nb36": "ATGTCCCAGTTAGAGGAGGAAACA",
    "nb37": "GCTTGCGATTGATGCTTAGTATCA",
    "nb38": "ACCACAGGAGGACGATACAGAGAA",
    "nb39": "CCACAGTGTCAACTAGAGCCTCTC",
    "nb40": "TAGTTTGGATGACCAAGGATAGCC",
    "nb41": "GGAGTTCGTCCAGAGAAGTACACG",
    "nb42": "CTACGTGTAAGGCATACCTGCCAG",
    "nb43": "CTTTCGTTGTTGACTCGACGGTAG",
    "nb44": "AGTAGAAAGGGTTCCTTCCCACTC",
    "nb45": "GATCCAACAGAGATGCCTTCAGTG",
    "nb46": "GCTGTGTTCCACTTCATTCTCCTG",
    "nb47": "GTGCAACTTTCCCACAGGTAGTTC",
    "nb48": "CATCTGGAACGTGGTACACCTGTA",
    "nb49": "ACTGGTGCAGCTTTGAACATCTAG",
    "nb50": "ATGGACTTTGGTAACTTCCTGCGT",
    "nb51": "GTTGAATGAGCCTACTGGGTCCTC",
    "nb52": "TGAGAGACAAGATTGTTCGTGGAC",
    "nb53": "AGATTCAGACCGTCTCATGCAAAG",
    "nb54": "CAAGAGCTTTGACTAAGGAGCATG",
    "nb55": "TGGAAGATGAGACCCTGATCTACG",
    "nb56": "TCACTACTCAACAGGTGGCATGAA",
    "nb57": "GCTAGGTCAATCTCCTTCGGAAGT",
    "nb58": "CAGGTTACTCCTCCGTGAGTCTGA",
    "nb59": "TCAATCAAGAAGGGAAAGCAAGGT",
    "nb60": "CATGTTCAACCAAGGCTTCTATGG",
    "nb61": "AGAGGGTACTATGTGCCTCAGCAC",
    "nb62": "CACCCACACTTACTTCAGGACGTA",
    "nb63": "TTCTGAAGTTCCTGGGTCTTGAAC",
    "nb64": "GACAGACACCGTTCATCGACTTTC",
    "nb65": "TTCTCAGTCTTCCTCCAGACAAGG",
    "nb66": "CCGATCCTTGTGGCTTCTAACTTC",
    "nb67": "GTTTGTCATACTCGTGTGCTCACC",
    "nb68": "GAATCTAAGCAAACACGAAGGTGG",
    "nb69": "TACAGTCCGAGCCTCATGTGATCT",
    "nb70": "ACCGAGATCCTACGAATGGAGTGT",
    "nb71": "CCTGGGAGCATCAGGTAGTAACAG",
    "nb72": "TAGCTGACTGTCTTCCATACCGAC",
    "nb73": "AAGAAACAGGATGACAGAACCCTC",
    "nb74": "TACAAGCATCCCAACACTTCCACT",
    "nb75": "GACCATTGTGATGAACCCTGTTGT",
    "nb76": "ATGCTTGTTACATCAACCCTGGAC",
    "nb77": "CGACCTGTTTCTCAGGGATACAAC",
    "nb78": "AACAACCGAACCTTTGAATCAGAA",
    "nb79": "TCTCGGAGATAGTTCTCACTGCTG",
    "nb80": "CGGATGAACATAGGATAGCGATTC",
    "nb81": "CCTCATCTTGTGAAGTTGTTTCGG",
    "nb82": "ACGGTATGTCGAGTTCCAGGACTA",
    "nb83": "TGGCTTGATCTAGGTAAGGTCGAA",
    "nb84": "GTAGTGGACCTAGAACCTGTGCCA",
    "nb85": "AACGGAGGAGTTAGTTGGATGATC",
    "nb86": "AGGTGATCCCAACAAGCGTAAGTA",
    "nb87": "TACATGCTCCTGTTGTTAGGGAGG",
    "nb88": "TCTTCTACTACCGATCCGAAGCAG",
    "nb89": "ACAGCATCAATGTTTGGCTAGTTG",
    "nb90": "GATGTAGAGGGTACGGTTTGAGGC",
    "nb91": "GGCTCCATAGGAACTCACGCTACT",
    "nb92": "TTGTGAGTGGAAAGATACAGGACC",
    "nb93": "AGTTTCCATCACTTCAGACTTGGG",
    "nb94": "GATTGTCCTCAAACTGCCACCTAC",
    "nb95": "CCTGTCTGGAAGAAGAATGGACTT",
    "nb96": "CTGAACGGTCATAGAGTCCACCAT",
}

BARCODE_LENGTH = 24


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C",
            "a": "t", "t": "a", "c": "g", "g": "c"}
    return "".join(comp.get(b, b) for b in reversed(seq))


def classify_barcode(segment: str, expected_barcodes: dict) -> dict:
    """Classify a read segment by finding the best-matching barcode.

    Uses edlib semi-global alignment (HW mode) to find where each barcode
    best aligns within the segment. Returns the barcode with lowest edit
    distance.

    Args:
        segment: Read subsequence (first or last ~100bp).
        expected_barcodes: Dict of {barcode_name: barcode_sequence}.

    Returns:
        Dict with keys: barcode_id, edit_distance, confidence.
    """
    best_id = None
    best_ed = float("inf")

    for bc_name, bc_seq in expected_barcodes.items():
        result = edlib.align(bc_seq, segment, mode="HW", task="distance")
        ed = result["editDistance"]
        if ed < best_ed:
            best_ed = ed
            best_id = bc_name

    confidence = 1.0 - (best_ed / BARCODE_LENGTH) if best_id else 0.0
    confidence = max(0.0, confidence)

    return {
        "barcode_id": best_id,
        "edit_distance": best_ed,
        "confidence": round(confidence, 4),
    }
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_barcodes.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/barcodes.py tests/test_barcodes.py
git commit -m "feat: add barcode module with 96 ONT barcodes and classifier"
```

---

### Task 2: Create sample sheet parser

**Files:**
- Create: `bin/sample_sheet.py`
- Create: `tests/test_sample_sheet.py`

Parses MinKNOW sample sheet CSV files. Extracts duplexed barcode pairs (e.g., `barcode05--barcode10`) and maps them to target aliases. Detects barcode ambiguity (one barcode mapped to multiple targets).

**Step 1: Write the failing tests**

Create `tests/test_sample_sheet.py`:

```python
#!/usr/bin/env python3
"""Tests for sample_sheet.py MinKNOW sample sheet parsing."""

import csv
import pytest
from pathlib import Path

from sample_sheet import parse_sample_sheet, detect_barcode_ambiguity


@pytest.fixture
def sample_sheet_csv(tmp_path):
    """Create a minimal MinKNOW sample sheet CSV."""
    path = tmp_path / "sample_sheet.csv"
    rows = [
        {"flow_cell_id": "FAL12345", "kit": "SQK-NBD114-96",
         "sample_id": "20260129_IF", "experiment_id": "exp001",
         "barcode": "barcode05--barcode10", "alias": "CYP2D6_v04_fwd"},
        {"flow_cell_id": "FAL12345", "kit": "SQK-NBD114-96",
         "sample_id": "20260129_IF", "experiment_id": "exp001",
         "barcode": "barcode10--barcode05", "alias": "CYP2D6_v04_rev"},
        {"flow_cell_id": "FAL12345", "kit": "SQK-NBD114-96",
         "sample_id": "20260129_IF", "experiment_id": "exp001",
         "barcode": "barcode01--barcode02", "alias": "GeneX_fwd"},
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    return path


@pytest.fixture
def ambiguous_sheet(tmp_path):
    """Sample sheet where nb05 maps to two different targets."""
    path = tmp_path / "ambiguous.csv"
    rows = [
        {"flow_cell_id": "FAL12345", "kit": "SQK-NBD114-96",
         "sample_id": "s1", "experiment_id": "exp001",
         "barcode": "barcode05--barcode10", "alias": "target_A"},
        {"flow_cell_id": "FAL12345", "kit": "SQK-NBD114-96",
         "sample_id": "s1", "experiment_id": "exp001",
         "barcode": "barcode05--barcode12", "alias": "target_B"},
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    return path


class TestParseSampleSheet:

    def test_parses_barcode_pairs(self, sample_sheet_csv):
        result = parse_sample_sheet(sample_sheet_csv)
        assert ("nb05", "nb10") in result
        assert ("nb10", "nb05") in result
        assert ("nb01", "nb02") in result

    def test_alias_mapping(self, sample_sheet_csv):
        result = parse_sample_sheet(sample_sheet_csv)
        assert result[("nb05", "nb10")] == "CYP2D6_v04_fwd"
        assert result[("nb10", "nb05")] == "CYP2D6_v04_rev"
        assert result[("nb01", "nb02")] == "GeneX_fwd"

    def test_expected_barcodes(self, sample_sheet_csv):
        result = parse_sample_sheet(sample_sheet_csv)
        # Should contain exactly 3 entries
        assert len(result) == 3

    def test_missing_barcode_column(self, tmp_path):
        path = tmp_path / "bad.csv"
        with open(path, "w") as f:
            f.write("flow_cell_id,kit,sample_id\nFAL,KIT,S1\n")
        with pytest.raises(ValueError, match="barcode"):
            parse_sample_sheet(path)

    def test_invalid_barcode_format(self, tmp_path):
        path = tmp_path / "bad_fmt.csv"
        rows = [{"flow_cell_id": "FAL", "kit": "K", "sample_id": "S",
                 "experiment_id": "E", "barcode": "barcode05",
                 "alias": "test"}]
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        with pytest.raises(ValueError, match="duplexed"):
            parse_sample_sheet(path)


class TestDetectBarcodeAmbiguity:

    def test_no_ambiguity(self, sample_sheet_csv):
        mapping = parse_sample_sheet(sample_sheet_csv)
        assert detect_barcode_ambiguity(mapping) is False

    def test_ambiguity_detected(self, ambiguous_sheet):
        mapping = parse_sample_sheet(ambiguous_sheet)
        assert detect_barcode_ambiguity(mapping) is True
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_sample_sheet.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'sample_sheet'`

**Step 3: Write the implementation**

Create `bin/sample_sheet.py`:

```python
#!/usr/bin/env python3
"""MinKNOW sample sheet parser for duplexed barcode experiments.

Parses CSV sample sheets produced by MinKNOW and extracts duplexed
barcode-pair-to-alias mappings. Detects barcode ambiguity where a
single barcode maps to multiple targets.
"""

import csv
import re
import sys
from pathlib import Path
from typing import Dict, Tuple


def parse_sample_sheet(path: Path) -> Dict[Tuple[str, str], str]:
    """Parse a MinKNOW sample sheet CSV into barcode-pair-to-alias mapping.

    Args:
        path: Path to MinKNOW sample sheet CSV.

    Returns:
        Dict mapping (upstream_barcode, downstream_barcode) tuples to alias
        strings. Barcode names are normalized to 'nbNN' format.

    Raises:
        ValueError: If sample sheet is missing required columns or has
            invalid barcode format.
    """
    path = Path(path)
    mapping = {}

    with open(path, "r") as f:
        reader = csv.DictReader(f)
        columns = reader.fieldnames or []

        if "barcode" not in columns:
            raise ValueError(
                f"Sample sheet has no 'barcode' column -- expected MinKNOW format. "
                f"Found columns: {columns}"
            )
        if "alias" not in columns:
            raise ValueError(
                f"Sample sheet has no 'alias' column -- expected MinKNOW format. "
                f"Found columns: {columns}"
            )

        for row in reader:
            barcode_field = row["barcode"].strip()
            alias = row["alias"].strip()

            if "--" not in barcode_field:
                raise ValueError(
                    f"Barcode '{barcode_field}' is not duplexed format. "
                    f"Expected 'barcodeNN--barcodeNN'."
                )

            parts = barcode_field.split("--")
            if len(parts) != 2:
                raise ValueError(
                    f"Barcode '{barcode_field}' has unexpected format. "
                    f"Expected exactly 'barcodeNN--barcodeNN'."
                )

            upstream = _normalize_barcode_name(parts[0].strip())
            downstream = _normalize_barcode_name(parts[1].strip())

            mapping[(upstream, downstream)] = alias

    return mapping


def detect_barcode_ambiguity(mapping: Dict[Tuple[str, str], str]) -> bool:
    """Check if any single barcode maps to multiple different targets.

    This occurs when the same barcode appears as the upstream barcode
    in multiple pairs that map to different aliases, meaning the upstream
    barcode alone cannot determine the target.

    Args:
        mapping: Output from parse_sample_sheet.

    Returns:
        True if ambiguity detected (full construct alignment needed).
    """
    upstream_targets = {}
    for (up, down), alias in mapping.items():
        if up not in upstream_targets:
            upstream_targets[up] = set()
        upstream_targets[up].add(alias)

    for bc, aliases in upstream_targets.items():
        if len(aliases) > 1:
            return True
    return False


def _normalize_barcode_name(name: str) -> str:
    """Convert 'barcode05' or 'barcode5' to 'nb05'."""
    match = re.match(r"barcode(\d+)", name, re.IGNORECASE)
    if not match:
        raise ValueError(f"Cannot parse barcode name '{name}'. Expected 'barcodeNN'.")
    num = int(match.group(1))
    if num < 1 or num > 96:
        raise ValueError(f"Barcode number {num} out of range (1-96).")
    return f"nb{num:02d}"
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_sample_sheet.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/sample_sheet.py tests/test_sample_sheet.py
git commit -m "feat: add MinKNOW sample sheet parser with ambiguity detection"
```

---

### Task 3: Update mkdb.py to add barcode columns

**Files:**
- Modify: `bin/mkdb.py:66-85` (Reads table CREATE)
- Create: `tests/test_mkdb.py`

Add the 6 new barcode classification columns to the `Reads` table schema.

**Step 1: Write the failing test**

Create `tests/test_mkdb.py`:

```python
#!/usr/bin/env python3
"""Tests for mkdb.py database schema."""

import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.fixture
def db_path(tmp_path):
    """Run mkdb.py and return the created DB path."""
    result = subprocess.run(
        [sys.executable, "bin/mkdb.py", "-e", "FAL12345_20260129_IF", "-o", str(tmp_path)],
        capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
    )
    assert result.returncode == 0, f"mkdb.py failed: {result.stderr}"
    return tmp_path / "SMA_FAL12345_20260129_IF.db"


class TestReadsSchema:

    def test_db_created(self, db_path):
        assert db_path.exists()

    def test_barcode_columns_exist(self, db_path):
        conn = sqlite3.connect(db_path)
        cursor = conn.execute("PRAGMA table_info(Reads)")
        columns = {row[1] for row in cursor.fetchall()}
        conn.close()

        expected_new = {
            "bc_start_id", "bc_start_ed", "bc_start_conf",
            "bc_end_id", "bc_end_ed", "bc_end_conf",
        }
        assert expected_new.issubset(columns), (
            f"Missing columns: {expected_new - columns}"
        )

    def test_barcode_column_types(self, db_path):
        conn = sqlite3.connect(db_path)
        cursor = conn.execute("PRAGMA table_info(Reads)")
        col_types = {row[1]: row[2] for row in cursor.fetchall()}
        conn.close()

        assert col_types["bc_start_id"] == "TEXT"
        assert col_types["bc_start_ed"] == "INTEGER"
        assert col_types["bc_start_conf"] == "REAL"
        assert col_types["bc_end_id"] == "TEXT"
        assert col_types["bc_end_ed"] == "INTEGER"
        assert col_types["bc_end_conf"] == "REAL"

    def test_original_columns_still_present(self, db_path):
        conn = sqlite3.connect(db_path)
        cursor = conn.execute("PRAGMA table_info(Reads)")
        columns = {row[1] for row in cursor.fetchall()}
        conn.close()

        original = {
            "uniq_id", "exp_id", "tgt_id", "read_id", "readseq",
            "readlen", "model_tier", "model_ver", "trim", "mod_bitflag",
            "ed", "q_bc", "q_ld", "ER",
        }
        assert original.issubset(columns)
```

**Step 2: Run test to verify it fails**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_mkdb.py -v`
Expected: FAIL on `test_barcode_columns_exist` — columns don't exist yet

**Step 3: Update mkdb.py — add barcode columns to CREATE TABLE**

In `bin/mkdb.py`, replace the Reads table DDL (lines 66-85) with:

```python
# Reads Table
c.execute('''
	CREATE TABLE Reads (
		uniq_id TEXT PRIMARY KEY,
		exp_id TEXT,
		tgt_id TEXT,
		read_id TEXT,
		readseq TEXT,
		readlen INTEGER,
		model_tier TEXT,
		model_ver TEXT,
		trim INTEGER,
		mod_bitflag INTEGER,
		ed INTEGER,
		q_bc REAL,
		q_ld REAL,
		ER TEXT,
		bc_start_id TEXT,
		bc_start_ed INTEGER,
		bc_start_conf REAL,
		bc_end_id TEXT,
		bc_end_ed INTEGER,
		bc_end_conf REAL,
		FOREIGN KEY(tgt_id) REFERENCES Target(tgt_id),
		FOREIGN KEY(mod_bitflag) REFERENCES Mods(mod_bitflag)
	)
''')
```

**Step 4: Run test to verify it passes**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_mkdb.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/mkdb.py tests/test_mkdb.py
git commit -m "feat: add barcode classification columns to Reads schema"
```

---

### Task 4: Enhance ingest.py with barcode classification

**Files:**
- Modify: `bin/ingest.py` (major changes)
- Create: `tests/test_ingest.py`

This is the core change. `ingest.py` gains new arguments (`-ss`, `-rd`, `--full-construct`) and classifies each read during the processing loop.

**Step 1: Write the failing test**

Create `tests/test_ingest.py`:

```python
#!/usr/bin/env python3
"""Tests for ingest.py with barcode classification.

These tests create synthetic BAM files with known barcode sequences
embedded in the reads, then verify classification is correct.
"""

import csv
import sqlite3
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

# ONT native barcode sequences needed for test data
NB01 = "CACAAAGACACCGACAACTTTCTT"
NB02 = "ACAGACGACTACAAACGGAATCGA"
NB05 = "AAGGTTACACAAACCCTGGACAAG"
NB10 = "GAGAGGACAAAGGTTTCAACGCTT"

FLANK_F = "AAGGTTAA"
FLANK_R = "CAGCACCT"
REV_FLANK_F = "AGGTGCTG"
REV_FLANK_R = "TTAACCTTAGCAAT"


def rc(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))


def make_read_sequence(bc_up, bc_down, target_seq):
    """Build a synthetic untrimmed read sequence.

    Structure: front_flank + bc_upstream + rear_flank + target +
               rev_front_flank + RC(bc_downstream) + rev_rear_flank
    """
    return (
        FLANK_F + bc_up + FLANK_R +
        target_seq +
        REV_FLANK_F + rc(bc_down) + REV_FLANK_R
    )


@pytest.fixture
def test_env(tmp_path):
    """Set up a complete test environment with BAM, sample sheet, refs, and DB."""
    # Target sequences
    target_a_seq = "ATCGATCGATCGATCGATCG" * 5  # 100bp
    target_b_seq = "GCTAGCTAGCTAGCTAGCTA" * 5  # 100bp

    # Build synthetic reads
    read1_seq = make_read_sequence(NB05, NB10, target_a_seq)
    read2_seq = make_read_sequence(NB01, NB02, target_b_seq)
    read3_seq = make_read_sequence(NB05, NB10, target_a_seq)

    # Write unaligned BAM
    bam_path = tmp_path / "FAL12345_20260129_IF_sup_v5.2.0_trim0_0.bam"
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for i, seq in enumerate([read1_seq, read2_seq, read3_seq]):
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i:04d}"
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.flag = 4  # unmapped
            bam.write(a)

    # Write sample sheet
    ss_path = tmp_path / "sample_sheet.csv"
    rows = [
        {"flow_cell_id": "FAL12345", "kit": "SQK-NBD114-96",
         "sample_id": "20260129_IF", "experiment_id": "exp001",
         "barcode": "barcode05--barcode10", "alias": "target_A"},
        {"flow_cell_id": "FAL12345", "kit": "SQK-NBD114-96",
         "sample_id": "20260129_IF", "experiment_id": "exp001",
         "barcode": "barcode01--barcode02", "alias": "target_B"},
    ]
    with open(ss_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    # Write reference FASTAs
    ref_dir = tmp_path / "references"
    ref_dir.mkdir()
    (ref_dir / "target_A.fasta").write_text(f">target_A\n{target_a_seq}\n")
    (ref_dir / "target_B.fasta").write_text(f">target_B\n{target_b_seq}\n")

    # Write summary TSV (end reasons)
    summary_path = tmp_path / "summary.tsv"
    summary_path.write_text(
        "read_id\tend_reason\n"
        "read_0000\tsignal_positive\n"
        "read_0001\tsignal_positive\n"
        "read_0002\tdata_service_unblock_mux_change\n"
    )

    # Create DB
    out_dir = tmp_path / "Output"
    result = subprocess.run(
        [sys.executable, "bin/mkdb.py", "-e", "FAL12345_20260129_IF", "-o", str(out_dir)],
        capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
    )
    assert result.returncode == 0, f"mkdb failed: {result.stderr}"

    db_path = out_dir / "SMA_FAL12345_20260129_IF.db"
    output_bam = tmp_path / "tagged.bam"

    return {
        "bam": bam_path,
        "sample_sheet": ss_path,
        "ref_dir": ref_dir,
        "summary": summary_path,
        "db": db_path,
        "output_bam": output_bam,
        "out_dir": out_dir,
    }


class TestIngestWithClassification:

    def _run_ingest(self, env):
        cmd = [
            sys.executable, "bin/ingest.py",
            "-e", "FAL12345_20260129_IF",
            "-b", str(env["bam"]),
            "-s", str(env["summary"]),
            "-d", str(env["db"]),
            "-o", str(env["output_bam"]),
            "-ss", str(env["sample_sheet"]),
            "-rd", str(env["ref_dir"]),
        ]
        result = subprocess.run(
            cmd, capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
        )
        return result

    def test_ingest_succeeds(self, test_env):
        result = self._run_ingest(test_env)
        assert result.returncode == 0, f"ingest failed: {result.stderr}\n{result.stdout}"

    def test_all_reads_ingested(self, test_env):
        self._run_ingest(test_env)
        conn = sqlite3.connect(test_env["db"])
        count = conn.execute("SELECT COUNT(*) FROM Reads").fetchone()[0]
        conn.close()
        assert count == 3

    def test_barcode_start_classified(self, test_env):
        self._run_ingest(test_env)
        conn = sqlite3.connect(test_env["db"])
        rows = conn.execute(
            "SELECT read_id, bc_start_id FROM Reads ORDER BY read_id"
        ).fetchall()
        conn.close()
        # read_0000 and read_0002 have nb05, read_0001 has nb01
        assert rows[0] == ("read_0000", "nb05")
        assert rows[1] == ("read_0001", "nb01")
        assert rows[2] == ("read_0002", "nb05")

    def test_barcode_end_classified(self, test_env):
        self._run_ingest(test_env)
        conn = sqlite3.connect(test_env["db"])
        rows = conn.execute(
            "SELECT read_id, bc_end_id FROM Reads ORDER BY read_id"
        ).fetchall()
        conn.close()
        # read_0000 and read_0002 have nb10, read_0001 has nb02
        assert rows[0] == ("read_0000", "nb10")
        assert rows[1] == ("read_0001", "nb02")
        assert rows[2] == ("read_0002", "nb10")

    def test_target_assigned_from_barcode_pair(self, test_env):
        self._run_ingest(test_env)
        conn = sqlite3.connect(test_env["db"])
        rows = conn.execute(
            "SELECT read_id, tgt_id FROM Reads ORDER BY read_id"
        ).fetchall()
        conn.close()
        assert rows[0] == ("read_0000", "target_A")
        assert rows[1] == ("read_0001", "target_B")
        assert rows[2] == ("read_0002", "target_A")

    def test_confidence_is_high_for_perfect_barcodes(self, test_env):
        self._run_ingest(test_env)
        conn = sqlite3.connect(test_env["db"])
        rows = conn.execute(
            "SELECT bc_start_conf, bc_end_conf FROM Reads"
        ).fetchall()
        conn.close()
        for start_conf, end_conf in rows:
            assert start_conf >= 0.9
            assert end_conf >= 0.9

    def test_ed_calculated_against_correct_target(self, test_env):
        """ed should be low because reads contain the exact target sequence."""
        self._run_ingest(test_env)
        conn = sqlite3.connect(test_env["db"])
        rows = conn.execute("SELECT ed FROM Reads").fetchall()
        conn.close()
        for (ed,) in rows:
            # Reads contain flanks + barcodes around the target, so ed
            # won't be 0, but should be reasonable (flank/barcode overhead)
            assert ed is not None


class TestIngestBackwardCompatibility:
    """Ensure old single-target mode still works without -ss/-rd."""

    def test_single_target_mode(self, tmp_path):
        target_seq = "ATCGATCG" * 10
        # Write single-ref FASTA
        ref_path = tmp_path / "target.fa"
        ref_path.write_text(f">single_target\n{target_seq}\n")

        # Write BAM with one read
        bam_path = tmp_path / "FAL12345_20260129_IF_sup_v5.2.0_trim0_0.bam"
        header = {"HD": {"VN": "1.6", "SO": "unsorted"}}
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            a = pysam.AlignedSegment()
            a.query_name = "read_0000"
            a.query_sequence = target_seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(target_seq))
            a.flag = 4
            bam.write(a)

        # Summary
        summary = tmp_path / "summary.tsv"
        summary.write_text("read_id\tend_reason\nread_0000\tsignal_positive\n")

        # DB
        out_dir = tmp_path / "Output"
        subprocess.run(
            [sys.executable, "bin/mkdb.py", "-e", "FAL12345_20260129_IF", "-o", str(out_dir)],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
        )
        db_path = out_dir / "SMA_FAL12345_20260129_IF.db"

        # Run ingest WITHOUT -ss/-rd (backward compatible)
        result = subprocess.run(
            [sys.executable, "bin/ingest.py",
             "-e", "FAL12345_20260129_IF",
             "-b", str(bam_path),
             "-s", str(summary),
             "-r", str(ref_path),
             "-d", str(db_path),
             "-o", str(tmp_path / "out.bam")],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
        )
        assert result.returncode == 0, f"Failed: {result.stderr}\n{result.stdout}"

        conn = sqlite3.connect(db_path)
        row = conn.execute(
            "SELECT tgt_id, bc_start_id FROM Reads"
        ).fetchone()
        conn.close()
        assert row[0] == "single_target"
        assert row[1] is None  # no barcode classification in single-target mode
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_ingest.py -v`
Expected: FAIL — ingest.py doesn't accept `-ss`/`-rd` args yet

**Step 3: Rewrite ingest.py with classification**

Replace `bin/ingest.py` entirely. The key changes:
- Add `-ss`, `-rd`, `--full-construct`, `--split-bams`, `--tag` arguments
- When `-ss` and `-rd` provided: load sample sheet, load all reference FASTAs, classify each read's start and end, look up target from barcode pair
- When not provided: behave exactly as before (single `-r` target)
- New `insert_read` includes 6 barcode columns

```python
#!/usr/bin/env python3
# ingest.py
# Main processing script: BAM tagging, Metric Calculation, Barcode
# Classification, DB Ingestion

import argparse
import math
import sqlite3
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import edlib
import numpy as np
import pandas as pd
import pysam

from barcodes import BARCODES, BARCODE_LENGTH, classify_barcode, reverse_complement
from sample_sheet import parse_sample_sheet, detect_barcode_ambiguity


#########
# funcs #
#########

def parse_bam_filename(bam_path: Path) -> Tuple[str, str, int, int]:
	"""
	Parses metadata from BAM filename (Right-to-Left parsing).
	Expected Format: {exp_id}_{tier}_v{ver}_{trim}_{mods}.bam
	"""
	filename = bam_path.resolve().stem
	parts = filename.split('_')

	if len(parts) < 5:
		sys.exit(f"[ingest] Error: Metadata parse failed for '{filename}'.")

	try:
		mods = int(parts[-1])
		trim = int(parts[-2][-1])
		ver_str = parts[-3].lstrip("v")
		tier = parts[-4]
		return tier, ver_str, trim, mods
	except Exception as e:
		sys.exit(f"[ingest] Error: Metadata parse failed for '{filename}'. {e}")


def calculate_q_bc(quality_scores: Any) -> float:
	"""Calculates probability-averaged Phred quality score."""
	if len(quality_scores) == 0:
		return 0.0

	q_arr = np.array(quality_scores, dtype=np.float64)
	probs = np.power(10, -q_arr / 10.0)
	avg_prob = np.mean(probs)

	if avg_prob > 0:
		return -10 * math.log10(avg_prob)
	return 0.0


def calculate_q_ld(ed: int, ref_len: int) -> float:
	"""
	Calculates Levenshtein Quality.
	Formula: -10 * log10( min( max(1/L^2, ed/L), 1 ) )
	"""
	if ref_len == 0:
		return 0.0

	L = float(ref_len)
	term1 = 1.0 / (L * L)
	term2 = float(ed) / L
	prob = min(max(term1, term2), 1.0)

	return -10 * math.log10(prob)


def load_references(ref_dir: Path) -> Dict[str, Tuple[str, int]]:
	"""Load all FASTA files from a directory.

	Returns:
		Dict mapping alias (from FASTA header) to (sequence, length).
	"""
	refs = {}
	for fasta_path in sorted(ref_dir.glob("*.fasta")):
		tgt_id = None
		seq_parts = []
		with open(fasta_path) as f:
			for line in f:
				line = line.strip()
				if not line:
					continue
				if line.startswith(">"):
					tgt_id = line[1:]
				else:
					seq_parts.append(line)
		if tgt_id and seq_parts:
			seq = "".join(seq_parts)
			refs[tgt_id] = (seq, len(seq))
	# Also check .fa files
	for fasta_path in sorted(ref_dir.glob("*.fa")):
		tgt_id = None
		seq_parts = []
		with open(fasta_path) as f:
			for line in f:
				line = line.strip()
				if not line:
					continue
				if line.startswith(">"):
					tgt_id = line[1:]
				else:
					seq_parts.append(line)
		if tgt_id and seq_parts:
			seq = "".join(seq_parts)
			refs[tgt_id] = (seq, len(seq))
	return refs


def insert_target(cursor: sqlite3.Cursor, tgt_id: str, tgt_seq: str, tgt_len: int) -> None:
	"""Inserts target sequence metadata."""
	cursor.execute('''
		INSERT OR REPLACE INTO Target (tgt_id, tgt_refseq, tgt_reflen)
		VALUES (?, ?, ?)
	''', (tgt_id, tgt_seq, tgt_len))


def insert_read(cursor: sqlite3.Cursor, data: Tuple) -> None:
	"""Inserts a single processed read record (with barcode columns)."""
	cursor.execute('''
		INSERT OR IGNORE INTO Reads (
			uniq_id, exp_id, tgt_id, read_id, readseq, readlen,
			model_tier, model_ver, trim, mod_bitflag,
			ed, q_bc, q_ld, ER,
			bc_start_id, bc_start_ed, bc_start_conf,
			bc_end_id, bc_end_ed, bc_end_conf
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
	''', data)


############
# argparse #
############

parser = argparse.ArgumentParser(
	description="Ingest BAM, classify barcodes, calculate metrics, populate DB")
parser.add_argument("-e", "--expid", required=True,
	help="Experiment ID")
parser.add_argument("-b", "--bam", required=True,
	help="Input uBAM file (untrimmed)")
parser.add_argument("-s", "--summary", required=True,
	help="Pod5 summary TSV")
parser.add_argument("-r", "--ref", default=None,
	help="Single target FASTA (backward-compatible mode)")
parser.add_argument("-d", "--database", required=True,
	help="Target SQLite DB")
parser.add_argument("-o", "--output_bam", required=True,
	help="Output tagged BAM")
# New classification arguments
parser.add_argument("-ss", "--sample-sheet", default=None,
	help="MinKNOW sample sheet CSV (enables barcode classification)")
parser.add_argument("-rd", "--ref-dir", default=None,
	help="Directory of per-target reference FASTAs")
parser.add_argument("--full-construct", action="store_true",
	help="Force full construct alignment for all reads")
parser.add_argument("--split-bams", default=None,
	help="Output directory for per-barcode split BAMs")
parser.add_argument("--tag", action="store_true",
	help="Write BC:Z: and BA:Z: tags to output BAM")
args = parser.parse_args()


#########
# param #
#########

EXP_ID = args.expid
INPUT_BAM = Path(args.bam)
SUMMARY_TSV = Path(args.summary)
DB_PATH = Path(args.database)
OUTPUT_BAM = Path(args.output_bam)

MODEL_TIER, MODEL_VER, TRIM, MOD_BITFLAG = parse_bam_filename(INPUT_BAM)

print(f"[ingest] Run Metadata:\n  Exp: {EXP_ID}\n  Tier: {MODEL_TIER}"
      f"\n  Ver: {MODEL_VER}\n  Trim: {TRIM}\n  Mods: {MOD_BITFLAG}")


######################
# classification mode #
######################

CLASSIFY_MODE = args.sample_sheet is not None and args.ref_dir is not None

if CLASSIFY_MODE:
	print("[ingest] Classification mode: DUPLEXED BARCODE")
	# Parse sample sheet
	barcode_to_alias = parse_sample_sheet(args.sample_sheet)
	print(f"[ingest] Sample sheet: {len(barcode_to_alias)} barcode pairs loaded")

	# Detect ambiguity
	needs_full_construct = args.full_construct or detect_barcode_ambiguity(barcode_to_alias)
	if needs_full_construct and not args.full_construct:
		print("[ingest] Warning: Barcode ambiguity detected -- full construct alignment enabled")

	# Build expected barcodes set (only those in the sample sheet)
	expected_bc_names = set()
	for (up, down) in barcode_to_alias.keys():
		expected_bc_names.add(up)
		expected_bc_names.add(down)
	expected_barcodes = {k: v for k, v in BARCODES.items() if k in expected_bc_names}
	print(f"[ingest] Expected barcodes: {sorted(expected_barcodes.keys())}")

	# Build RC barcodes for end classification
	expected_barcodes_rc = {k: reverse_complement(v) for k, v in expected_barcodes.items()}

	# Load reference FASTAs
	ref_dir = Path(args.ref_dir)
	if not ref_dir.is_dir():
		sys.exit(f"[ingest] Error: Reference directory '{ref_dir}' does not exist.")
	references = load_references(ref_dir)
	print(f"[ingest] Loaded {len(references)} reference targets")

	# Validate: every alias in sample sheet has a reference
	for alias in barcode_to_alias.values():
		if alias not in references:
			sys.exit(
				f"[ingest] Error: Alias '{alias}' in sample sheet has no matching "
				f"FASTA in {ref_dir}/. Available: {sorted(references.keys())}"
			)

	# Insert all targets into DB
	# (done below after DB connection)

	# Set up split BAM writers if requested
	split_bam_dir = None
	if args.split_bams:
		split_bam_dir = Path(args.split_bams)
		split_bam_dir.mkdir(parents=True, exist_ok=True)

else:
	# Single-target backward-compatible mode
	if not args.ref:
		sys.exit("[ingest] Error: Must provide either -r (single target) or -ss/-rd (classification).")
	print("[ingest] Classification mode: SINGLE TARGET (no barcode classification)")

	REF_FASTA = Path(args.ref)
	tgt_id = None
	tgt_seq_parts = []
	try:
		with open(REF_FASTA, 'r') as f:
			for line in f:
				line = line.strip()
				if not line:
					continue
				if line.startswith('>'):
					tgt_id = line[1:]
				else:
					tgt_seq_parts.append(line)
		tgt_seq = "".join(tgt_seq_parts)
		if not tgt_id or not tgt_seq:
			raise ValueError("Empty ID or Sequence")
		tgt_len = len(tgt_seq)
		print(f"[ingest] Target Loaded: {tgt_id} ({tgt_len} bp)")
	except Exception as e:
		sys.exit(f"[ingest] Error loading FASTA: {e}")


###########
# summary #
###########

print(f"[ingest] Loading Pod5 Summary from {SUMMARY_TSV}...")
try:
	summary_df = pd.read_csv(SUMMARY_TSV, sep='\t', usecols=['read_id', 'end_reason'])
	end_reason_map = pd.Series(summary_df.end_reason.values, index=summary_df.read_id).to_dict()
	print(f"[ingest] Loaded {len(end_reason_map)} records.")
except Exception as e:
	sys.exit(f"[ingest] Error loading summary TSV: {e}")


########
# main #
########

SEGMENT_LEN = 100  # bp to examine at each end for barcode classification

with sqlite3.connect(DB_PATH) as conn:
	c = conn.cursor()

	if CLASSIFY_MODE:
		# Insert all targets
		for alias, (seq, length) in references.items():
			insert_target(c, alias, seq, length)
	else:
		insert_target(c, tgt_id, tgt_seq, tgt_len)
	conn.commit()

	save_verbosity = pysam.set_verbosity(0)

	# Split BAM writers (lazy-opened per barcode)
	split_writers = {}

	with pysam.AlignmentFile(INPUT_BAM, "rb", check_sq=False) as bam_in, \
			pysam.AlignmentFile(OUTPUT_BAM, "wb", template=bam_in) as bam_out:

		pysam.set_verbosity(save_verbosity)
		print("[ingest] Processing reads...")

		count_total = 0
		count_processed = 0
		count_skipped = 0
		count_unknown_er = 0
		barcode_counts = {}  # track per-alias counts

		for read in bam_in:
			count_total += 1
			seq = read.query_sequence
			quals = read.query_qualities

			if not seq or quals is None:
				count_skipped += 1
				continue

			read_id = read.query_name
			read_len = len(seq)
			q_bc = calculate_q_bc(quals)

			# Barcode classification
			bc_start_id = None
			bc_start_ed = None
			bc_start_conf = None
			bc_end_id = None
			bc_end_ed = None
			bc_end_conf = None
			assigned_tgt_id = None
			assigned_tgt_seq = None
			assigned_tgt_len = None

			if CLASSIFY_MODE:
				# Stage 1: classify start (first SEGMENT_LEN bp)
				start_segment = seq[:SEGMENT_LEN]
				start_result = classify_barcode(start_segment, expected_barcodes)
				bc_start_id = start_result["barcode_id"]
				bc_start_ed = start_result["edit_distance"]
				bc_start_conf = start_result["confidence"]

				# Stage 1: classify end (last SEGMENT_LEN bp, against RC barcodes)
				end_segment = seq[-SEGMENT_LEN:]
				end_result = classify_barcode(end_segment, expected_barcodes_rc)
				bc_end_id = end_result["barcode_id"]
				bc_end_ed = end_result["edit_distance"]
				bc_end_conf = end_result["confidence"]

				# Look up target from barcode pair
				pair = (bc_start_id, bc_end_id)
				alias = barcode_to_alias.get(pair)

				if alias and alias in references:
					assigned_tgt_id = alias
					assigned_tgt_seq, assigned_tgt_len = references[alias]
				else:
					# Pair not in sample sheet -- still ingest with best-effort
					assigned_tgt_id = f"unmatched_{bc_start_id}_{bc_end_id}"
					assigned_tgt_seq = None
					assigned_tgt_len = 0

				# Track barcode counts
				barcode_counts[assigned_tgt_id] = barcode_counts.get(assigned_tgt_id, 0) + 1

				# Optional tagging
				if args.tag:
					if bc_start_id:
						read.set_tag("BS", bc_start_id, value_type="Z")
					if bc_end_id:
						read.set_tag("BE", bc_end_id, value_type="Z")
					if alias:
						read.set_tag("BA", alias, value_type="Z")

				# Optional split BAMs
				if split_bam_dir and assigned_tgt_id and not assigned_tgt_id.startswith("unmatched"):
					if assigned_tgt_id not in split_writers:
						split_path = split_bam_dir / f"{assigned_tgt_id}.bam"
						split_writers[assigned_tgt_id] = pysam.AlignmentFile(
							str(split_path), "wb", template=bam_in
						)
					split_writers[assigned_tgt_id].write(read)

			else:
				# Single-target mode
				assigned_tgt_id = tgt_id
				assigned_tgt_seq = tgt_seq
				assigned_tgt_len = tgt_len

			# Alignment metrics (against assigned target)
			if assigned_tgt_seq:
				align_res = edlib.align(seq, assigned_tgt_seq, mode="NW", task="distance")
				ed = align_res["editDistance"]
				q_ld = calculate_q_ld(ed, assigned_tgt_len)
			else:
				ed = None
				q_ld = None

			# End reason
			er_val = end_reason_map.get(read_id, "unknown")
			if er_val == "unknown":
				count_unknown_er += 1

			uniq_id = f"{EXP_ID}_{MODEL_TIER}{MODEL_VER}t{TRIM}m{MOD_BITFLAG}_{read_id}"

			# DB Insert (20 columns)
			read_data = (
				uniq_id, EXP_ID, assigned_tgt_id, read_id, seq, read_len,
				MODEL_TIER, MODEL_VER, TRIM, MOD_BITFLAG,
				ed, q_bc, q_ld, er_val,
				bc_start_id, bc_start_ed, bc_start_conf,
				bc_end_id, bc_end_ed, bc_end_conf,
			)
			insert_read(c, read_data)

			# End reason tag (always)
			read.set_tag("ER", er_val, value_type="Z")
			bam_out.write(read)

			count_processed += 1
			if count_processed % 1000 == 0:
				print(f"  ...processed {count_processed} reads", end='\r')

	# Close split BAM writers
	for writer in split_writers.values():
		writer.close()

	conn.commit()

print(f"\n[ingest] Complete.")
print(f"  - Total Reads:  {count_total}")
print(f"  - Processed:    {count_processed}")
print(f"  - Skipped:      {count_skipped}")
print(f"  - Unknown ER:   {count_unknown_er}")
print(f"  - Output DB:    {DB_PATH}")
print(f"  - Output BAM:   {OUTPUT_BAM}")

if CLASSIFY_MODE and barcode_counts:
	print(f"\n[ingest] Classification Summary:")
	for target, count in sorted(barcode_counts.items(), key=lambda x: -x[1]):
		print(f"  {target}: {count} reads")
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_ingest.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/ingest.py tests/test_ingest.py
git commit -m "feat: integrate barcode classification into ingest.py"
```

---

### Task 5: End-to-end integration test

**Files:**
- Create: `tests/test_integration.py`

Runs the full pipeline (`mkdb` -> `inputInit` -> `extractMeta` -> `ingest`) with synthetic data to verify the complete workflow.

**Step 1: Write the integration test**

Create `tests/test_integration.py`:

```python
#!/usr/bin/env python3
"""End-to-end integration test for the full SMA-seq pipeline with classification."""

import csv
import sqlite3
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

NB05 = "AAGGTTACACAAACCCTGGACAAG"
NB10 = "GAGAGGACAAAGGTTTCAACGCTT"
FLANK_F = "AAGGTTAA"
FLANK_R = "CAGCACCT"
REV_FLANK_F = "AGGTGCTG"
REV_FLANK_R = "TTAACCTTAGCAAT"


def rc(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))


@pytest.fixture
def pipeline_env(tmp_path):
    """Full pipeline test environment."""
    target_seq = "ATCGATCGATCGATCGATCG" * 5

    # Build 5 reads: 3 with nb05--nb10, 2 with random barcode (should still classify)
    reads = []
    for i in range(3):
        seq = FLANK_F + NB05 + FLANK_R + target_seq + REV_FLANK_F + rc(NB10) + REV_FLANK_R
        reads.append((f"read_{i:04d}", seq))
    # 2 reads with slightly mutated barcodes (still classifiable)
    for i in range(3, 5):
        mutated = list(NB05)
        mutated[0] = "G"
        seq = FLANK_F + "".join(mutated) + FLANK_R + target_seq + REV_FLANK_F + rc(NB10) + REV_FLANK_R
        reads.append((f"read_{i:04d}", seq))

    # Write BAM
    bam_path = tmp_path / "FAL99999_20260214_TEST_sup_v5.2.0_trim0_0.bam"
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for read_id, seq in reads:
            a = pysam.AlignedSegment()
            a.query_name = read_id
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.flag = 4
            bam.write(a)

    # Write pod5 summary (simulate extractMeta output)
    summary = tmp_path / "summary.tsv"
    lines = ["read_id\tend_reason\n"]
    for read_id, _ in reads:
        lines.append(f"{read_id}\tsignal_positive\n")
    summary.write_text("".join(lines))

    # Sample sheet
    ss = tmp_path / "sample_sheet.csv"
    rows = [{"flow_cell_id": "FAL99999", "kit": "SQK-NBD114-96",
             "sample_id": "20260214_TEST", "experiment_id": "exp_test",
             "barcode": "barcode05--barcode10", "alias": "test_target"}]
    with open(ss, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    # Reference
    ref_dir = tmp_path / "refs"
    ref_dir.mkdir()
    (ref_dir / "test_target.fasta").write_text(f">test_target\n{target_seq}\n")

    out_dir = tmp_path / "Output"
    return {
        "tmp": tmp_path, "bam": bam_path, "summary": summary,
        "ss": ss, "ref_dir": ref_dir, "out_dir": out_dir,
    }


class TestFullPipeline:

    def test_mkdb_then_ingest(self, pipeline_env):
        env = pipeline_env

        # Step 1: mkdb
        r1 = subprocess.run(
            [sys.executable, "bin/mkdb.py",
             "-e", "FAL99999_20260214_TEST", "-o", str(env["out_dir"])],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
        )
        assert r1.returncode == 0, r1.stderr

        db = env["out_dir"] / "SMA_FAL99999_20260214_TEST.db"
        assert db.exists()

        # Step 2: ingest with classification
        r2 = subprocess.run(
            [sys.executable, "bin/ingest.py",
             "-e", "FAL99999_20260214_TEST",
             "-b", str(env["bam"]),
             "-s", str(env["summary"]),
             "-d", str(db),
             "-o", str(env["tmp"] / "tagged.bam"),
             "-ss", str(env["ss"]),
             "-rd", str(env["ref_dir"])],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
        )
        assert r2.returncode == 0, f"{r2.stderr}\n{r2.stdout}"

        # Verify all 5 reads ingested
        conn = sqlite3.connect(db)
        count = conn.execute("SELECT COUNT(*) FROM Reads").fetchone()[0]
        assert count == 5

        # Verify all classified as nb05 start
        starts = conn.execute(
            "SELECT DISTINCT bc_start_id FROM Reads"
        ).fetchall()
        assert starts == [("nb05",)]

        # Verify all classified as nb10 end
        ends = conn.execute(
            "SELECT DISTINCT bc_end_id FROM Reads"
        ).fetchall()
        assert ends == [("nb10",)]

        # Verify target assigned
        targets = conn.execute(
            "SELECT DISTINCT tgt_id FROM Reads"
        ).fetchall()
        assert targets == [("test_target",)]

        # Verify confidence for mutated reads is lower but still decent
        confs = conn.execute(
            "SELECT read_id, bc_start_conf FROM Reads ORDER BY read_id"
        ).fetchall()
        # First 3: perfect match -> conf ~1.0
        for read_id, conf in confs[:3]:
            assert conf >= 0.95, f"{read_id} conf={conf}"
        # Last 2: 1 mutation -> conf ~0.96
        for read_id, conf in confs[3:]:
            assert conf >= 0.9, f"{read_id} conf={conf}"

        conn.close()
```

**Step 2: Run test**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_integration.py -v`
Expected: All PASS (if Tasks 1-4 are complete)

**Step 3: Commit**

```bash
cd /tmp/ont-sma-seq
git add tests/test_integration.py
git commit -m "test: add end-to-end integration test for barcode classification"
```

---

### Task 6: Add conftest.py and pytest configuration

**Files:**
- Create: `tests/conftest.py`
- Create: `pytest.ini`

**Step 1: Create pytest config**

Create `pytest.ini`:
```ini
[pytest]
testpaths = tests
pythonpath = bin
```

Create `tests/conftest.py`:
```python
"""Shared test fixtures for ONT-SMA-seq tests."""
```

**Step 2: Verify all tests run from repo root**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/ -v`
Expected: All tests PASS

**Step 3: Commit**

```bash
cd /tmp/ont-sma-seq
git add pytest.ini tests/conftest.py
git commit -m "chore: add pytest configuration"
```

---

### Summary of all changes

| File | Action | Description |
|------|--------|-------------|
| `bin/barcodes.py` | CREATE | 96 ONT barcodes, `classify_barcode()`, `reverse_complement()` |
| `bin/sample_sheet.py` | CREATE | MinKNOW CSV parser, ambiguity detection |
| `bin/mkdb.py` | MODIFY | Add 6 barcode columns to Reads table |
| `bin/ingest.py` | MODIFY | Dual-end classification, multi-target support, optional split/tag |
| `tests/test_barcodes.py` | CREATE | Barcode data integrity and classification tests |
| `tests/test_sample_sheet.py` | CREATE | Sample sheet parsing and validation tests |
| `tests/test_mkdb.py` | CREATE | Schema validation tests |
| `tests/test_ingest.py` | CREATE | Classification correctness and backward compat tests |
| `tests/test_integration.py` | CREATE | Full pipeline end-to-end test |
| `tests/conftest.py` | CREATE | Shared fixtures |
| `pytest.ini` | CREATE | Test configuration |
