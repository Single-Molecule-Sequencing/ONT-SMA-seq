# SMA-seq Prepare Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build `bin/prepare.py` (experiment consolidation), `bin/align.py` (all-vs-all alignment metrics), and `bin/qc.py` (diagnostic visualizations) for the SMA-seq pipeline.

**Architecture:** Staged pipeline with 5 sequential stages (discover -> plan -> merge -> init -> align+qc). Builds on existing `sma_merge` package for discovery/validation. New alignment engine uses edlib NW mode (forward-strand only) with CIGAR parsing for segmented metrics. QC generates 10 publication-quality matplotlib plots.

**Tech Stack:** Python 3.13, edlib, pysam, matplotlib, scipy, numpy, sqlite3, pod5. Test runner: pytest. Existing infra: `bin/sma_merge/` package (discover, validate, tag, basecall, subsample).

**Design doc:** `docs/plans/2026-02-14-sma-prepare-design.md`

---

### Task 1: CIGAR Parsing Utilities

The alignment engine needs a CIGAR parser that walks edlib's CIGAR string and maps alignment columns to reference coordinates. This is the foundational building block for all metrics.

**Files:**
- Create: `bin/align.py`
- Test: `tests/test_align.py`

**Step 1: Write failing tests for CIGAR parsing**

```python
# tests/test_align.py
"""Tests for alignment engine."""
from __future__ import annotations

import pytest


class TestParseCigar:
    """Test CIGAR string parsing into (op, length) tuples."""

    def test_simple_match(self):
        from align import parse_cigar
        assert parse_cigar("10=") == [("=", 10)]

    def test_mixed_ops(self):
        from align import parse_cigar
        result = parse_cigar("5=1X3=2I4=1D2=")
        assert result == [("=", 5), ("X", 1), ("=", 3), ("I", 2), ("=", 4), ("D", 1), ("=", 2)]

    def test_empty_string(self):
        from align import parse_cigar
        assert parse_cigar("") == []


class TestCigarMetrics:
    """Test whole-alignment metrics computed from CIGAR."""

    def test_perfect_match(self):
        from align import compute_cigar_metrics
        m = compute_cigar_metrics("10=", ref_len=10, read_len=10)
        assert m["matches"] == 10
        assert m["identity"] == 1.0
        assert m["ref_coverage"] == 1.0
        assert m["max_ins"] == 0
        assert m["max_del"] == 0
        assert m["n_sig_indels"] == 0

    def test_with_mismatch(self):
        from align import compute_cigar_metrics
        m = compute_cigar_metrics("4=1X5=", ref_len=10, read_len=10)
        assert m["matches"] == 9
        assert m["mismatches"] == 1
        assert m["identity"] == pytest.approx(0.9)

    def test_with_insertion(self):
        from align import compute_cigar_metrics
        # 5 match, 3 insert, 5 match = read_len 13, ref consumed 10
        m = compute_cigar_metrics("5=3I5=", ref_len=10, read_len=13)
        assert m["max_ins"] == 3
        assert m["n_sig_indels"] == 0  # 3 < 5 threshold

    def test_significant_indel(self):
        from align import compute_cigar_metrics
        m = compute_cigar_metrics("5=6I5=", ref_len=10, read_len=16)
        assert m["max_ins"] == 6
        assert m["n_sig_indels"] == 1

    def test_with_deletion(self):
        from align import compute_cigar_metrics
        # 5 match, 2 del, 3 match = read_len 8, ref consumed 10
        m = compute_cigar_metrics("5=2D3=", ref_len=10, read_len=8)
        assert m["max_del"] == 2
        assert m["ref_coverage"] == pytest.approx(0.8)  # 8 ref bases covered by alignment matches
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_align.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'align'`

**Step 3: Implement CIGAR parsing and metrics**

```python
# bin/align.py
"""Alignment engine: all-vs-all forward-strand alignment with segmented metrics."""
from __future__ import annotations

import re
from dataclasses import dataclass, field

CIGAR_RE = re.compile(r"(\d+)([=XIDM])")
SIG_INDEL_THRESHOLD = 5


def parse_cigar(cigar: str) -> list[tuple[str, int]]:
    """Parse an edlib extended CIGAR string into (op, length) tuples."""
    return [(op, int(length)) for length, op in CIGAR_RE.findall(cigar)]


def compute_cigar_metrics(cigar: str, ref_len: int, read_len: int) -> dict:
    """Compute whole-alignment metrics from a CIGAR string.

    Returns dict with: matches, mismatches, insertions, deletions,
    identity, ref_coverage, max_ins, max_del, n_sig_indels, alignment_length.
    """
    ops = parse_cigar(cigar)
    matches = 0
    mismatches = 0
    ins_total = 0
    del_total = 0
    max_ins = 0
    max_del = 0
    n_sig_indels = 0
    ref_bases_covered = 0

    for op, length in ops:
        if op == "=" or op == "M":
            matches += length
            ref_bases_covered += length
        elif op == "X":
            mismatches += length
            ref_bases_covered += length
        elif op == "I":
            ins_total += length
            max_ins = max(max_ins, length)
            if length >= SIG_INDEL_THRESHOLD:
                n_sig_indels += 1
        elif op == "D":
            del_total += length
            ref_bases_covered += length
            max_del = max(max_del, length)
            if length >= SIG_INDEL_THRESHOLD:
                n_sig_indels += 1

    alignment_length = matches + mismatches + ins_total + del_total
    identity = matches / alignment_length if alignment_length > 0 else 0.0
    ref_coverage = (matches + mismatches) / ref_len if ref_len > 0 else 0.0

    return {
        "matches": matches,
        "mismatches": mismatches,
        "insertions": ins_total,
        "deletions": del_total,
        "alignment_length": alignment_length,
        "identity": identity,
        "ref_coverage": ref_coverage,
        "max_ins": max_ins,
        "max_del": max_del,
        "n_sig_indels": n_sig_indels,
    }
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_align.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/align.py tests/test_align.py
git commit -m "feat(align): add CIGAR parsing and whole-alignment metrics"
```

---

### Task 2: Segmented Alignment Metrics

Walk the CIGAR to compute identity and contiguity within the 5', middle, and 3' thirds of the reference.

**Files:**
- Modify: `bin/align.py`
- Modify: `tests/test_align.py`

**Step 1: Write failing tests for segmented metrics**

Add to `tests/test_align.py`:

```python
class TestSegmentedMetrics:
    """Test metrics computed within 5'/mid/3' reference segments."""

    def test_perfect_match_all_segments_identical(self):
        from align import compute_segmented_metrics
        # 30bp perfect match, ref_len=30 -> segments at 10, 20
        segs = compute_segmented_metrics("30=", ref_len=30)
        for seg in ["five_prime", "middle", "three_prime"]:
            assert segs[f"{seg}_identity"] == 1.0
            assert segs[f"{seg}_contiguity"] == 10
            assert segs[f"{seg}_gap_count"] == 0

    def test_mismatch_in_five_prime(self):
        from align import compute_segmented_metrics
        # First 10bp of ref: 3=1X6= -> 5' segment has identity 9/10
        segs = compute_segmented_metrics("3=1X6=10=10=", ref_len=30)
        assert segs["five_prime_identity"] == pytest.approx(0.9)
        assert segs["middle_identity"] == 1.0
        assert segs["three_prime_identity"] == 1.0

    def test_insertion_in_middle(self):
        from align import compute_segmented_metrics
        # ref positions 0-9 perfect, 10-19 has a 3bp insertion, 20-29 perfect
        segs = compute_segmented_metrics("10=5=3I5=10=", ref_len=30)
        assert segs["five_prime_gap_count"] == 0
        assert segs["middle_gap_count"] == 1
        assert segs["three_prime_gap_count"] == 0

    def test_contiguity_broken_by_deletion(self):
        from align import compute_segmented_metrics
        # 5' segment: 4=2D4= -> contiguity = 4 (longest run without gap)
        segs = compute_segmented_metrics("4=2D4=10=10=", ref_len=30)
        assert segs["five_prime_contiguity"] == 4
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_align.py::TestSegmentedMetrics -v`
Expected: FAIL — `cannot import name 'compute_segmented_metrics'`

**Step 3: Implement segmented metrics**

Add to `bin/align.py`:

```python
def compute_segmented_metrics(cigar: str, ref_len: int) -> dict:
    """Compute identity/contiguity/gap_count within 5', middle, 3' reference thirds.

    Walks the CIGAR, tracking the current reference position to assign each
    operation to the appropriate segment.
    """
    seg1_end = ref_len // 3
    seg2_end = 2 * ref_len // 3
    seg_len = [seg1_end, seg2_end - seg1_end, ref_len - seg2_end]

    # Per-segment accumulators
    seg_matches = [0, 0, 0]
    seg_mismatches = [0, 0, 0]
    seg_gap_count = [0, 0, 0]
    seg_contiguity = [0, 0, 0]
    seg_current_run = [0, 0, 0]

    def _seg_index(ref_pos: int) -> int:
        if ref_pos < seg1_end:
            return 0
        elif ref_pos < seg2_end:
            return 1
        else:
            return 2

    ref_pos = 0
    ops = parse_cigar(cigar)

    for op, length in ops:
        if op in ("=", "M"):
            for _ in range(length):
                if ref_pos < ref_len:
                    s = _seg_index(ref_pos)
                    seg_matches[s] += 1
                    seg_current_run[s] += 1
                    seg_contiguity[s] = max(seg_contiguity[s], seg_current_run[s])
                ref_pos += 1
        elif op == "X":
            for _ in range(length):
                if ref_pos < ref_len:
                    s = _seg_index(ref_pos)
                    seg_mismatches[s] += 1
                    seg_current_run[s] = 0
                ref_pos += 1
        elif op == "I":
            # Insertion: ref_pos doesn't advance. Assign to current segment.
            if ref_pos < ref_len:
                s = _seg_index(ref_pos)
            else:
                s = 2
            seg_gap_count[s] += 1
            seg_current_run[s] = 0
        elif op == "D":
            for _ in range(length):
                if ref_pos < ref_len:
                    s = _seg_index(ref_pos)
                    seg_gap_count[s] += 1 if _ == 0 else 0  # count deletion event once
                    seg_current_run[s] = 0
                ref_pos += 1

    names = ["five_prime", "middle", "three_prime"]
    result = {}
    for i, name in enumerate(names):
        total = seg_matches[i] + seg_mismatches[i]
        result[f"{name}_identity"] = seg_matches[i] / seg_len[i] if seg_len[i] > 0 else 0.0
        result[f"{name}_contiguity"] = seg_contiguity[i]
        result[f"{name}_gap_count"] = seg_gap_count[i]
    return result
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_align.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/align.py tests/test_align.py
git commit -m "feat(align): add segmented 5'/mid/3' identity and contiguity metrics"
```

---

### Task 3: 5'/3' Alignment Quality and Indel Position Metrics

Add terminal alignment quality checks and indel position tracking.

**Files:**
- Modify: `bin/align.py`
- Modify: `tests/test_align.py`

**Step 1: Write failing tests**

Add to `tests/test_align.py`:

```python
class TestTerminalMetrics:
    """Test 5' and 3' alignment quality metrics."""

    def test_perfect_alignment_no_offset(self):
        from align import compute_terminal_metrics
        m = compute_terminal_metrics("20=", ref_len=20)
        assert m["five_prime_offset"] == 0
        assert m["three_prime_offset"] == 0
        assert m["five_prime_identity_20"] == 1.0
        assert m["three_prime_identity_20"] == 1.0
        assert m["five_prime_first_match"] == 0

    def test_five_prime_deletion_offset(self):
        from align import compute_terminal_metrics
        # 3bp deletion at start, then 17 matches
        m = compute_terminal_metrics("3D17=", ref_len=20)
        assert m["five_prime_offset"] == 3
        assert m["five_prime_first_match"] == 3

    def test_three_prime_deletion_offset(self):
        from align import compute_terminal_metrics
        m = compute_terminal_metrics("18=2D", ref_len=20)
        assert m["three_prime_offset"] == 2

    def test_short_ref_identity_window(self):
        from align import compute_terminal_metrics
        # ref_len=10 < 20bp window -> use full ref
        m = compute_terminal_metrics("10=", ref_len=10)
        assert m["five_prime_identity_20"] == 1.0
        assert m["three_prime_identity_20"] == 1.0


class TestIndelPositions:
    """Test significant indel position tracking."""

    def test_no_indels(self):
        from align import compute_indel_positions
        result = compute_indel_positions("20=", min_size=3)
        assert result == []

    def test_small_indel_ignored(self):
        from align import compute_indel_positions
        result = compute_indel_positions("10=2I8=", min_size=3)
        assert result == []

    def test_significant_insertion_tracked(self):
        from align import compute_indel_positions
        result = compute_indel_positions("10=5I10=", min_size=3)
        assert len(result) == 1
        assert result[0] == ("I", 5, 10)  # type, size, ref_position

    def test_significant_deletion_tracked(self):
        from align import compute_indel_positions
        result = compute_indel_positions("5=4D15=", min_size=3)
        assert len(result) == 1
        assert result[0] == ("D", 4, 5)
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_align.py::TestTerminalMetrics tests/test_align.py::TestIndelPositions -v`
Expected: FAIL

**Step 3: Implement terminal metrics and indel positions**

Add to `bin/align.py`:

```python
def compute_terminal_metrics(cigar: str, ref_len: int, window: int = 20) -> dict:
    """Compute 5' and 3' alignment quality metrics.

    five_prime_offset: number of ref bases at 5' end covered by deletions
    five_prime_identity_20: identity within first `window` ref bases
    five_prime_first_match: ref position of first match/mismatch
    three_prime_offset: ref bases at 3' end covered by deletions
    three_prime_identity_20: identity within last `window` ref bases
    """
    ops = parse_cigar(cigar)
    win = min(window, ref_len)

    # Build per-ref-position array: 1=match, 0=mismatch, -1=deletion
    ref_status = [-1] * ref_len  # default: uncovered (deletion)
    ref_pos = 0
    for op, length in ops:
        if op in ("=", "M"):
            for _ in range(length):
                if ref_pos < ref_len:
                    ref_status[ref_pos] = 1
                ref_pos += 1
        elif op == "X":
            for _ in range(length):
                if ref_pos < ref_len:
                    ref_status[ref_pos] = 0
                ref_pos += 1
        elif op == "I":
            pass  # no ref consumed
        elif op == "D":
            for _ in range(length):
                if ref_pos < ref_len:
                    ref_status[ref_pos] = -1
                ref_pos += 1

    # 5' offset: leading deletions
    five_prime_offset = 0
    for s in ref_status:
        if s == -1:
            five_prime_offset += 1
        else:
            break

    # 3' offset: trailing deletions
    three_prime_offset = 0
    for s in reversed(ref_status):
        if s == -1:
            three_prime_offset += 1
        else:
            break

    # 5' first match position
    five_prime_first_match = 0
    for i, s in enumerate(ref_status):
        if s >= 0:
            five_prime_first_match = i
            break
    else:
        five_prime_first_match = ref_len

    # 5' identity within window
    five_window = ref_status[:win]
    five_prime_identity = sum(1 for s in five_window if s == 1) / len(five_window) if five_window else 0.0

    # 3' identity within window
    three_window = ref_status[-win:] if win > 0 else []
    three_prime_identity = sum(1 for s in three_window if s == 1) / len(three_window) if three_window else 0.0

    return {
        "five_prime_offset": five_prime_offset,
        "three_prime_offset": three_prime_offset,
        "five_prime_identity_20": five_prime_identity,
        "three_prime_identity_20": three_prime_identity,
        "five_prime_first_match": five_prime_first_match,
    }


def compute_indel_positions(cigar: str, min_size: int = 3) -> list[tuple[str, int, int]]:
    """Return list of (type, size, ref_position) for indels >= min_size."""
    ops = parse_cigar(cigar)
    result = []
    ref_pos = 0

    for op, length in ops:
        if op in ("=", "M", "X"):
            ref_pos += length
        elif op == "I":
            if length >= min_size:
                result.append(("I", length, ref_pos))
        elif op == "D":
            if length >= min_size:
                result.append(("D", length, ref_pos))
            ref_pos += length

    return result
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_align.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/align.py tests/test_align.py
git commit -m "feat(align): add 5'/3' terminal quality and indel position tracking"
```

---

### Task 4: All-vs-All Alignment and Classification

Wire up edlib to align each read against every reference, compute all metrics, rank targets, and classify reads.

**Files:**
- Modify: `bin/align.py`
- Modify: `tests/test_align.py`

**Step 1: Write failing tests**

Add to `tests/test_align.py`:

```python
class TestAlignReadToRef:
    """Test single read-vs-ref alignment with full metric computation."""

    def test_identical_sequences(self):
        from align import align_read_to_ref
        ref = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32bp
        m = align_read_to_ref(ref, ref, "target1")
        assert m["ref_id"] == "target1"
        assert m["ed"] == 0
        assert m["ned"] == 0.0
        assert m["identity"] == 1.0

    def test_with_errors(self):
        from align import align_read_to_ref
        ref  = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32bp
        read = "ACGTACTTACGTACGTACGTACGTACGTACGT"  # 1 mismatch at pos 6
        m = align_read_to_ref(read, ref, "target1")
        assert m["ed"] == 1
        assert m["ned"] == pytest.approx(1 / 32)

    def test_returns_all_metric_keys(self):
        from align import align_read_to_ref, METRIC_COLUMNS
        ref = "ACGTACGTAC"
        m = align_read_to_ref(ref, ref, "t1")
        for col in METRIC_COLUMNS:
            assert col in m, f"Missing metric: {col}"


class TestClassifyRead:
    """Test read classification against multiple references."""

    def test_assigns_best_match(self):
        from align import classify_read
        refs = {
            "targetA": "ACGTACGTACGTACGT",
            "targetB": "TTTTTTTTTTTTTTTT",
        }
        read = "ACGTACGTACGTACGT"
        result = classify_read(read, "read1", refs)
        assert result["assigned_ref"] == "targetA"
        assert result["best_ned"] == 0.0
        assert result["margin"] > 0

    def test_poor_match_flagged(self):
        from align import classify_read
        refs = {"targetA": "ACGTACGTACGTACGT"}
        read = "NNNNNNNNNNNNNNNN"
        result = classify_read(read, "read1", refs)
        assert result["confidence_flag"] == "POOR"

    def test_returns_all_pairwise_metrics(self):
        from align import classify_read
        refs = {"t1": "AAAA", "t2": "CCCC"}
        result = classify_read("AAAA", "r1", refs)
        assert len(result["pairwise"]) == 2
        assert result["pairwise"][0]["rank"] == 1
        assert result["pairwise"][1]["rank"] == 2
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_align.py::TestAlignReadToRef tests/test_align.py::TestClassifyRead -v`
Expected: FAIL

**Step 3: Implement alignment and classification**

Add to `bin/align.py`:

```python
import edlib

METRIC_COLUMNS = [
    "read_id", "ref_id", "read_len", "ref_len",
    "ed", "ned", "identity", "ref_coverage", "read_to_ref_ratio",
    "seg5_identity", "seg5_contiguity",
    "segM_identity", "segM_contiguity",
    "seg3_identity", "seg3_contiguity",
    "max_ins", "max_del", "n_sig_indels",
    "five_prime_offset", "five_prime_identity_20",
    "three_prime_offset", "three_prime_identity_20",
    "rank", "margin",
]


def align_read_to_ref(read_seq: str, ref_seq: str, ref_id: str) -> dict:
    """Align a single read to a single reference (forward-strand, NW mode).

    Returns a dict with all metric columns except rank/margin (added during classification).
    """
    result = edlib.align(read_seq, ref_seq, mode="NW", task="path")
    cigar = result.get("cigar", "")
    ed = result.get("editDistance", len(ref_seq))

    ref_len = len(ref_seq)
    read_len = len(read_seq)

    cigar_m = compute_cigar_metrics(cigar, ref_len, read_len)
    seg_m = compute_segmented_metrics(cigar, ref_len)
    term_m = compute_terminal_metrics(cigar, ref_len)

    return {
        "ref_id": ref_id,
        "read_len": read_len,
        "ref_len": ref_len,
        "ed": ed,
        "ned": ed / ref_len if ref_len > 0 else float("inf"),
        "identity": cigar_m["identity"],
        "ref_coverage": cigar_m["ref_coverage"],
        "read_to_ref_ratio": read_len / ref_len if ref_len > 0 else 0.0,
        "seg5_identity": seg_m["five_prime_identity"],
        "seg5_contiguity": seg_m["five_prime_contiguity"],
        "segM_identity": seg_m["middle_identity"],
        "segM_contiguity": seg_m["middle_contiguity"],
        "seg3_identity": seg_m["three_prime_identity"],
        "seg3_contiguity": seg_m["three_prime_contiguity"],
        "max_ins": cigar_m["max_ins"],
        "max_del": cigar_m["max_del"],
        "n_sig_indels": cigar_m["n_sig_indels"],
        "five_prime_offset": term_m["five_prime_offset"],
        "five_prime_identity_20": term_m["five_prime_identity_20"],
        "three_prime_offset": term_m["three_prime_offset"],
        "three_prime_identity_20": term_m["three_prime_identity_20"],
        "rank": 0,
        "margin": 0.0,
    }


def classify_read(
    read_seq: str, read_id: str, refs: dict[str, str],
) -> dict:
    """Align read against all references and classify by best NED.

    Returns dict with: read_id, read_len, assigned_ref, best_ned, second_ned,
    margin, confidence_flag, pairwise (list of per-ref metric dicts).
    """
    pairwise = []
    for ref_id, ref_seq in refs.items():
        m = align_read_to_ref(read_seq, ref_seq, ref_id)
        m["read_id"] = read_id
        pairwise.append(m)

    pairwise.sort(key=lambda x: x["ned"])
    for i, m in enumerate(pairwise):
        m["rank"] = i + 1

    best_ned = pairwise[0]["ned"]
    second_ned = pairwise[1]["ned"] if len(pairwise) > 1 else float("inf")
    margin = second_ned - best_ned

    for m in pairwise:
        m["margin"] = margin

    if best_ned > 0.5:
        flag = "POOR"
    elif margin <= 0.1:
        flag = "LOW"
    else:
        flag = "HIGH"

    return {
        "read_id": read_id,
        "read_len": len(read_seq),
        "assigned_ref": pairwise[0]["ref_id"],
        "best_ned": best_ned,
        "second_ned": second_ned,
        "margin": margin,
        "confidence_flag": flag,
        "pairwise": pairwise,
    }
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_align.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/align.py tests/test_align.py
git commit -m "feat(align): add all-vs-all alignment engine with read classification"
```

---

### Task 5: BAM Streaming and TSV Output

Add functions to stream reads from a BAM, run all-vs-all alignment, and write `alignments.tsv` + `classification.tsv`.

**Files:**
- Modify: `bin/align.py`
- Modify: `tests/test_align.py`

**Step 1: Write failing tests**

Add to `tests/test_align.py`:

```python
import tempfile
from pathlib import Path


def _make_test_bam(tmp_path: Path, sequences: dict[str, str]) -> Path:
    """Create a minimal unaligned BAM with given {read_id: sequence} pairs."""
    import pysam
    bam_path = tmp_path / "test.bam"
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "unknown"}})
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
        for rid, seq in sequences.items():
            a = pysam.AlignedSegment(header)
            a.query_name = rid
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.flag = 4  # unmapped
            f.write(a)
    return bam_path


class TestProcessBam:
    """Test BAM streaming + TSV output."""

    def test_writes_alignments_tsv(self, tmp_path):
        from align import process_bam
        bam = _make_test_bam(tmp_path, {"r1": "ACGTACGT", "r2": "TTTTTTTT"})
        refs = {"t1": "ACGTACGT", "t2": "TTTTTTTT"}
        align_tsv = tmp_path / "alignments.tsv"
        class_tsv = tmp_path / "classification.tsv"
        process_bam(bam, refs, align_tsv, class_tsv)
        assert align_tsv.exists()
        lines = align_tsv.read_text().strip().split("\n")
        assert len(lines) == 5  # header + 2 reads * 2 refs

    def test_writes_classification_tsv(self, tmp_path):
        from align import process_bam
        bam = _make_test_bam(tmp_path, {"r1": "ACGTACGT"})
        refs = {"t1": "ACGTACGT", "t2": "TTTTTTTT"}
        align_tsv = tmp_path / "alignments.tsv"
        class_tsv = tmp_path / "classification.tsv"
        process_bam(bam, refs, align_tsv, class_tsv)
        assert class_tsv.exists()
        lines = class_tsv.read_text().strip().split("\n")
        assert len(lines) == 2  # header + 1 read
        assert "t1" in lines[1]  # assigned to perfect match

    def test_end_reason_included_when_tagged(self, tmp_path):
        from align import process_bam
        import pysam
        bam_path = tmp_path / "tagged.bam"
        header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
            a = pysam.AlignedSegment(header)
            a.query_name = "r1"
            a.query_sequence = "ACGT"
            a.query_qualities = pysam.qualitystring_to_array("IIII")
            a.flag = 4
            a.set_tag("er", "signal_positive", "Z")
            f.write(a)
        refs = {"t1": "ACGT"}
        class_tsv = tmp_path / "class.tsv"
        process_bam(bam_path, refs, tmp_path / "align.tsv", class_tsv)
        content = class_tsv.read_text()
        assert "signal_positive" in content
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_align.py::TestProcessBam -v`
Expected: FAIL

**Step 3: Implement BAM processing**

Add to `bin/align.py`:

```python
from pathlib import Path
import csv

import pysam


def parse_fasta(fasta_path: Path) -> dict[str, str]:
    """Parse a FASTA file into {header: sequence} dict."""
    refs = {}
    current_id = None
    current_seq = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    refs[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
    if current_id is not None:
        refs[current_id] = "".join(current_seq)
    return refs


def process_bam(
    bam_path: Path,
    refs: dict[str, str],
    alignments_tsv: Path,
    classification_tsv: Path,
) -> int:
    """Stream reads from BAM, align all-vs-all, write TSVs.

    Returns the number of reads processed.
    """
    alignments_tsv.parent.mkdir(parents=True, exist_ok=True)
    classification_tsv.parent.mkdir(parents=True, exist_ok=True)

    class_fields = [
        "read_id", "read_len", "assigned_ref",
        "best_ned", "second_ned", "margin", "confidence_flag", "end_reason",
    ]

    n_processed = 0

    with (
        open(alignments_tsv, "w", newline="") as af,
        open(classification_tsv, "w", newline="") as cf,
        pysam.AlignmentFile(str(bam_path), check_sq=False) as bam,
    ):
        aw = csv.DictWriter(af, fieldnames=METRIC_COLUMNS, delimiter="\t", extrasaction="ignore")
        aw.writeheader()
        cw = csv.DictWriter(cf, fieldnames=class_fields, delimiter="\t")
        cw.writeheader()

        for read in bam:
            seq = read.query_sequence
            if seq is None:
                continue

            rid = read.query_name
            result = classify_read(seq, rid, refs)

            for pw in result["pairwise"]:
                aw.writerow(pw)

            er = ""
            try:
                er = read.get_tag("er")
            except KeyError:
                pass

            cw.writerow({
                "read_id": result["read_id"],
                "read_len": result["read_len"],
                "assigned_ref": result["assigned_ref"],
                "best_ned": f"{result['best_ned']:.6f}",
                "second_ned": f"{result['second_ned']:.6f}",
                "margin": f"{result['margin']:.6f}",
                "confidence_flag": result["confidence_flag"],
                "end_reason": er,
            })
            n_processed += 1

    return n_processed
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_align.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/align.py tests/test_align.py
git commit -m "feat(align): add BAM streaming and TSV output for alignments + classification"
```

---

### Task 6: QC Plotting Module — Core Utilities + Margin KDE

Build `bin/qc.py` with the shared plotting setup and the first diagnostic plot (margin KDE).

**Files:**
- Create: `bin/qc.py`
- Create: `tests/test_qc.py`

**Step 1: Write failing tests**

```python
# tests/test_qc.py
"""Tests for QC visualization module."""
from __future__ import annotations

import csv
from pathlib import Path

import pytest


def _make_classification_tsv(tmp_path: Path, rows: list[dict]) -> Path:
    """Write a classification.tsv for testing."""
    tsv = tmp_path / "classification.tsv"
    fields = ["read_id", "read_len", "assigned_ref", "best_ned", "second_ned", "margin", "confidence_flag", "end_reason"]
    with open(tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return tsv


def _make_alignments_tsv(tmp_path: Path, rows: list[dict]) -> Path:
    """Write an alignments.tsv for testing."""
    from align import METRIC_COLUMNS
    tsv = tmp_path / "alignments.tsv"
    with open(tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=METRIC_COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return tsv


class TestLoadData:
    """Test TSV loading utilities."""

    def test_load_classification(self, tmp_path):
        from qc import load_classification
        tsv = _make_classification_tsv(tmp_path, [
            {"read_id": "r1", "read_len": "100", "assigned_ref": "t1",
             "best_ned": "0.05", "second_ned": "0.4", "margin": "0.35",
             "confidence_flag": "HIGH", "end_reason": "signal_positive"},
        ])
        df = load_classification(tsv)
        assert len(df) == 1
        assert df["margin"].dtype == float

    def test_load_alignments(self, tmp_path):
        from qc import load_alignments
        from align import METRIC_COLUMNS
        row = {col: "0" for col in METRIC_COLUMNS}
        row["read_id"] = "r1"
        row["ref_id"] = "t1"
        tsv = _make_alignments_tsv(tmp_path, [row])
        df = load_alignments(tsv)
        assert len(df) == 1


class TestMarginKDE:
    """Test margin KDE plot generation."""

    def test_creates_png(self, tmp_path):
        from qc import plot_margin_kde
        tsv = _make_classification_tsv(tmp_path, [
            {"read_id": f"r{i}", "read_len": "100", "assigned_ref": "t1",
             "best_ned": "0.05", "second_ned": "0.4", "margin": f"{0.1 + i * 0.01}",
             "confidence_flag": "HIGH", "end_reason": "signal_positive"}
            for i in range(50)
        ])
        from qc import load_classification
        df = load_classification(tsv)
        out = tmp_path / "qc"
        out.mkdir()
        plot_margin_kde(df, out)
        assert (out / "margin_kde.png").exists()
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_qc.py -v`
Expected: FAIL

**Step 3: Implement QC core + margin KDE**

```python
# bin/qc.py
"""QC visualization module for SMA-seq alignment classification."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

# Plotting defaults
STYLE = "seaborn-v0_8-whitegrid"
FIGSIZE = (8, 5)
DPI = 300
SAVGOL_WINDOW = 51
SAVGOL_POLY = 3

FLOAT_COLS_CLASS = ["best_ned", "second_ned", "margin"]
INT_COLS_CLASS = ["read_len"]
FLOAT_COLS_ALIGN = [
    "ned", "identity", "ref_coverage", "read_to_ref_ratio",
    "seg5_identity", "seg5_contiguity", "segM_identity", "segM_contiguity",
    "seg3_identity", "seg3_contiguity", "five_prime_identity_20", "three_prime_identity_20",
    "margin",
]
INT_COLS_ALIGN = [
    "read_len", "ref_len", "ed", "max_ins", "max_del", "n_sig_indels",
    "five_prime_offset", "three_prime_offset", "rank",
]


def load_classification(tsv: Path) -> pd.DataFrame:
    """Load classification.tsv into a DataFrame with correct dtypes."""
    df = pd.read_csv(tsv, sep="\t")
    for col in FLOAT_COLS_CLASS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in INT_COLS_CLASS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    return df


def load_alignments(tsv: Path) -> pd.DataFrame:
    """Load alignments.tsv into a DataFrame with correct dtypes."""
    df = pd.read_csv(tsv, sep="\t")
    for col in FLOAT_COLS_ALIGN:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in INT_COLS_ALIGN:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    return df


def _kde_smooth(values: np.ndarray, bw: float = 0.02, n_points: int = 500) -> tuple[np.ndarray, np.ndarray]:
    """Compute Gaussian KDE with Savitzky-Golay smoothing."""
    from scipy.stats import gaussian_kde
    vmin, vmax = values.min(), values.max()
    pad = (vmax - vmin) * 0.1
    x = np.linspace(vmin - pad, vmax + pad, n_points)
    try:
        kde = gaussian_kde(values, bw_method=bw)
        y = kde(x)
    except (np.linalg.LinAlgError, ValueError):
        y = np.zeros_like(x)
    window = min(SAVGOL_WINDOW, len(y) - 1)
    if window % 2 == 0:
        window -= 1
    if window >= SAVGOL_POLY + 2:
        y = savgol_filter(y, window, SAVGOL_POLY)
        y = np.maximum(y, 0)
    return x, y


def _save(fig: plt.Figure, outdir: Path, name: str) -> None:
    """Save figure as PNG and PDF."""
    fig.savefig(outdir / f"{name}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(outdir / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_margin_kde(df: pd.DataFrame, outdir: Path) -> None:
    """Plot 2: Score margin distribution KDE."""
    plt.style.use(STYLE)
    fig, ax = plt.subplots(figsize=FIGSIZE)

    margins = df["margin"].dropna().values
    x, y = _kde_smooth(margins)
    ax.fill_between(x, y, alpha=0.3)
    ax.plot(x, y, linewidth=1.5)

    # Threshold candidates
    for thresh in [0.05, 0.1, 0.2]:
        ax.axvline(thresh, color="red", linestyle="--", alpha=0.5, label=f"threshold={thresh}")

    ax.set_xlabel("Classification Margin (NED_2nd - NED_best)")
    ax.set_ylabel("Density")
    ax.set_title("Score Margin Distribution")
    ax.legend(fontsize=8)
    _save(fig, outdir, "margin_kde")
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_qc.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/qc.py tests/test_qc.py
git commit -m "feat(qc): add data loading utilities and margin KDE plot"
```

---

### Task 7: QC Plots — NED by Target, Length vs NED, Segmented Violins

Add the three core classification diagnostic plots.

**Files:**
- Modify: `bin/qc.py`
- Modify: `tests/test_qc.py`

**Step 1: Write failing tests**

Add to `tests/test_qc.py`:

```python
def _sample_classification_df(n=100):
    """Generate a sample classification DataFrame for plotting tests."""
    import pandas as pd
    import numpy as np
    rng = np.random.default_rng(42)
    targets = ["t1", "t2"]
    ers = ["signal_positive", "data_service_unblock_mux_change"]
    return pd.DataFrame({
        "read_id": [f"r{i}" for i in range(n)],
        "read_len": rng.integers(100, 300, n),
        "assigned_ref": rng.choice(targets, n),
        "best_ned": rng.uniform(0.01, 0.3, n),
        "second_ned": rng.uniform(0.3, 0.8, n),
        "margin": rng.uniform(0.05, 0.5, n),
        "confidence_flag": rng.choice(["HIGH", "LOW", "POOR"], n, p=[0.8, 0.15, 0.05]),
        "end_reason": rng.choice(ers, n),
    })


def _sample_alignments_df(n_reads=50, n_refs=2):
    """Generate a sample alignments DataFrame."""
    import pandas as pd
    import numpy as np
    rng = np.random.default_rng(42)
    rows = []
    for i in range(n_reads):
        for j in range(n_refs):
            rows.append({
                "read_id": f"r{i}", "ref_id": f"t{j}",
                "read_len": 200, "ref_len": 200,
                "ed": int(rng.integers(1, 50)),
                "ned": float(rng.uniform(0.01, 0.3)),
                "identity": float(rng.uniform(0.7, 1.0)),
                "ref_coverage": float(rng.uniform(0.8, 1.0)),
                "read_to_ref_ratio": float(rng.uniform(0.8, 1.2)),
                "seg5_identity": float(rng.uniform(0.7, 1.0)),
                "seg5_contiguity": int(rng.integers(10, 60)),
                "segM_identity": float(rng.uniform(0.7, 1.0)),
                "segM_contiguity": int(rng.integers(10, 60)),
                "seg3_identity": float(rng.uniform(0.7, 1.0)),
                "seg3_contiguity": int(rng.integers(10, 60)),
                "max_ins": int(rng.integers(0, 10)),
                "max_del": int(rng.integers(0, 10)),
                "n_sig_indels": int(rng.integers(0, 3)),
                "five_prime_offset": int(rng.integers(0, 5)),
                "five_prime_identity_20": float(rng.uniform(0.7, 1.0)),
                "three_prime_offset": int(rng.integers(0, 5)),
                "three_prime_identity_20": float(rng.uniform(0.7, 1.0)),
                "rank": j + 1,
                "margin": float(rng.uniform(0.05, 0.5)),
            })
    return pd.DataFrame(rows)


class TestNedByTarget:
    def test_creates_png(self, tmp_path):
        from qc import plot_ned_by_target
        df_align = _sample_alignments_df()
        df_class = _sample_classification_df()
        out = tmp_path / "qc"
        out.mkdir()
        plot_ned_by_target(df_align, df_class, out)
        assert (out / "ned_by_target.png").exists()


class TestLengthVsNed:
    def test_creates_png(self, tmp_path):
        from qc import plot_length_vs_ned
        df = _sample_classification_df()
        out = tmp_path / "qc"
        out.mkdir()
        plot_length_vs_ned(df, out)
        assert (out / "length_vs_ned_scatter.png").exists()


class TestSegmentedViolins:
    def test_creates_png(self, tmp_path):
        from qc import plot_segmented_violins
        df_align = _sample_alignments_df()
        df_class = _sample_classification_df(50)
        out = tmp_path / "qc"
        out.mkdir()
        plot_segmented_violins(df_align, df_class, out)
        assert (out / "segmented_identity_violins.png").exists()
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_qc.py::TestNedByTarget tests/test_qc.py::TestLengthVsNed tests/test_qc.py::TestSegmentedViolins -v`
Expected: FAIL

**Step 3: Implement the three plots**

Add to `bin/qc.py`:

```python
def plot_ned_by_target(df_align: pd.DataFrame, df_class: pd.DataFrame, outdir: Path) -> None:
    """Plot 3: NED distribution by assigned target (best + second-best)."""
    plt.style.use(STYLE)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    # Merge classification assignment into alignments
    assign = df_class[["read_id", "assigned_ref"]].drop_duplicates()
    merged = df_align.merge(assign, on="read_id", how="left")

    targets = sorted(merged["assigned_ref"].dropna().unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(targets), 1)))

    for i, tgt in enumerate(targets):
        # Best (rank 1)
        best = merged[(merged["assigned_ref"] == tgt) & (merged["rank"] == 1)]
        if len(best) > 2:
            x, y = _kde_smooth(best["ned"].values)
            ax1.plot(x, y, color=colors[i], label=tgt, linewidth=1.5)
            ax1.fill_between(x, y, alpha=0.2, color=colors[i])

        # Second best (rank 2)
        second = merged[(merged["assigned_ref"] == tgt) & (merged["rank"] == 2)]
        if len(second) > 2:
            x, y = _kde_smooth(second["ned"].values)
            ax2.plot(x, y, color=colors[i], label=tgt, linewidth=1.5)
            ax2.fill_between(x, y, alpha=0.2, color=colors[i])

    ax1.set_ylabel("Density")
    ax1.set_title("NED to Assigned Target (Rank 1)")
    ax1.legend(fontsize=8)
    ax2.set_xlabel("Normalized Edit Distance")
    ax2.set_ylabel("Density")
    ax2.set_title("NED to Second-Best Target (Rank 2)")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    _save(fig, outdir, "ned_by_target")


def plot_length_vs_ned(df: pd.DataFrame, outdir: Path, ref_lengths: dict[str, int] | None = None) -> None:
    """Plot 5: Read length vs NED scatter."""
    plt.style.use(STYLE)
    fig, ax = plt.subplots(figsize=FIGSIZE)

    targets = sorted(df["assigned_ref"].dropna().unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(targets), 1)))

    for i, tgt in enumerate(targets):
        sub = df[df["assigned_ref"] == tgt]
        ax.scatter(sub["read_len"], sub["best_ned"], s=5, alpha=0.4, color=colors[i], label=tgt)

    # Threshold lines
    for thresh in [0.1, 0.3, 0.5]:
        ax.axhline(thresh, color="gray", linestyle="--", alpha=0.5, linewidth=0.8)

    if ref_lengths:
        for name, length in ref_lengths.items():
            ax.axvline(length, color="red", linestyle=":", alpha=0.5, linewidth=0.8)

    ax.set_xlabel("Read Length (bp)")
    ax.set_ylabel("NED to Assigned Target")
    ax.set_title("Read Length vs. Normalized Edit Distance")
    ax.legend(fontsize=8, markerscale=3)
    _save(fig, outdir, "length_vs_ned_scatter")


def plot_segmented_violins(df_align: pd.DataFrame, df_class: pd.DataFrame, outdir: Path) -> None:
    """Plot 4: Segmented identity violin plots per target."""
    plt.style.use(STYLE)

    assign = df_class[["read_id", "assigned_ref"]].drop_duplicates()
    best = df_align[df_align["rank"] == 1].merge(assign, on="read_id", how="left")

    targets = sorted(best["assigned_ref"].dropna().unique())
    n_targets = max(len(targets), 1)
    fig, axes = plt.subplots(1, n_targets, figsize=(4 * n_targets, 5), sharey=True, squeeze=False)

    segments = ["seg5_identity", "segM_identity", "seg3_identity"]
    labels = ["5'", "Mid", "3'"]

    for idx, tgt in enumerate(targets):
        ax = axes[0, idx]
        sub = best[best["assigned_ref"] == tgt]
        data = [sub[seg].dropna().values for seg in segments]
        if all(len(d) > 0 for d in data):
            vp = ax.violinplot(data, showmedians=True)
            ax.set_xticks([1, 2, 3])
            ax.set_xticklabels(labels)
        ax.set_title(tgt, fontsize=10)
        if idx == 0:
            ax.set_ylabel("Identity")

    fig.suptitle("Segmented Alignment Identity (5' / Mid / 3')")
    fig.tight_layout()
    _save(fig, outdir, "segmented_identity_violins")
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_qc.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/qc.py tests/test_qc.py
git commit -m "feat(qc): add NED-by-target, length-vs-NED, and segmented violin plots"
```

---

### Task 8: QC Plots — 5' Quality, Indel Profile, Threshold Sweep

Add the remaining diagnostic plots focused on alignment quality and classification threshold selection.

**Files:**
- Modify: `bin/qc.py`
- Modify: `tests/test_qc.py`

**Step 1: Write failing tests**

Add to `tests/test_qc.py`:

```python
class TestFivePrimeQuality:
    def test_creates_png(self, tmp_path):
        from qc import plot_five_prime_quality
        df_align = _sample_alignments_df()
        df_class = _sample_classification_df(50)
        out = tmp_path / "qc"
        out.mkdir()
        plot_five_prime_quality(df_align, df_class, out)
        assert (out / "five_prime_quality.png").exists()


class TestIndelProfile:
    def test_creates_png(self, tmp_path):
        from qc import plot_indel_profile
        df_align = _sample_alignments_df()
        df_class = _sample_classification_df(50)
        out = tmp_path / "qc"
        out.mkdir()
        plot_indel_profile(df_align, df_class, out)
        assert (out / "indel_profile.png").exists()


class TestThresholdSweep:
    def test_creates_png(self, tmp_path):
        from qc import plot_threshold_sweep
        df = _sample_classification_df()
        out = tmp_path / "qc"
        out.mkdir()
        plot_threshold_sweep(df, out)
        assert (out / "threshold_sweep.png").exists()
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_qc.py::TestFivePrimeQuality tests/test_qc.py::TestIndelProfile tests/test_qc.py::TestThresholdSweep -v`
Expected: FAIL

**Step 3: Implement the three plots**

Add to `bin/qc.py`:

```python
def plot_five_prime_quality(df_align: pd.DataFrame, df_class: pd.DataFrame, outdir: Path) -> None:
    """Plot 6: 5' alignment quality (offset + identity KDEs)."""
    plt.style.use(STYLE)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

    assign = df_class[["read_id", "assigned_ref"]].drop_duplicates()
    best = df_align[df_align["rank"] == 1].merge(assign, on="read_id", how="left")
    targets = sorted(best["assigned_ref"].dropna().unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(targets), 1)))

    for i, tgt in enumerate(targets):
        sub = best[best["assigned_ref"] == tgt]

        offsets = sub["five_prime_offset"].dropna().values.astype(float)
        if len(offsets) > 2:
            x, y = _kde_smooth(offsets, bw=0.5)
            ax1.plot(x, y, color=colors[i], label=tgt)

        ident = sub["five_prime_identity_20"].dropna().values
        if len(ident) > 2:
            x, y = _kde_smooth(ident)
            ax2.plot(x, y, color=colors[i], label=tgt)

    ax1.set_xlabel("5' Offset (ref bases uncovered)")
    ax1.set_ylabel("Density")
    ax1.set_title("5' Alignment Offset")
    ax1.legend(fontsize=8)

    ax2.set_xlabel("5' Identity (first 20bp)")
    ax2.set_ylabel("Density")
    ax2.set_title("5' Alignment Identity")
    ax2.legend(fontsize=8)

    fig.tight_layout()
    _save(fig, outdir, "five_prime_quality")


def plot_indel_profile(df_align: pd.DataFrame, df_class: pd.DataFrame, outdir: Path) -> None:
    """Plot 7: Significant indel profile (size KDE + position histogram)."""
    plt.style.use(STYLE)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

    assign = df_class[["read_id", "assigned_ref"]].drop_duplicates()
    best = df_align[df_align["rank"] == 1].merge(assign, on="read_id", how="left")
    targets = sorted(best["assigned_ref"].dropna().unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(targets), 1)))

    for i, tgt in enumerate(targets):
        sub = best[best["assigned_ref"] == tgt]
        ins = sub["max_ins"].dropna().values.astype(float)
        dels = sub["max_del"].dropna().values.astype(float)
        if len(ins) > 2:
            x, y = _kde_smooth(ins, bw=0.5)
            ax1.plot(x, y, color=colors[i], linestyle="-", label=f"{tgt} ins")
        if len(dels) > 2:
            x, y = _kde_smooth(dels, bw=0.5)
            ax1.plot(x, y, color=colors[i], linestyle="--", label=f"{tgt} del")

    ax1.set_xlabel("Indel Size (bp)")
    ax1.set_ylabel("Density")
    ax1.set_title("Largest Indel Size Distribution")
    ax1.legend(fontsize=7)

    # Position histogram: bin indel counts by 5'/mid/3' segment
    seg_counts = {"5'": 0, "Mid": 0, "3'": 0}
    for _, row in best.iterrows():
        ref_len = row.get("ref_len", 0)
        if ref_len <= 0:
            continue
        offset = row.get("five_prime_offset", 0)
        if offset > ref_len / 3:
            seg_counts["5'"] += 1
        three_offset = row.get("three_prime_offset", 0)
        if three_offset > ref_len / 3:
            seg_counts["3'"] += 1

    n_sig = best["n_sig_indels"].dropna()
    ax2.bar(["5'", "Mid", "3'"], [seg_counts["5'"], len(n_sig[n_sig > 0]) - seg_counts["5'"] - seg_counts["3'"], seg_counts["3'"]])
    ax2.set_xlabel("Reference Segment")
    ax2.set_ylabel("Reads with Significant Indels")
    ax2.set_title("Indel Hotspot by Segment")

    fig.tight_layout()
    _save(fig, outdir, "indel_profile")


def plot_threshold_sweep(df: pd.DataFrame, outdir: Path) -> None:
    """Plot 8: Misclassification threshold sweep."""
    plt.style.use(STYLE)
    fig, ax1 = plt.subplots(figsize=FIGSIZE)
    ax2 = ax1.twinx()

    thresholds = np.arange(0.1, 0.52, 0.02)
    frac_classified = []
    misclass_rate = []
    n_total = len(df)

    for t in thresholds:
        classified = df[df["best_ned"] <= t]
        frac = len(classified) / n_total if n_total > 0 else 0
        frac_classified.append(frac)

        if len(classified) > 0:
            ambiguous = classified[classified["margin"] <= 0.05]
            misclass_rate.append(len(ambiguous) / len(classified))
        else:
            misclass_rate.append(0)

    ax1.plot(thresholds, frac_classified, "b-", linewidth=2, label="Fraction Classified")
    ax2.plot(thresholds, misclass_rate, "r--", linewidth=2, label="Est. Misclassification Rate")

    ax1.set_xlabel("NED Threshold")
    ax1.set_ylabel("Fraction of Reads Classified", color="blue")
    ax2.set_ylabel("Est. Misclassification Rate", color="red")
    ax1.set_title("Classification Threshold Sweep")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=8)

    _save(fig, outdir, "threshold_sweep")
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_qc.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/qc.py tests/test_qc.py
git commit -m "feat(qc): add 5' quality, indel profile, and threshold sweep plots"
```

---

### Task 9: QC Plots — Heatmap, Per-Read Profiles, End Reason Facets

Add the remaining three plot types.

**Files:**
- Modify: `bin/qc.py`
- Modify: `tests/test_qc.py`

**Step 1: Write failing tests**

Add to `tests/test_qc.py`:

```python
class TestClassificationHeatmap:
    def test_creates_png(self, tmp_path):
        from qc import plot_classification_heatmap
        df_align = _sample_alignments_df()
        df_class = _sample_classification_df(50)
        out = tmp_path / "qc"
        out.mkdir()
        plot_classification_heatmap(df_align, df_class, out)
        assert (out / "classification_heatmap.png").exists()


class TestPerReadProfiles:
    def test_creates_png(self, tmp_path):
        from qc import plot_per_read_profiles
        df_align = _sample_alignments_df()
        out = tmp_path / "qc"
        out.mkdir()
        plot_per_read_profiles(df_align, out, n_sample=10)
        assert (out / "per_read_profiles_sample.png").exists()


class TestEndReasonFacets:
    def test_creates_faceted_plots(self, tmp_path):
        from qc import plot_end_reason_facets
        df_align = _sample_alignments_df()
        df_class = _sample_classification_df(50)
        out = tmp_path / "qc"
        out.mkdir()
        plot_end_reason_facets(df_align, df_class, out)
        assert (out / "margin_by_endreason.png").exists()
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_qc.py::TestClassificationHeatmap tests/test_qc.py::TestPerReadProfiles tests/test_qc.py::TestEndReasonFacets -v`
Expected: FAIL

**Step 3: Implement the three plots**

Add to `bin/qc.py`:

```python
def plot_classification_heatmap(df_align: pd.DataFrame, df_class: pd.DataFrame, outdir: Path, max_reads: int = 5000) -> None:
    """Plot 1: Classification landscape heatmap (reads x refs, color=NED)."""
    plt.style.use(STYLE)

    assign = df_class[["read_id", "assigned_ref"]].drop_duplicates()

    # Subsample if too many reads
    read_ids = df_align["read_id"].unique()
    if len(read_ids) > max_reads:
        rng = np.random.default_rng(42)
        read_ids = rng.choice(read_ids, max_reads, replace=False)

    sub = df_align[df_align["read_id"].isin(read_ids)].merge(assign, on="read_id", how="left")
    refs = sorted(sub["ref_id"].unique())

    pivot = sub.pivot_table(index="read_id", columns="ref_id", values="ned", aggfunc="first")
    pivot = pivot.reindex(columns=refs)

    # Sort by assigned target, then by NED to that target
    pivot = pivot.merge(assign.set_index("read_id"), left_index=True, right_index=True)
    pivot = pivot.sort_values(["assigned_ref", refs[0] if refs else ""])
    pivot = pivot.drop(columns=["assigned_ref"])

    fig, ax = plt.subplots(figsize=(max(6, len(refs) * 2), 8))
    im = ax.imshow(pivot.values, aspect="auto", cmap="viridis_r", vmin=0, vmax=0.8)
    ax.set_xticks(range(len(refs)))
    ax.set_xticklabels(refs, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(f"Reads (n={len(pivot)})")
    ax.set_yticks([])
    ax.set_title("Classification Landscape (color = NED)")
    fig.colorbar(im, ax=ax, label="Normalized Edit Distance")
    fig.tight_layout()
    _save(fig, outdir, "classification_heatmap")


def plot_per_read_profiles(df_align: pd.DataFrame, outdir: Path, n_sample: int = 50) -> None:
    """Plot 9: Per-read bar charts showing NED to all targets (sampled)."""
    plt.style.use(STYLE)
    read_ids = df_align["read_id"].unique()
    rng = np.random.default_rng(42)
    sampled = rng.choice(read_ids, min(n_sample, len(read_ids)), replace=False)

    n = len(sampled)
    ncols = min(5, n)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 2.5 * nrows), squeeze=False)

    for idx, rid in enumerate(sampled):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        sub = df_align[df_align["read_id"] == rid].sort_values("ned")
        colors = ["green" if r == 1 else "gray" for r in sub["rank"]]
        ax.barh(sub["ref_id"], sub["ned"], color=colors)
        ax.set_xlim(0, 1)
        ax.set_title(rid[:8], fontsize=7)
        ax.tick_params(labelsize=6)

    # Hide empty subplots
    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    fig.suptitle(f"Per-Read NED Profiles (n={n} sampled)", fontsize=11)
    fig.tight_layout()
    _save(fig, outdir, "per_read_profiles_sample")


def plot_end_reason_facets(df_align: pd.DataFrame, df_class: pd.DataFrame, outdir: Path) -> None:
    """Plot 10: Margin KDE faceted by end reason."""
    plt.style.use(STYLE)

    ers = df_class["end_reason"].dropna().unique()
    if len(ers) < 2:
        # Not enough end reasons to facet — skip
        fig, ax = plt.subplots(figsize=FIGSIZE)
        ax.text(0.5, 0.5, "Insufficient end reason diversity for faceting",
                ha="center", va="center", transform=ax.transAxes)
        _save(fig, outdir, "margin_by_endreason")
        return

    fig, axes = plt.subplots(1, len(ers), figsize=(5 * len(ers), 5), sharey=True, squeeze=False)

    for idx, er in enumerate(sorted(ers)):
        ax = axes[0, idx]
        sub = df_class[df_class["end_reason"] == er]
        margins = sub["margin"].dropna().values
        if len(margins) > 2:
            x, y = _kde_smooth(margins)
            ax.fill_between(x, y, alpha=0.3)
            ax.plot(x, y, linewidth=1.5)
        ax.set_title(er.replace("_", "\n"), fontsize=8)
        ax.set_xlabel("Margin")
        if idx == 0:
            ax.set_ylabel("Density")

    fig.suptitle("Classification Margin by End Reason")
    fig.tight_layout()
    _save(fig, outdir, "margin_by_endreason")
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_qc.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/qc.py tests/test_qc.py
git commit -m "feat(qc): add heatmap, per-read profiles, and end-reason faceted plots"
```

---

### Task 10: QC Orchestrator Function

Add a `run_qc()` function that generates all 10 plots from the TSV files in one call.

**Files:**
- Modify: `bin/qc.py`
- Modify: `tests/test_qc.py`

**Step 1: Write failing test**

Add to `tests/test_qc.py`:

```python
class TestRunQC:
    def test_generates_all_plots(self, tmp_path):
        from qc import run_qc, load_classification, load_alignments
        import pandas as pd

        df_class = _sample_classification_df(100)
        df_align = _sample_alignments_df(100, 2)

        class_tsv = tmp_path / "classification.tsv"
        align_tsv = tmp_path / "alignments.tsv"
        df_class.to_csv(class_tsv, sep="\t", index=False)
        df_align.to_csv(align_tsv, sep="\t", index=False)

        qc_dir = tmp_path / "qc"
        run_qc(align_tsv, class_tsv, qc_dir)

        expected = [
            "classification_heatmap.png",
            "margin_kde.png",
            "ned_by_target.png",
            "segmented_identity_violins.png",
            "length_vs_ned_scatter.png",
            "five_prime_quality.png",
            "indel_profile.png",
            "threshold_sweep.png",
            "per_read_profiles_sample.png",
            "margin_by_endreason.png",
        ]
        for name in expected:
            assert (qc_dir / name).exists(), f"Missing: {name}"
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_qc.py::TestRunQC -v`
Expected: FAIL

**Step 3: Implement run_qc**

Add to `bin/qc.py`:

```python
def run_qc(
    alignments_tsv: Path,
    classification_tsv: Path,
    outdir: Path,
    ref_lengths: dict[str, int] | None = None,
) -> list[Path]:
    """Generate all QC plots from alignment and classification TSVs.

    Returns list of generated plot file paths.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    df_align = load_alignments(alignments_tsv)
    df_class = load_classification(classification_tsv)

    plots = []

    plot_classification_heatmap(df_align, df_class, outdir)
    plots.append(outdir / "classification_heatmap.png")

    plot_margin_kde(df_class, outdir)
    plots.append(outdir / "margin_kde.png")

    plot_ned_by_target(df_align, df_class, outdir)
    plots.append(outdir / "ned_by_target.png")

    plot_segmented_violins(df_align, df_class, outdir)
    plots.append(outdir / "segmented_identity_violins.png")

    plot_length_vs_ned(df_class, outdir, ref_lengths=ref_lengths)
    plots.append(outdir / "length_vs_ned_scatter.png")

    plot_five_prime_quality(df_align, df_class, outdir)
    plots.append(outdir / "five_prime_quality.png")

    plot_indel_profile(df_align, df_class, outdir)
    plots.append(outdir / "indel_profile.png")

    plot_threshold_sweep(df_class, outdir)
    plots.append(outdir / "threshold_sweep.png")

    plot_per_read_profiles(df_align, outdir)
    plots.append(outdir / "per_read_profiles_sample.png")

    plot_end_reason_facets(df_align, df_class, outdir)
    plots.append(outdir / "margin_by_endreason.png")

    return plots
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_qc.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/qc.py tests/test_qc.py
git commit -m "feat(qc): add run_qc orchestrator for all 10 diagnostic plots"
```

---

### Task 11: Manifest and Database Initialization

Build the manifest (checkpoint/provenance) and DB init logic that `prepare.py` will use.

**Files:**
- Create: `bin/prepare.py`
- Create: `tests/test_prepare.py`

**Step 1: Write failing tests**

```python
# tests/test_prepare.py
"""Tests for prepare.py — manifest and DB initialization."""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest


class TestManifest:

    def test_create_and_save(self, tmp_path):
        from prepare import Manifest
        m = Manifest(exp_id="TEST_EXP", expdir="/fake/path")
        m.mark_stage("discover")
        path = tmp_path / "manifest.json"
        m.save(path)
        assert path.exists()
        data = json.loads(path.read_text())
        assert data["exp_id"] == "TEST_EXP"
        assert "discover" in data["stages_completed"]

    def test_load_existing(self, tmp_path):
        from prepare import Manifest
        m = Manifest(exp_id="TEST_EXP", expdir="/fake/path")
        m.mark_stage("discover")
        m.mark_stage("plan")
        path = tmp_path / "manifest.json"
        m.save(path)

        loaded = Manifest.load(path)
        assert loaded.exp_id == "TEST_EXP"
        assert "plan" in loaded.stages_completed

    def test_load_nonexistent_returns_none(self, tmp_path):
        from prepare import Manifest
        assert Manifest.load(tmp_path / "nope.json") is None


class TestInitDB:

    def test_creates_database_with_schema(self, tmp_path):
        from prepare import init_database
        db_path = tmp_path / "test.db"
        init_database(db_path, "TEST_EXP", "FC001", "S001", "alias1")
        assert db_path.exists()

        conn = sqlite3.connect(str(db_path))
        c = conn.cursor()
        # Check tables exist
        c.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = {row[0] for row in c.fetchall()}
        assert "Reads" in tables
        assert "Mods" in tables
        assert "Exp" in tables
        assert "Target" in tables
        assert "RunMetadata" in tables

        # Check Mods populated
        c.execute("SELECT COUNT(*) FROM Mods")
        assert c.fetchone()[0] == 10

        # Check Exp populated
        c.execute("SELECT exp_id FROM Exp")
        assert c.fetchone()[0] == "TEST_EXP"

        conn.close()

    def test_populate_run_metadata(self, tmp_path):
        from prepare import init_database, insert_run_metadata
        db_path = tmp_path / "test.db"
        init_database(db_path, "TEST_EXP", "FC001", "S001", "alias1")

        insert_run_metadata(db_path, {
            "run_id": "run1",
            "flow_cell_id": "FC001",
            "device_id": "MD-100",
            "sample_id": "S001",
            "experiment_id": "TEST_EXP",
            "kit": "SQK-NBD114-24",
            "basecall_model": "sup@v5.2.0",
            "source_bam_count": 100,
        })

        conn = sqlite3.connect(str(db_path))
        c = conn.cursor()
        c.execute("SELECT run_id, basecall_model FROM RunMetadata")
        row = c.fetchone()
        assert row[0] == "run1"
        assert row[1] == "sup@v5.2.0"
        conn.close()
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_prepare.py -v`
Expected: FAIL

**Step 3: Implement manifest and DB init**

```python
# bin/prepare.py
"""SMA-seq experiment preparation: discover, plan, merge, init, align+qc."""
from __future__ import annotations

import json
import os
import sqlite3
import time
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Manifest:
    """Checkpoint / provenance manifest for the prepare pipeline."""

    exp_id: str
    expdir: str
    timestamp: str = ""
    runs: list[dict] = field(default_factory=list)
    merge_groups: list[list[str]] = field(default_factory=list)
    output_bam: str = ""
    output_bam_reads: int = 0
    database_path: str = ""
    ref_path: str = ""
    ref_targets: list[str] = field(default_factory=list)
    stages_completed: list[str] = field(default_factory=list)
    qc_plots: list[str] = field(default_factory=list)

    def __post_init__(self):
        if not self.timestamp:
            self.timestamp = time.strftime("%Y-%m-%dT%H:%M:%S")

    def mark_stage(self, stage: str) -> None:
        if stage not in self.stages_completed:
            self.stages_completed.append(stage)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.__dict__, f, indent=2, default=str)

    @classmethod
    def load(cls, path: Path) -> "Manifest | None":
        if not path.exists():
            return None
        with open(path) as f:
            data = json.load(f)
        m = cls(exp_id=data["exp_id"], expdir=data["expdir"])
        for k, v in data.items():
            if hasattr(m, k):
                setattr(m, k, v)
        return m


def init_database(
    db_path: Path, exp_id: str, flow_cell_id: str, sample_id: str, alias: str,
) -> None:
    """Create SMA-seq SQLite database with full schema."""
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    c = conn.cursor()

    c.execute('''CREATE TABLE Reads (
        uniq_id TEXT PRIMARY KEY, exp_id TEXT, tgt_id TEXT, read_id TEXT,
        readseq TEXT, readlen INTEGER, model_tier TEXT, model_ver TEXT,
        trim INTEGER, mod_bitflag INTEGER, ed INTEGER, q_bc REAL, q_ld REAL,
        ER TEXT, bc_start_id TEXT, bc_start_ed INTEGER, bc_start_conf REAL,
        bc_end_id TEXT, bc_end_ed INTEGER, bc_end_conf REAL,
        trunc_level TEXT, signal_duration_s REAL, mean_qscore REAL,
        FOREIGN KEY(tgt_id) REFERENCES Target(tgt_id),
        FOREIGN KEY(mod_bitflag) REFERENCES Mods(mod_bitflag)
    )''')

    c.execute('''CREATE TABLE Mods (
        mod_bitflag INTEGER PRIMARY KEY, mods TEXT
    )''')

    c.execute('''CREATE TABLE Exp (
        exp_id TEXT PRIMARY KEY, flow_cell_id TEXT,
        sample_id TEXT, alias TEXT, exp_desc TEXT
    )''')

    c.execute('''CREATE TABLE Target (
        tgt_id TEXT PRIMARY KEY, tgt_refseq TEXT, tgt_reflen INTEGER
    )''')

    c.execute('''CREATE TABLE IF NOT EXISTS RunMetadata (
        run_id TEXT PRIMARY KEY, flow_cell_id TEXT, device_id TEXT,
        sample_id TEXT, experiment_id TEXT, kit TEXT,
        protocol_run_id TEXT, start_time TEXT, basecall_model TEXT,
        source_bam_count INTEGER, source_bam_paths TEXT, merge_timestamp TEXT
    )''')

    c.execute('''CREATE TABLE IF NOT EXISTS ReadRun (
        read_id TEXT PRIMARY KEY, run_id TEXT,
        FOREIGN KEY (run_id) REFERENCES RunMetadata(run_id)
    )''')

    mods_data = [
        (0, "non"), (1, "6mA"), (2, "5mCG_5hmCG"), (4, "5mC_5hmC"),
        (8, "4mC_5mC"), (16, "5mC"), (3, "6mA,5mCG_5hmCG"),
        (5, "6mA,5mC_5hmC"), (9, "6mA,4mC_5mC"), (17, "6mA,5mC"),
    ]
    c.executemany("INSERT OR IGNORE INTO Mods (mod_bitflag, mods) VALUES (?, ?)", mods_data)

    c.execute(
        "INSERT INTO Exp (exp_id, flow_cell_id, sample_id, alias, exp_desc) VALUES (?, ?, ?, ?, ?)",
        (exp_id, flow_cell_id, sample_id, alias, "Initialized via prepare"),
    )

    conn.commit()
    conn.close()


def insert_run_metadata(db_path: Path, run_meta: dict) -> None:
    """Insert a single RunMetadata row."""
    conn = sqlite3.connect(str(db_path))
    c = conn.cursor()
    c.execute(
        """INSERT OR REPLACE INTO RunMetadata
        (run_id, flow_cell_id, device_id, sample_id, experiment_id, kit,
         protocol_run_id, start_time, basecall_model, source_bam_count,
         source_bam_paths, merge_timestamp)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (
            run_meta.get("run_id", ""),
            run_meta.get("flow_cell_id", ""),
            run_meta.get("device_id", ""),
            run_meta.get("sample_id", ""),
            run_meta.get("experiment_id", ""),
            run_meta.get("kit", ""),
            run_meta.get("protocol_run_id", ""),
            run_meta.get("start_time", ""),
            run_meta.get("basecall_model", ""),
            run_meta.get("source_bam_count", 0),
            run_meta.get("source_bam_paths", ""),
            time.strftime("%Y-%m-%dT%H:%M:%S"),
        ),
    )
    conn.commit()
    conn.close()
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_prepare.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/prepare.py tests/test_prepare.py
git commit -m "feat(prepare): add manifest checkpoint and database initialization"
```

---

### Task 12: Prepare CLI — POD5 Symlinking and BAM Merge

Add POD5 consolidation (symlinks) and BAM merge functions, plus the merge plan display.

**Files:**
- Modify: `bin/prepare.py`
- Modify: `tests/test_prepare.py`

**Step 1: Write failing tests**

Add to `tests/test_prepare.py`:

```python
class TestSymlinkPod5:

    def test_creates_symlinks(self, tmp_path):
        from prepare import symlink_pod5s
        from sma_merge.models import RunInfo

        # Create fake POD5 files
        run_dir = tmp_path / "run1" / "pod5_pass"
        run_dir.mkdir(parents=True)
        (run_dir / "chunk1.pod5").write_text("fake")
        (run_dir / "chunk2.pod5").write_text("fake")

        run = RunInfo(
            run_dir=run_dir.parent, flow_cell_id="FC1", device_id="D1",
            protocol_group_id="p1", basecall_model="sup", sample_id="",
            run_id="run1_abc", sample_rate=5000, pod5_dir=run_dir,
            pod5_count=2, mod_base_models="",
        )

        out = tmp_path / "output" / "pod5"
        symlink_pod5s([run], out)
        links = list(out.rglob("*.pod5"))
        assert len(links) == 2
        assert all(link.is_symlink() for link in links)


class TestCollectBams:

    def test_finds_all_bams(self, tmp_path):
        from prepare import collect_bam_files
        from sma_merge.models import RunInfo

        run_dir = tmp_path / "run1"
        bam_dir = run_dir / "bam_pass" / "barcode01"
        bam_dir.mkdir(parents=True)
        (bam_dir / "chunk1.bam").write_text("fake")
        (bam_dir / "chunk2.bam").write_text("fake")
        unclass = run_dir / "bam_pass" / "unclassified"
        unclass.mkdir(parents=True)
        (unclass / "chunk3.bam").write_text("fake")

        run = RunInfo(
            run_dir=run_dir, flow_cell_id="FC1", device_id="D1",
            protocol_group_id="p1", basecall_model="sup", sample_id="",
            run_id="run1", sample_rate=5000, pod5_dir=run_dir / "pod5_pass",
            pod5_count=0, mod_base_models="",
        )

        bams = collect_bam_files([run])
        assert len(bams) == 3


class TestFormatMergePlan:

    def test_formats_plan(self):
        from prepare import format_merge_plan
        from sma_merge.models import RunGroup, RunInfo
        from pathlib import Path

        run = RunInfo(
            run_dir=Path("/tmp"), flow_cell_id="FC1", device_id="D1",
            protocol_group_id="p1", basecall_model="sup@v5.2.0",
            sample_id="", run_id="run1", sample_rate=5000,
            pod5_dir=Path("/tmp/pod5"), pod5_count=10, mod_base_models="",
        )
        group = RunGroup(flow_cell_id="FC1", runs=[run],
                        basecall_model="sup@v5.2.0", is_consistent=True)

        text = format_merge_plan([group], "TEST_EXP")
        assert "FC1" in text
        assert "REUSE" in text or "sup" in text
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_prepare.py::TestSymlinkPod5 tests/test_prepare.py::TestCollectBams tests/test_prepare.py::TestFormatMergePlan -v`
Expected: FAIL

**Step 3: Implement POD5 symlinking, BAM collection, and merge plan**

Add to `bin/prepare.py`:

```python
from sma_merge.models import RunInfo, RunGroup
from sma_merge.validate import validate_runs


def symlink_pod5s(runs: list[RunInfo], output_dir: Path) -> None:
    """Create organized symlinks to POD5 files from each run."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for run in runs:
        run_subdir = output_dir / run.run_id
        run_subdir.mkdir(exist_ok=True)
        for pod5_file in sorted(run.pod5_dir.rglob("*.pod5")):
            link = run_subdir / pod5_file.name
            if link.exists() or link.is_symlink():
                link.unlink()
            link.symlink_to(pod5_file.resolve())


def collect_bam_files(runs: list[RunInfo]) -> list[Path]:
    """Collect all BAM files from bam_pass/ across runs."""
    bams: list[Path] = []
    for run in runs:
        bam_pass = run.run_dir / "bam_pass"
        if bam_pass.is_dir():
            bams.extend(sorted(bam_pass.rglob("*.bam")))
    return bams


def format_merge_plan(groups: list[RunGroup], exp_id: str) -> str:
    """Format the merge plan for user display."""
    lines = [
        f"[prepare] Merge Plan for {exp_id}",
        "-" * (28 + len(exp_id)),
    ]
    for g in groups:
        for i, run in enumerate(g.runs, 1):
            pod5_size_gb = 0
            for f in run.pod5_dir.rglob("*.pod5"):
                try:
                    pod5_size_gb += f.stat().st_size
                except OSError:
                    pass
            pod5_size_gb /= 1e9

            model_short = run.basecall_model.split("@")[-1] if "@" in run.basecall_model else run.basecall_model
            lines.append(f"Run {i}: {run.run_id}")
            lines.append(f"  Flow Cell: {g.flow_cell_id} | POD5: {run.pod5_count} files ({pod5_size_gb:.0f} GB) | Model: {model_short}")
            action = "REUSE existing BAMs" if g.is_consistent else "REBASECALL from POD5"
            lines.append(f"  Action: {action}")
            lines.append("")

        if len(g.runs) > 1:
            lines.append(f"Merge group: [{' + '.join(f'Run {i+1}' for i in range(len(g.runs)))}] -> same flow cell")
        for issue in g.issues:
            lines.append(f"  NOTE: {issue}")
        lines.append("")

    return "\n".join(lines)
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_prepare.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/prepare.py tests/test_prepare.py
git commit -m "feat(prepare): add POD5 symlinking, BAM collection, and merge plan display"
```

---

### Task 13: Prepare CLI Entry Point

Wire up the argparse CLI and the main `run_prepare()` function that orchestrates all stages.

**Files:**
- Modify: `bin/prepare.py`
- Modify: `tests/test_prepare.py`

**Step 1: Write failing tests**

Add to `tests/test_prepare.py`:

```python
class TestBuildParser:

    def test_all_required_args(self):
        from prepare import build_parser
        parser = build_parser()
        # Should fail without required args
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parses_valid_args(self):
        from prepare import build_parser
        parser = build_parser()
        args = parser.parse_args([
            "-d", "/fake/exp", "-e", "FC1_S1_A", "-r", "/fake/ref.fa",
        ])
        assert args.expdir == Path("/fake/exp")
        assert args.expid == "FC1_S1_A"
        assert args.ref == Path("/fake/ref.fa")
        assert args.outdir == Path("Output")  # default

    def test_optional_flags(self):
        from prepare import build_parser
        parser = build_parser()
        args = parser.parse_args([
            "-d", "/x", "-e", "E1", "-r", "/r.fa",
            "--force-rebasecall", "--dry-run", "-o", "/out",
        ])
        assert args.force_rebasecall is True
        assert args.dry_run is True
        assert args.outdir == Path("/out")
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_prepare.py::TestBuildParser -v`
Expected: FAIL

**Step 3: Implement CLI entry point**

Add to `bin/prepare.py`:

```python
import argparse
import sys


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="prepare",
        description="SMA-seq experiment preparation: discover, plan, merge, init, align+qc.",
    )
    parser.add_argument("-d", "--expdir", type=Path, required=True,
                        help="Path to experiment directory")
    parser.add_argument("-e", "--expid", required=True,
                        help="Experiment ID (FlowCell_Sample_Alias)")
    parser.add_argument("-r", "--ref", type=Path, required=True,
                        help="Multi-sequence FASTA with target sequences")
    parser.add_argument("-o", "--outdir", type=Path, default=Path("Output"),
                        help="Output directory (default: Output)")
    parser.add_argument("--force-rebasecall", action="store_true",
                        help="Skip smart detection, always re-basecall from POD5")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print merge plan without executing")
    return parser


def run_prepare(args: argparse.Namespace) -> None:
    """Execute the full prepare pipeline."""
    from sma_merge.discover import discover_runs
    from align import parse_fasta, process_bam
    from qc import run_qc

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / "prepare_manifest.json"

    # Check for existing manifest (resume)
    existing = Manifest.load(manifest_path)
    if existing and existing.exp_id == args.expid:
        print(f"[prepare] Found existing manifest. Stages completed: {existing.stages_completed}")
        manifest = existing
    else:
        manifest = Manifest(exp_id=args.expid, expdir=str(args.expdir))

    # --- Stage 1: DISCOVER ---
    if "discover" not in manifest.stages_completed:
        print(f"[prepare] Stage 1: Discovering runs in {args.expdir}...")
        runs = discover_runs(args.expdir)
        if not runs:
            sys.exit("[prepare] No MinKNOW runs found.")
        print(f"[prepare] Found {len(runs)} run(s).")
        manifest.runs = [
            {"run_id": r.run_id, "flow_cell_id": r.flow_cell_id,
             "device_id": r.device_id, "pod5_count": r.pod5_count,
             "basecall_model": r.basecall_model}
            for r in runs
        ]
        manifest.mark_stage("discover")
        manifest.save(manifest_path)
    else:
        print("[prepare] Stage 1: DISCOVER (cached)")
        runs = discover_runs(args.expdir)

    # --- Stage 2: PLAN ---
    if "plan" not in manifest.stages_completed:
        groups = validate_runs(runs)
        plan_text = format_merge_plan(groups, args.expid)
        print(plan_text)

        if args.dry_run:
            print("[prepare] Dry run — stopping after plan.")
            return

        response = input("Proceed? [Y/n] ").strip().lower()
        if response and response != "y":
            print("[prepare] Aborted by user.")
            return

        manifest.mark_stage("plan")
        manifest.save(manifest_path)
    else:
        print("[prepare] Stage 2: PLAN (cached)")
        groups = validate_runs(runs)

    # --- Stage 3: MERGE ---
    if "merge" not in manifest.stages_completed:
        print("[prepare] Stage 3: Merging...")
        for g in groups:
            # 3a: POD5 symlinks
            pod5_out = outdir / "pod5"
            symlink_pod5s(g.runs, pod5_out)

            # 3b: BAM merge
            if g.is_consistent and not args.force_rebasecall:
                import pysam
                bams = collect_bam_files(g.runs)
                if not bams:
                    print(f"[prepare] Warning: No BAMs found for {g.flow_cell_id}. Need re-basecall.")
                    continue
                merged_bam = outdir / "merged.bam"
                print(f"[prepare] Merging {len(bams)} BAM files...")
                pysam.merge("-f", "-n", str(merged_bam), *[str(b) for b in bams])
                print(f"[prepare] Sorting...")
                sorted_bam = outdir / "merged_sorted.bam"
                pysam.sort("-o", str(sorted_bam), str(merged_bam))
                merged_bam.unlink()
                sorted_bam.rename(merged_bam)
                print(f"[prepare] Indexing...")
                pysam.index(str(merged_bam))
                manifest.output_bam = str(merged_bam)
            else:
                from sma_merge.basecall import merge_pod5s, basecall
                print(f"[prepare] Re-basecalling from POD5...")
                merged_pod5 = outdir / f"{g.flow_cell_id}_merged.pod5"
                merge_pod5s([r.pod5_dir for r in g.runs], merged_pod5)
                merged_bam = outdir / "merged.bam"
                basecall(merged_pod5, merged_bam, g.basecall_model)
                manifest.output_bam = str(merged_bam)

        manifest.mark_stage("merge")
        manifest.save(manifest_path)
    else:
        print("[prepare] Stage 3: MERGE (cached)")

    # --- Stage 4: INIT ---
    if "init" not in manifest.stages_completed:
        print("[prepare] Stage 4: Initializing database...")
        exp_parts = args.expid.split("_")
        fc = exp_parts[0] if len(exp_parts) >= 1 else ""
        sid = exp_parts[1] if len(exp_parts) >= 2 else ""
        alias = exp_parts[-1] if len(exp_parts) >= 3 else ""

        db_path = outdir / f"SMA_{args.expid}.db"
        init_database(db_path, args.expid, fc, sid, alias)

        for g in groups:
            for run in g.runs:
                insert_run_metadata(db_path, {
                    "run_id": run.run_id,
                    "flow_cell_id": run.flow_cell_id,
                    "device_id": run.device_id,
                    "sample_id": run.sample_id,
                    "experiment_id": args.expid,
                    "basecall_model": run.basecall_model,
                    "source_bam_count": run.pod5_count,
                })

        manifest.database_path = str(db_path)
        manifest.mark_stage("init")
        manifest.save(manifest_path)
    else:
        print("[prepare] Stage 4: INIT (cached)")

    # --- Stage 5: ALIGN + QC ---
    if "align" not in manifest.stages_completed:
        print("[prepare] Stage 5a: Aligning all reads against all targets...")
        refs = parse_fasta(args.ref)
        manifest.ref_path = str(args.ref)
        manifest.ref_targets = list(refs.keys())

        merged_bam = Path(manifest.output_bam)
        align_tsv = outdir / "alignments.tsv"
        class_tsv = outdir / "classification.tsv"
        n = process_bam(merged_bam, refs, align_tsv, class_tsv)
        print(f"[prepare] Processed {n} reads.")
        manifest.output_bam_reads = n
        manifest.mark_stage("align")
        manifest.save(manifest_path)
    else:
        print("[prepare] Stage 5a: ALIGN (cached)")

    if "qc" not in manifest.stages_completed:
        print("[prepare] Stage 5b: Generating QC plots...")
        qc_dir = outdir / "qc"
        align_tsv = outdir / "alignments.tsv"
        class_tsv = outdir / "classification.tsv"
        ref_lengths = {k: len(v) for k, v in parse_fasta(args.ref).items()} if args.ref.exists() else None
        plots = run_qc(align_tsv, class_tsv, qc_dir, ref_lengths=ref_lengths)
        manifest.qc_plots = [str(p) for p in plots]
        manifest.mark_stage("qc")
        manifest.save(manifest_path)
    else:
        print("[prepare] Stage 5b: QC (cached)")

    print(f"\n[prepare] Complete. Outputs in {outdir}/")
    print(f"  - Merged BAM: {manifest.output_bam}")
    print(f"  - Database:   {manifest.database_path}")
    print(f"  - Alignments: {outdir}/alignments.tsv")
    print(f"  - Classified: {outdir}/classification.tsv")
    print(f"  - QC Plots:   {outdir}/qc/")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    run_prepare(args)


if __name__ == "__main__":
    main()
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_prepare.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/prepare.py tests/test_prepare.py
git commit -m "feat(prepare): add CLI entry point and full pipeline orchestration"
```

---

### Task 14: Integration Test

End-to-end test with synthetic data: create a fake experiment dir, fake BAMs, fake reference, and run the full pipeline.

**Files:**
- Create: `tests/test_prepare_integration.py`

**Step 1: Write the integration test**

```python
# tests/test_prepare_integration.py
"""Integration test for the full prepare pipeline."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pysam
import pytest


def _make_fake_experiment(tmp_path: Path) -> tuple[Path, Path]:
    """Create a minimal fake experiment directory with BAMs and a reference."""
    exp_dir = tmp_path / "experiment"
    run_dir = exp_dir / "no_sample" / "20260101_1200_MD-100_FC001_abc123"
    pod5_dir = run_dir / "pod5_pass"
    bam_dir = run_dir / "bam_pass" / "unclassified"
    pod5_dir.mkdir(parents=True)
    bam_dir.mkdir(parents=True)

    # Fake POD5 (just a marker file — discover will need to be mocked)
    (pod5_dir / "chunk.pod5").write_bytes(b"")

    # Create real BAMs
    seq_a = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32bp
    seq_b = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
    bam_path = bam_dir / "chunk1.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
        for i, seq in enumerate([seq_a, seq_b, seq_a, seq_b]):
            a = pysam.AlignedSegment(header)
            a.query_name = f"read_{i}"
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.flag = 4
            a.set_tag("er", "signal_positive", "Z")
            f.write(a)

    # Reference FASTA
    ref_path = tmp_path / "targets.fa"
    ref_path.write_text(">targetA\nACGTACGTACGTACGTACGTACGTACGTACGT\n>targetB\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n")

    return exp_dir, ref_path


def test_full_pipeline_produces_outputs(tmp_path, monkeypatch):
    """Test that prepare produces all expected output files."""
    from prepare import Manifest, init_database, insert_run_metadata, collect_bam_files
    from prepare import symlink_pod5s, format_merge_plan
    from align import parse_fasta, process_bam
    from qc import run_qc
    from sma_merge.models import RunInfo, RunGroup

    exp_dir, ref_path = _make_fake_experiment(tmp_path)
    outdir = tmp_path / "output"
    outdir.mkdir()

    # Simulate discovery (we can't use real POD5 discovery in tests)
    run_dir = list((exp_dir / "no_sample").iterdir())[0]
    run = RunInfo(
        run_dir=run_dir, flow_cell_id="FC001", device_id="MD-100",
        protocol_group_id="test", basecall_model="sup@v5.2.0",
        sample_id="", run_id="20260101_1200_MD-100_FC001_abc123",
        sample_rate=5000, pod5_dir=run_dir / "pod5_pass",
        pod5_count=1, mod_base_models="",
    )
    group = RunGroup(flow_cell_id="FC001", runs=[run],
                     basecall_model="sup@v5.2.0", is_consistent=True)

    # Stage 3: Merge BAMs (REUSE path — just merge existing BAMs)
    bams = collect_bam_files([run])
    assert len(bams) == 1

    merged_bam = outdir / "merged.bam"
    pysam.merge("-f", str(merged_bam), *[str(b) for b in bams])

    # Stage 4: Init DB
    db_path = outdir / "SMA_TEST.db"
    init_database(db_path, "TEST_EXP", "FC001", "S001", "alias")

    # Stage 5a: Align
    refs = parse_fasta(ref_path)
    align_tsv = outdir / "alignments.tsv"
    class_tsv = outdir / "classification.tsv"
    n = process_bam(merged_bam, refs, align_tsv, class_tsv)
    assert n == 4

    # Check classification
    import csv
    with open(class_tsv) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == 4
    # Reads 0,2 should match targetA, reads 1,3 should match targetB
    for row in rows:
        if row["read_id"] in ("read_0", "read_2"):
            assert row["assigned_ref"] == "targetA"
        else:
            assert row["assigned_ref"] == "targetB"

    # Stage 5b: QC
    qc_dir = outdir / "qc"
    plots = run_qc(align_tsv, class_tsv, qc_dir)
    assert len(plots) >= 10
    for p in plots:
        assert p.exists()

    # Check DB
    conn = sqlite3.connect(str(db_path))
    c = conn.cursor()
    c.execute("SELECT COUNT(*) FROM Mods")
    assert c.fetchone()[0] == 10
    conn.close()
```

**Step 2: Run the integration test**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_prepare_integration.py -v`
Expected: PASS

**Step 3: Commit**

```bash
cd /tmp/ont-sma-seq
git add tests/test_prepare_integration.py
git commit -m "test(prepare): add end-to-end integration test with synthetic data"
```

---

### Task 15: Final Cleanup and Full Test Suite

Run the full test suite, fix any issues, update README.

**Files:**
- Modify: `README.md`

**Step 1: Run full test suite**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/ -v --tb=short`
Expected: All tests PASS (including existing 361 tests + new tests)

**Step 2: Update README with prepare.py usage**

Add `prepare.py` section to `README.md` between `extractMeta.py` and `ingest.py` sections. Include:

```markdown
### `prepare.py`

Discovers MinKNOW runs, merges POD5/BAM data, initializes the database, aligns all reads against target references, and generates diagnostic QC plots.

* `-d`, `--expdir`: Path to experiment directory (Required).
* `-e`, `--expid`: Experiment ID (Required).
* `-r`, `--ref`: Multi-sequence FASTA with target sequences (Required).
* `-o`, `--outdir`: Output directory (Default: `Output`).
* `--force-rebasecall`: Always re-basecall from POD5.
* `--dry-run`: Print merge plan without executing.

\```bash
python3 bin/prepare.py -d /path/to/experiment -e EXP_ID -r targets.fa -o Output/
\```
```

**Step 3: Commit**

```bash
cd /tmp/ont-sma-seq
git add README.md
git commit -m "docs: add prepare.py usage to README"
```

---

## Task Summary

| Task | Description | New Files | Test File |
|------|-------------|-----------|-----------|
| 1 | CIGAR parsing + whole-alignment metrics | `bin/align.py` | `tests/test_align.py` |
| 2 | Segmented 5'/mid/3' metrics | (modify align.py) | (modify test_align.py) |
| 3 | Terminal quality + indel positions | (modify align.py) | (modify test_align.py) |
| 4 | All-vs-all alignment + classification | (modify align.py) | (modify test_align.py) |
| 5 | BAM streaming + TSV output | (modify align.py) | (modify test_align.py) |
| 6 | QC core utilities + margin KDE | `bin/qc.py` | `tests/test_qc.py` |
| 7 | NED-by-target, length-vs-NED, violins | (modify qc.py) | (modify test_qc.py) |
| 8 | 5' quality, indels, threshold sweep | (modify qc.py) | (modify test_qc.py) |
| 9 | Heatmap, per-read profiles, ER facets | (modify qc.py) | (modify test_qc.py) |
| 10 | QC orchestrator run_qc() | (modify qc.py) | (modify test_qc.py) |
| 11 | Manifest + DB initialization | `bin/prepare.py` | `tests/test_prepare.py` |
| 12 | POD5 symlinking + BAM merge | (modify prepare.py) | (modify test_prepare.py) |
| 13 | CLI entry point + pipeline orchestration | (modify prepare.py) | (modify test_prepare.py) |
| 14 | Integration test | | `tests/test_prepare_integration.py` |
| 15 | Full test suite + README update | (modify README.md) | |
