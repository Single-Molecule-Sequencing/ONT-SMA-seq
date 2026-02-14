# ONT-Native Configuration & Truncated Molecule Detection Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace ad-hoc CLI flags with Dorado-compatible TOML construct config, add 5-level truncated molecule detection, auto-generate truncated references, and provide an interactive construct wizard.

**Architecture:** Extended Dorado TOML arrangement format with `[sma]` section for barcode pairing, confidence thresholds, and truncation rules. A `construct.py` module parses and validates the TOML. `mkrefs.py` generates truncated reference FASTAs. `ingest.py` uses structural element detection + confidence gating to classify each read into one of five truncation levels. An interactive `construct_wizard.py` guides users through TOML creation.

**Tech Stack:** Python 3.11+, tomllib (stdlib), tomli_w, edlib, pysam, sqlite3, pytest

**Design doc:** `docs/plans/2026-02-14-ont-native-config-design.md`

---

## Task 1: Construct TOML Parser (`bin/construct.py`)

Parse and validate the extended Dorado TOML arrangement format. This is the foundation module that all other tasks depend on.

**Files:**
- Create: `bin/construct.py`
- Create: `tests/test_construct.py`

**Step 1: Write the failing tests**

```python
# tests/test_construct.py
"""Tests for the construct TOML parser."""

import textwrap
from pathlib import Path

import pytest

from construct import ConstructConfig, parse_construct_toml, ValidationError


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

MINIMAL_DUAL_TOML = textwrap.dedent("""\
    [arrangement]
    name = "test_dual"
    kit = "SQK-NBD114-96"
    mask1_front = "AAGGTTAA"
    mask1_rear  = "CAGCACCT"
    mask2_front = "AGGTGCTG"
    mask2_rear  = "TTAACCTT"
    barcode1_pattern = "NB%02i"
    barcode2_pattern = "NB%02i"
    first_index = 1
    last_index = 96

    [scoring]
    max_barcode_penalty = 11
    min_barcode_penalty_dist = 3
    front_barcode_window = 100
    rear_barcode_window = 100

    [sma]
    mode = "dual_independent"

    [sma.confidence]
    full_length_threshold = 0.75
    start_barcode_min = 0.6
    flank_max_error_rate = 0.5

    [sma.truncation]
    auto_generate_refs = true
    min_target_length = 20

    [[sma.targets]]
    barcode1 = "NB05"
    barcode2 = "NB10"
    alias = "target_A"
    reference = "references/target_A.fasta"
""")

MINIMAL_START_ONLY_TOML = textwrap.dedent("""\
    [arrangement]
    name = "test_start_only"
    kit = "SQK-NBD114-96"
    mask1_front = "AAGGTTAA"
    mask1_rear  = "CAGCACCT"
    barcode1_pattern = "NB%02i"
    first_index = 1
    last_index = 96

    [sma]
    mode = "start_only"

    [[sma.targets]]
    barcode1 = "NB05"
    alias = "target_A"
    reference = "references/target_A.fasta"
""")


@pytest.fixture
def dual_toml(tmp_path):
    p = tmp_path / "construct.toml"
    p.write_text(MINIMAL_DUAL_TOML)
    return p


@pytest.fixture
def start_only_toml(tmp_path):
    p = tmp_path / "construct.toml"
    p.write_text(MINIMAL_START_ONLY_TOML)
    return p


# ---------------------------------------------------------------------------
# Tests: parse_construct_toml
# ---------------------------------------------------------------------------

class TestParseConstructToml:

    def test_parses_dual_mode(self, dual_toml):
        cfg = parse_construct_toml(dual_toml)
        assert cfg.arrangement.name == "test_dual"
        assert cfg.arrangement.kit == "SQK-NBD114-96"
        assert cfg.sma.mode == "dual_independent"
        assert len(cfg.sma.targets) == 1
        assert cfg.sma.targets[0].barcode1 == "NB05"
        assert cfg.sma.targets[0].barcode2 == "NB10"

    def test_parses_start_only_mode(self, start_only_toml):
        cfg = parse_construct_toml(start_only_toml)
        assert cfg.sma.mode == "start_only"
        assert cfg.sma.targets[0].barcode2 is None

    def test_arrangement_flanks(self, dual_toml):
        cfg = parse_construct_toml(dual_toml)
        assert cfg.arrangement.mask1_front == "AAGGTTAA"
        assert cfg.arrangement.mask1_rear == "CAGCACCT"
        assert cfg.arrangement.mask2_front == "AGGTGCTG"
        assert cfg.arrangement.mask2_rear == "TTAACCTT"

    def test_scoring_defaults(self, start_only_toml):
        """Scoring section is optional — defaults should be applied."""
        cfg = parse_construct_toml(start_only_toml)
        assert cfg.scoring.front_barcode_window == 100
        assert cfg.scoring.max_barcode_penalty == 11

    def test_confidence_defaults(self, start_only_toml):
        cfg = parse_construct_toml(start_only_toml)
        assert cfg.sma.confidence.full_length_threshold == 0.75
        assert cfg.sma.confidence.start_barcode_min == 0.6

    def test_truncation_defaults(self, start_only_toml):
        cfg = parse_construct_toml(start_only_toml)
        assert cfg.sma.truncation.auto_generate_refs is True
        assert cfg.sma.truncation.min_target_length == 20

    def test_barcode_pair_to_alias(self, dual_toml):
        cfg = parse_construct_toml(dual_toml)
        mapping = cfg.barcode_pair_to_alias()
        assert mapping == {("nb05", "nb10"): "target_A"}

    def test_used_barcodes(self, dual_toml):
        cfg = parse_construct_toml(dual_toml)
        used = cfg.used_barcode_ids()
        assert used == {"nb05", "nb10"}

    def test_flank_sequences(self, dual_toml):
        cfg = parse_construct_toml(dual_toml)
        assert cfg.flank_front == "AAGGTTAA"
        assert cfg.flank_rear == "CAGCACCT"


class TestConstructValidation:

    def test_rejects_missing_arrangement(self, tmp_path):
        p = tmp_path / "bad.toml"
        p.write_text('[sma]\nmode = "dual_independent"\n')
        with pytest.raises(ValidationError):
            parse_construct_toml(p)

    def test_rejects_missing_sma_section(self, tmp_path):
        p = tmp_path / "bad.toml"
        p.write_text(textwrap.dedent("""\
            [arrangement]
            name = "x"
            kit = "SQK-NBD114-96"
            mask1_front = "AAGG"
            mask1_rear = "CCTT"
            barcode1_pattern = "NB%02i"
            first_index = 1
            last_index = 96
        """))
        with pytest.raises(ValidationError):
            parse_construct_toml(p)

    def test_rejects_invalid_mode(self, tmp_path):
        toml_text = MINIMAL_DUAL_TOML.replace(
            'mode = "dual_independent"', 'mode = "triple_barcode"'
        )
        p = tmp_path / "bad.toml"
        p.write_text(toml_text)
        with pytest.raises(ValidationError):
            parse_construct_toml(p)

    def test_rejects_dual_target_without_barcode2(self, tmp_path):
        toml_text = MINIMAL_DUAL_TOML.replace(
            'barcode2 = "NB10"\n', ''
        )
        p = tmp_path / "bad.toml"
        p.write_text(toml_text)
        with pytest.raises(ValidationError):
            parse_construct_toml(p)

    def test_rejects_invalid_barcode_id(self, tmp_path):
        toml_text = MINIMAL_DUAL_TOML.replace(
            'barcode1 = "NB05"', 'barcode1 = "NB999"'
        )
        p = tmp_path / "bad.toml"
        p.write_text(toml_text)
        with pytest.raises(ValidationError):
            parse_construct_toml(p)

    def test_rejects_non_dna_flank(self, tmp_path):
        toml_text = MINIMAL_DUAL_TOML.replace(
            'mask1_front = "AAGGTTAA"', 'mask1_front = "AAGGXXAA"'
        )
        p = tmp_path / "bad.toml"
        p.write_text(toml_text)
        with pytest.raises(ValidationError):
            parse_construct_toml(p)
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_construct.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'construct'`

**Step 3: Write the implementation**

```python
# bin/construct.py
"""Parse and validate extended Dorado TOML construct configuration.

The construct TOML uses Dorado's standard [arrangement] and [scoring]
sections, extended with an [sma] section for barcode pairing, confidence
thresholds, and truncation rules.  Dorado ignores the [sma] section.
"""

import re
import tomllib
from dataclasses import dataclass, field
from pathlib import Path


_DNA_RE = re.compile(r"^[ACGTacgt]+$")
_BARCODE_ID_RE = re.compile(r"^[Nn][Bb](\d+)$")
_VALID_MODES = {"dual_independent", "start_only"}


class ValidationError(Exception):
    """Raised when construct TOML fails validation."""


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ArrangementConfig:
    name: str
    kit: str
    mask1_front: str
    mask1_rear: str
    barcode1_pattern: str
    first_index: int
    last_index: int
    mask2_front: str | None = None
    mask2_rear: str | None = None
    barcode2_pattern: str | None = None


@dataclass
class ScoringConfig:
    max_barcode_penalty: int = 11
    min_barcode_penalty_dist: int = 3
    front_barcode_window: int = 100
    rear_barcode_window: int = 100
    min_flank_score: float = 0.5


@dataclass
class ConfidenceConfig:
    full_length_threshold: float = 0.75
    start_barcode_min: float = 0.6
    flank_max_error_rate: float = 0.5


@dataclass
class TruncationConfig:
    auto_generate_refs: bool = True
    min_target_length: int = 20


@dataclass
class TargetEntry:
    barcode1: str
    alias: str
    reference: str
    barcode2: str | None = None


@dataclass
class SmaConfig:
    mode: str
    targets: list[TargetEntry]
    confidence: ConfidenceConfig = field(default_factory=ConfidenceConfig)
    truncation: TruncationConfig = field(default_factory=TruncationConfig)
    barcode_fasta: str | None = None


@dataclass
class ConstructConfig:
    arrangement: ArrangementConfig
    scoring: ScoringConfig
    sma: SmaConfig

    @property
    def flank_front(self) -> str:
        return self.arrangement.mask1_front

    @property
    def flank_rear(self) -> str:
        return self.arrangement.mask1_rear

    def barcode_pair_to_alias(self) -> dict[tuple[str, str], str]:
        """Return mapping compatible with sample_sheet.parse_sample_sheet()."""
        result = {}
        for t in self.sma.targets:
            bc1 = _normalize_bc(t.barcode1)
            if t.barcode2:
                bc2 = _normalize_bc(t.barcode2)
                result[(bc1, bc2)] = t.alias
            else:
                result[(bc1, "")] = t.alias
        return result

    def used_barcode_ids(self) -> set[str]:
        """Return set of all barcode IDs referenced in targets."""
        ids = set()
        for t in self.sma.targets:
            ids.add(_normalize_bc(t.barcode1))
            if t.barcode2:
                ids.add(_normalize_bc(t.barcode2))
        return ids


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _normalize_bc(bc_id: str) -> str:
    """Normalize 'NB05' or 'nb05' to 'nb05'."""
    m = _BARCODE_ID_RE.match(bc_id)
    if not m:
        raise ValidationError(f"Invalid barcode ID: {bc_id!r}")
    num = int(m.group(1))
    if num < 1 or num > 96:
        raise ValidationError(f"Barcode number out of range 1-96: {bc_id!r}")
    return f"nb{num:02d}"


def _validate_dna(seq: str, field_name: str) -> None:
    if not _DNA_RE.match(seq):
        raise ValidationError(
            f"{field_name} contains non-DNA characters: {seq!r}"
        )


def _require(data: dict, key: str, context: str) -> object:
    if key not in data:
        raise ValidationError(f"Missing required field '{key}' in {context}")
    return data[key]


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def parse_construct_toml(path: Path | str) -> ConstructConfig:
    """Parse and validate a construct TOML file.

    Args:
        path: Path to the TOML file.

    Returns:
        Validated ConstructConfig.

    Raises:
        ValidationError: If the TOML is invalid or missing required fields.
    """
    path = Path(path)
    with open(path, "rb") as f:
        raw = tomllib.load(f)

    # --- [arrangement] (required) ---
    if "arrangement" not in raw:
        raise ValidationError("Missing required [arrangement] section")
    arr_raw = raw["arrangement"]

    for fld in ("name", "kit", "mask1_front", "mask1_rear",
                "barcode1_pattern", "first_index", "last_index"):
        _require(arr_raw, fld, "[arrangement]")

    for fld in ("mask1_front", "mask1_rear"):
        _validate_dna(arr_raw[fld], f"[arrangement].{fld}")
    for fld in ("mask2_front", "mask2_rear"):
        if fld in arr_raw:
            _validate_dna(arr_raw[fld], f"[arrangement].{fld}")

    arrangement = ArrangementConfig(
        name=arr_raw["name"],
        kit=arr_raw["kit"],
        mask1_front=arr_raw["mask1_front"].upper(),
        mask1_rear=arr_raw["mask1_rear"].upper(),
        barcode1_pattern=arr_raw["barcode1_pattern"],
        first_index=arr_raw["first_index"],
        last_index=arr_raw["last_index"],
        mask2_front=arr_raw.get("mask2_front", "").upper() or None,
        mask2_rear=arr_raw.get("mask2_rear", "").upper() or None,
        barcode2_pattern=arr_raw.get("barcode2_pattern"),
    )

    # --- [scoring] (optional, all defaults) ---
    scr_raw = raw.get("scoring", {})
    scoring = ScoringConfig(**{
        k: scr_raw[k] for k in (
            "max_barcode_penalty", "min_barcode_penalty_dist",
            "front_barcode_window", "rear_barcode_window", "min_flank_score",
        ) if k in scr_raw
    })

    # --- [sma] (required) ---
    if "sma" not in raw:
        raise ValidationError("Missing required [sma] section")
    sma_raw = raw["sma"]

    mode = _require(sma_raw, "mode", "[sma]")
    if mode not in _VALID_MODES:
        raise ValidationError(
            f"Invalid [sma].mode: {mode!r}. Must be one of {_VALID_MODES}"
        )

    # Confidence
    conf_raw = sma_raw.get("confidence", {})
    confidence = ConfidenceConfig(**{
        k: conf_raw[k] for k in (
            "full_length_threshold", "start_barcode_min", "flank_max_error_rate",
        ) if k in conf_raw
    })

    # Truncation
    trunc_raw = sma_raw.get("truncation", {})
    truncation = TruncationConfig(**{
        k: trunc_raw[k] for k in (
            "auto_generate_refs", "min_target_length",
        ) if k in trunc_raw
    })

    # Targets
    targets_raw = sma_raw.get("targets", [])
    if not targets_raw:
        raise ValidationError("No [[sma.targets]] entries defined")

    targets = []
    for i, t in enumerate(targets_raw):
        bc1 = _require(t, "barcode1", f"[[sma.targets]][{i}]")
        _normalize_bc(bc1)  # validate
        alias = _require(t, "alias", f"[[sma.targets]][{i}]")
        ref = _require(t, "reference", f"[[sma.targets]][{i}]")
        bc2 = t.get("barcode2")
        if mode == "dual_independent" and bc2 is None:
            raise ValidationError(
                f"[[sma.targets]][{i}]: barcode2 is required in dual_independent mode"
            )
        if bc2:
            _normalize_bc(bc2)  # validate
        targets.append(TargetEntry(
            barcode1=bc1, barcode2=bc2, alias=alias, reference=ref,
        ))

    sma = SmaConfig(
        mode=mode,
        targets=targets,
        confidence=confidence,
        truncation=truncation,
        barcode_fasta=sma_raw.get("barcode_fasta"),
    )

    return ConstructConfig(arrangement=arrangement, scoring=scoring, sma=sma)
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_construct.py -v`
Expected: All 15 tests PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/construct.py tests/test_construct.py
git commit -m "feat: add construct TOML parser with validation"
```

---

## Task 2: Barcode FASTA Loader (`bin/barcodes.py` modification)

Add ability to load barcode sequences from an external FASTA file (as specified in `barcode_fasta` field of the construct TOML), while keeping the built-in 96 barcodes as default.

**Files:**
- Modify: `bin/barcodes.py` (add `load_barcodes_from_fasta` function)
- Modify: `tests/test_barcodes.py` (add tests for FASTA loading)

**Step 1: Write the failing tests**

Add to `tests/test_barcodes.py`:

```python
# Add at the top:
from barcodes import load_barcodes_from_fasta

# Add new test class:
class TestLoadBarcodesFromFasta:

    def test_loads_sequences(self, tmp_path):
        fasta = tmp_path / "barcodes.fasta"
        fasta.write_text(">NB01\nACGTACGTACGTACGTACGTACGT\n>NB02\nTGCATGCATGCATGCATGCATGCA\n")
        result = load_barcodes_from_fasta(fasta)
        assert len(result) == 2
        assert result["nb01"] == "ACGTACGTACGTACGTACGTACGT"
        assert result["nb02"] == "TGCATGCATGCATGCATGCATGCA"

    def test_normalizes_ids(self, tmp_path):
        fasta = tmp_path / "barcodes.fasta"
        fasta.write_text(">NB05\nACGTACGTACGTACGTACGTACGT\n")
        result = load_barcodes_from_fasta(fasta)
        assert "nb05" in result

    def test_uppercase_sequences(self, tmp_path):
        fasta = tmp_path / "barcodes.fasta"
        fasta.write_text(">NB01\nacgtacgtacgtacgtacgtacgt\n")
        result = load_barcodes_from_fasta(fasta)
        assert result["nb01"] == "ACGTACGTACGTACGTACGTACGT"

    def test_raises_on_empty(self, tmp_path):
        fasta = tmp_path / "barcodes.fasta"
        fasta.write_text("")
        with pytest.raises(ValueError):
            load_barcodes_from_fasta(fasta)
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_barcodes.py::TestLoadBarcodesFromFasta -v`
Expected: FAIL with `ImportError: cannot import name 'load_barcodes_from_fasta'`

**Step 3: Write the implementation**

Add to `bin/barcodes.py` (after the existing `BARCODES` dict and before `reverse_complement`):

```python
def load_barcodes_from_fasta(path: Path | str) -> dict[str, str]:
    """Load barcode sequences from a FASTA file.

    Header format must match NB\\d+ (case-insensitive).
    Returns dict mapping normalized IDs (nb01, nb02, ...) to uppercase sequences.
    """
    import re
    path = Path(path)
    bc_re = re.compile(r"^>?\s*([Nn][Bb]\d+)")
    barcodes = {}
    current_id = None
    current_seq = []

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id and current_seq:
                barcodes[current_id] = "".join(current_seq).upper()
            m = bc_re.match(line)
            if not m:
                raise ValueError(f"Invalid barcode header: {line!r}")
            num = int(re.search(r"\d+", m.group(1)).group())
            current_id = f"nb{num:02d}"
            current_seq = []
        else:
            current_seq.append(line)

    if current_id and current_seq:
        barcodes[current_id] = "".join(current_seq).upper()

    if not barcodes:
        raise ValueError(f"No barcode sequences found in {path}")

    return barcodes
```

Also add `from pathlib import Path` to the imports if not already present.

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_barcodes.py -v`
Expected: All existing + 4 new tests PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/barcodes.py tests/test_barcodes.py
git commit -m "feat: add barcode FASTA loader for custom barcode sets"
```

---

## Task 3: Truncated Reference Generator (`bin/mkrefs.py`)

CLI tool that reads the construct TOML and user-provided reference FASTAs, then generates truncated reference variants.

**Files:**
- Create: `bin/mkrefs.py`
- Create: `tests/test_mkrefs.py`

**Step 1: Write the failing tests**

```python
# tests/test_mkrefs.py
"""Tests for truncated reference generator."""

from pathlib import Path

import pytest

from mkrefs import generate_truncated_refs, write_manifest


NB05 = "CACAAAGACACCGACAACTTTCTT"
NB10 = "GAGAGGACAAAGGTTTCAACGCTT"
FLANK_F = "AAGGTTAA"
FLANK_R = "CAGCACCT"
REV_FLANK_F = "AGGTGCTG"
REV_FLANK_R = "TTAACCTT"
TARGET = "ATCGATCG" * 10  # 80bp


def rc(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))


@pytest.fixture
def ref_dir(tmp_path):
    d = tmp_path / "references"
    d.mkdir()
    (d / "target_A.fasta").write_text(f">target_A\n{TARGET}\n")
    return d


class TestGenerateTruncatedRefs:

    def test_generates_full_construct(self, ref_dir):
        refs = generate_truncated_refs(
            alias="target_A",
            target_seq=TARGET,
            bc1_seq=NB05,
            bc2_seq=NB10,
            mask1_front=FLANK_F,
            mask1_rear=FLANK_R,
            mask2_front=REV_FLANK_F,
            mask2_rear=REV_FLANK_R,
        )
        full = refs["full"]
        expected = FLANK_F + NB05 + FLANK_R + TARGET + REV_FLANK_F + rc(NB10) + REV_FLANK_R
        assert full == expected

    def test_generates_bc1_target(self, ref_dir):
        refs = generate_truncated_refs(
            alias="target_A",
            target_seq=TARGET,
            bc1_seq=NB05,
            bc2_seq=NB10,
            mask1_front=FLANK_F,
            mask1_rear=FLANK_R,
            mask2_front=REV_FLANK_F,
            mask2_rear=REV_FLANK_R,
        )
        bc1_target = refs["bc1_target"]
        expected = FLANK_F + NB05 + FLANK_R + TARGET
        assert bc1_target == expected

    def test_generates_bc1_only(self, ref_dir):
        refs = generate_truncated_refs(
            alias="target_A",
            target_seq=TARGET,
            bc1_seq=NB05,
            bc2_seq=NB10,
            mask1_front=FLANK_F,
            mask1_rear=FLANK_R,
            mask2_front=REV_FLANK_F,
            mask2_rear=REV_FLANK_R,
        )
        bc1_only = refs["bc1_only"]
        expected = FLANK_F + NB05 + FLANK_R
        assert bc1_only == expected

    def test_returns_three_levels(self, ref_dir):
        refs = generate_truncated_refs(
            alias="target_A",
            target_seq=TARGET,
            bc1_seq=NB05,
            bc2_seq=NB10,
            mask1_front=FLANK_F,
            mask1_rear=FLANK_R,
            mask2_front=REV_FLANK_F,
            mask2_rear=REV_FLANK_R,
        )
        assert set(refs.keys()) == {"full", "bc1_target", "bc1_only"}

    def test_start_only_mode(self, ref_dir):
        refs = generate_truncated_refs(
            alias="target_A",
            target_seq=TARGET,
            bc1_seq=NB05,
            bc2_seq=None,
            mask1_front=FLANK_F,
            mask1_rear=FLANK_R,
            mask2_front=None,
            mask2_rear=None,
        )
        # No bc2, so only bc1_target and bc1_only
        assert "full" not in refs
        assert set(refs.keys()) == {"bc1_target", "bc1_only"}


class TestWriteManifest:

    def test_writes_tsv(self, tmp_path):
        entries = [
            {"alias": "target_A", "level": "full", "length": 200, "path": "target_A_full.fasta"},
            {"alias": "target_A", "level": "bc1_target", "length": 150, "path": "target_A_bc1_target.fasta"},
        ]
        manifest = tmp_path / "manifest.tsv"
        write_manifest(entries, manifest)
        lines = manifest.read_text().strip().split("\n")
        assert len(lines) == 3  # header + 2 entries
        assert lines[0].startswith("alias\t")
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_mkrefs.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'mkrefs'`

**Step 3: Write the implementation**

```python
#!/usr/bin/env python3
"""Generate truncated reference FASTAs from construct TOML.

For each target defined in the construct TOML, generates reference sequences
for each truncation level: full construct, bc1+target, and bc1-only.  These
are used by ingest.py for alignment-based validation of truncated reads.

Usage:
    python bin/mkrefs.py -c construct.toml -o references/
"""

import argparse
import csv
import sys
from pathlib import Path

from barcodes import BARCODES, load_barcodes_from_fasta, reverse_complement
from construct import parse_construct_toml


# ---------------------------------------------------------------------------
# Core logic (no I/O)
# ---------------------------------------------------------------------------

def generate_truncated_refs(
    alias: str,
    target_seq: str,
    bc1_seq: str,
    bc2_seq: str | None,
    mask1_front: str,
    mask1_rear: str,
    mask2_front: str | None,
    mask2_rear: str | None,
) -> dict[str, str]:
    """Generate truncated reference sequences for a target.

    Returns dict mapping truncation level to reference sequence.
    """
    refs = {}

    # bc1_only: flank_F + BC1 + flank_R
    refs["bc1_only"] = mask1_front + bc1_seq + mask1_rear

    # bc1_target: flank_F + BC1 + flank_R + target
    refs["bc1_target"] = mask1_front + bc1_seq + mask1_rear + target_seq

    # full: flank_F + BC1 + flank_R + target + revF + RC(BC2) + revR
    if bc2_seq and mask2_front and mask2_rear:
        refs["full"] = (
            mask1_front + bc1_seq + mask1_rear
            + target_seq
            + mask2_front + reverse_complement(bc2_seq) + mask2_rear
        )

    return refs


def write_manifest(
    entries: list[dict[str, object]], path: Path,
) -> None:
    """Write a TSV manifest of all generated references."""
    fieldnames = ["alias", "level", "length", "path"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for entry in entries:
            writer.writerow(entry)


def _read_first_sequence(fasta_path: Path) -> str:
    """Read the first sequence from a FASTA file."""
    seq_lines = []
    started = False
    for line in fasta_path.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if started:
                break
            started = True
            continue
        if started:
            seq_lines.append(line)
    if not seq_lines:
        raise ValueError(f"No sequence found in {fasta_path}")
    return "".join(seq_lines).upper()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate truncated reference FASTAs from construct TOML"
    )
    parser.add_argument("-c", "--construct", required=True,
        help="Construct TOML file")
    parser.add_argument("-o", "--outdir", required=True,
        help="Output directory for truncated references")
    args = parser.parse_args()

    cfg = parse_construct_toml(args.construct)
    outdir = Path(args.outdir) / "truncated"
    outdir.mkdir(parents=True, exist_ok=True)

    # Load barcode sequences
    if cfg.sma.barcode_fasta:
        bc_seqs = load_barcodes_from_fasta(cfg.sma.barcode_fasta)
    else:
        bc_seqs = BARCODES

    arr = cfg.arrangement
    manifest_entries = []

    for target in cfg.sma.targets:
        bc1_id = target.barcode1.lower()
        bc1_seq = bc_seqs.get(bc1_id)
        if not bc1_seq:
            sys.exit(f"[mkrefs] Error: Barcode {target.barcode1} not found")

        bc2_seq = None
        if target.barcode2:
            bc2_id = target.barcode2.lower()
            bc2_seq = bc_seqs.get(bc2_id)
            if not bc2_seq:
                sys.exit(f"[mkrefs] Error: Barcode {target.barcode2} not found")

        # Read target reference sequence
        ref_path = Path(target.reference)
        if not ref_path.is_absolute():
            ref_path = Path(args.construct).parent / ref_path
        if not ref_path.exists():
            sys.exit(f"[mkrefs] Error: Reference not found: {ref_path}")
        target_seq = _read_first_sequence(ref_path)

        # Generate truncated references
        refs = generate_truncated_refs(
            alias=target.alias,
            target_seq=target_seq,
            bc1_seq=bc1_seq,
            bc2_seq=bc2_seq,
            mask1_front=arr.mask1_front,
            mask1_rear=arr.mask1_rear,
            mask2_front=arr.mask2_front,
            mask2_rear=arr.mask2_rear,
        )

        for level, seq in refs.items():
            fname = f"{target.alias}_{level}.fasta"
            fpath = outdir / fname
            fpath.write_text(f">{target.alias}_{level}\n{seq}\n")
            manifest_entries.append({
                "alias": target.alias,
                "level": level,
                "length": len(seq),
                "path": fname,
            })
            print(f"  {fname} ({len(seq)} bp)")

    write_manifest(manifest_entries, outdir / "manifest.tsv")
    print(f"[mkrefs] Generated {len(manifest_entries)} references in {outdir}")


if __name__ == "__main__":
    main()
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_mkrefs.py -v`
Expected: All 6 tests PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/mkrefs.py tests/test_mkrefs.py
git commit -m "feat: add truncated reference generator (mkrefs.py)"
```

---

## Task 4: Database Schema Extension (`bin/mkdb.py`)

Add `trunc_level` column to the Reads table.

**Files:**
- Modify: `bin/mkdb.py` (add column to CREATE TABLE)
- Modify: `tests/test_mkdb.py` (add test for new column)

**Step 1: Write the failing test**

Add to `tests/test_mkdb.py` in the `TestReadsSchema` class:

```python
    def test_trunc_level_column_exists(self, db_path):
        conn = sqlite3.connect(str(db_path))
        cursor = conn.execute("PRAGMA table_info(Reads)")
        columns = {row[1]: row[2] for row in cursor.fetchall()}
        conn.close()
        assert "trunc_level" in columns
        assert columns["trunc_level"] == "TEXT"
```

**Step 2: Run test to verify it fails**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_mkdb.py::TestReadsSchema::test_trunc_level_column_exists -v`
Expected: FAIL with `AssertionError: assert 'trunc_level' in columns`

**Step 3: Write the implementation**

In `bin/mkdb.py`, find the `CREATE TABLE Reads` statement and add the new column after `bc_end_conf`:

```sql
bc_end_conf REAL,
trunc_level TEXT,
```

The column goes after `bc_end_conf REAL` and before the `FOREIGN KEY` constraints.

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_mkdb.py -v`
Expected: All tests PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/mkdb.py tests/test_mkdb.py
git commit -m "feat: add trunc_level column to Reads schema"
```

---

## Task 5: Truncation Detection in Ingest (`bin/ingest.py`)

Add construct TOML support and 5-level truncation classification to the ingest pipeline.

**Files:**
- Modify: `bin/ingest.py` (add `-c` flag, truncation logic, 21-column insert)
- Modify: `tests/test_ingest.py` (add truncation tests)

**Step 1: Write the failing tests**

Add to `tests/test_ingest.py`:

```python
class TestTruncationClassification:
    """Test truncation detection when using construct TOML."""

    @pytest.fixture
    def construct_env(self, tmp_path):
        """Build environment with construct TOML instead of sample sheet."""
        target_seq = "ATCGATCGATCGATCGATCG" * 5  # 100bp

        # Write construct TOML
        construct = tmp_path / "construct.toml"
        construct.write_text(textwrap.dedent(f"""\
            [arrangement]
            name = "test"
            kit = "SQK-NBD114-96"
            mask1_front = "{FLANK_F}"
            mask1_rear  = "{FLANK_R}"
            mask2_front = "{REV_FLANK_F}"
            mask2_rear  = "{REV_FLANK_R}"
            barcode1_pattern = "NB%02i"
            barcode2_pattern = "NB%02i"
            first_index = 1
            last_index = 96

            [sma]
            mode = "dual_independent"

            [sma.confidence]
            full_length_threshold = 0.75
            start_barcode_min = 0.6

            [sma.truncation]
            auto_generate_refs = false
            min_target_length = 20

            [[sma.targets]]
            barcode1 = "NB05"
            barcode2 = "NB10"
            alias = "target_A"
            reference = "references/target_A.fasta"
        """))

        # Write reference
        ref_dir = tmp_path / "references"
        ref_dir.mkdir()
        (ref_dir / "target_A.fasta").write_text(f">target_A\n{target_seq}\n")

        # Build reads with different truncation levels
        bc_up = BARCODES["nb05"]
        bc_down = BARCODES["nb10"]

        reads = []
        # Full-length read
        full_seq = FLANK_F + bc_up + FLANK_R + target_seq + REV_FLANK_F + rc(bc_down) + REV_FLANK_R
        reads.append(("read_full", full_seq))
        # Truncated: bc1 + target only (no bc2)
        trunc_seq = FLANK_F + bc_up + FLANK_R + target_seq
        reads.append(("read_trunc", trunc_seq))
        # Truncated: bc1 only (very short)
        bc_only_seq = FLANK_F + bc_up + FLANK_R
        reads.append(("read_bc_only", bc_only_seq))

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

        # Summary TSV
        summary = tmp_path / "summary.tsv"
        lines = ["read_id\tend_reason\n"]
        for read_id, _ in reads:
            lines.append(f"{read_id}\tsignal_positive\n")
        summary.write_text("".join(lines))

        # Create DB
        out_dir = tmp_path / "Output"
        r1 = subprocess.run(
            [sys.executable, "bin/mkdb.py", "-e", "FAL99999_20260214_TEST", "-o", str(out_dir)],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
        )
        assert r1.returncode == 0, r1.stderr
        db = out_dir / "SMA_FAL99999_20260214_TEST.db"

        return {
            "db": db, "bam": bam_path, "summary": summary,
            "construct": construct, "ref_dir": ref_dir,
            "out_dir": out_dir, "tmp": tmp_path,
        }

    def test_full_length_classified(self, construct_env):
        env = construct_env
        r = subprocess.run(
            [sys.executable, "bin/ingest.py",
             "-e", "FAL99999_20260214_TEST",
             "-b", str(env["bam"]), "-s", str(env["summary"]),
             "-d", str(env["db"]), "-o", str(env["tmp"] / "tagged.bam"),
             "-c", str(env["construct"]), "-rd", str(env["ref_dir"])],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
            env={**os.environ, "PYTHONPATH": "bin"},
        )
        assert r.returncode == 0, f"{r.stderr}\n{r.stdout}"

        conn = sqlite3.connect(str(env["db"]))
        conn.row_factory = sqlite3.Row
        row = conn.execute(
            "SELECT trunc_level FROM Reads WHERE read_id = 'read_full'"
        ).fetchone()
        conn.close()
        assert row is not None
        assert row["trunc_level"] == "full_length"

    def test_bc1_target_classified(self, construct_env):
        env = construct_env
        subprocess.run(
            [sys.executable, "bin/ingest.py",
             "-e", "FAL99999_20260214_TEST",
             "-b", str(env["bam"]), "-s", str(env["summary"]),
             "-d", str(env["db"]), "-o", str(env["tmp"] / "tagged.bam"),
             "-c", str(env["construct"]), "-rd", str(env["ref_dir"])],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
            env={**os.environ, "PYTHONPATH": "bin"},
        )
        conn = sqlite3.connect(str(env["db"]))
        conn.row_factory = sqlite3.Row
        row = conn.execute(
            "SELECT trunc_level FROM Reads WHERE read_id = 'read_trunc'"
        ).fetchone()
        conn.close()
        assert row is not None
        assert row["trunc_level"] == "bc1_target"

    def test_bc1_only_classified(self, construct_env):
        env = construct_env
        subprocess.run(
            [sys.executable, "bin/ingest.py",
             "-e", "FAL99999_20260214_TEST",
             "-b", str(env["bam"]), "-s", str(env["summary"]),
             "-d", str(env["db"]), "-o", str(env["tmp"] / "tagged.bam"),
             "-c", str(env["construct"]), "-rd", str(env["ref_dir"])],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
            env={**os.environ, "PYTHONPATH": "bin"},
        )
        conn = sqlite3.connect(str(env["db"]))
        conn.row_factory = sqlite3.Row
        row = conn.execute(
            "SELECT trunc_level FROM Reads WHERE read_id = 'read_bc_only'"
        ).fetchone()
        conn.close()
        assert row is not None
        assert row["trunc_level"] == "bc1_only"

    def test_backward_compat_no_construct(self, construct_env):
        """Without -c flag, trunc_level should be NULL (uses -ss mode)."""
        env = construct_env
        # Write a sample sheet for backward-compat mode
        ss = env["tmp"] / "sample_sheet.csv"
        ss.write_text(
            "flow_cell_id,kit,sample_id,experiment_id,barcode,alias\n"
            "FAL99999,SQK-NBD114-96,TEST,exp,barcode05--barcode10,target_A\n"
        )
        subprocess.run(
            [sys.executable, "bin/ingest.py",
             "-e", "FAL99999_20260214_TEST",
             "-b", str(env["bam"]), "-s", str(env["summary"]),
             "-d", str(env["db"]), "-o", str(env["tmp"] / "tagged.bam"),
             "-ss", str(ss), "-rd", str(env["ref_dir"])],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
            env={**os.environ, "PYTHONPATH": "bin"},
        )
        conn = sqlite3.connect(str(env["db"]))
        conn.row_factory = sqlite3.Row
        row = conn.execute(
            "SELECT trunc_level FROM Reads WHERE read_id = 'read_full'"
        ).fetchone()
        conn.close()
        assert row is not None
        assert row["trunc_level"] is None
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_ingest.py::TestTruncationClassification::test_full_length_classified -v`
Expected: FAIL (no `-c` flag, no truncation logic)

**Step 3: Write the implementation**

Modify `bin/ingest.py` with these changes:

1. **Add imports** at the top:
```python
from construct import parse_construct_toml
```

2. **Add `-c` CLI argument** in the argparse section:
```python
parser.add_argument("-c", "--construct",
    help="Construct TOML file (enables truncation detection)")
```

3. **Add truncation classification function** (new function before main loop):
```python
def classify_truncation(
    read_seq: str,
    bc_start_conf: float,
    bc_end_conf: float,
    bc_start_ed: int,
    bc_end_ed: int,
    flank_front: str,
    flank_rear: str,
    flank_rev_front: str | None,
    segment_len: int,
    confidence: object,
    min_target_length: int,
) -> str:
    """Classify truncation level using structural detection + confidence gating.

    Returns one of: 'full_length', 'bc1_target_bc2', 'bc1_target',
    'bc1_only', 'adapter_only'.
    """
    from report_analysis import find_flank_position

    read_len = len(read_seq)

    # Gate 1: Is start barcode confident?
    if bc_start_conf < confidence.start_barcode_min:
        return "adapter_only"

    # Structural: find front flank (mask1_rear, after barcode)
    front_flank = find_flank_position(
        read_seq, flank_rear, 0, min(segment_len + len(flank_rear), read_len)
    )
    front_flank_ok = (
        front_flank is not None
        and front_flank["edit_distance"] / max(len(flank_rear), 1)
        <= confidence.flank_max_error_rate
    )

    # How long is the region after front flank?
    if front_flank_ok:
        region_after_flank = read_len - front_flank["end"]
    else:
        # Estimate: after barcode window
        region_after_flank = read_len - segment_len

    # Gate 2: Is read long enough for target?
    if region_after_flank < min_target_length:
        return "bc1_only"

    # Structural: find rear flank (mask2_front) near read end
    rear_flank_ok = False
    if flank_rev_front:
        rear_flank = find_flank_position(
            read_seq, flank_rev_front,
            max(0, read_len - segment_len - len(flank_rev_front)),
            read_len,
        )
        rear_flank_ok = (
            rear_flank is not None
            and rear_flank["edit_distance"] / max(len(flank_rev_front), 1)
            <= confidence.flank_max_error_rate
        )

    # Gate 3: Is end barcode confident?
    if bc_end_conf >= confidence.full_length_threshold and rear_flank_ok:
        return "full_length"

    if bc_end_conf >= 0.3 and bc_end_conf < confidence.full_length_threshold:
        return "bc1_target_bc2"

    return "bc1_target"
```

4. **Add construct mode** to the main logic. After the existing classification mode block (where `-ss` and `-rd` are checked), add a new branch:

```python
# If construct TOML provided, use it for configuration
construct_cfg = None
if args.construct:
    construct_cfg = parse_construct_toml(args.construct)
    # Use TOML targets instead of sample sheet
    if not barcode_pair_to_alias:
        barcode_pair_to_alias = construct_cfg.barcode_pair_to_alias()
    # Use TOML flanks
    if not flank_front:
        flank_front = construct_cfg.flank_front
    if not flank_rear:
        flank_rear = construct_cfg.flank_rear
```

5. **In the per-read classification loop**, after barcode classification and target assignment, add truncation detection:

```python
# Truncation classification (only with construct TOML)
trunc_level = None
if construct_cfg:
    trunc_level = classify_truncation(
        read_seq=seq,
        bc_start_conf=start_result["confidence"],
        bc_end_conf=end_result["confidence"],
        bc_start_ed=start_result["edit_distance"],
        bc_end_ed=end_result["edit_distance"],
        flank_front=construct_cfg.arrangement.mask1_front,
        flank_rear=construct_cfg.arrangement.mask1_rear,
        flank_rev_front=construct_cfg.arrangement.mask2_front,
        segment_len=SEGMENT_LEN,
        confidence=construct_cfg.sma.confidence,
        min_target_length=construct_cfg.sma.truncation.min_target_length,
    )
```

6. **Update insert_read** to accept 21 columns (add `trunc_level` at position 21):

Change the INSERT statement to include `trunc_level` as the 21st column. In both classification mode and single-target mode, pass `trunc_level` (or `None`) as the 21st element of the data tuple.

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_ingest.py -v`
Expected: All existing + 4 new tests PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/ingest.py tests/test_ingest.py
git commit -m "feat: add construct TOML support and truncation detection to ingest"
```

---

## Task 6: Report Truncation Statistics

Update the report pipeline to display truncation level breakdowns.

**Files:**
- Modify: `bin/report.py` (add `-c` flag, pass truncation data)
- Modify: `bin/report_analysis.py` (add truncation level grouping)
- Modify: `bin/report_template.py` (add truncation breakdown to summary tab)
- Modify: `tests/test_report.py` (add truncation tests)

**Step 1: Write the failing tests**

Add to `tests/test_report.py`:

```python
class TestTruncationAnalysis:

    def test_truncation_counts(self):
        pair_to_alias = {("nb05", "nb10"): "target_A"}
        full_seq = make_full_length_read()
        trunc_seq = make_truncated_read()

        reads = [
            {**make_db_row("r1", full_seq), "trunc_level": "full_length"},
            {**make_db_row("r2", full_seq), "trunc_level": "full_length"},
            {**make_db_row("r3", trunc_seq, tgt_id="nb05_target_only",
                           bc_end="nb09", bc_end_ed=10, bc_end_conf=0.3,
                           ed=None, q_ld=None), "trunc_level": "bc1_target"},
            {**make_db_row("r4", trunc_seq[:40], tgt_id="nb05_bc_only",
                           bc_end="nb09", bc_end_ed=10, bc_end_conf=0.1,
                           ed=None, q_ld=None), "trunc_level": "bc1_only"},
        ]
        result = analyze_classification(reads, pair_to_alias)
        trunc = result.get("truncation_counts")
        assert trunc is not None
        assert trunc["full_length"] == 2
        assert trunc["bc1_target"] == 1
        assert trunc["bc1_only"] == 1

    def test_truncation_absent_when_no_levels(self):
        """When no reads have trunc_level, truncation_counts should be empty."""
        pair_to_alias = {("nb05", "nb10"): "target_A"}
        full_seq = make_full_length_read()
        reads = [make_db_row("r1", full_seq)]
        result = analyze_classification(reads, pair_to_alias)
        trunc = result.get("truncation_counts", {})
        assert trunc == {}
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_report.py::TestTruncationAnalysis -v`
Expected: FAIL with `AssertionError`

**Step 3: Write the implementation**

In `bin/report_analysis.py`, in the `analyze_classification` function, after the existing summary computation, add:

```python
    # Truncation level counts
    truncation_counts = {}
    for r in db_reads:
        level = r.get("trunc_level")
        if level:
            truncation_counts[level] = truncation_counts.get(level, 0) + 1
```

Add `"truncation_counts": truncation_counts` to the returned dict.

In `bin/report.py`, update the SQL query to include `trunc_level`:

```python
"SELECT read_id, readseq, readlen, tgt_id, ed, q_bc, q_ld, ER, "
"bc_start_id, bc_start_ed, bc_start_conf, "
"bc_end_id, bc_end_ed, bc_end_conf, trunc_level FROM Reads"
```

Add `-c` / `--construct` argument and use it to pass flank info:

```python
parser.add_argument("-c", "--construct",
    help="Construct TOML file (for flank and truncation info)")
```

In `bin/report_template.py`, add a truncation breakdown section to the Summary tab. After the per-target stats table, add a truncation level table if `truncation_counts` is non-empty:

```python
    trunc = analysis.get("truncation_counts", {})
    if trunc:
        html += '<h3>Truncation Levels</h3>'
        html += '<table class="data-table"><thead><tr>'
        html += '<th>Level</th><th>Count</th><th>%</th></tr></thead><tbody>'
        total = sum(trunc.values())
        level_order = ["full_length", "bc1_target_bc2", "bc1_target", "bc1_only", "adapter_only"]
        for level in level_order:
            count = trunc.get(level, 0)
            if count > 0:
                pct = 100.0 * count / total if total else 0
                html += f'<tr><td>{level}</td><td>{count}</td><td>{pct:.1f}%</td></tr>'
        html += '</tbody></table>'
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_report.py -v`
Expected: All tests PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/report.py bin/report_analysis.py bin/report_template.py tests/test_report.py
git commit -m "feat: add truncation level statistics to report"
```

---

## Task 7: Interactive Construct Wizard (`bin/construct_wizard.py`)

Interactive CLI that guides users through construct TOML creation with step-by-step explanations and diagrams.

**Files:**
- Create: `bin/construct_wizard.py`
- Create: `tests/test_construct_wizard.py`

**Step 1: Write the failing tests**

```python
# tests/test_construct_wizard.py
"""Tests for the construct wizard (non-interactive parts)."""

import textwrap
from pathlib import Path

import pytest

from construct_wizard import (
    build_construct_toml,
    discover_references,
    format_construct_diagram,
    KNOWN_KITS,
)


@pytest.fixture
def ref_dir(tmp_path):
    d = tmp_path / "references"
    d.mkdir()
    (d / "target_A.fasta").write_text(">target_A\nATCGATCG\n")
    (d / "target_B.fasta").write_text(">target_B\nGCTAGCTA\n")
    (d / "README.txt").write_text("Not a FASTA\n")
    return d


class TestDiscoverReferences:

    def test_finds_fasta_files(self, ref_dir):
        refs = discover_references(ref_dir)
        assert len(refs) == 2
        assert any(r["name"] == "target_A" for r in refs)
        assert any(r["name"] == "target_B" for r in refs)

    def test_reports_lengths(self, ref_dir):
        refs = discover_references(ref_dir)
        for r in refs:
            assert "length" in r
            assert r["length"] == 8

    def test_ignores_non_fasta(self, ref_dir):
        refs = discover_references(ref_dir)
        names = [r["name"] for r in refs]
        assert "README" not in names


class TestFormatConstructDiagram:

    def test_dual_mode(self):
        diagram = format_construct_diagram(
            mode="dual_independent",
            mask1_front="AAGG",
            mask1_rear="CCTT",
            mask2_front="AATT",
            mask2_rear="GGCC",
        )
        assert "BC1" in diagram
        assert "RC(BC2)" in diagram
        assert "TARGET" in diagram

    def test_start_only_mode(self):
        diagram = format_construct_diagram(
            mode="start_only",
            mask1_front="AAGG",
            mask1_rear="CCTT",
        )
        assert "BC1" in diagram
        assert "TARGET" in diagram
        assert "BC2" not in diagram


class TestBuildConstructToml:

    def test_produces_valid_toml(self, tmp_path):
        output = tmp_path / "construct.toml"
        build_construct_toml(
            output_path=output,
            name="test",
            kit="SQK-NBD114-96",
            mode="dual_independent",
            mask1_front="AAGGTTAA",
            mask1_rear="CAGCACCT",
            mask2_front="AGGTGCTG",
            mask2_rear="TTAACCTT",
            targets=[
                {"barcode1": "NB05", "barcode2": "NB10",
                 "alias": "target_A", "reference": "refs/target_A.fasta"},
            ],
        )
        assert output.exists()
        text = output.read_text()
        assert "[arrangement]" in text
        assert "[sma]" in text
        assert "[[sma.targets]]" in text

    def test_roundtrips_through_parser(self, tmp_path):
        from construct import parse_construct_toml
        output = tmp_path / "construct.toml"
        build_construct_toml(
            output_path=output,
            name="test",
            kit="SQK-NBD114-96",
            mode="dual_independent",
            mask1_front="AAGGTTAA",
            mask1_rear="CAGCACCT",
            mask2_front="AGGTGCTG",
            mask2_rear="TTAACCTT",
            targets=[
                {"barcode1": "NB05", "barcode2": "NB10",
                 "alias": "target_A", "reference": "refs/target_A.fasta"},
            ],
        )
        cfg = parse_construct_toml(output)
        assert cfg.arrangement.name == "test"
        assert cfg.sma.mode == "dual_independent"
        assert len(cfg.sma.targets) == 1


class TestKnownKits:

    def test_has_common_kits(self):
        assert "SQK-NBD114-96" in KNOWN_KITS
        assert KNOWN_KITS["SQK-NBD114-96"]["count"] == 96
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_construct_wizard.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'construct_wizard'`

**Step 3: Write the implementation**

```python
#!/usr/bin/env python3
"""Interactive construct TOML builder for SMA-seq experiments.

Guides users through a 7-step process to create a construct TOML file
with visual diagrams and input validation at each step.

Usage:
    python bin/construct_wizard.py -o construct.toml [-rd references/]
    python bin/construct_wizard.py -i existing.toml -o updated.toml
"""

import argparse
import sys
from pathlib import Path

import tomli_w

from construct import parse_construct_toml, ValidationError


# ---------------------------------------------------------------------------
# Known ONT barcode kits
# ---------------------------------------------------------------------------

KNOWN_KITS = {
    "SQK-NBD114-96": {"count": 96, "length": 24, "desc": "Native Barcoding Kit 96"},
    "SQK-NBD114-24": {"count": 24, "length": 24, "desc": "Native Barcoding Kit 24"},
    "SQK-RBK114-96": {"count": 96, "length": 24, "desc": "Rapid Barcoding Kit 96"},
    "SQK-RBK114-24": {"count": 24, "length": 24, "desc": "Rapid Barcoding Kit 24"},
}


# ---------------------------------------------------------------------------
# Non-interactive helpers (testable)
# ---------------------------------------------------------------------------

def discover_references(ref_dir: Path) -> list[dict]:
    """Find FASTA files in a directory and report their names and lengths."""
    refs = []
    for ext in ("*.fasta", "*.fa"):
        for p in sorted(ref_dir.glob(ext)):
            seq_lines = []
            for line in p.read_text().splitlines():
                if not line.startswith(">"):
                    seq_lines.append(line.strip())
            seq = "".join(seq_lines)
            refs.append({
                "name": p.stem,
                "path": str(p),
                "length": len(seq),
            })
    return refs


def format_construct_diagram(
    mode: str,
    mask1_front: str,
    mask1_rear: str,
    mask2_front: str | None = None,
    mask2_rear: str | None = None,
) -> str:
    """Format an ASCII construct diagram."""
    if mode == "dual_independent":
        return (
            f"5'-adapter-[{mask1_front}]-[BC1]-[{mask1_rear}]"
            f"-TARGET-"
            f"[{mask2_front}]-[RC(BC2)]-[{mask2_rear}]-adapter-3'"
        )
    else:
        return (
            f"5'-adapter-[{mask1_front}]-[BC1]-[{mask1_rear}]"
            f"-TARGET-3'"
        )


def build_construct_toml(
    output_path: Path,
    name: str,
    kit: str,
    mode: str,
    mask1_front: str,
    mask1_rear: str,
    targets: list[dict],
    mask2_front: str | None = None,
    mask2_rear: str | None = None,
    barcode_fasta: str | None = None,
    confidence: dict | None = None,
    truncation: dict | None = None,
) -> None:
    """Build and write a construct TOML file."""
    doc = {
        "arrangement": {
            "name": name,
            "kit": kit,
            "mask1_front": mask1_front,
            "mask1_rear": mask1_rear,
            "barcode1_pattern": "NB%02i",
            "first_index": 1,
            "last_index": KNOWN_KITS.get(kit, {}).get("count", 96),
        },
        "scoring": {
            "max_barcode_penalty": 11,
            "min_barcode_penalty_dist": 3,
            "front_barcode_window": 100,
            "rear_barcode_window": 100,
        },
        "sma": {
            "mode": mode,
            "targets": targets,
        },
    }

    if mode == "dual_independent" and mask2_front and mask2_rear:
        doc["arrangement"]["mask2_front"] = mask2_front
        doc["arrangement"]["mask2_rear"] = mask2_rear
        doc["arrangement"]["barcode2_pattern"] = "NB%02i"

    if barcode_fasta:
        doc["sma"]["barcode_fasta"] = barcode_fasta

    if confidence:
        doc["sma"]["confidence"] = confidence
    if truncation:
        doc["sma"]["truncation"] = truncation

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "wb") as f:
        tomli_w.dump(doc, f)


# ---------------------------------------------------------------------------
# Interactive wizard
# ---------------------------------------------------------------------------

def _prompt(msg: str, default: str | None = None) -> str:
    """Prompt user for input with optional default."""
    if default:
        val = input(f"  {msg} [{default}]: ").strip()
        return val if val else default
    return input(f"  {msg}: ").strip()


def _prompt_choice(msg: str, options: list[str]) -> int:
    """Prompt user to select from numbered options. Returns 0-based index."""
    print(f"\n  {msg}")
    for i, opt in enumerate(options, 1):
        print(f"    [{i}] {opt}")
    while True:
        try:
            choice = int(input("  > ")) - 1
            if 0 <= choice < len(options):
                return choice
        except (ValueError, EOFError):
            pass
        print(f"  Please enter 1-{len(options)}")


def run_wizard(output_path: Path, ref_dir: Path | None = None,
               input_toml: Path | None = None) -> None:
    """Run the interactive 7-step construct wizard."""
    print()
    print("=" * 62)
    print("         SMA-seq Construct Definition Wizard")
    print("=" * 62)
    print()

    # Load existing config if editing
    existing = None
    if input_toml:
        try:
            existing = parse_construct_toml(input_toml)
            print(f"  Loaded existing config: {input_toml}")
        except ValidationError as e:
            print(f"  Warning: Could not parse {input_toml}: {e}")

    # Step 1: Mode
    print("\nStep 1/7: Library Design Mode")
    print("-" * 30)
    print("\n  Your SMA-seq library has barcodes at one or both ends:")
    print()
    print("  Dual independent:")
    print("  5'-adapter-[BC1]-flank-TARGET-flank-[RC(BC2)]-adapter-3'")
    print("  (BC1 and BC2 can be different barcodes)")
    print()
    print("  Start-only:")
    print("  5'-adapter-[BC1]-flank-TARGET-3'")
    print("  (Only one barcode, target may be truncated at 3' end)")

    mode_idx = _prompt_choice("Select mode:", [
        "Dual independent barcodes (recommended for SMA-seq)",
        "Start-only barcode",
    ])
    mode = "dual_independent" if mode_idx == 0 else "start_only"
    print(f"\n  Mode: {mode}")

    # Step 2: Kit
    print("\nStep 2/7: Barcode Kit")
    print("-" * 30)
    kit_options = []
    kit_keys = list(KNOWN_KITS.keys())
    for k in kit_keys:
        info = KNOWN_KITS[k]
        kit_options.append(f"{k} ({info['desc']}, {info['count']} barcodes)")
    kit_options.append("Custom (provide barcode FASTA)")
    kit_idx = _prompt_choice("Select barcode kit:", kit_options)

    barcode_fasta = None
    if kit_idx < len(kit_keys):
        kit = kit_keys[kit_idx]
        info = KNOWN_KITS[kit]
        print(f"\n  Using {info['count']} barcodes ({info['length']}bp each)")
    else:
        kit = _prompt("Kit name")
        barcode_fasta = _prompt("Path to barcode FASTA")

    # Step 3: Flanks
    print("\nStep 3/7: Flanking Sequences")
    print("-" * 30)
    print("\n  These sequences flank the barcode in your construct:")
    print(f"  5'-...[mask1_front][BARCODE1][mask1_rear]-TARGET-...-3'")

    mask1_front = _prompt("mask1_front (before BC1)").upper()
    mask1_rear = _prompt("mask1_rear (after BC1)").upper()

    mask2_front = None
    mask2_rear = None
    if mode == "dual_independent":
        print(f"\n  3' end: ...-TARGET-[mask2_front][RC(BC2)][mask2_rear]-...-3'")
        mask2_front = _prompt("mask2_front (before RC(BC2))").upper()
        mask2_rear = _prompt("mask2_rear (after RC(BC2))").upper()

    diagram = format_construct_diagram(mode, mask1_front, mask1_rear,
                                       mask2_front, mask2_rear)
    print(f"\n  Construct: {diagram}")

    # Step 4: Targets
    print("\nStep 4/7: Barcode Pair -> Target Mapping")
    print("-" * 30)

    available_refs = []
    if ref_dir:
        available_refs = discover_references(ref_dir)
        if available_refs:
            print("\n  Available references:")
            for i, r in enumerate(available_refs, 1):
                print(f"    [{i}] {r['name']} ({r['length']} bp)")

    targets = []
    while True:
        print(f"\n  Pair {len(targets) + 1}:")
        bc1_num = _prompt("BC1 (upstream barcode number)")
        bc1 = f"NB{int(bc1_num):02d}"

        bc2 = None
        if mode == "dual_independent":
            bc2_num = _prompt("BC2 (downstream barcode number)")
            bc2 = f"NB{int(bc2_num):02d}"

        if available_refs:
            ref_idx = _prompt_choice("Target reference:", [
                f"{r['name']} ({r['length']}bp)" for r in available_refs
            ])
            ref_path = available_refs[ref_idx]["path"]
            default_alias = available_refs[ref_idx]["name"]
        else:
            ref_path = _prompt("Path to reference FASTA")
            default_alias = Path(ref_path).stem

        alias = _prompt("Alias", default_alias)

        entry = {"barcode1": bc1, "alias": alias, "reference": ref_path}
        if bc2:
            entry["barcode2"] = bc2

        targets.append(entry)
        bc_desc = f"{bc1}/{bc2}" if bc2 else bc1
        print(f"\n  Added: {bc_desc} -> {alias}")

        more = _prompt("Add another pair? [Y/n]", "y")
        if more.lower().startswith("n"):
            break

    # Step 5: Confidence
    print("\nStep 5/7: Confidence Thresholds")
    print("-" * 30)
    fl_thresh = float(_prompt("Full-length threshold (bc_end_conf >=)", "0.75"))
    start_min = float(_prompt("Start barcode minimum (bc_start_conf >=)", "0.6"))
    flank_err = float(_prompt("Flank max error rate", "0.5"))
    confidence = {
        "full_length_threshold": fl_thresh,
        "start_barcode_min": start_min,
        "flank_max_error_rate": flank_err,
    }

    # Step 6: Truncation
    print("\nStep 6/7: Truncation Settings")
    print("-" * 30)
    auto_gen = _prompt("Auto-generate truncated references? [Y/n]", "y")
    auto_gen_bool = not auto_gen.lower().startswith("n")
    min_tgt = int(_prompt("Minimum target length (bp)", "20"))
    truncation = {
        "auto_generate_refs": auto_gen_bool,
        "min_target_length": min_tgt,
    }

    # Step 7: Review & Save
    print("\nStep 7/7: Review & Save")
    print("-" * 30)
    name = _prompt("Construct name", f"SMA_{kit.replace('-', '_')}")

    print()
    print("  " + "=" * 50)
    print(f"  Construct: {name}")
    print(f"  Mode:      {mode}")
    print(f"  Kit:       {kit}")
    print(f"  Structure: {diagram}")
    print(f"  Targets:   {len(targets)} pair(s)")
    for t in targets:
        bc_desc = f"{t['barcode1']}/{t.get('barcode2', '-')}"
        print(f"    {bc_desc} -> {t['alias']}")
    print(f"  Truncation: {'auto-generate' if auto_gen_bool else 'manual'}")
    print("  " + "=" * 50)

    confirm = _prompt("Save? [Y/n]", "y")
    if confirm.lower().startswith("n"):
        print("  Cancelled.")
        return

    build_construct_toml(
        output_path=output_path,
        name=name,
        kit=kit,
        mode=mode,
        mask1_front=mask1_front,
        mask1_rear=mask1_rear,
        mask2_front=mask2_front,
        mask2_rear=mask2_rear,
        targets=targets,
        barcode_fasta=barcode_fasta,
        confidence=confidence,
        truncation=truncation,
    )
    print(f"\n  Written to {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Interactive SMA-seq construct TOML builder"
    )
    parser.add_argument("-o", "--output", required=True,
        help="Output TOML file path")
    parser.add_argument("-rd", "--ref-dir",
        help="Directory with reference FASTAs (for auto-discovery)")
    parser.add_argument("-i", "--input",
        help="Existing TOML to edit")
    args = parser.parse_args()

    ref_dir = Path(args.ref_dir) if args.ref_dir else None
    input_toml = Path(args.input) if args.input else None

    run_wizard(
        output_path=Path(args.output),
        ref_dir=ref_dir,
        input_toml=input_toml,
    )


if __name__ == "__main__":
    main()
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_construct_wizard.py -v`
Expected: All 7 tests PASS

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add bin/construct_wizard.py tests/test_construct_wizard.py
git commit -m "feat: add interactive construct wizard for TOML generation"
```

---

## Task 8: Integration Test

End-to-end test: wizard output → mkrefs → ingest → report.

**Files:**
- Modify: `tests/test_construct.py` (add integration test)

**Step 1: Write the integration test**

Add to `tests/test_construct.py`:

```python
class TestIntegrationPipeline:
    """End-to-end: TOML → mkrefs → ingest → report."""

    def test_full_pipeline_with_construct(self, tmp_path):
        import csv
        import subprocess

        target_seq = "ATCGATCGATCGATCGATCG" * 5  # 100bp

        # 1. Write construct TOML
        construct = tmp_path / "construct.toml"
        construct.write_text(MINIMAL_DUAL_TOML.replace(
            'reference = "references/target_A.fasta"',
            f'reference = "{tmp_path / "references" / "target_A.fasta"}"',
        ))

        # 2. Write reference FASTA
        ref_dir = tmp_path / "references"
        ref_dir.mkdir()
        (ref_dir / "target_A.fasta").write_text(f">target_A\n{target_seq}\n")

        # 3. Run mkrefs
        r1 = subprocess.run(
            [sys.executable, "bin/mkrefs.py",
             "-c", str(construct), "-o", str(ref_dir)],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
            env={**os.environ, "PYTHONPATH": "bin"},
        )
        assert r1.returncode == 0, r1.stderr
        assert (ref_dir / "truncated" / "manifest.tsv").exists()

        # 4. Build BAM with test reads
        from barcodes import BARCODES, reverse_complement
        bc_up = BARCODES["nb05"]
        bc_down = BARCODES["nb10"]
        flank_f = "AAGGTTAA"
        flank_r = "CAGCACCT"
        rev_f = "AGGTGCTG"
        rev_r = "TTAACCTT"

        full_seq = flank_f + bc_up + flank_r + target_seq + rev_f + reverse_complement(bc_down) + rev_r
        trunc_seq = flank_f + bc_up + flank_r + target_seq

        bam_path = tmp_path / "FAL99999_20260214_TEST_sup_v5.2.0_trim0_0.bam"
        header = {"HD": {"VN": "1.6", "SO": "unsorted"}}
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            for rid, seq in [("r_full", full_seq), ("r_trunc", trunc_seq)]:
                a = pysam.AlignedSegment()
                a.query_name = rid
                a.query_sequence = seq
                a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
                a.flag = 4
                bam.write(a)

        summary = tmp_path / "summary.tsv"
        summary.write_text("read_id\tend_reason\nr_full\tsignal_positive\nr_trunc\tsignal_positive\n")

        # 5. mkdb
        out_dir = tmp_path / "Output"
        r2 = subprocess.run(
            [sys.executable, "bin/mkdb.py",
             "-e", "FAL99999_20260214_TEST", "-o", str(out_dir)],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
        )
        assert r2.returncode == 0, r2.stderr
        db = out_dir / "SMA_FAL99999_20260214_TEST.db"

        # 6. ingest with construct TOML
        r3 = subprocess.run(
            [sys.executable, "bin/ingest.py",
             "-e", "FAL99999_20260214_TEST",
             "-b", str(bam_path), "-s", str(summary),
             "-d", str(db), "-o", str(tmp_path / "tagged.bam"),
             "-c", str(construct), "-rd", str(ref_dir)],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
            env={**os.environ, "PYTHONPATH": "bin"},
        )
        assert r3.returncode == 0, f"{r3.stderr}\n{r3.stdout}"

        # 7. Verify truncation levels
        import sqlite3
        conn = sqlite3.connect(str(db))
        conn.row_factory = sqlite3.Row
        rows = {
            r["read_id"]: r["trunc_level"]
            for r in conn.execute("SELECT read_id, trunc_level FROM Reads")
        }
        conn.close()
        assert rows["r_full"] == "full_length"
        assert rows["r_trunc"] in ("bc1_target", "bc1_target_bc2")
```

Add needed imports to the top of `tests/test_construct.py`:

```python
import os
import sqlite3
import subprocess
import sys
import pysam
```

**Step 2: Run test to verify it fails**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_construct.py::TestIntegrationPipeline -v`
Expected: FAIL (until all prior tasks are complete)

**Step 3: Run test after all prior tasks**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_construct.py -v`
Expected: All tests PASS

**Step 4: Run full test suite**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/ -v`
Expected: All tests PASS (existing 69 + new ~30)

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq
git add tests/test_construct.py
git commit -m "test: add end-to-end integration test for construct pipeline"
```

---

## Verification

After all tasks are complete:

```bash
# 1. Run full test suite
cd /tmp/ont-sma-seq && python -m pytest tests/ -v

# 2. Test with demo data
python bin/construct_wizard.py -o /tmp/sma-demo/construct.toml -rd /tmp/sma-demo/references/

# 3. Generate truncated references
python bin/mkrefs.py -c /tmp/sma-demo/construct.toml -o /tmp/sma-demo/references/

# 4. Re-ingest with construct TOML
python bin/ingest.py \
  -e FBD66244_20251230_IF \
  -b /path/to/reads.bam -s /path/to/summary.tsv \
  -d /tmp/sma-demo/Output/SMA_FBD66244_20251230_IF.db \
  -o /tmp/sma-demo/tagged.bam \
  -c /tmp/sma-demo/construct.toml \
  -rd /tmp/sma-demo/references/

# 5. Generate report with truncation data
python bin/report.py \
  -d /tmp/sma-demo/Output/SMA_FBD66244_20251230_IF.db \
  -c /tmp/sma-demo/construct.toml \
  -o /tmp/sma-demo/report_v4.html
```
