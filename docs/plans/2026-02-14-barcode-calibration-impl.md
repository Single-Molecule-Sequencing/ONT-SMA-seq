# Barcode Calibration & Threshold Visualization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a CLI pipeline + browser visualization tool for data-driven barcode threshold calibration from MinKNOW output directories.

**Architecture:** CLI (`sma_calibrate.py`) auto-discovers MinKNOW runs, validates basecalling, merges BAMs, and ingests with SMA-seq demultiplexing. Browser app (`calibrate_viz`) provides interactive KDE distributions, draggable thresholds, confusion matrices, and cross-experiment comparison. Both share an extended SQLite schema.

**Tech Stack:** Python 3.11+, FastAPI, HTMX, D3.js, edlib, pysam, pandas, scipy (KDE + Savitzky-Golay), SQLite

**Design doc:** `docs/plans/2026-02-14-barcode-calibration-design.md`

---

## Phase 1: CLI Pipeline — Discovery & Merge

### Task 1: Scaffolding & Discovery Module

**Files:**
- Create: `bin/calibrate/__init__.py`
- Create: `bin/calibrate/discover.py`
- Create: `tests/test_discover.py`

**Context:** The discovery module walks a root directory looking for MinKNOW output directories. It identifies them by the presence of `final_summary_*.txt` files and parses run metadata from those files plus `sequencing_summary_*.txt` and `sample_sheet_*.csv`.

MinKNOW output directory naming: `{start_time}_{device_id}_{flow_cell_id}_{short_protocol_run_id}/`
Contains: `final_summary_*.txt`, `sequencing_summary_*.txt`, `sample_sheet_*.csv`, `bam_pass/`, `pod5_pass/`

**Step 1: Write the failing tests**

```python
"""Tests for MinKNOW output directory discovery."""
from __future__ import annotations

from pathlib import Path

import pytest

from calibrate.discover import (
    discover_runs,
    parse_final_summary,
    RunInfo,
)


@pytest.fixture
def minknow_tree(tmp_path: Path) -> Path:
    """Create a minimal MinKNOW output directory tree with two runs."""
    root = tmp_path / "experiment"
    root.mkdir()

    # Run 1
    run1 = root / "no_sample_id" / "20251228_2219_MD-100098_FBD69411_34fa833d"
    run1.mkdir(parents=True)
    (run1 / "final_summary_FBD69411_34fa833d_abcd1234.txt").write_text(
        "protocol_run_id=run-uuid-1\n"
        "acquisition_run_id=acq-uuid-1\n"
        "started=2025-12-28T22:19:00Z\n"
        "acquisition_stopped=2025-12-29T03:20:00Z\n"
        "flow_cell_id=FBD69411\n"
        "device_id=MD-100098\n"
        "sample_id=no_sample_id\n"
        "experiment_id=my_experiment\n"
        "protocol=sequencing/sequencing_MIN114_DNA_e8_2_400K:FLO-MIN114:SQK-NBD114-96\n"
    )
    seq_summary = run1 / "sequencing_summary_FBD69411_34fa833d_abcd1234.txt"
    seq_summary.write_text(
        "read_id\tduration\tend_reason\tsequence_length_template\tmean_qscore_template\tbarcode_arrangement\n"
        "read-001\t1.5\tsignal_positive\t500\t12.3\tbarcode05\n"
        "read-002\t0.8\tdata_service_unblock_mux_change\t200\t9.1\tunclassified\n"
    )
    bam_dir = run1 / "bam_pass"
    bam_dir.mkdir()

    # Run 2 — same flow cell, later start (should merge with Run 1)
    run2 = root / "no_sample_id" / "20251229_1055_MD-100098_FBD69411_5b5c57a9"
    run2.mkdir(parents=True)
    (run2 / "final_summary_FBD69411_5b5c57a9_efgh5678.txt").write_text(
        "protocol_run_id=run-uuid-2\n"
        "acquisition_run_id=acq-uuid-2\n"
        "started=2025-12-29T10:55:00Z\n"
        "acquisition_stopped=2025-12-29T14:00:00Z\n"
        "flow_cell_id=FBD69411\n"
        "device_id=MD-100098\n"
        "sample_id=no_sample_id\n"
        "experiment_id=my_experiment\n"
        "protocol=sequencing/sequencing_MIN114_DNA_e8_2_400K:FLO-MIN114:SQK-NBD114-96\n"
    )
    (run2 / "bam_pass").mkdir()

    return root


class TestParsingFinalSummary:
    def test_parses_key_fields(self, minknow_tree: Path):
        fs_path = next(minknow_tree.rglob("final_summary_FBD69411_34fa833d_*.txt"))
        info = parse_final_summary(fs_path)
        assert info.flow_cell_id == "FBD69411"
        assert info.device_id == "MD-100098"
        assert info.sample_id == "no_sample_id"
        assert info.experiment_id == "my_experiment"

    def test_extracts_start_time(self, minknow_tree: Path):
        fs_path = next(minknow_tree.rglob("final_summary_FBD69411_34fa833d_*.txt"))
        info = parse_final_summary(fs_path)
        assert info.started == "2025-12-28T22:19:00Z"


class TestDiscoverRuns:
    def test_finds_two_runs(self, minknow_tree: Path):
        runs = discover_runs(minknow_tree)
        assert len(runs) == 2

    def test_groups_by_flow_cell(self, minknow_tree: Path):
        runs = discover_runs(minknow_tree)
        flow_cells = {r.flow_cell_id for r in runs}
        assert flow_cells == {"FBD69411"}

    def test_returns_run_dirs(self, minknow_tree: Path):
        runs = discover_runs(minknow_tree)
        for r in runs:
            assert r.run_dir.is_dir()
            assert (r.run_dir / "bam_pass").is_dir()
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_discover.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'calibrate'`

**Step 3: Implement discovery module**

```python
"""Auto-discover MinKNOW output directories and parse run metadata."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class RunInfo:
    """Parsed metadata from a single MinKNOW sequencing run."""

    run_dir: Path
    flow_cell_id: str
    device_id: str
    sample_id: str
    experiment_id: str
    started: str
    protocol_run_id: str
    final_summary_path: Path
    sequencing_summary_path: Path | None
    sample_sheet_path: Path | None


def parse_final_summary(path: Path) -> RunInfo:
    """Parse a final_summary_*.txt key=value file into RunInfo.

    Parameters
    ----------
    path : Path
        Path to a final_summary_*.txt file.

    Returns
    -------
    RunInfo
        Parsed run metadata.
    """
    kv: dict[str, str] = {}
    for line in path.read_text().strip().splitlines():
        if "=" in line:
            key, _, value = line.partition("=")
            kv[key.strip()] = value.strip()

    run_dir = path.parent

    # Find companion files in same directory
    seq_summaries = list(run_dir.glob("sequencing_summary_*.txt"))
    sample_sheets = list(run_dir.glob("sample_sheet_*.csv"))

    return RunInfo(
        run_dir=run_dir,
        flow_cell_id=kv.get("flow_cell_id", ""),
        device_id=kv.get("device_id", ""),
        sample_id=kv.get("sample_id", ""),
        experiment_id=kv.get("experiment_id", ""),
        started=kv.get("started", ""),
        protocol_run_id=kv.get("protocol_run_id", ""),
        final_summary_path=path,
        sequencing_summary_path=seq_summaries[0] if seq_summaries else None,
        sample_sheet_path=sample_sheets[0] if sample_sheets else None,
    )


def discover_runs(root: Path) -> list[RunInfo]:
    """Walk directory tree and find all MinKNOW output directories.

    Identifies MinKNOW output dirs by presence of final_summary_*.txt.

    Parameters
    ----------
    root : Path
        Root directory to search recursively.

    Returns
    -------
    list[RunInfo]
        All discovered runs sorted by start time.
    """
    runs: list[RunInfo] = []
    for fs_path in sorted(root.rglob("final_summary_*.txt")):
        runs.append(parse_final_summary(fs_path))
    runs.sort(key=lambda r: r.started)
    return runs
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_discover.py -v`
Expected: PASS (5 tests)

**Step 5: Commit**

```bash
git add bin/calibrate/__init__.py bin/calibrate/discover.py tests/test_discover.py
git commit -m "feat: add MinKNOW output directory discovery module"
```

---

### Task 2: Run Grouping & Merge Logic

**Files:**
- Create: `bin/calibrate/merge.py`
- Create: `tests/test_merge.py`

**Context:** Group discovered runs by `(flow_cell_id, sample_id)`. Runs within a group that started within 24 hours of each other are merged. The merge module collects all BAM files from per-barcode subdirectories across grouped runs and concatenates them with `pysam.merge`. It also validates that all BAMs share the same basecalling model (from `@PG` header lines).

The existing `parse_bam_filename()` in `bin/ingest.py:128-147` extracts model tier and version from filenames, but we also need to check BAM `@PG` headers for basecall model info and trim/demux flags. MinKNOW BAM files have `@PG` lines like:
```
@PG  ID:basecaller  PN:dorado  VN:0.8.0  CL:dorado basecaller ...model_name... --no-trim
```

**Step 1: Write the failing tests**

```python
"""Tests for run grouping and BAM merge logic."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pysam
import pytest

from calibrate.discover import RunInfo
from calibrate.merge import (
    group_runs,
    find_bam_files,
    validate_basecalling,
    merge_bams,
    MergeGroup,
    ValidationError,
)


def _make_run_info(
    tmp_path: Path,
    name: str,
    flow_cell_id: str = "FBD69411",
    sample_id: str = "sample_A",
    started: str = "2025-12-28T22:19:00Z",
) -> RunInfo:
    """Helper to create a RunInfo with a real directory."""
    d = tmp_path / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "bam_pass").mkdir(exist_ok=True)
    fs = d / "final_summary.txt"
    fs.write_text(f"flow_cell_id={flow_cell_id}\nstarted={started}\n")
    return RunInfo(
        run_dir=d,
        flow_cell_id=flow_cell_id,
        device_id="MD-100098",
        sample_id=sample_id,
        experiment_id="exp1",
        started=started,
        protocol_run_id="proto-1",
        final_summary_path=fs,
        sequencing_summary_path=None,
        sample_sheet_path=None,
    )


class TestGroupRuns:
    def test_groups_same_flow_cell_and_sample(self, tmp_path: Path):
        r1 = _make_run_info(tmp_path, "run1", started="2025-12-28T22:19:00Z")
        r2 = _make_run_info(tmp_path, "run2", started="2025-12-29T10:55:00Z")
        groups = group_runs([r1, r2])
        assert len(groups) == 1
        assert len(groups[0].runs) == 2

    def test_separates_different_flow_cells(self, tmp_path: Path):
        r1 = _make_run_info(tmp_path, "run1", flow_cell_id="FBD69411")
        r2 = _make_run_info(tmp_path, "run2", flow_cell_id="FBD99999")
        groups = group_runs([r1, r2])
        assert len(groups) == 2

    def test_separates_runs_more_than_24h_apart(self, tmp_path: Path):
        r1 = _make_run_info(tmp_path, "run1", started="2025-12-28T00:00:00Z")
        r2 = _make_run_info(tmp_path, "run2", started="2025-12-30T12:00:00Z")
        groups = group_runs([r1, r2])
        assert len(groups) == 2

    def test_separates_different_samples(self, tmp_path: Path):
        r1 = _make_run_info(tmp_path, "run1", sample_id="sample_A")
        r2 = _make_run_info(tmp_path, "run2", sample_id="sample_B")
        groups = group_runs([r1, r2])
        assert len(groups) == 2


class TestFindBamFiles:
    def test_finds_flat_bams(self, tmp_path: Path):
        bam_dir = tmp_path / "bam_pass"
        bam_dir.mkdir()
        (bam_dir / "reads_0.bam").touch()
        (bam_dir / "reads_1.bam").touch()
        bams = find_bam_files(tmp_path)
        assert len(bams) == 2

    def test_finds_per_barcode_bams(self, tmp_path: Path):
        for alias in ["sample_A", "sample_B", "unclassified"]:
            d = tmp_path / "bam_pass" / alias
            d.mkdir(parents=True)
            (d / f"{alias}_0.bam").touch()
        bams = find_bam_files(tmp_path)
        assert len(bams) == 3

    def test_returns_empty_for_no_bams(self, tmp_path: Path):
        (tmp_path / "bam_pass").mkdir()
        bams = find_bam_files(tmp_path)
        assert len(bams) == 0
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_merge.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'calibrate.merge'`

**Step 3: Implement merge module**

```python
"""Group MinKNOW runs and merge BAM files for SMA-seq ingestion."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pysam

from calibrate.discover import RunInfo


class ValidationError(Exception):
    """Raised when BAM validation fails."""


@dataclass
class MergeGroup:
    """A group of runs to merge into a single BAM."""

    flow_cell_id: str
    sample_id: str
    runs: list[RunInfo] = field(default_factory=list)


def _parse_iso_time(s: str) -> datetime:
    """Parse an ISO 8601 timestamp."""
    s = s.replace("Z", "+00:00")
    return datetime.fromisoformat(s)


def group_runs(
    runs: list[RunInfo],
    max_gap_hours: float = 24.0,
) -> list[MergeGroup]:
    """Group runs by (flow_cell_id, sample_id) and time proximity.

    Runs on the same flow cell with the same sample that started within
    ``max_gap_hours`` of each other are grouped together.
    """
    # Sort by start time
    sorted_runs = sorted(runs, key=lambda r: r.started)

    groups: list[MergeGroup] = []
    for run in sorted_runs:
        placed = False
        for group in groups:
            if (
                group.flow_cell_id == run.flow_cell_id
                and group.sample_id == run.sample_id
            ):
                # Check time gap against last run in group
                last_start = _parse_iso_time(group.runs[-1].started)
                this_start = _parse_iso_time(run.started)
                if this_start - last_start <= timedelta(hours=max_gap_hours):
                    group.runs.append(run)
                    placed = True
                    break
        if not placed:
            groups.append(
                MergeGroup(
                    flow_cell_id=run.flow_cell_id,
                    sample_id=run.sample_id,
                    runs=[run],
                )
            )
    return groups


def find_bam_files(run_dir: Path) -> list[Path]:
    """Find all BAM files in a run directory.

    Searches ``bam_pass/`` for BAM files, both flat and in per-alias
    subdirectories.
    """
    bam_pass = run_dir / "bam_pass"
    if not bam_pass.is_dir():
        return []
    return sorted(bam_pass.rglob("*.bam"))


def validate_basecalling(bam_paths: list[Path]) -> str:
    """Validate that all BAM files share the same basecalling model.

    Extracts the basecall model from the @PG header of each BAM file.
    Raises ValidationError if models differ or trimming/demuxing detected.

    Returns
    -------
    str
        The basecalling model string common to all files.
    """
    if not bam_paths:
        raise ValidationError("No BAM files to validate")

    models: set[str] = set()
    for bam_path in bam_paths:
        save = pysam.set_verbosity(0)
        try:
            with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bf:
                header = bf.header.to_dict()
        finally:
            pysam.set_verbosity(save)

        pg_lines = header.get("PG", [])
        model = "unknown"
        for pg in pg_lines:
            cl = pg.get("CL", "")
            if "basecaller" in cl.lower() or "dorado" in cl.lower():
                model = cl
                if "--trim" in cl and "--no-trim" not in cl:
                    raise ValidationError(
                        f"BAM {bam_path.name} was basecalled with trimming enabled"
                    )
                break
        models.add(model)

    if len(models) > 1:
        raise ValidationError(
            f"Inconsistent basecalling models across BAMs: {models}"
        )
    return models.pop()


def merge_bams(
    group: MergeGroup,
    output_path: Path,
) -> dict:
    """Merge all BAM files from a run group into a single BAM.

    Returns
    -------
    dict
        Merge provenance metadata.
    """
    all_bams: list[Path] = []
    for run in group.runs:
        all_bams.extend(find_bam_files(run.run_dir))

    if not all_bams:
        raise ValidationError(
            f"No BAM files found in group {group.flow_cell_id}/{group.sample_id}"
        )

    bam_strs = [str(b) for b in all_bams]
    pysam.merge("-f", "-o", str(output_path), *bam_strs)

    return {
        "flow_cell_id": group.flow_cell_id,
        "sample_id": group.sample_id,
        "source_bam_count": len(all_bams),
        "source_bam_paths": [str(p) for p in all_bams],
        "run_ids": [r.protocol_run_id for r in group.runs],
    }
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_merge.py -v`
Expected: PASS (7 tests)

**Step 5: Commit**

```bash
git add bin/calibrate/merge.py tests/test_merge.py
git commit -m "feat: add run grouping and BAM merge with basecalling validation"
```

---

### Task 3: Schema Extensions & Enhanced Ingestion

**Files:**
- Modify: `bin/mkdb.py` (add RunMetadata + ReadRun tables, add signal_duration_s + mean_qscore to Reads)
- Modify: `bin/ingest.py` (load full sequencing summary, populate new columns)
- Create: `bin/calibrate/signal.py` (load sequencing summary for signal duration)
- Create: `tests/test_signal.py`
- Modify: `tests/test_mkdb.py` (verify new tables exist)

**Context:** Extend the existing database schema to capture signal duration (from sequencing summary `duration` column), mean Q-score, and run provenance. The existing `bin/ingest.py` currently loads only `read_id` and `end_reason` from the summary (line 367-368). We need to also load `duration`, `mean_qscore_template`, and `sequence_length_template`.

The existing Reads INSERT uses 21 columns (line 192-201 in ingest.py). We're adding 2 new columns: `signal_duration_s` and `mean_qscore`, bringing the total to 23.

**Step 1: Write the failing tests**

`tests/test_signal.py`:
```python
"""Tests for sequencing summary parsing."""
from __future__ import annotations

from pathlib import Path

import pytest

from calibrate.signal import load_sequencing_summary


@pytest.fixture
def summary_tsv(tmp_path: Path) -> Path:
    p = tmp_path / "sequencing_summary.txt"
    p.write_text(
        "read_id\tduration\tend_reason\tsequence_length_template\t"
        "mean_qscore_template\tbarcode_arrangement\n"
        "read-001\t1.5\tsignal_positive\t500\t12.3\tbarcode05\n"
        "read-002\t0.8\tdata_service_unblock_mux_change\t200\t9.1\tunclassified\n"
        "read-003\t2.1\tsignal_positive\t750\t15.0\tbarcode10\n"
    )
    return p


class TestLoadSequencingSummary:
    def test_returns_dict_keyed_by_read_id(self, summary_tsv: Path):
        data = load_sequencing_summary(summary_tsv)
        assert "read-001" in data
        assert "read-002" in data
        assert "read-003" in data

    def test_contains_duration(self, summary_tsv: Path):
        data = load_sequencing_summary(summary_tsv)
        assert data["read-001"]["duration"] == pytest.approx(1.5)

    def test_contains_end_reason(self, summary_tsv: Path):
        data = load_sequencing_summary(summary_tsv)
        assert data["read-001"]["end_reason"] == "signal_positive"

    def test_contains_mean_qscore(self, summary_tsv: Path):
        data = load_sequencing_summary(summary_tsv)
        assert data["read-001"]["mean_qscore"] == pytest.approx(12.3)

    def test_contains_sequence_length(self, summary_tsv: Path):
        data = load_sequencing_summary(summary_tsv)
        assert data["read-001"]["sequence_length"] == 500
```

Add to `tests/test_mkdb.py`:
```python
def test_reads_table_has_signal_columns(self, db_path):
    """Reads table must have signal_duration_s and mean_qscore columns."""
    conn = sqlite3.connect(db_path)
    cols = [row[1] for row in conn.execute("PRAGMA table_info(Reads)").fetchall()]
    conn.close()
    assert "signal_duration_s" in cols
    assert "mean_qscore" in cols

def test_run_metadata_table_exists(self, db_path):
    """RunMetadata table must exist."""
    conn = sqlite3.connect(db_path)
    tables = [row[0] for row in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'"
    ).fetchall()]
    conn.close()
    assert "RunMetadata" in tables
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_signal.py tests/test_mkdb.py -v`
Expected: FAIL

**Step 3: Implement**

`bin/calibrate/signal.py`:
```python
"""Load per-read metadata from MinKNOW sequencing summary files."""
from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_sequencing_summary(path: Path) -> dict[str, dict]:
    """Load sequencing summary TSV into a dict keyed by read_id.

    Returns
    -------
    dict[str, dict]
        Mapping of read_id to dict with keys: duration, end_reason,
        mean_qscore, sequence_length.
    """
    cols = [
        "read_id",
        "duration",
        "end_reason",
        "sequence_length_template",
        "mean_qscore_template",
    ]
    df = pd.read_csv(path, sep="\t", usecols=cols)
    result: dict[str, dict] = {}
    for _, row in df.iterrows():
        result[row["read_id"]] = {
            "duration": float(row["duration"]),
            "end_reason": str(row["end_reason"]),
            "mean_qscore": float(row["mean_qscore_template"]),
            "sequence_length": int(row["sequence_length_template"]),
        }
    return result
```

Modify `bin/mkdb.py`: Add `signal_duration_s REAL` and `mean_qscore REAL` after `trunc_level TEXT` in the Reads CREATE TABLE. Add new tables:

```sql
CREATE TABLE IF NOT EXISTS RunMetadata (
    run_id TEXT PRIMARY KEY,
    flow_cell_id TEXT,
    device_id TEXT,
    sample_id TEXT,
    experiment_id TEXT,
    kit TEXT,
    protocol_run_id TEXT,
    start_time TEXT,
    basecall_model TEXT,
    source_bam_count INTEGER,
    source_bam_paths TEXT,
    merge_timestamp TEXT
);

CREATE TABLE IF NOT EXISTS ReadRun (
    read_id TEXT PRIMARY KEY,
    run_id TEXT,
    FOREIGN KEY (run_id) REFERENCES RunMetadata(run_id)
);
```

Modify `bin/ingest.py`: Update `insert_read()` SQL to include 23 columns. Update both `read_data` tuples to pass `signal_duration_s` and `mean_qscore` from the enhanced summary data.

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_signal.py tests/test_mkdb.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add bin/calibrate/signal.py bin/mkdb.py bin/ingest.py tests/test_signal.py tests/test_mkdb.py
git commit -m "feat: extend schema with signal duration, mean qscore, and run provenance"
```

---

### Task 4: CLI Entry Point & Config Generation

**Files:**
- Create: `bin/calibrate/cli.py`
- Create: `bin/calibrate/config_gen.py`
- Create: `tests/test_cli.py`
- Create: `tests/test_config_gen.py`

**Context:** The CLI ties together discover, validate, merge, and ingest. The `config_gen` module generates a construct TOML and sample sheet from user-provided barcode-target mappings. CLI subcommands: `discover`, `validate`, `ingest`, `viz`, `export`.

**Step 1: Write the failing tests**

`tests/test_config_gen.py`:
```python
"""Tests for SMA-seq config file generation."""
from __future__ import annotations

from pathlib import Path

import pytest

from calibrate.config_gen import generate_construct_toml, generate_sample_sheet


class TestGenerateConstructToml:
    def test_writes_valid_toml(self, tmp_path: Path):
        targets = [
            {"barcode1": "NB05", "barcode2": "NB10", "alias": "V04_2_fwd"},
            {"barcode1": "NB10", "barcode2": "NB05", "alias": "V04_2_rev"},
        ]
        out = tmp_path / "construct.toml"
        generate_construct_toml(
            targets=targets,
            kit="SQK-NBD114-96",
            mode="dual_independent",
            output_path=out,
            mask1_front="AAGGTTAA",
            mask1_rear="CAGCACCT",
            mask2_front="AGGTGCTG",
            mask2_rear="TTAACCTT",
        )
        assert out.exists()
        import tomllib
        with out.open("rb") as f:
            data = tomllib.load(f)
        assert data["arrangement"]["kit"] == "SQK-NBD114-96"
        assert data["sma"]["mode"] == "dual_independent"
        assert len(data["sma"]["targets"]) == 2

    def test_start_only_mode_omits_mask2(self, tmp_path: Path):
        targets = [{"barcode1": "NB05", "alias": "V04_2_fwd"}]
        out = tmp_path / "construct.toml"
        generate_construct_toml(
            targets=targets,
            kit="SQK-NBD114-96",
            mode="start_only",
            output_path=out,
            mask1_front="AAGGTTAA",
            mask1_rear="CAGCACCT",
        )
        import tomllib
        with out.open("rb") as f:
            data = tomllib.load(f)
        assert "mask2_front" not in data["arrangement"]


class TestGenerateSampleSheet:
    def test_writes_csv(self, tmp_path: Path):
        entries = [
            {"barcode": "barcode05--barcode10", "alias": "V04_2_fwd", "type": "test_sample"},
        ]
        out = tmp_path / "sample_sheet.csv"
        generate_sample_sheet(
            entries=entries,
            flow_cell_id="FBD69411",
            kit="SQK-NBD114-96",
            experiment_id="my_exp",
            output_path=out,
        )
        assert out.exists()
        lines = out.read_text().strip().splitlines()
        assert len(lines) == 2  # header + 1 entry
        assert "barcode05--barcode10" in lines[1]
```

`tests/test_cli.py`:
```python
"""Tests for the calibrate CLI entry point."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


class TestCLIDiscoverSubcommand:
    def test_discover_prints_runs(self, tmp_path: Path):
        """discover subcommand should list found runs."""
        run_dir = tmp_path / "sample" / "20251228_2219_MD_FBD69411_abc123"
        run_dir.mkdir(parents=True)
        (run_dir / "final_summary_FBD69411_abc123_xyz.txt").write_text(
            "flow_cell_id=FBD69411\n"
            "device_id=MD-100098\n"
            "sample_id=sample\n"
            "experiment_id=exp1\n"
            "started=2025-12-28T22:19:00Z\n"
            "protocol_run_id=proto-1\n"
        )
        (run_dir / "bam_pass").mkdir()

        result = subprocess.run(
            [sys.executable, "-m", "calibrate.cli", "discover", str(tmp_path)],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq/bin",
        )
        assert result.returncode == 0
        assert "FBD69411" in result.stdout
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_config_gen.py tests/test_cli.py -v`
Expected: FAIL

**Step 3: Implement config generation and CLI**

`bin/calibrate/config_gen.py` — generates construct TOML using `tomli_w` and sample sheet CSV matching MinKNOW format. Functions: `generate_construct_toml()`, `generate_sample_sheet()`.

`bin/calibrate/cli.py` — argparse with subcommands:
- `discover <root_dir>` — lists discovered runs
- `validate <root_dir>` — checks basecalling consistency
- `ingest <root_dir> -c <construct.toml> -o <output_dir>` — full pipeline (discover → validate → merge → ingest)
- `viz <db_path> [db_path...]` — launches browser visualization
- `export <db_path> [--output <dir>] [--compare]` — static HTML export

Add `__main__.py` to `bin/calibrate/` so `python -m calibrate.cli` works.

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_config_gen.py tests/test_cli.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add bin/calibrate/cli.py bin/calibrate/config_gen.py bin/calibrate/__main__.py tests/test_config_gen.py tests/test_cli.py
git commit -m "feat: add CLI entry point and config generation"
```

---

## Phase 2: Visualization App — Core

### Task 5: Viz App Scaffolding & Distribution Computation

**Files:**
- Create: `bin/calibrate_viz/__init__.py`
- Create: `bin/calibrate_viz/app.py`
- Create: `bin/calibrate_viz/api.py`
- Create: `bin/calibrate_viz/distributions.py`
- Create: `tests/test_distributions.py`
- Create: `bin/calibrate_viz/templates/base.html`
- Create: `bin/calibrate_viz/static/style.css`

**Context:** The distribution computation module is the analytical core. It loads data from SQLite, computes KDE distributions using `scipy.stats.gaussian_kde`, applies Savitzky-Golay smoothing for peak detection via `scipy.signal.savgol_filter` and `scipy.signal.find_peaks`. Must support grouping by barcode, target, end_reason. The FastAPI app follows the same pattern as `bin/viz/app.py`.

**Step 1: Write the failing tests**

`tests/test_distributions.py`:
```python
"""Tests for KDE distribution computation."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
import pytest

from calibrate_viz.distributions import (
    compute_kde,
    load_distribution_data,
    find_kde_peaks,
)


@pytest.fixture
def sample_db(tmp_path: Path) -> Path:
    """Create a minimal SQLite DB with reads for distribution testing."""
    db_path = tmp_path / "test.db"
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE Reads (
            read_id TEXT PRIMARY KEY,
            readlen INTEGER,
            signal_duration_s REAL,
            mean_qscore REAL,
            ER TEXT,
            bc_start_id TEXT,
            bc_end_id TEXT,
            bc_start_conf REAL,
            bc_end_conf REAL,
            tgt_id TEXT,
            trunc_level TEXT
        )
    """)
    # Insert reads with known distributions
    reads = [
        ("r1", 500, 1.5, 12.0, "signal_positive", "nb05", "nb10", 0.9, 0.8, "V04_2_fwd", "full_length"),
        ("r2", 480, 1.4, 11.5, "signal_positive", "nb05", "nb10", 0.85, 0.75, "V04_2_fwd", "full_length"),
        ("r3", 510, 1.6, 12.5, "signal_positive", "nb05", "nb10", 0.95, 0.9, "V04_2_fwd", "full_length"),
        ("r4", 200, 0.7, 9.0, "data_service_unblock_mux_change", "nb05", None, 0.7, 0.1, "V04_2_fwd", "bc1_target"),
        ("r5", 300, 0.9, 10.0, "signal_positive", "nb10", "nb05", 0.88, 0.82, "V04_2_rev", "full_length"),
    ]
    conn.executemany(
        "INSERT INTO Reads VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        reads,
    )
    conn.commit()
    conn.close()
    return db_path


class TestComputeKDE:
    def test_returns_x_and_y_arrays(self):
        data = np.array([100, 200, 300, 400, 500] * 10)
        x, y = compute_kde(data)
        assert len(x) == len(y)
        assert len(x) > 0

    def test_y_values_are_nonnegative(self):
        data = np.array([100, 200, 300, 400, 500] * 10)
        _, y = compute_kde(data)
        assert np.all(y >= 0)


class TestLoadDistributionData:
    def test_loads_readlen(self, sample_db: Path):
        data = load_distribution_data(sample_db, column="readlen")
        assert len(data) == 5

    def test_loads_signal_duration(self, sample_db: Path):
        data = load_distribution_data(sample_db, column="signal_duration_s")
        assert len(data) == 5

    def test_groups_by_barcode(self, sample_db: Path):
        groups = load_distribution_data(
            sample_db, column="readlen", group_by="bc_start_id"
        )
        assert "nb05" in groups
        assert "nb10" in groups

    def test_groups_by_end_reason(self, sample_db: Path):
        groups = load_distribution_data(
            sample_db, column="readlen", group_by="ER"
        )
        assert "signal_positive" in groups


class TestFindKDEPeaks:
    def test_finds_single_peak(self):
        x = np.linspace(0, 10, 500)
        y = np.exp(-((x - 5) ** 2) / 0.5)
        peaks = find_kde_peaks(x, y)
        assert len(peaks) >= 1
        assert any(abs(p - 5.0) < 0.5 for p in peaks)
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_distributions.py -v`
Expected: FAIL

**Step 3: Implement**

`bin/calibrate_viz/distributions.py`:
```python
"""KDE distribution computation with peak detection for barcode calibration."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks, savgol_filter
from scipy.stats import gaussian_kde


def compute_kde(
    data: np.ndarray,
    n_points: int = 512,
    bw_method: str | float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute a KDE over the data range.

    Returns (x, y) arrays suitable for plotting.
    """
    if len(data) < 2:
        return np.array([]), np.array([])
    kde = gaussian_kde(data, bw_method=bw_method)
    x_min, x_max = data.min(), data.max()
    margin = (x_max - x_min) * 0.1
    x = np.linspace(x_min - margin, x_max + margin, n_points)
    y = kde(x)
    return x, y


def find_kde_peaks(
    x: np.ndarray,
    y: np.ndarray,
    window_length: int = 31,
    polyorder: int = 3,
    prominence: float = 0.01,
) -> list[float]:
    """Find peaks in a KDE curve using Savitzky-Golay smoothing.

    Returns x-values of detected peaks.
    """
    if len(y) < window_length:
        return []
    smoothed = savgol_filter(y, window_length, polyorder)
    peak_indices, _ = find_peaks(smoothed, prominence=prominence)
    return [float(x[i]) for i in peak_indices]


def load_distribution_data(
    db_path: Path,
    column: str,
    group_by: str | None = None,
    where: str | None = None,
) -> dict[str, np.ndarray] | np.ndarray:
    """Load a numeric column from the Reads table.

    If ``group_by`` is provided, returns a dict mapping group values
    to arrays. Otherwise returns a single array.
    """
    conn = sqlite3.connect(db_path)

    if group_by:
        query = f"SELECT {group_by}, {column} FROM Reads WHERE {column} IS NOT NULL"
        if where:
            query += f" AND {where}"
        rows = conn.execute(query).fetchall()
        conn.close()

        groups: dict[str, list[float]] = {}
        for grp, val in rows:
            key = str(grp) if grp is not None else "unknown"
            groups.setdefault(key, []).append(float(val))
        return {k: np.array(v) for k, v in groups.items()}
    else:
        query = f"SELECT {column} FROM Reads WHERE {column} IS NOT NULL"
        if where:
            query += f" AND {where}"
        rows = conn.execute(query).fetchall()
        conn.close()
        return np.array([float(r[0]) for r in rows])
```

`bin/calibrate_viz/app.py` — FastAPI app with multi-DB support. Accepts list of DB paths on startup, stores in global state. Serves templates and static files.

`bin/calibrate_viz/api.py` — REST endpoints:
- `GET /api/overview` — summary stats per loaded DB
- `GET /api/distributions?column=readlen&group_by=bc_start_id&db=0` — KDE data as JSON
- `GET /api/distributions/peaks?column=readlen&db=0` — peak positions

`bin/calibrate_viz/templates/base.html` — same sidebar pattern as config visualizer, with DB selector dropdown when multiple databases loaded.

`bin/calibrate_viz/static/style.css` — reuse from config visualizer with additions for KDE plots.

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_distributions.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add bin/calibrate_viz/ tests/test_distributions.py
git commit -m "feat: add calibrate viz app scaffolding and KDE distribution computation"
```

---

### Task 6: D3 Distribution Page with Draggable Thresholds

**Files:**
- Create: `bin/calibrate_viz/static/distributions.js`
- Create: `bin/calibrate_viz/templates/distributions.html`
- Modify: `bin/calibrate_viz/api.py` (add distribution-related endpoints)
- Modify: `bin/calibrate_viz/app.py` (add page route)

**Context:** The distribution page renders two paired KDE plots (signal duration + read length) with D3.js. Each plot supports:
- Multiple curves (one per group) with a legend
- Draggable vertical threshold lines that fire HTMX requests on release
- Reference markers showing expected product sizes
- End reason subplot below each main plot

The D3 rendering follows the pattern in `bin/viz/static/construct.js` — a `renderKDE()` function that takes data from a JSON script block or HTMX endpoint.

The draggable threshold interaction:
1. User drags a vertical line on the KDE plot
2. On drag end, the threshold value is sent to `/api/threshold-impact?column=readlen&threshold=350`
3. Server computes how many reads fall above/below the threshold per group
4. Returns an HTMX fragment updating the stats table below the plot

**Step 1: Write the D3 visualization code**

`bin/calibrate_viz/static/distributions.js`:
```javascript
// KDE distribution plots with draggable thresholds
// Depends on D3.js v7 (loaded from CDN in base.html)

function renderKDE(containerId, data, options) {
    // data: { groups: { "nb05": { x: [...], y: [...] }, ... }, peaks: { "nb05": [123, 456] } }
    // options: { xlabel, ylabel, thresholds: [{ value, label, color }], expectedSizes: [{ value, label }] }

    const container = d3.select(containerId);
    container.selectAll("*").remove();

    const margin = { top: 30, right: 120, bottom: 50, left: 60 };
    const width = 800 - margin.left - margin.right;
    const height = 350 - margin.top - margin.bottom;

    const svg = container.append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);

    const groups = Object.entries(data.groups);
    if (groups.length === 0) return;

    // Scales
    let allX = [], allY = [];
    groups.forEach(([_, g]) => { allX.push(...g.x); allY.push(...g.y); });

    const xScale = d3.scaleLinear()
        .domain([d3.min(allX), d3.max(allX)])
        .range([0, width]);
    const yScale = d3.scaleLinear()
        .domain([0, d3.max(allY) * 1.1])
        .range([height, 0]);
    const color = d3.scaleOrdinal(d3.schemeTableau10);

    // Axes
    svg.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale));
    svg.append("g").call(d3.axisLeft(yScale).ticks(5));

    // Axis labels
    svg.append("text").attr("x", width / 2).attr("y", height + 40)
        .attr("text-anchor", "middle").text(options.xlabel || "");
    svg.append("text").attr("transform", "rotate(-90)")
        .attr("x", -height / 2).attr("y", -45)
        .attr("text-anchor", "middle").text(options.ylabel || "Density");

    // KDE curves
    groups.forEach(([name, g], i) => {
        const line = d3.line().x((_, j) => xScale(g.x[j])).y((_, j) => yScale(g.y[j]));
        svg.append("path")
            .datum(g.y)
            .attr("d", line)
            .attr("fill", "none")
            .attr("stroke", color(i))
            .attr("stroke-width", 2);
    });

    // Expected size markers
    if (options.expectedSizes) {
        options.expectedSizes.forEach(s => {
            svg.append("line")
                .attr("x1", xScale(s.value)).attr("x2", xScale(s.value))
                .attr("y1", 0).attr("y2", height)
                .attr("stroke", "#999").attr("stroke-dasharray", "4,4");
            svg.append("text")
                .attr("x", xScale(s.value) + 4).attr("y", 12)
                .attr("font-size", "11px").attr("fill", "#666")
                .text(s.label);
        });
    }

    // Draggable thresholds
    if (options.thresholds) {
        options.thresholds.forEach(t => {
            const threshLine = svg.append("line")
                .attr("x1", xScale(t.value)).attr("x2", xScale(t.value))
                .attr("y1", 0).attr("y2", height)
                .attr("stroke", t.color || "red").attr("stroke-width", 2)
                .attr("cursor", "ew-resize");

            const threshLabel = svg.append("text")
                .attr("x", xScale(t.value) + 4).attr("y", -5)
                .attr("font-size", "12px").attr("fill", t.color || "red")
                .text(`${t.label}: ${t.value.toFixed(1)}`);

            const drag = d3.drag()
                .on("drag", function(event) {
                    const newX = Math.max(0, Math.min(width, event.x));
                    const newVal = xScale.invert(newX);
                    threshLine.attr("x1", newX).attr("x2", newX);
                    threshLabel.attr("x", newX + 4)
                        .text(`${t.label}: ${newVal.toFixed(1)}`);
                })
                .on("end", function(event) {
                    const newVal = xScale.invert(Math.max(0, Math.min(width, event.x)));
                    // Fire HTMX request to update stats
                    const statsDiv = document.getElementById(t.statsTarget || "threshold-stats");
                    if (statsDiv) {
                        htmx.ajax("GET",
                            `/api/threshold-impact?column=${options.column}&threshold=${newVal}&group_by=${options.groupBy || ""}`,
                            {target: statsDiv});
                    }
                });

            threshLine.call(drag);
        });
    }

    // Legend
    const legend = svg.append("g").attr("transform", `translate(${width + 10}, 0)`);
    groups.forEach(([name, _], i) => {
        const g = legend.append("g").attr("transform", `translate(0, ${i * 20})`);
        g.append("rect").attr("width", 12).attr("height", 12).attr("fill", color(i));
        g.append("text").attr("x", 16).attr("y", 10).attr("font-size", "12px").text(name);
    });

    return svg;
}


function renderEndReasonSubplot(containerId, data, xScale) {
    // Stacked area chart of end reasons using shared x-axis from main KDE
    const container = d3.select(containerId);
    container.selectAll("*").remove();

    const margin = { top: 5, right: 120, bottom: 30, left: 60 };
    const width = 800 - margin.left - margin.right;
    const height = 100 - margin.top - margin.bottom;

    const svg = container.append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);

    // End reason colors
    const erColors = {
        "signal_positive": "#2ecc71",
        "data_service_unblock_mux_change": "#e74c3c",
        "unblock_mux_change": "#f39c12",
        "signal_negative": "#9b59b6",
        "mux_change": "#3498db",
    };

    const groups = Object.entries(data.groups);
    const yScale = d3.scaleLinear().domain([0, 1]).range([height, 0]);

    groups.forEach(([name, g], i) => {
        const barWidth = width / g.counts.length;
        g.counts.forEach((c, j) => {
            let yPos = 0;
            Object.entries(c).forEach(([er, frac]) => {
                svg.append("rect")
                    .attr("x", j * barWidth)
                    .attr("y", yScale(yPos + frac))
                    .attr("width", barWidth)
                    .attr("height", yScale(yPos) - yScale(yPos + frac))
                    .attr("fill", erColors[er] || "#bdc3c7")
                    .attr("opacity", 0.7);
                yPos += frac;
            });
        });
    });

    svg.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale));
}
```

**Step 2: Create the template**

`bin/calibrate_viz/templates/distributions.html` extends `base.html`. Contains:
- Grouping mode tabs (All / Per-barcode / Per-target / Per-end_reason)
- Two D3 container divs (signal duration + read length)
- End reason subplot divs below each
- Threshold stats table div (updated via HTMX)
- JS block that fetches `/api/distributions` and calls `renderKDE()`

**Step 3: Wire up API endpoints**

Add to `bin/calibrate_viz/api.py`:
- `GET /api/distributions` — returns JSON with grouped KDE data
- `GET /api/distributions/peaks` — returns peak positions
- `GET /api/threshold-impact` — returns HTML fragment with classification counts at a given threshold

**Step 4: Verify page renders**

Run app manually: `cd /tmp/ont-sma-seq && python -c "from calibrate_viz.app import app; print('OK')"`
Expected: App imports successfully

**Step 5: Commit**

```bash
git add bin/calibrate_viz/
git commit -m "feat: add D3 distribution page with draggable thresholds and end reason subplots"
```

---

### Task 7: Barcode Confusion Matrix & Confidence Page

**Files:**
- Create: `bin/calibrate_viz/confusion.py`
- Create: `bin/calibrate_viz/static/confidence.js`
- Create: `bin/calibrate_viz/templates/confidence.html`
- Create: `tests/test_confusion.py`
- Modify: `bin/calibrate_viz/api.py`

**Context:** The confusion matrix computes true vs. classified barcode assignments. Ground truth is established by checking whether a read aligns well to its assigned target (low edit distance relative to target length). The confidence page shows KDE distributions of `bc_start_conf` and `bc_end_conf`, a scatter plot of start vs. end confidence colored by trunc_level, and the confusion matrix — all updating as thresholds are dragged.

**Step 1: Write the failing tests**

`tests/test_confusion.py`:
```python
"""Tests for barcode confusion matrix computation."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from calibrate_viz.confusion import (
    compute_confusion_matrix,
    compute_threshold_impact,
)


@pytest.fixture
def classified_db(tmp_path: Path) -> Path:
    """DB with reads that have known correct/incorrect barcode assignments."""
    db_path = tmp_path / "test.db"
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE Reads (
            read_id TEXT PRIMARY KEY,
            readlen INTEGER,
            ed INTEGER,
            bc_start_id TEXT,
            bc_start_conf REAL,
            bc_end_id TEXT,
            bc_end_conf REAL,
            tgt_id TEXT,
            trunc_level TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE Target (
            tgt_id TEXT PRIMARY KEY,
            tgt_reflen INTEGER
        )
    """)
    conn.execute("INSERT INTO Target VALUES ('V04_2_fwd', 500)")
    conn.execute("INSERT INTO Target VALUES ('V04_2_rev', 500)")
    # Correctly classified reads (low ed relative to target length)
    conn.execute("INSERT INTO Reads VALUES ('r1', 500, 5, 'nb05', 0.9, 'nb10', 0.8, 'V04_2_fwd', 'full_length')")
    conn.execute("INSERT INTO Reads VALUES ('r2', 490, 8, 'nb05', 0.85, 'nb10', 0.7, 'V04_2_fwd', 'full_length')")
    # Misclassified read (high ed = wrong barcode assignment likely)
    conn.execute("INSERT INTO Reads VALUES ('r3', 480, 200, 'nb05', 0.4, 'nb10', 0.3, 'V04_2_fwd', 'bc1_target')")
    # Different target
    conn.execute("INSERT INTO Reads VALUES ('r4', 500, 3, 'nb10', 0.92, 'nb05', 0.85, 'V04_2_rev', 'full_length')")
    conn.commit()
    conn.close()
    return db_path


class TestComputeConfusionMatrix:
    def test_returns_matrix_dict(self, classified_db: Path):
        matrix = compute_confusion_matrix(classified_db, ed_threshold=0.1)
        assert isinstance(matrix, dict)
        assert "labels" in matrix
        assert "counts" in matrix

    def test_correct_assignments_on_diagonal(self, classified_db: Path):
        matrix = compute_confusion_matrix(classified_db, ed_threshold=0.1)
        # Most reads should be on the diagonal
        total_diagonal = sum(
            matrix["counts"][i][i] for i in range(len(matrix["labels"]))
        )
        assert total_diagonal >= 3


class TestComputeThresholdImpact:
    def test_returns_counts_per_trunc_level(self, classified_db: Path):
        result = compute_threshold_impact(
            classified_db,
            start_barcode_min=0.6,
            full_length_threshold=0.75,
        )
        assert "full_length" in result
        assert "bc1_target" in result
        assert isinstance(result["full_length"], int)
```

**Step 2-5: Implement, test, commit**

Implement `confusion.py` with:
- `compute_confusion_matrix(db_path, ed_threshold)` — queries Reads + Target, classifies correct assignment as `ed / tgt_reflen < ed_threshold`, builds an NxN matrix
- `compute_threshold_impact(db_path, start_barcode_min, full_length_threshold)` — re-classifies all reads at given thresholds, returns dict of trunc_level → count
- `get_affected_reads(db_path, old_thresholds, new_thresholds)` — returns list of reads that change classification

Implement `confidence.js` with D3 scatter plot (`bc_start_conf` vs `bc_end_conf`), confidence KDEs, and draggable threshold lines that update the confusion matrix and classification counts via HTMX.

Implement `confidence.html` template.

```bash
git add bin/calibrate_viz/confusion.py bin/calibrate_viz/static/confidence.js bin/calibrate_viz/templates/confidence.html tests/test_confusion.py
git commit -m "feat: add barcode confusion matrix and confidence page with interactive thresholds"
```

---

### Task 8: Barcode Separation Heatmap

**Files:**
- Create: `bin/calibrate_viz/separation.py`
- Create: `bin/calibrate_viz/static/separation.js`
- Create: `bin/calibrate_viz/templates/separation.html`
- Create: `tests/test_separation.py`
- Modify: `bin/calibrate_viz/api.py`

**Context:** Computes pairwise edit distances between all barcodes used in an experiment, then analyzes empirical separation: for each barcode, what is the edit distance distribution of true matches vs. best-alternative matches? The heatmap is a D3 matrix with clickable cells. The separation table shows per-barcode metrics sortable by separation gap.

The existing `bin/barcodes.py` has the 96 native barcodes and `classify_barcode()` with edlib HW alignment. We reuse these for computing pairwise distances. The key computation: for each read assigned to barcode X, compute the edit distance to barcode X (true match) AND to the next-best barcode Y (best alternative). The gap between these two distributions is the "separation" metric.

**Step 1: Write the failing tests**

`tests/test_separation.py`:
```python
"""Tests for barcode separation analysis."""
from __future__ import annotations

import pytest

from calibrate_viz.separation import (
    compute_pairwise_distances,
    compute_separation_metrics,
)


class TestPairwiseDistances:
    def test_returns_symmetric_matrix(self):
        barcodes = {"nb05": "CACAAAGACACCGACAACTTTCTT", "nb10": "GAGAGGACAAAGGTTTCAACGCTT"}
        matrix = compute_pairwise_distances(barcodes)
        assert matrix["nb05"]["nb10"] == matrix["nb10"]["nb05"]

    def test_diagonal_is_zero(self):
        barcodes = {"nb05": "CACAAAGACACCGACAACTTTCTT", "nb10": "GAGAGGACAAAGGTTTCAACGCTT"}
        matrix = compute_pairwise_distances(barcodes)
        assert matrix["nb05"]["nb05"] == 0
        assert matrix["nb10"]["nb10"] == 0

    def test_different_barcodes_have_nonzero_distance(self):
        barcodes = {"nb05": "CACAAAGACACCGACAACTTTCTT", "nb10": "GAGAGGACAAAGGTTTCAACGCTT"}
        matrix = compute_pairwise_distances(barcodes)
        assert matrix["nb05"]["nb10"] > 0


class TestSeparationMetrics:
    def test_returns_per_barcode_metrics(self):
        # Synthetic data: barcode nb05 has reads with low ED to itself, high to nb10
        read_eds = {
            "nb05": {"true_match": [2, 3, 1, 2], "next_best": [15, 16, 14, 15]},
            "nb10": {"true_match": [1, 2, 3, 1], "next_best": [14, 13, 15, 14]},
        }
        metrics = compute_separation_metrics(read_eds)
        assert "nb05" in metrics
        assert metrics["nb05"]["mean_true_ed"] < metrics["nb05"]["mean_next_best_ed"]
        assert metrics["nb05"]["separation_gap"] > 10
```

**Step 2-5: Implement, test, commit**

```bash
git add bin/calibrate_viz/separation.py bin/calibrate_viz/static/separation.js bin/calibrate_viz/templates/separation.html tests/test_separation.py
git commit -m "feat: add barcode separation heatmap and per-barcode metrics"
```

---

## Phase 3: Threshold Optimization, Comparison & Export

### Task 9: Threshold Optimization (ROC Curves)

**Files:**
- Create: `bin/calibrate_viz/thresholds.py`
- Create: `bin/calibrate_viz/static/thresholds.js`
- Create: `bin/calibrate_viz/templates/thresholds.html`
- Create: `tests/test_thresholds.py`
- Modify: `bin/calibrate_viz/api.py`

**Context:** Computes ROC-style curves for each threshold parameter. For `bc_start_conf`: at each candidate threshold value, compute true positive rate (correctly classified reads above threshold) and false positive rate (incorrectly classified reads above threshold). True/false is determined by alignment quality to assigned target.

Also provides a "recommend" function: given a target false-positive rate, find the optimal threshold and report sensitivity, specificity, F1 score. Selected thresholds can be exported to a construct TOML `[sma.confidence]` section.

**Step 1: Write the failing tests**

`tests/test_thresholds.py`:
```python
"""Tests for threshold optimization and ROC computation."""
from __future__ import annotations

import numpy as np
import pytest

from calibrate_viz.thresholds import (
    compute_roc,
    recommend_threshold,
    export_thresholds_to_toml,
)


class TestComputeROC:
    def test_returns_fpr_tpr_thresholds(self):
        # True positives have high confidence, false positives have low
        true_confs = np.array([0.9, 0.85, 0.8, 0.95, 0.88])
        false_confs = np.array([0.3, 0.2, 0.4, 0.15, 0.35])
        fpr, tpr, thresholds = compute_roc(true_confs, false_confs)
        assert len(fpr) == len(tpr) == len(thresholds)
        assert fpr[0] <= fpr[-1]  # sorted

    def test_perfect_separation_gives_auc_near_1(self):
        true_confs = np.array([0.9, 0.95, 0.85, 0.92])
        false_confs = np.array([0.1, 0.05, 0.15, 0.08])
        fpr, tpr, _ = compute_roc(true_confs, false_confs)
        # AUC should be very high
        auc = np.trapz(tpr, fpr)
        assert auc > 0.95


class TestRecommendThreshold:
    def test_returns_threshold_and_metrics(self):
        true_confs = np.array([0.9, 0.85, 0.8, 0.95])
        false_confs = np.array([0.3, 0.2, 0.4, 0.15])
        result = recommend_threshold(true_confs, false_confs, target_fpr=0.05)
        assert "threshold" in result
        assert "sensitivity" in result
        assert "specificity" in result
        assert "f1" in result
        assert result["threshold"] > 0.4


class TestExportToToml:
    def test_writes_valid_toml(self, tmp_path):
        out = tmp_path / "construct.toml"
        export_thresholds_to_toml(
            output_path=out,
            start_barcode_min=0.65,
            full_length_threshold=0.78,
            flank_max_error_rate=0.5,
        )
        assert out.exists()
        import tomllib
        with out.open("rb") as f:
            data = tomllib.load(f)
        assert data["sma"]["confidence"]["start_barcode_min"] == 0.65
```

**Step 2-5: Implement, test, commit**

```bash
git add bin/calibrate_viz/thresholds.py bin/calibrate_viz/static/thresholds.js bin/calibrate_viz/templates/thresholds.html tests/test_thresholds.py
git commit -m "feat: add ROC-based threshold optimization with TOML export"
```

---

### Task 10: Cross-Experiment Comparison

**Files:**
- Create: `bin/calibrate_viz/comparison.py`
- Create: `bin/calibrate_viz/templates/compare.html`
- Create: `tests/test_comparison.py`
- Modify: `bin/calibrate_viz/api.py`
- Modify: `bin/calibrate_viz/app.py`

**Context:** When multiple databases are loaded, the comparison page shows overlay distributions (multiple experiments on same axes) and side-by-side panels. The comparison table shows per-experiment summary metrics. The app needs to support loading N databases at startup and switching between them.

**Step 1: Write the failing tests**

`tests/test_comparison.py`:
```python
"""Tests for cross-experiment comparison."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from calibrate_viz.comparison import (
    compute_experiment_summary,
    compute_overlay_distributions,
)


def _make_db(tmp_path: Path, name: str, reads: list[tuple]) -> Path:
    db_path = tmp_path / f"{name}.db"
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE Reads (
            read_id TEXT PRIMARY KEY, readlen INTEGER,
            signal_duration_s REAL, mean_qscore REAL, ER TEXT,
            bc_start_id TEXT, bc_start_conf REAL,
            bc_end_id TEXT, bc_end_conf REAL,
            tgt_id TEXT, trunc_level TEXT, ed INTEGER
        )
    """)
    conn.executemany(
        "INSERT INTO Reads VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        reads,
    )
    conn.commit()
    conn.close()
    return db_path


class TestExperimentSummary:
    def test_returns_key_metrics(self, tmp_path: Path):
        db = _make_db(tmp_path, "exp1", [
            ("r1", 500, 1.5, 12.0, "signal_positive", "nb05", 0.9, "nb10", 0.8, "t1", "full_length", 5),
            ("r2", 300, 0.9, 10.0, "signal_positive", "nb05", 0.7, None, 0.1, "t1", "bc1_target", 20),
        ])
        summary = compute_experiment_summary(db)
        assert summary["total_reads"] == 2
        assert summary["mean_readlen"] == 400.0
        assert "trunc_proportions" in summary
        assert "end_reason_proportions" in summary


class TestOverlayDistributions:
    def test_overlay_returns_per_db_kde(self, tmp_path: Path):
        db1 = _make_db(tmp_path, "exp1", [
            ("r1", 500, 1.5, 12.0, "sp", "nb05", 0.9, "nb10", 0.8, "t1", "fl", 5),
            ("r2", 480, 1.4, 11.5, "sp", "nb05", 0.85, "nb10", 0.7, "t1", "fl", 8),
        ])
        db2 = _make_db(tmp_path, "exp2", [
            ("r1", 300, 0.9, 10.0, "sp", "nb05", 0.7, "nb10", 0.5, "t1", "bt", 20),
            ("r2", 290, 0.85, 9.5, "sp", "nb05", 0.65, "nb10", 0.4, "t1", "bt", 22),
        ])
        result = compute_overlay_distributions([db1, db2], column="readlen")
        assert len(result) == 2
```

**Step 2-5: Implement, test, commit**

```bash
git add bin/calibrate_viz/comparison.py bin/calibrate_viz/templates/compare.html tests/test_comparison.py
git commit -m "feat: add cross-experiment comparison with overlay and panel modes"
```

---

### Task 11: Static HTML Export

**Files:**
- Create: `bin/calibrate_viz/export.py`
- Create: `tests/test_calibrate_export.py`
- Modify: `bin/calibrate_viz/api.py` (add export endpoint)
- Modify: `bin/calibrate/cli.py` (add export subcommand)

**Context:** Follows the same self-contained export pattern as `bin/viz/export.py`. Each page is rendered as a standalone HTML file with inline CSS, D3.js bundled, and data embedded as JSON script blocks. SVG figures are also exported for manuscript use. The export reads from the SQLite database, computes all distributions/matrices/ROC curves, and bakes the results into the HTML.

**Step 1: Write the failing tests**

`tests/test_calibrate_export.py`:
```python
"""Tests for static HTML export of calibration results."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from calibrate_viz.export import export_calibration_report


@pytest.fixture
def populated_db(tmp_path: Path) -> Path:
    """Create a DB with enough data for a full export."""
    db_path = tmp_path / "test.db"
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE Reads (
            read_id TEXT PRIMARY KEY, readlen INTEGER,
            signal_duration_s REAL, mean_qscore REAL, ER TEXT,
            bc_start_id TEXT, bc_start_conf REAL, bc_start_ed INTEGER,
            bc_end_id TEXT, bc_end_conf REAL, bc_end_ed INTEGER,
            tgt_id TEXT, trunc_level TEXT, ed INTEGER
        )
    """)
    conn.execute("""
        CREATE TABLE Target (tgt_id TEXT PRIMARY KEY, tgt_reflen INTEGER)
    """)
    conn.execute("INSERT INTO Target VALUES ('V04_2_fwd', 500)")
    for i in range(50):
        conn.execute(
            "INSERT INTO Reads VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (f"r{i}", 400 + i * 5, 1.0 + i * 0.05, 10.0 + i * 0.1,
             "signal_positive", "nb05", 0.7 + i * 0.005, 3,
             "nb10", 0.6 + i * 0.005, 5,
             "V04_2_fwd", "full_length", 5 + i),
        )
    conn.commit()
    conn.close()
    return db_path


class TestExportCalibrationReport:
    def test_creates_output_directory(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        assert out.is_dir()

    def test_generates_all_html_pages(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        expected = ["index.html", "distributions.html", "confidence.html",
                     "separation.html", "thresholds.html"]
        for name in expected:
            assert (out / name).exists(), f"Missing {name}"

    def test_html_is_self_contained(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        html = (out / "index.html").read_text()
        assert "<style>" in html
        assert "d3.min.js" in html or "d3.v7" in html or "var DATA" in html

    def test_generates_svg_figures(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        figures = out / "figures"
        assert figures.is_dir()
        svgs = list(figures.glob("*.svg"))
        assert len(svgs) >= 2  # at least signal + readlen KDEs
```

**Step 2-5: Implement, test, commit**

```bash
git add bin/calibrate_viz/export.py tests/test_calibrate_export.py
git commit -m "feat: add static HTML and SVG export for calibration reports"
```

---

### Task 12: Integration Test & Final Polish

**Files:**
- Create: `tests/test_calibrate_integration.py`
- Modify: `bin/calibrate/cli.py` (any final wiring)
- Modify: `env/env.yml` (add scipy dependency)

**Context:** End-to-end test: create a minimal MinKNOW output directory with synthetic BAMs, run the full pipeline (discover → validate → merge → ingest), then verify the database has the expected data and the visualization app can load it.

**Step 1: Write the integration test**

`tests/test_calibrate_integration.py`:
```python
"""Integration test for the full calibrate pipeline."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pysam
import pytest

from calibrate.discover import discover_runs
from calibrate.merge import group_runs, find_bam_files
from calibrate.signal import load_sequencing_summary


@pytest.fixture
def full_experiment(tmp_path: Path) -> Path:
    """Create a complete minimal experiment with BAMs and summaries."""
    root = tmp_path / "experiment" / "sample" / "20251228_2219_MD_FBD69411_abc123"
    root.mkdir(parents=True)

    # Final summary
    (root / "final_summary_FBD69411_abc123_xyz.txt").write_text(
        "flow_cell_id=FBD69411\ndevice_id=MD-100098\n"
        "sample_id=sample\nexperiment_id=exp1\n"
        "started=2025-12-28T22:19:00Z\nprotocol_run_id=proto-1\n"
    )

    # Sequencing summary
    (root / "sequencing_summary_FBD69411_abc123_xyz.txt").write_text(
        "read_id\tduration\tend_reason\tsequence_length_template\t"
        "mean_qscore_template\tbarcode_arrangement\n"
        "read-001\t1.5\tsignal_positive\t500\t12.3\tbarcode05\n"
        "read-002\t0.8\tsignal_positive\t300\t10.1\tbarcode10\n"
    )

    # BAM with 2 unmapped reads
    bam_dir = root / "bam_pass"
    bam_dir.mkdir()
    bam_path = bam_dir / "reads_0.bam"
    header = {"HD": {"VN": "1.6", "SO": "unknown"}}
    save = pysam.set_verbosity(0)
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_out:
        for rid, seq in [("read-001", "ACGT" * 125), ("read-002", "ACGT" * 75)]:
            a = pysam.AlignedSegment()
            a.query_name = rid
            a.query_sequence = seq
            a.flag = 4
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            bam_out.write(a)
    pysam.set_verbosity(save)

    return tmp_path / "experiment"


class TestFullPipeline:
    def test_discover_finds_experiment(self, full_experiment: Path):
        runs = discover_runs(full_experiment)
        assert len(runs) == 1
        assert runs[0].flow_cell_id == "FBD69411"

    def test_groups_into_single_merge_group(self, full_experiment: Path):
        runs = discover_runs(full_experiment)
        groups = group_runs(runs)
        assert len(groups) == 1

    def test_finds_bam_files(self, full_experiment: Path):
        runs = discover_runs(full_experiment)
        bams = find_bam_files(runs[0].run_dir)
        assert len(bams) == 1

    def test_loads_sequencing_summary(self, full_experiment: Path):
        runs = discover_runs(full_experiment)
        assert runs[0].sequencing_summary_path is not None
        data = load_sequencing_summary(runs[0].sequencing_summary_path)
        assert len(data) == 2
        assert data["read-001"]["duration"] == pytest.approx(1.5)
        assert data["read-001"]["end_reason"] == "signal_positive"
```

**Step 2: Verify env.yml has scipy**

Add to `env/env.yml` pip section:
```yaml
    - scipy
```

**Step 3: Run full test suite**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/ -v --tb=short`
Expected: All tests pass

**Step 4: Commit**

```bash
git add tests/test_calibrate_integration.py env/env.yml
git commit -m "test: add integration test for full calibrate pipeline"
```

---

## Task Dependency Graph

```
Task 1 (Discovery) ──→ Task 2 (Merge) ──→ Task 3 (Schema + Ingest) ──→ Task 4 (CLI + Config Gen)
                                                      │
                                                      ↓
Task 5 (Viz Scaffolding + Distributions) ──→ Task 6 (D3 Distribution Page)
         │                                            │
         ↓                                            ↓
Task 7 (Confusion + Confidence) ──→ Task 8 (Separation Heatmap)
         │                                    │
         ↓                                    ↓
Task 9 (Threshold Optimization) ──→ Task 10 (Cross-Experiment Comparison)
                                              │
                                              ↓
                                   Task 11 (Static Export)
                                              │
                                              ↓
                                   Task 12 (Integration Test)
```

**Parallelizable groups:**
- Tasks 1-4 are sequential (CLI pipeline)
- Tasks 5-6 can start after Task 3 (need schema)
- Tasks 7-8 can run in parallel (independent pages)
- Tasks 9-10 depend on 7-8
- Task 11 depends on all viz tasks
- Task 12 is final
