# SMA-seq Experiment Preparation Tool — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a three-stage pipeline (scan, basecall, init) that prepares SMA-seq experiment raw data into a single untrimmed BAM and QC report.

**Architecture:** Three independent CLI tools sharing a JSON manifest. `sma_scan.py` discovers MinKNOW runs and writes a manifest. `sma_basecall.py` re-basecalls from pod5 (or merges existing BAMs) using the manifest. `sma_init.py` creates the SMA-seq database and QC report from the merged BAM.

**Tech Stack:** Python 3.13, pysam, pod5, pathlib, subprocess (for dorado/samtools), sqlite3, json. HTML report reuses patterns from `bin/report_template.py`.

**Codebase:** `/tmp/ont-sma-seq/` — existing SMA-seq pipeline with `bin/` scripts and `tests/` using pytest. Design doc: `docs/plans/2026-02-14-experiment-prep-design.md`.

---

## Task 1: `sma_scan.py` — Parse final_summary.txt

Core parser for the key-value `final_summary*.txt` files that MinKNOW writes per run.

**Files:**
- Create: `bin/sma_scan.py`
- Test: `tests/test_sma_scan.py`

**Step 1: Write failing test for parse_final_summary**

```python
# tests/test_sma_scan.py
"""Tests for sma_scan.py experiment scanner."""

from pathlib import Path
import json
import pytest

# Module under test lives in bin/ — add to path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))

from sma_scan import parse_final_summary


class TestParseFinalSummary:
    def test_extracts_flow_cell_id(self, tmp_path):
        fs = tmp_path / "final_summary_FBD69411_abc_123.txt"
        fs.write_text(
            "instrument=MD-100098\n"
            "flow_cell_id=FBD69411\n"
            "protocol_group_id=my_experiment\n"
            "protocol=sequencing/sequencing_MIN114_DNA_e8_2_400K:FLO-MIN114:SQK-NBD114-24:400\n"
            "protocol_run_id=abc-def\n"
            "started=2025-12-29T10:56:12.341976-05:00\n"
            "acquisition_stopped=2026-01-01T10:56:12.413307-05:00\n"
            "pod5_files_in_final_dest=72\n"
            "bam_files_in_final_dest=2165\n"
        )
        result = parse_final_summary(fs)
        assert result["flow_cell_id"] == "FBD69411"
        assert result["instrument"] == "MD-100098"
        assert result["protocol_group_id"] == "my_experiment"
        assert result["protocol_run_id"] == "abc-def"
        assert result["started"] == "2025-12-29T10:56:12.341976-05:00"
        assert result["acquisition_stopped"] == "2026-01-01T10:56:12.413307-05:00"
        assert result["pod5_files_in_final_dest"] == 72
        assert result["bam_files_in_final_dest"] == 2165

    def test_missing_keys_return_none(self, tmp_path):
        fs = tmp_path / "final_summary_X_Y_Z.txt"
        fs.write_text("instrument=MD-100098\nflow_cell_id=FBD69411\n")
        result = parse_final_summary(fs)
        assert result["flow_cell_id"] == "FBD69411"
        assert result.get("started") is None
```

**Step 2: Run test — expect FAIL (module not found)**

```bash
cd /tmp/ont-sma-seq && python -m pytest tests/test_sma_scan.py::TestParseFinalSummary -v
```

**Step 3: Implement parse_final_summary**

```python
#!/usr/bin/env python3
# sma_scan.py — Scan experiment directory, discover MinKNOW runs, output manifest

from pathlib import Path
from typing import Optional


def parse_final_summary(path: Path) -> dict:
    """Parse a MinKNOW final_summary*.txt file (key=value format).

    Integer-valued keys (file counts) are converted to int.
    All other values remain strings. Missing keys return None via dict.get().
    """
    INT_KEYS = {
        "pod5_files_in_final_dest", "pod5_files_in_fallback",
        "bam_files_in_final_dest", "bam_files_in_fallback",
        "fast5_files_in_final_dest", "fastq_files_in_final_dest",
        "bai_files_in_final_dest",
    }
    data: dict = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if "=" not in line:
                continue
            key, _, val = line.partition("=")
            if key in INT_KEYS:
                data[key] = int(val)
            else:
                data[key] = val
    return data
```

**Step 4: Run test — expect PASS**

```bash
cd /tmp/ont-sma-seq && python -m pytest tests/test_sma_scan.py::TestParseFinalSummary -v
```

**Step 5: Commit**

```bash
cd /tmp/ont-sma-seq && git add bin/sma_scan.py tests/test_sma_scan.py
git commit -m "feat(scan): add final_summary.txt parser"
```

---

## Task 2: Parse MinKNOW sample_sheet.csv

Extract experiment_id, kit, flow_cell_product_code, and barcode→alias mapping from the MinKNOW-format sample sheet.

**Files:**
- Modify: `bin/sma_scan.py`
- Modify: `tests/test_sma_scan.py`

**Step 1: Write failing test**

```python
class TestParseMinkNowSampleSheet:
    def test_extracts_experiment_and_kit(self, tmp_path):
        ss = tmp_path / "sample_sheet_FBD69411_20251228_2101_abc.csv"
        ss.write_text(
            "protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,"
            "flow_cell_product_code,kit,barcode,alias,type,rowNumber\n"
            "abc-def,MD-100098,FBD69411,,my_experiment,FLO-MIN114,"
            "SQK-NBD114-24,barcode02,V04_2,test_sample,1\n"
            "abc-def,MD-100098,FBD69411,,my_experiment,FLO-MIN114,"
            "SQK-NBD114-24,barcode04,V04_4,test_sample,2\n"
        )
        result = parse_minknow_sample_sheet(ss)
        assert result["experiment_id"] == "my_experiment"
        assert result["kit"] == "SQK-NBD114-24"
        assert result["flow_cell_product_code"] == "FLO-MIN114"
        assert result["barcode_aliases"] == {
            "barcode02": "V04_2",
            "barcode04": "V04_4",
        }

    def test_empty_sample_sheet_returns_empty_aliases(self, tmp_path):
        ss = tmp_path / "sample_sheet.csv"
        ss.write_text(
            "protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,"
            "flow_cell_product_code,kit,barcode,alias,type,rowNumber\n"
        )
        result = parse_minknow_sample_sheet(ss)
        assert result["barcode_aliases"] == {}
```

**Step 2: Run test — expect FAIL**

**Step 3: Implement parse_minknow_sample_sheet**

```python
import csv

def parse_minknow_sample_sheet(path: Path) -> dict:
    """Parse MinKNOW-format sample_sheet*.csv.

    Returns dict with experiment_id, kit, flow_cell_product_code,
    and barcode_aliases mapping.
    """
    result: dict = {
        "experiment_id": None,
        "kit": None,
        "flow_cell_product_code": None,
        "barcode_aliases": {},
    }
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if result["experiment_id"] is None:
                result["experiment_id"] = row.get("experiment_id") or None
                result["kit"] = row.get("kit") or None
                result["flow_cell_product_code"] = row.get("flow_cell_product_code") or None
            bc = row.get("barcode", "")
            alias = row.get("alias", "")
            if bc and alias:
                result["barcode_aliases"][bc] = alias
    return result
```

**Step 4: Run test — expect PASS**

**Step 5: Commit**

```bash
git commit -m "feat(scan): add MinKNOW sample sheet parser"
```

---

## Task 3: Parse report.json for software versions

Extract basecaller version, MinKNOW version, and bream version from the MinKNOW report JSON.

**Files:**
- Modify: `bin/sma_scan.py`
- Modify: `tests/test_sma_scan.py`

**Step 1: Write failing test**

```python
class TestParseReportJson:
    def test_extracts_software_versions(self, tmp_path):
        report = tmp_path / "report_FBD69411_20251228_abc.json"
        report.write_text(json.dumps({
            "protocol_run_info": {
                "software_versions": {
                    "basecaller_build_version": "7.11.0+5d1db4a52",
                    "bream": "8.8.3",
                    "distribution_version": "25.09.16",
                    "minknow": {"full": "25.09.5"},
                },
                "flow_cell": {
                    "flow_cell_id": "FBD69411",
                    "channel_count": 512,
                    "product_code": "FLO-MIN114",
                },
                "device": {
                    "device_type": "MINION_MK1D",
                    "device_id": "MD-100098",
                },
            }
        }))
        result = parse_report_json(report)
        assert result["basecaller_version"] == "7.11.0+5d1db4a52"
        assert result["bream_version"] == "8.8.3"
        assert result["distribution_version"] == "25.09.16"
        assert result["channel_count"] == 512
        assert result["device_type"] == "MINION_MK1D"

    def test_missing_keys_return_none(self, tmp_path):
        report = tmp_path / "report.json"
        report.write_text(json.dumps({"protocol_run_info": {}}))
        result = parse_report_json(report)
        assert result["basecaller_version"] is None
```

**Step 2: Run test — expect FAIL**

**Step 3: Implement parse_report_json**

```python
import json as _json

def parse_report_json(path: Path) -> dict:
    """Parse MinKNOW report*.json for software versions and device info.

    Navigates nested JSON safely, returning None for missing keys.
    """
    with open(path) as fh:
        data = _json.load(fh)
    pri = data.get("protocol_run_info", {})
    sw = pri.get("software_versions", {})
    fc = pri.get("flow_cell", {})
    dev = pri.get("device", {})
    mk = sw.get("minknow", {})
    return {
        "basecaller_version": sw.get("basecaller_build_version"),
        "bream_version": sw.get("bream"),
        "distribution_version": sw.get("distribution_version"),
        "minknow_version": mk.get("full") if isinstance(mk, dict) else None,
        "channel_count": fc.get("channel_count"),
        "device_type": dev.get("device_type"),
    }
```

**Step 4: Run test — expect PASS**

**Step 5: Commit**

```bash
git commit -m "feat(scan): add report.json parser for software versions"
```

---

## Task 4: Discover MinKNOW run folders

Recursively scan an experiment directory to find MinKNOW run folders, identified by their naming pattern and contents.

**Files:**
- Modify: `bin/sma_scan.py`
- Modify: `tests/test_sma_scan.py`

**Step 1: Write failing test**

```python
class TestDiscoverRuns:
    def _make_run_dir(self, base: Path, name: str,
                      flow_cell_id: str = "FBD69411",
                      instrument: str = "MD-100098",
                      run_id: str = "abc12345",
                      with_final_summary: bool = True,
                      with_pod5: bool = True,
                      with_bam: bool = True,
                      bam_demuxed: bool = False,
                      pod5_count: int = 3,
                      bam_count: int = 5) -> Path:
        """Helper to create a realistic MinKNOW run directory."""
        run_dir = base / name
        run_dir.mkdir(parents=True)
        if with_final_summary:
            fs = run_dir / f"final_summary_{flow_cell_id}_{run_id[:8]}_deadbeef.txt"
            fs.write_text(
                f"instrument={instrument}\n"
                f"flow_cell_id={flow_cell_id}\n"
                f"protocol_group_id=test_experiment\n"
                f"protocol=sequencing/seq:FLO-MIN114:SQK-NBD114-24:400\n"
                f"protocol_run_id={run_id}\n"
                f"started=2025-12-29T10:00:00-05:00\n"
                f"acquisition_stopped=2025-12-30T10:00:00-05:00\n"
                f"pod5_files_in_final_dest={pod5_count}\n"
                f"bam_files_in_final_dest={bam_count}\n"
            )
        if with_pod5:
            pod_dir = run_dir / "pod5_pass"
            pod_dir.mkdir()
            for i in range(pod5_count):
                (pod_dir / f"file_{i}.pod5").write_bytes(b"")
        if with_bam:
            bam_dir = run_dir / "bam_pass"
            bam_dir.mkdir()
            if bam_demuxed:
                for sub in ["unclassified", "V04_2", "V04_4"]:
                    (bam_dir / sub).mkdir()
                    (bam_dir / sub / "reads.bam").write_bytes(b"")
            else:
                for i in range(bam_count):
                    (bam_dir / f"reads_{i}.bam").write_bytes(b"")
        # Sample sheet
        ss = run_dir / f"sample_sheet_{flow_cell_id}_20251228_2101_{run_id[:8]}.csv"
        ss.write_text(
            "protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,"
            "flow_cell_product_code,kit,barcode,alias,type,rowNumber\n"
            f"{run_id},{instrument},{flow_cell_id},,test_experiment,FLO-MIN114,"
            "SQK-NBD114-24,barcode02,V04_2,test_sample,1\n"
        )
        return run_dir

    def test_finds_single_run(self, tmp_path):
        exp_dir = tmp_path / "experiment"
        exp_dir.mkdir()
        self._make_run_dir(
            exp_dir / "no_sample_id",
            "20251229_1055_MD-100098_FBD69411_abc12345",
        )
        runs = discover_runs(exp_dir)
        assert len(runs) == 1
        assert runs[0]["flow_cell_id"] == "FBD69411"

    def test_finds_multiple_runs(self, tmp_path):
        exp_dir = tmp_path / "experiment"
        exp_dir.mkdir()
        parent = exp_dir / "no_sample_id"
        self._make_run_dir(parent, "20251228_2219_MD-100098_FBD69411_run1",
                          run_id="run1xxxx")
        self._make_run_dir(parent, "20251229_1055_MD-100098_FBD69411_run2",
                          run_id="run2xxxx")
        runs = discover_runs(exp_dir)
        assert len(runs) == 2

    def test_detects_demuxed_bams(self, tmp_path):
        exp_dir = tmp_path / "experiment"
        exp_dir.mkdir()
        self._make_run_dir(
            exp_dir / "no_sample_id",
            "20251229_1055_MD-100098_FBD69411_abc12345",
            bam_demuxed=True,
        )
        runs = discover_runs(exp_dir)
        assert runs[0]["bam_demuxed"] is True

    def test_run_without_final_summary(self, tmp_path):
        exp_dir = tmp_path / "experiment"
        exp_dir.mkdir()
        self._make_run_dir(
            exp_dir / "no_sample_id",
            "20251228_2219_MD-100098_FBD69411_abc12345",
            with_final_summary=False,
        )
        runs = discover_runs(exp_dir)
        assert len(runs) == 1
        assert runs[0]["final_summary"] is None
```

**Step 2: Run test — expect FAIL**

**Step 3: Implement discover_runs**

```python
def discover_runs(experiment_dir: Path) -> list[dict]:
    """Recursively discover MinKNOW run folders within an experiment directory.

    A run folder is identified by containing pod5*/ or pod5_pass/ directories,
    or final_summary*.txt files. Returns list of run info dicts.
    """
    experiment_dir = Path(experiment_dir)
    runs = []

    # MinKNOW run folders follow pattern: YYYYMMDD_HHMM_DEVICE_FLOWCELL_RUNID
    for candidate in sorted(experiment_dir.rglob("*")):
        if not candidate.is_dir():
            continue
        # Check if this looks like a run folder
        has_pod5 = any(
            d.is_dir() and d.name in ("pod5", "pod5_pass")
            for d in candidate.iterdir()
        )
        has_final_summary = any(
            f.name.startswith("final_summary") and f.name.endswith(".txt")
            for f in candidate.iterdir() if f.is_file()
        )
        has_bam = any(
            d.is_dir() and d.name in ("bam_pass",)
            for d in candidate.iterdir()
        )
        if not (has_pod5 or has_final_summary):
            continue

        run_info = _extract_run_info(candidate)
        runs.append(run_info)

    return runs


def _extract_run_info(run_dir: Path) -> dict:
    """Extract metadata from a single MinKNOW run directory."""
    info: dict = {
        "run_dir": str(run_dir),
        "run_dir_name": run_dir.name,
        "flow_cell_id": None,
        "instrument": None,
        "run_id": None,
        "started": None,
        "stopped": None,
        "protocol": None,
        "protocol_run_id": None,
        "experiment_id": None,
        "kit": None,
        "flow_cell_product_code": None,
        "pod5_dir": None,
        "pod5_count": 0,
        "bam_dir": None,
        "bam_count": 0,
        "bam_demuxed": False,
        "sequencing_summary": None,
        "final_summary": None,
        "sample_sheet": None,
        "report_json": None,
        "barcode_aliases": {},
        "software": {},
    }

    # --- Parse run dir name for basic info ---
    parts = run_dir.name.split("_")
    if len(parts) >= 5:
        info["run_id"] = parts[-1]

    # --- final_summary ---
    for f in sorted(run_dir.iterdir()):
        if f.name.startswith("final_summary") and f.name.endswith(".txt"):
            info["final_summary"] = str(f.relative_to(run_dir.parent.parent))
            fs_data = parse_final_summary(f)
            info["flow_cell_id"] = fs_data.get("flow_cell_id")
            info["instrument"] = fs_data.get("instrument")
            info["started"] = fs_data.get("started")
            info["stopped"] = fs_data.get("acquisition_stopped")
            info["protocol"] = fs_data.get("protocol")
            info["protocol_run_id"] = fs_data.get("protocol_run_id")
            info["pod5_count"] = fs_data.get("pod5_files_in_final_dest", 0)
            info["bam_count"] = fs_data.get("bam_files_in_final_dest", 0)
            break

    # --- sample_sheet ---
    for f in sorted(run_dir.iterdir()):
        if f.name.startswith("sample_sheet") and f.name.endswith(".csv"):
            info["sample_sheet"] = str(f.relative_to(run_dir.parent.parent))
            ss_data = parse_minknow_sample_sheet(f)
            info["experiment_id"] = ss_data.get("experiment_id")
            info["kit"] = ss_data.get("kit")
            info["flow_cell_product_code"] = ss_data.get("flow_cell_product_code")
            info["barcode_aliases"] = ss_data.get("barcode_aliases", {})
            break

    # --- report.json ---
    for f in sorted(run_dir.iterdir()):
        if f.name.startswith("report_") and f.name.endswith(".json"):
            info["report_json"] = str(f.relative_to(run_dir.parent.parent))
            try:
                info["software"] = parse_report_json(f)
            except Exception:
                pass
            break

    # --- pod5 directory ---
    for name in ("pod5_pass", "pod5"):
        pod_dir = run_dir / name
        if pod_dir.is_dir():
            info["pod5_dir"] = name
            if info["pod5_count"] == 0:
                info["pod5_count"] = sum(1 for f in pod_dir.iterdir() if f.suffix == ".pod5")
            break

    # --- bam directory ---
    bam_dir = run_dir / "bam_pass"
    if bam_dir.is_dir():
        info["bam_dir"] = "bam_pass"
        subdirs = [d for d in bam_dir.iterdir() if d.is_dir()]
        info["bam_demuxed"] = len(subdirs) > 0
        if info["bam_count"] == 0:
            info["bam_count"] = sum(
                1 for f in bam_dir.rglob("*.bam")
            )

    # --- sequencing_summary ---
    for f in sorted(run_dir.iterdir()):
        if f.name.startswith("sequencing_summary") and f.name.endswith(".txt"):
            info["sequencing_summary"] = str(f.relative_to(run_dir.parent.parent))
            break

    # --- Extract flow_cell_id from dir name if not in final_summary ---
    if info["flow_cell_id"] is None and len(parts) >= 4:
        info["flow_cell_id"] = parts[3]
    if info["instrument"] is None and len(parts) >= 3:
        info["instrument"] = parts[2]

    return info
```

**Step 4: Run test — expect PASS**

**Step 5: Commit**

```bash
git commit -m "feat(scan): add MinKNOW run folder discovery"
```

---

## Task 5: Merge decision logic

Group runs by flow_cell_id and decide whether to merge based on protocol, kit, experiment_id, and time gap.

**Files:**
- Modify: `bin/sma_scan.py`
- Modify: `tests/test_sma_scan.py`

**Step 1: Write failing test**

```python
from sma_scan import decide_merges


class TestDecideMerges:
    def test_auto_merge_same_flow_cell(self):
        runs = [
            {"flow_cell_id": "FBD69411", "instrument": "MD-100098",
             "experiment_id": "exp1", "kit": "SQK-NBD114-24",
             "protocol": "seq:FLO:SQK:400",
             "started": "2025-12-28T22:00:00-05:00",
             "stopped": "2025-12-29T08:00:00-05:00"},
            {"flow_cell_id": "FBD69411", "instrument": "MD-100098",
             "experiment_id": "exp1", "kit": "SQK-NBD114-24",
             "protocol": "seq:FLO:SQK:400",
             "started": "2025-12-29T10:00:00-05:00",
             "stopped": "2025-12-30T10:00:00-05:00"},
        ]
        groups = decide_merges(runs, max_gap_hours=24)
        assert len(groups) == 1
        assert groups[0]["merge_decision"] == "auto"
        assert len(groups[0]["runs"]) == 2
        assert groups[0]["time_gap_hours"] == pytest.approx(2.0, abs=0.1)

    def test_flag_when_gap_exceeds_threshold(self):
        runs = [
            {"flow_cell_id": "FBD69411", "instrument": "MD-100098",
             "experiment_id": "exp1", "kit": "SQK-NBD114-24",
             "protocol": "seq:FLO:SQK:400",
             "started": "2025-12-28T10:00:00-05:00",
             "stopped": "2025-12-28T20:00:00-05:00"},
            {"flow_cell_id": "FBD69411", "instrument": "MD-100098",
             "experiment_id": "exp1", "kit": "SQK-NBD114-24",
             "protocol": "seq:FLO:SQK:400",
             "started": "2025-12-31T10:00:00-05:00",
             "stopped": None},
        ]
        groups = decide_merges(runs, max_gap_hours=24)
        assert groups[0]["merge_decision"] == "review"
        assert "gap" in groups[0]["warnings"][0].lower()

    def test_separate_flow_cells(self):
        runs = [
            {"flow_cell_id": "FBD69411", "instrument": "MD-100098",
             "experiment_id": "exp1", "kit": "K1", "protocol": "P1",
             "started": "2025-12-28T10:00:00", "stopped": None},
            {"flow_cell_id": "FBD70599", "instrument": "MD-101527",
             "experiment_id": "exp1", "kit": "K1", "protocol": "P1",
             "started": "2025-12-28T10:00:00", "stopped": None},
        ]
        groups = decide_merges(runs, max_gap_hours=24)
        assert len(groups) == 2

    def test_single_run_no_merge(self):
        runs = [
            {"flow_cell_id": "FBD66244", "instrument": "MD-101527",
             "experiment_id": "exp1", "kit": "K1", "protocol": "P1",
             "started": "2025-12-30T17:10:00", "stopped": None},
        ]
        groups = decide_merges(runs, max_gap_hours=24)
        assert len(groups) == 1
        assert groups[0]["merge_decision"] == "single"
        assert len(groups[0]["runs"]) == 1
```

**Step 2: Run test — expect FAIL**

**Step 3: Implement decide_merges**

```python
from datetime import datetime, timezone
from collections import defaultdict


def _parse_timestamp(ts: str | None) -> datetime | None:
    """Parse an ISO timestamp string, handling timezone offsets."""
    if not ts:
        return None
    # Handle various formats MinKNOW uses
    for fmt in (
        "%Y-%m-%dT%H:%M:%S.%f%z",
        "%Y-%m-%dT%H:%M:%S%z",
        "%Y-%m-%dT%H:%M:%S.%f",
        "%Y-%m-%dT%H:%M:%S",
    ):
        try:
            return datetime.strptime(ts, fmt)
        except ValueError:
            continue
    return None


def decide_merges(runs: list[dict], max_gap_hours: float = 24.0) -> list[dict]:
    """Group runs by flow_cell_id and decide merge eligibility.

    Returns list of flow cell groups, each with merge_decision and warnings.
    """
    by_fc: dict[str, list[dict]] = defaultdict(list)
    for run in runs:
        fc_id = run.get("flow_cell_id") or "unknown"
        by_fc[fc_id].append(run)

    groups = []
    for fc_id, fc_runs in sorted(by_fc.items()):
        # Sort by start time
        fc_runs.sort(key=lambda r: r.get("started") or "")

        group: dict = {
            "flow_cell_id": fc_id,
            "instrument": fc_runs[0].get("instrument"),
            "runs": fc_runs,
            "time_gap_hours": None,
            "merge_decision": "single",
            "warnings": [],
        }

        if len(fc_runs) == 1:
            groups.append(group)
            continue

        # Compute time gap between consecutive runs
        max_observed_gap = 0.0
        warnings = []
        for i in range(1, len(fc_runs)):
            prev_stop = _parse_timestamp(fc_runs[i - 1].get("stopped"))
            curr_start = _parse_timestamp(fc_runs[i].get("started"))
            if prev_stop and curr_start:
                gap_hours = (curr_start - prev_stop).total_seconds() / 3600
                max_observed_gap = max(max_observed_gap, gap_hours)

            # Check consistency
            if fc_runs[i].get("experiment_id") != fc_runs[i - 1].get("experiment_id"):
                warnings.append(f"Different experiment_id between runs {i-1} and {i}")
            if fc_runs[i].get("kit") != fc_runs[i - 1].get("kit"):
                warnings.append(f"Different kit between runs {i-1} and {i}")

        group["time_gap_hours"] = round(max_observed_gap, 1)

        if warnings:
            group["merge_decision"] = "review"
            group["warnings"] = warnings
        elif max_observed_gap > max_gap_hours:
            group["merge_decision"] = "review"
            group["warnings"] = [f"Time gap {max_observed_gap:.1f}h exceeds {max_gap_hours}h threshold"]
        else:
            group["merge_decision"] = "auto"

        groups.append(group)

    return groups
```

**Step 4: Run test — expect PASS**

**Step 5: Commit**

```bash
git commit -m "feat(scan): add merge decision logic for multi-run experiments"
```

---

## Task 6: Build manifest and CLI

Assemble the full manifest JSON from discovered runs and merge decisions. Add argparse CLI.

**Files:**
- Modify: `bin/sma_scan.py`
- Modify: `tests/test_sma_scan.py`

**Step 1: Write failing test for build_manifest**

```python
class TestBuildManifest:
    def test_manifest_structure(self, tmp_path):
        runs = [
            {
                "run_dir": str(tmp_path / "run1"),
                "flow_cell_id": "FBD69411",
                "instrument": "MD-100098",
                "experiment_id": "my_exp",
                "kit": "SQK-NBD114-24",
                "flow_cell_product_code": "FLO-MIN114",
                "pod5_dir": "pod5_pass",
                "pod5_count": 12,
                "bam_dir": "bam_pass",
                "bam_count": 5,
                "bam_demuxed": True,
                "started": "2025-12-28T22:00:00",
                "stopped": None,
                "protocol": "seq:FLO:SQK:400",
                "sequencing_summary": None,
                "final_summary": None,
                "sample_sheet": None,
                "report_json": None,
                "barcode_aliases": {"barcode02": "V04_2"},
                "software": {"basecaller_version": "7.11.0"},
            },
        ]
        manifest = build_manifest(tmp_path, runs, max_gap_hours=24)
        assert manifest["experiment_id"] == "my_exp"
        assert manifest["kit"] == "SQK-NBD114-24"
        assert len(manifest["flow_cells"]) == 1
        assert manifest["action"]["needs_rebasecall"] is True  # demuxed BAMs

    def test_needs_rebasecall_when_no_bams(self, tmp_path):
        runs = [
            {
                "run_dir": str(tmp_path / "run1"),
                "flow_cell_id": "FBD44097",
                "instrument": "MD-101527",
                "experiment_id": "exp",
                "kit": "K1",
                "flow_cell_product_code": "FLO-MIN114",
                "pod5_dir": "pod5_pass",
                "pod5_count": 51,
                "bam_dir": None,
                "bam_count": 0,
                "bam_demuxed": False,
                "started": "2026-02-09T18:17:00",
                "stopped": None,
                "protocol": "P1",
                "sequencing_summary": "seq_summary.txt",
                "final_summary": "final_summary.txt",
                "sample_sheet": None,
                "report_json": None,
                "barcode_aliases": {},
                "software": {},
            },
        ]
        manifest = build_manifest(tmp_path, runs, max_gap_hours=24)
        assert manifest["action"]["needs_rebasecall"] is True
```

**Step 2: Run test — expect FAIL**

**Step 3: Implement build_manifest + CLI**

```python
def build_manifest(experiment_dir: Path, runs: list[dict],
                   max_gap_hours: float = 24.0) -> dict:
    """Build the complete manifest JSON from discovered runs."""
    groups = decide_merges(runs, max_gap_hours)

    # Aggregate experiment-level info from first run
    first_run = runs[0] if runs else {}
    barcode_aliases = {}
    seq_summary_paths = []
    pod5_sources = []
    total_pod5 = 0
    total_bam = 0
    any_demuxed = False

    for run in runs:
        barcode_aliases.update(run.get("barcode_aliases", {}))
        if run.get("sequencing_summary"):
            seq_summary_paths.append(run["sequencing_summary"])
        if run.get("pod5_dir"):
            pod5_sources.append(str(Path(run["run_dir"]) / run["pod5_dir"]))
        total_pod5 += run.get("pod5_count", 0)
        total_bam += run.get("bam_count", 0)
        if run.get("bam_demuxed"):
            any_demuxed = True

    needs_rebasecall = any_demuxed or total_bam == 0

    manifest = {
        "experiment_id": first_run.get("experiment_id"),
        "experiment_dir": str(experiment_dir),
        "scan_timestamp": datetime.now(timezone.utc).isoformat(),
        "kit": first_run.get("kit"),
        "flow_cell_product": first_run.get("flow_cell_product_code"),
        "software": first_run.get("software", {}),
        "barcode_aliases": barcode_aliases,
        "flow_cells": groups,
        "action": {
            "needs_rebasecall": needs_rebasecall,
            "pod5_sources": pod5_sources,
            "sequencing_summary_paths": seq_summary_paths,
            "total_pod5_files": total_pod5,
            "total_bam_files": total_bam,
        },
    }
    return manifest


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Scan SMA-seq experiment directory and generate manifest"
    )
    parser.add_argument("experiment_dir", type=Path,
                        help="Path to experiment directory")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Output manifest JSON path")
    parser.add_argument("--max-gap", type=float, default=24.0,
                        help="Max time gap (hours) for auto-merge [%(default)s]")
    parser.add_argument("--interactive", action="store_true",
                        help="Always prompt before merging")
    args = parser.parse_args()

    exp_dir = args.experiment_dir.resolve()
    if not exp_dir.is_dir():
        sys.exit(f"[sma_scan] Error: {exp_dir} is not a directory")

    print(f"[sma_scan] Scanning {exp_dir}")
    runs = discover_runs(exp_dir)
    print(f"[sma_scan] Found {len(runs)} run(s)")

    for run in runs:
        fc = run.get("flow_cell_id", "?")
        pod5 = run.get("pod5_count", 0)
        bam = run.get("bam_count", 0)
        demux = " (demuxed)" if run.get("bam_demuxed") else ""
        print(f"  - {run['run_dir_name']}: FC={fc} pod5={pod5} bam={bam}{demux}")

    manifest = build_manifest(exp_dir, runs, args.max_gap)

    # Print merge decisions
    for group in manifest["flow_cells"]:
        fc = group["flow_cell_id"]
        n_runs = len(group["runs"])
        decision = group["merge_decision"]
        gap = group.get("time_gap_hours")
        if n_runs > 1:
            print(f"[sma_scan] Flow cell {fc}: {n_runs} runs, "
                  f"merge={decision}, gap={gap}h")
            for w in group.get("warnings", []):
                print(f"  WARNING: {w}")

    if manifest["action"]["needs_rebasecall"]:
        print("[sma_scan] BAMs need re-basecalling (demuxed or missing)")

    # Write manifest
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as fh:
        _json.dump(manifest, fh, indent=2)
    print(f"[sma_scan] Manifest written to {args.output}")


if __name__ == "__main__":
    import sys
    main()
```

**Step 4: Run test — expect PASS**

**Step 5: Write integration test (subprocess)**

```python
class TestScanCLI:
    def test_cli_produces_manifest(self, tmp_path):
        """Integration test: run sma_scan.py via subprocess."""
        # Create a fake experiment directory
        exp_dir = tmp_path / "experiment"
        exp_dir.mkdir()
        run_dir = exp_dir / "no_sample_id" / "20251229_1055_MD-100098_FBD69411_abc12345"
        run_dir.mkdir(parents=True)

        # Final summary
        (run_dir / "final_summary_FBD69411_abc_123.txt").write_text(
            "instrument=MD-100098\nflow_cell_id=FBD69411\n"
            "protocol_group_id=test_exp\n"
            "protocol=seq:FLO:SQK:400\nprotocol_run_id=abc-def\n"
            "started=2025-12-29T10:00:00\n"
            "pod5_files_in_final_dest=5\nbam_files_in_final_dest=10\n"
        )
        # Pod5 dir
        pod5_dir = run_dir / "pod5_pass"
        pod5_dir.mkdir()
        for i in range(5):
            (pod5_dir / f"f_{i}.pod5").write_bytes(b"")
        # Sample sheet
        (run_dir / "sample_sheet_FBD69411_20251228_abc.csv").write_text(
            "protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,"
            "flow_cell_product_code,kit,barcode,alias,type,rowNumber\n"
            "abc-def,MD-100098,FBD69411,,test_exp,FLO-MIN114,"
            "SQK-NBD114-24,barcode02,V04_2,test_sample,1\n"
        )

        out = tmp_path / "manifest.json"
        result = subprocess.run(
            [sys.executable, "bin/sma_scan.py", str(exp_dir), "-o", str(out)],
            capture_output=True, text=True,
            cwd="/tmp/ont-sma-seq",
        )
        assert result.returncode == 0, result.stderr
        assert out.exists()
        manifest = _json.loads(out.read_text())
        assert manifest["experiment_id"] == "test_exp"
        assert manifest["kit"] == "SQK-NBD114-24"
        assert len(manifest["flow_cells"]) == 1
```

**Step 6: Run all scan tests**

```bash
cd /tmp/ont-sma-seq && python -m pytest tests/test_sma_scan.py -v
```

**Step 7: Commit**

```bash
git commit -m "feat(scan): add manifest builder and CLI entry point"
```

---

## Task 7: Run scanner on real experiments

Validate `sma_scan.py` against the actual experiment directories on D: drive.

**Files:** None (validation only)

**Step 1: Run on no_trim (2 runs, merge candidate)**

```bash
cd /tmp/ont-sma-seq && python bin/sma_scan.py \
  "/mnt/d/12282025_IF_DoubleBC_SMA_seq_no_trim" \
  -o /tmp/sma-manifests/no_trim_manifest.json
```

Expected: 2 runs found, FBD69411, auto-merge, gap ~12h.

**Step 2: Run on one_nick (single run)**

```bash
python bin/sma_scan.py \
  "/mnt/d/12302025_IF_DoubleBC_SMA_seq_one_nick" \
  -o /tmp/sma-manifests/one_nick_manifest.json
```

Expected: 1 run found, FBD66244, single.

**Step 3: Run on extended (single run under GG_Part5_Odd)**

```bash
python bin/sma_scan.py \
  "/mnt/d/20260129_IF_GG_Part5_Odd/20260129_IF_SMAseq_Ggallphos_extended" \
  -o /tmp/sma-manifests/extended_manifest.json
```

**Step 4: Run on level1 (no BAMs)**

```bash
python bin/sma_scan.py \
  "/mnt/d/20260209_IF_GG_Part5_Level1" \
  -o /tmp/sma-manifests/level1_manifest.json
```

Expected: needs_rebasecall=true because bam_count=0.

**Step 5: Inspect manifests and fix any issues**

**Step 6: Commit any fixes**

---

## Task 8: `sma_basecall.py` — BAM merge from existing files

Start with the simpler `--from-bam` mode that merges existing demuxed BAM files.

**Files:**
- Create: `bin/sma_basecall.py`
- Create: `tests/test_sma_basecall.py`

**Step 1: Write failing test for collect_bam_files**

```python
# tests/test_sma_basecall.py
"""Tests for sma_basecall.py BAM merging and basecalling."""

from pathlib import Path
import json
import subprocess
import sys
import pytest
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))

from sma_basecall import collect_bam_files


class TestCollectBamFiles:
    def test_collects_from_demuxed_dirs(self, tmp_path):
        """Collect BAMs from demuxed subdirectories."""
        bam_dir = tmp_path / "bam_pass"
        bam_dir.mkdir()
        for sub in ["unclassified", "V04_2", "V04_4"]:
            d = bam_dir / sub
            d.mkdir()
            # Create minimal valid BAM
            header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
            with pysam.AlignmentFile(str(d / "reads.bam"), "wb", header=header) as f:
                pass
        bams = collect_bam_files([str(bam_dir)])
        assert len(bams) == 3

    def test_collects_from_flat_dir(self, tmp_path):
        """Collect BAMs from flat directory."""
        bam_dir = tmp_path / "bam_pass"
        bam_dir.mkdir()
        for i in range(3):
            header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
            with pysam.AlignmentFile(str(bam_dir / f"reads_{i}.bam"), "wb", header=header) as f:
                pass
        bams = collect_bam_files([str(bam_dir)])
        assert len(bams) == 3
```

**Step 2: Run test — expect FAIL**

**Step 3: Implement collect_bam_files and merge logic**

```python
#!/usr/bin/env python3
# sma_basecall.py — Basecall from pod5 or merge existing BAMs

import argparse
import json
import subprocess
import sys
from pathlib import Path
from datetime import datetime, timezone


def collect_bam_files(bam_dirs: list[str]) -> list[Path]:
    """Collect all .bam files from directories, including subdirectories."""
    bams = []
    for d in bam_dirs:
        dp = Path(d)
        if not dp.is_dir():
            continue
        for f in sorted(dp.rglob("*.bam")):
            if f.is_file() and f.stat().st_size > 0:
                bams.append(f)
    return bams


def merge_bams(bam_files: list[Path], output: Path) -> None:
    """Merge multiple BAM files into one sorted, indexed BAM."""
    if not bam_files:
        raise ValueError("No BAM files to merge")

    if len(bam_files) == 1:
        # Sort single file
        subprocess.run(
            ["samtools", "sort", "-o", str(output), str(bam_files[0])],
            check=True,
        )
    else:
        # Merge then sort
        merged_tmp = output.with_suffix(".unsorted.bam")
        subprocess.run(
            ["samtools", "merge", "-f", str(merged_tmp)] + [str(b) for b in bam_files],
            check=True,
        )
        subprocess.run(
            ["samtools", "sort", "-o", str(output), str(merged_tmp)],
            check=True,
        )
        merged_tmp.unlink()

    # Index
    subprocess.run(["samtools", "index", str(output)], check=True)
```

**Step 4: Run test — expect PASS**

**Step 5: Commit**

```bash
git commit -m "feat(basecall): add BAM collection and merge utilities"
```

---

## Task 9: `sma_basecall.py` — Dorado basecalling from pod5

Add the `--from-pod5` mode that runs Dorado basecaller.

**Files:**
- Modify: `bin/sma_basecall.py`
- Modify: `tests/test_sma_basecall.py`

**Step 1: Write failing test for build_dorado_command**

```python
from sma_basecall import build_dorado_command


class TestBuildDoradoCommand:
    def test_basic_command(self):
        cmd = build_dorado_command(
            model="sup",
            pod5_dir="/data/pod5",
            dorado_bin="dorado",
        )
        assert cmd == [
            "dorado", "basecaller", "sup", "/data/pod5",
            "--no-trim", "--emit-moves",
        ]

    def test_custom_model(self):
        cmd = build_dorado_command(
            model="dna_r10.4.1_e8.2_400bps_hac@v5.0.0",
            pod5_dir="/data/pod5",
            dorado_bin="/home/user/dorado/bin/dorado",
        )
        assert cmd[0] == "/home/user/dorado/bin/dorado"
        assert "dna_r10.4.1_e8.2_400bps_hac@v5.0.0" in cmd

    def test_subsample_adds_read_ids(self, tmp_path):
        ids_file = tmp_path / "ids.txt"
        ids_file.write_text("read1\nread2\n")
        cmd = build_dorado_command(
            model="sup", pod5_dir="/data/pod5",
            dorado_bin="dorado",
            read_ids_file=str(ids_file),
        )
        assert "--read-ids" in cmd
        assert str(ids_file) in cmd
```

**Step 2: Run test — expect FAIL**

**Step 3: Implement build_dorado_command**

```python
DORADO_BIN = Path.home() / "dorado" / "bin" / "dorado"


def build_dorado_command(
    model: str,
    pod5_dir: str,
    dorado_bin: str = str(DORADO_BIN),
    read_ids_file: str | None = None,
) -> list[str]:
    """Build Dorado basecaller command line."""
    cmd = [
        dorado_bin, "basecaller", model, pod5_dir,
        "--no-trim", "--emit-moves",
    ]
    if read_ids_file:
        cmd.extend(["--read-ids", read_ids_file])
    return cmd


def run_dorado(
    model: str,
    pod5_dir: str,
    output_bam: Path,
    dorado_bin: str = str(DORADO_BIN),
    read_ids_file: str | None = None,
) -> None:
    """Run Dorado basecaller and write output BAM."""
    cmd = build_dorado_command(model, pod5_dir, dorado_bin, read_ids_file)
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    with open(output_bam, "wb") as out_fh:
        proc = subprocess.run(cmd, stdout=out_fh, check=True)
```

**Step 4: Run test — expect PASS**

**Step 5: Commit**

```bash
git commit -m "feat(basecall): add Dorado command builder and runner"
```

---

## Task 10: `sma_basecall.py` — Pod5 subsampling

Add `--subsample N` flag that extracts a random subset of reads from pod5 files.

**Files:**
- Modify: `bin/sma_basecall.py`
- Modify: `tests/test_sma_basecall.py`

**Step 1: Write failing test**

```python
from sma_basecall import subsample_read_ids


class TestSubsampleReadIds:
    def test_subsamples_correct_count(self):
        all_ids = [f"read_{i}" for i in range(1000)]
        sampled = subsample_read_ids(all_ids, n=100, seed=42)
        assert len(sampled) == 100
        assert len(set(sampled)) == 100  # unique

    def test_returns_all_when_n_exceeds_total(self):
        all_ids = [f"read_{i}" for i in range(50)]
        sampled = subsample_read_ids(all_ids, n=100, seed=42)
        assert len(sampled) == 50

    def test_deterministic_with_seed(self):
        all_ids = [f"read_{i}" for i in range(1000)]
        s1 = subsample_read_ids(all_ids, n=100, seed=42)
        s2 = subsample_read_ids(all_ids, n=100, seed=42)
        assert s1 == s2
```

**Step 2: Run test — expect FAIL**

**Step 3: Implement subsample_read_ids**

```python
import random


def subsample_read_ids(all_ids: list[str], n: int, seed: int = 42) -> list[str]:
    """Randomly subsample n read IDs from a list."""
    if n >= len(all_ids):
        return list(all_ids)
    rng = random.Random(seed)
    return rng.sample(all_ids, n)


def get_pod5_read_ids(pod5_dir: str) -> list[str]:
    """Get all read IDs from pod5 files in a directory."""
    import pod5
    ids = []
    pod5_path = Path(pod5_dir)
    for pod5_file in sorted(pod5_path.glob("*.pod5")):
        with pod5.Reader(pod5_file) as reader:
            for read in reader.reads():
                ids.append(str(read.read_id))
    return ids
```

**Step 4: Run test — expect PASS**

**Step 5: Commit**

```bash
git commit -m "feat(basecall): add pod5 read ID extraction and subsampling"
```

---

## Task 11: `sma_basecall.py` — Sequencing summary merge

Merge sequencing_summary.txt files from multiple runs, preserving run_id.

**Files:**
- Modify: `bin/sma_basecall.py`
- Modify: `tests/test_sma_basecall.py`

**Step 1: Write failing test**

```python
from sma_basecall import merge_sequencing_summaries


class TestMergeSequencingSummaries:
    def test_merges_two_files(self, tmp_path):
        header = "read_id\trun_id\tend_reason\tduration\n"
        (tmp_path / "s1.txt").write_text(header + "r1\trun1\tsignal_positive\t1.5\n")
        (tmp_path / "s2.txt").write_text(header + "r2\trun2\tunblock\t0.3\n")
        out = tmp_path / "merged.tsv"
        merge_sequencing_summaries(
            [str(tmp_path / "s1.txt"), str(tmp_path / "s2.txt")],
            str(out),
        )
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 3  # header + 2 data lines
        assert "r1" in lines[1]
        assert "r2" in lines[2]

    def test_single_file_copies(self, tmp_path):
        header = "read_id\trun_id\tend_reason\n"
        (tmp_path / "s1.txt").write_text(header + "r1\trun1\tsignal_positive\n")
        out = tmp_path / "merged.tsv"
        merge_sequencing_summaries([str(tmp_path / "s1.txt")], str(out))
        assert out.exists()
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 2
```

**Step 2: Run test — expect FAIL**

**Step 3: Implement merge_sequencing_summaries**

```python
def merge_sequencing_summaries(summary_paths: list[str], output: str) -> None:
    """Merge multiple sequencing_summary.txt files into one TSV."""
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    header_written = False
    with open(out_path, "w") as out_fh:
        for sp in summary_paths:
            with open(sp) as in_fh:
                for i, line in enumerate(in_fh):
                    if i == 0:
                        if not header_written:
                            out_fh.write(line)
                            header_written = True
                        continue
                    out_fh.write(line)
```

**Step 4: Run test — expect PASS**

**Step 5: Commit**

```bash
git commit -m "feat(basecall): add sequencing summary merge"
```

---

## Task 12: `sma_basecall.py` — CLI and provenance

Wire up the CLI entry point with argparse and provenance JSON output.

**Files:**
- Modify: `bin/sma_basecall.py`
- Modify: `tests/test_sma_basecall.py`

**Step 1: Write integration test**

```python
class TestBasecallCLI:
    def test_from_bam_merges_and_indexes(self, tmp_path):
        """Integration: --from-bam collects, merges, sorts, indexes."""
        # Create a manifest pointing to our test BAMs
        bam_dir = tmp_path / "run1" / "bam_pass"
        bam_dir.mkdir(parents=True)

        # Create 2 small BAMs with reads
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "unsorted"},
            "SQ": [{"SN": "unmapped", "LN": 1000}],
        })
        for i in range(2):
            with pysam.AlignmentFile(str(bam_dir / f"reads_{i}.bam"), "wb", header=header) as f:
                seg = pysam.AlignedSegment(header)
                seg.query_name = f"read_{i}"
                seg.query_sequence = "ACGT" * 25
                seg.query_qualities = pysam.qualitystring_to_array("I" * 100)
                seg.flag = 4  # unmapped
                f.write(seg)

        manifest = {
            "experiment_id": "test_exp",
            "experiment_dir": str(tmp_path),
            "flow_cells": [
                {"runs": [{"run_dir": str(tmp_path / "run1"), "bam_dir": "bam_pass"}]}
            ],
            "action": {
                "pod5_sources": [],
                "sequencing_summary_paths": [],
            },
        }
        manifest_path = tmp_path / "manifest.json"
        manifest_path.write_text(json.dumps(manifest))

        out_dir = tmp_path / "output"
        result = subprocess.run(
            [sys.executable, "bin/sma_basecall.py",
             str(manifest_path), "--from-bam", "-o", str(out_dir)],
            capture_output=True, text=True,
            cwd="/tmp/ont-sma-seq",
        )
        assert result.returncode == 0, result.stderr
        merged_bam = out_dir / "test_exp_merged.bam"
        assert merged_bam.exists()
        assert (out_dir / "test_exp_merged.bam.bai").exists()

        # Verify read count
        with pysam.AlignmentFile(str(merged_bam)) as f:
            reads = list(f)
        assert len(reads) == 2
```

**Step 2: Run test — expect FAIL**

**Step 3: Implement main() CLI**

Add `main()` function to `sma_basecall.py` with argparse:
- `manifest` positional arg
- `--from-pod5` / `--from-bam` mutually exclusive
- `-o` / `--output` required
- `--model` (default "sup")
- `--subsample N`
- `--seed` (default 42)
- `--dorado-bin`

The `--from-bam` path: read manifest → collect BAM dirs → collect BAMs → merge → index → write provenance.
The `--from-pod5` path: read manifest → get pod5 dirs → optionally subsample → run dorado per dir → merge outputs → index → write provenance.

**Step 4: Run test — expect PASS**

**Step 5: Commit**

```bash
git commit -m "feat(basecall): add CLI with --from-bam and --from-pod5 modes"
```

---

## Task 13: Validate on real experiment (from-bam mode)

Test the `--from-bam` merge on the one_nick experiment.

**Step 1: Run scan**

```bash
cd /tmp/ont-sma-seq
python bin/sma_scan.py "/mnt/d/12302025_IF_DoubleBC_SMA_seq_one_nick" \
  -o /tmp/sma-manifests/one_nick_manifest.json
```

**Step 2: Run from-bam merge**

```bash
python bin/sma_basecall.py /tmp/sma-manifests/one_nick_manifest.json \
  --from-bam -o /tmp/sma-one-nick/
```

**Step 3: Verify output**

```bash
samtools flagstat /tmp/sma-one-nick/one_nick_merged.bam
samtools view -c /tmp/sma-one-nick/one_nick_merged.bam
```

**Step 4: Commit any fixes**

---

## Task 14: `sma_init.py` — Database creation from merged BAM

Create the SMA-seq database from a merged BAM file, without classification.

**Files:**
- Create: `bin/sma_init.py`
- Create: `tests/test_sma_init.py`

**Step 1: Write failing test**

```python
# tests/test_sma_init.py
"""Tests for sma_init.py database initialization and QC."""

from pathlib import Path
import json
import sqlite3
import subprocess
import sys
import pytest
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))

from sma_init import ingest_bam_to_db


class TestIngestBamToDb:
    def _make_bam(self, path: Path, reads: list[tuple[str, str]]) -> None:
        """Create a BAM file with given (read_id, sequence) pairs."""
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "unsorted"},
            "PG": [{"ID": "basecaller", "PN": "dorado",
                     "VN": "0.8.4", "CL": "dorado basecaller sup --no-trim"}],
        })
        with pysam.AlignmentFile(str(path), "wb", header=header) as f:
            for read_id, seq in reads:
                seg = pysam.AlignedSegment(header)
                seg.query_name = read_id
                seg.query_sequence = seq
                seg.query_qualities = pysam.qualitystring_to_array(
                    "I" * len(seq)
                )
                seg.flag = 4
                f.write(seg)

    def test_creates_database_with_reads(self, tmp_path):
        bam = tmp_path / "merged.bam"
        self._make_bam(bam, [
            ("read_1", "ACGTACGT" * 25),
            ("read_2", "TGCATGCA" * 30),
        ])
        db_path = tmp_path / "test.db"
        ingest_bam_to_db(
            bam_path=bam,
            db_path=db_path,
            exp_id="test_exp",
            summary_map={},
        )
        conn = sqlite3.connect(str(db_path))
        count = conn.execute("SELECT COUNT(*) FROM Reads").fetchone()[0]
        assert count == 2
        # Check fields
        row = conn.execute(
            "SELECT read_id, readlen, q_bc FROM Reads WHERE read_id='read_1'"
        ).fetchone()
        assert row[0] == "read_1"
        assert row[1] == 200  # 8*25
        assert row[2] > 0  # some Q score
        conn.close()

    def test_attaches_end_reason_from_summary(self, tmp_path):
        bam = tmp_path / "merged.bam"
        self._make_bam(bam, [("read_1", "ACGT" * 50)])
        db_path = tmp_path / "test.db"
        summary_map = {"read_1": {"end_reason": "signal_positive", "duration": 1.5}}
        ingest_bam_to_db(bam, db_path, "exp", summary_map)
        conn = sqlite3.connect(str(db_path))
        row = conn.execute(
            "SELECT ER, signal_duration_s FROM Reads WHERE read_id='read_1'"
        ).fetchone()
        assert row[0] == "signal_positive"
        assert row[1] == 1.5
        conn.close()
```

**Step 2: Run test — expect FAIL**

**Step 3: Implement sma_init.py with ingest_bam_to_db**

This function:
1. Creates DB using mkdb.py schema (or calls mkdb directly)
2. Iterates BAM reads
3. Computes q_bc from quality scores
4. Looks up end_reason/duration from summary map
5. Inserts into Reads table (no tgt_id, no barcode fields — those come from ingest.py later)

**Step 4: Run test — expect PASS**

**Step 5: Commit**

```bash
git commit -m "feat(init): add BAM-to-database ingestion without classification"
```

---

## Task 15: `sma_init.py` — QC metrics computation

Compute pre-classification QC metrics from the database.

**Files:**
- Modify: `bin/sma_init.py`
- Modify: `tests/test_sma_init.py`

**Step 1: Write failing test for compute_qc_metrics**

```python
from sma_init import compute_qc_metrics


class TestComputeQcMetrics:
    def test_computes_length_distribution(self, tmp_path):
        db_path = tmp_path / "test.db"
        # ... (create DB with known reads)
        metrics = compute_qc_metrics(db_path)
        assert "read_lengths" in metrics
        assert len(metrics["read_lengths"]) > 0

    def test_computes_end_reason_counts(self, tmp_path):
        # ... (create DB with known end reasons)
        metrics = compute_qc_metrics(db_path)
        assert "end_reason_counts" in metrics
        assert metrics["end_reason_counts"]["signal_positive"] > 0

    def test_computes_qscore_stats(self, tmp_path):
        # ...
        metrics = compute_qc_metrics(db_path)
        assert "mean_qbc" in metrics
        assert metrics["mean_qbc"] > 0
```

**Step 2-5: Implement, test, commit**

---

## Task 16: `sma_init.py` — QC HTML report template

Create the interactive HTML QC report with 6 tabs.

**Files:**
- Create: `bin/sma_init_template.py`
- Modify: `tests/test_sma_init.py`

This follows the exact same pattern as `report_template.py` — a single `generate_qc_html(metrics, manifest)` function that returns an HTML string with embedded CSS/JS.

**Tabs:** Overview, Read Lengths, Signal & Duration, End Reasons, Quality, Channel Activity.

**Step 1: Write test that HTML contains all tab IDs**

**Step 2-5: Implement, test, commit**

---

## Task 17: `sma_init.py` — CLI and integration test

Wire up CLI and test end-to-end.

**Step 1: Integration test**

```python
class TestInitCLI:
    def test_full_pipeline(self, tmp_path):
        """End-to-end: BAM → DB + QC report."""
        # Create test BAM, summary, manifest
        # Run sma_init.py via subprocess
        # Assert DB exists, QC report exists, report > 10KB
```

**Step 2-5: Implement, test, commit**

---

## Task 18: End-to-end validation on real data

Run the full scan → basecall → init pipeline on a real experiment.

**Step 1: Run on one_nick experiment (simplest: single run, 437 pod5)**

```bash
# Scan
python bin/sma_scan.py "/mnt/d/12302025_IF_DoubleBC_SMA_seq_one_nick" \
  -o /tmp/sma-one-nick/manifest.json

# Merge existing BAMs (quick test)
python bin/sma_basecall.py /tmp/sma-one-nick/manifest.json \
  --from-bam -o /tmp/sma-one-nick/

# Init DB + QC report
python bin/sma_init.py \
  --bam /tmp/sma-one-nick/one_nick_merged.bam \
  --sequencing-summary /tmp/sma-one-nick/one_nick_summary.tsv \
  --manifest /tmp/sma-one-nick/manifest.json \
  -o /tmp/sma-one-nick/
```

**Step 2: Verify QC report in browser**

**Step 3: Run on no_trim (merge test: 2 runs)**

**Step 4: Commit any fixes**

---

## Verification

```bash
# All tests pass
cd /tmp/ont-sma-seq && python -m pytest tests/ -v

# Manifests generated for all experiments
ls /tmp/sma-manifests/*.json

# QC report viewable in browser
# Copy to Google Drive for verification
```
