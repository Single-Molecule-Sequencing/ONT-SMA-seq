# SMA-seq Config Visualizer Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build an interactive browser-based tool for inspecting, creating, and managing SMA-seq experiment configurations with static HTML export.

**Architecture:** FastAPI backend serving Jinja2 templates with HTMX for interactivity and D3.js for construct diagrams. ONT-native file formats (MinKNOW CSV, Dorado TOML/FASTA) as source of truth, with a minimal SMA-seq overlay TOML for pipeline-specific config.

**Tech Stack:** Python 3, FastAPI, uvicorn, Jinja2, HTMX (CDN), D3.js (CDN), Pico CSS (CDN), tomli/tomli_w, edlib

**Design Doc:** `docs/plans/2026-02-14-sma-seq-config-visualizer-design.md`

---

### Task 1: Project Scaffolding & Dependencies

**Files:**
- Create: `bin/viz/__init__.py`
- Create: `bin/viz/app.py`
- Create: `bin/viz/templates/.gitkeep`
- Create: `bin/viz/static/.gitkeep`
- Modify: `env/env.yml`

**Step 1: Update conda env with new dependencies**

Add to `env/env.yml`:
```yaml
dependencies:
  - edlib
  - matplotlib
  - pandas
  - pod5
  - pysam
  - seaborn
  - pip:
    - fastapi
    - uvicorn[standard]
    - jinja2
    - python-multipart
    - tomli
    - tomli_w
```

**Step 2: Create package structure**

Create `bin/viz/__init__.py`:
```python
"""SMA-seq experiment configuration visualizer."""
```

**Step 3: Create minimal FastAPI app**

Create `bin/viz/app.py`:
```python
"""SMA-seq Config Visualizer - FastAPI application."""

from __future__ import annotations

import sys
import webbrowser
from pathlib import Path

import uvicorn
from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

APP_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = APP_DIR / "templates"
STATIC_DIR = APP_DIR / "static"

app = FastAPI(title="SMA-seq Config Visualizer")
templates = Jinja2Templates(directory=str(TEMPLATES_DIR))
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

# Experiment directory - set at startup
experiment_dir: Path | None = None


@app.get("/")
async def dashboard(request: Request):
    return templates.TemplateResponse("dashboard.html", {
        "request": request,
        "experiment_dir": str(experiment_dir) if experiment_dir else "Not set",
    })


def main():
    global experiment_dir

    if len(sys.argv) < 2:
        print("Usage: python -m bin.viz.app <experiment_dir>")
        print("  experiment_dir: path to experiment directory containing config files")
        sys.exit(1)

    experiment_dir = Path(sys.argv[1]).resolve()
    if not experiment_dir.is_dir():
        print(f"Error: {experiment_dir} is not a directory")
        sys.exit(1)

    host = "127.0.0.1"
    port = 8050

    print(f"[viz] Serving experiment: {experiment_dir}")
    print(f"[viz] Open http://{host}:{port}")
    webbrowser.open(f"http://{host}:{port}")
    uvicorn.run(app, host=host, port=port, log_level="warning")


if __name__ == "__main__":
    main()
```

**Step 4: Create minimal dashboard template**

Create `bin/viz/templates/dashboard.html`:
```html
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>SMA-seq Config Visualizer</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@picocss/pico@2/css/pico.min.css">
    <script src="https://unpkg.com/htmx.org@2.0.4"></script>
</head>
<body>
  <main class="container">
    <h1>SMA-seq Config Visualizer</h1>
    <p>Experiment: {{ experiment_dir }}</p>
  </main>
</body>
</html>
```

**Step 5: Verify app launches**

Run: `cd /tmp/ont-sma-seq && mkdir -p /tmp/test-experiment && python -m bin.viz.app /tmp/test-experiment &`
Expected: Server starts on port 8050, browser opens.
Kill: `kill %1`

**Step 6: Commit**

```bash
git add bin/viz/ env/env.yml
git commit -m "feat: scaffold viz app with FastAPI + HTMX + Pico CSS"
```

---

### Task 2: Pydantic Models

**Files:**
- Create: `bin/viz/models.py`
- Create: `tests/test_viz_models.py`

**Step 1: Write tests for Pydantic models**

Create `tests/test_viz_models.py`:
```python
"""Tests for viz Pydantic models."""

from __future__ import annotations

import pytest

from viz.models import (
    Assumption,
    BarcodeAssignment,
    ConstructConfig,
    DemultiplexingConfig,
    ExperimentConfig,
    SampleSheetEntry,
    ScoringParams,
    ArrangementConfig,
    TargetRef,
    TruncationConfig,
    TruncationRules,
)


class TestExperimentConfig:
    def test_minimal_config(self):
        cfg = ExperimentConfig(description="Test experiment")
        assert cfg.description == "Test experiment"

    def test_full_config_roundtrip(self):
        cfg = ExperimentConfig(
            description="CYP2D6 variant library",
            construct=ConstructConfig(
                adapter_5prime="AATGTACTTCGTTCAGTTACGTATTGCT",
                adapter_3prime="AATGTACTTCGTTCAGTTACGTATTGCT",
                insert_type="amplicon",
            ),
            demultiplexing=DemultiplexingConfig(
                mode="dual",
                start_barcode=BarcodeAssignment(
                    assignments={"nb05": "library_A"}
                ),
                end_barcode=BarcodeAssignment(
                    assignments={"nb10": "forward"}
                ),
            ),
        )
        assert cfg.demultiplexing.mode == "dual"
        assert cfg.construct.insert_type == "amplicon"


class TestTruncationRules:
    def test_defaults(self):
        rules = TruncationRules()
        assert rules.full == "assign_to_target"
        assert rules.adapter_only == "unclassified"
        assert rules.chimeric == "flag_for_review"

    def test_valid_rules(self):
        rules = TruncationRules(
            trunc_target="assign_by_end_barcode",
        )
        assert rules.trunc_target == "assign_by_end_barcode"


class TestSampleSheetEntry:
    def test_duplexed_barcode(self):
        entry = SampleSheetEntry(
            flow_cell_id="FAL12345",
            kit="SQK-NBD114-96",
            barcode="barcode05--barcode10",
            alias="CYP2D6_v04_fwd",
            type="test_sample",
        )
        assert entry.upstream_barcode == "nb05"
        assert entry.downstream_barcode == "nb10"

    def test_single_barcode_rejected(self):
        with pytest.raises(ValueError, match="duplexed"):
            SampleSheetEntry(
                flow_cell_id="FAL12345",
                kit="SQK-NBD114-96",
                barcode="barcode05",
                alias="target_A",
            )


class TestArrangementConfig:
    def test_defaults(self):
        arr = ArrangementConfig(
            name="custom",
            barcode1_pattern="BC%02i",
            first_index=1,
            last_index=96,
        )
        assert arr.name == "custom"
        assert arr.scoring is not None
        assert arr.scoring.max_barcode_penalty == 11


class TestTargetRef:
    def test_gc_content(self):
        t = TargetRef(tgt_id="test", sequence="AATTCCGG", length=8)
        assert t.gc_content == 50.0

    def test_gc_empty(self):
        t = TargetRef(tgt_id="test", sequence="", length=0)
        assert t.gc_content == 0.0
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && python -m pytest tests/test_viz_models.py -v`
Expected: ImportError - `viz.models` not found

**Step 3: Implement models**

Create `bin/viz/models.py`:
```python
"""Pydantic models for SMA-seq experiment configuration.

Maps to ONT-native file formats (MinKNOW CSV, Dorado TOML/FASTA)
plus an SMA-seq overlay TOML for pipeline-specific settings.
"""

from __future__ import annotations

import re
from typing import Any

from pydantic import BaseModel, Field, field_validator, computed_field


# ---------------------------------------------------------------------------
# Barcode normalization (reused from sample_sheet.py)
# ---------------------------------------------------------------------------

_BARCODE_RE = re.compile(r"barcode(\d+)")
_DUPLEXED_RE = re.compile(r"^(barcode\d+)--(barcode\d+)$")


def _normalize_barcode(name: str) -> str:
    m = _BARCODE_RE.fullmatch(name)
    if not m:
        raise ValueError(f"Invalid barcode name '{name}'")
    num = int(m.group(1))
    if num < 1 or num > 96:
        raise ValueError(f"Barcode number {num} out of range 1-96")
    return f"nb{num:02d}"


# ---------------------------------------------------------------------------
# Scoring parameters (Dorado arrangement.toml [scoring] section)
# ---------------------------------------------------------------------------


class ScoringParams(BaseModel):
    """Dorado barcode scoring parameters."""
    max_barcode_penalty: int = 11
    min_barcode_penalty_dist: int = 3
    min_separation_only_dist: int = 6
    barcode_end_proximity: int = 75
    flank_left_pad: int = 5
    flank_right_pad: int = 10
    front_barcode_window: int = 175
    rear_barcode_window: int = 175
    midstrand_flank_score: float = 0.95


# ---------------------------------------------------------------------------
# Arrangement config (Dorado arrangement.toml [arrangement] section)
# ---------------------------------------------------------------------------


class ArrangementConfig(BaseModel):
    """Dorado custom barcode arrangement configuration."""
    name: str
    kit: str = ""
    mask1_front: str = ""
    mask1_rear: str = ""
    mask2_front: str = ""
    mask2_rear: str = ""
    barcode1_pattern: str = "BC%02i"
    barcode2_pattern: str = "BC%02i"
    first_index: int = 1
    last_index: int = 96
    rear_only_barcodes: bool = False
    scoring: ScoringParams = Field(default_factory=ScoringParams)


# ---------------------------------------------------------------------------
# Sample sheet entry (MinKNOW CSV row)
# ---------------------------------------------------------------------------


class SampleSheetEntry(BaseModel):
    """One row of a MinKNOW sample sheet CSV."""
    flow_cell_id: str = ""
    kit: str = ""
    barcode: str
    alias: str
    experiment_id: str = ""
    sample_id: str = ""
    type: str = "test_sample"

    @field_validator("barcode")
    @classmethod
    def validate_duplexed(cls, v: str) -> str:
        if not _DUPLEXED_RE.match(v):
            raise ValueError(
                f"Barcode '{v}' must be in duplexed 'barcodeNN--barcodeNN' format"
            )
        return v

    @computed_field
    @property
    def upstream_barcode(self) -> str:
        m = _DUPLEXED_RE.match(self.barcode)
        return _normalize_barcode(m.group(1)) if m else ""

    @computed_field
    @property
    def downstream_barcode(self) -> str:
        m = _DUPLEXED_RE.match(self.barcode)
        return _normalize_barcode(m.group(2)) if m else ""


# ---------------------------------------------------------------------------
# Target reference
# ---------------------------------------------------------------------------


class TargetRef(BaseModel):
    """A target reference sequence."""
    tgt_id: str
    sequence: str
    length: int

    @computed_field
    @property
    def gc_content(self) -> float:
        if not self.sequence:
            return 0.0
        gc = sum(1 for c in self.sequence.upper() if c in "GC")
        return round(gc / len(self.sequence) * 100, 1)


# ---------------------------------------------------------------------------
# Demultiplexing config (SMA-seq overlay)
# ---------------------------------------------------------------------------


class BarcodeAssignment(BaseModel):
    """Independent barcode-to-label assignments for one end."""
    assignments: dict[str, str] = Field(default_factory=dict)


class DemuxPairEntry(BaseModel):
    """One dual-barcode pair assignment."""
    start: str
    end: str
    alias: str


class DemultiplexingConfig(BaseModel):
    """Demultiplexing strategy configuration."""
    mode: str = "dual"  # start_only | end_only | dual
    start_barcode: BarcodeAssignment = Field(default_factory=BarcodeAssignment)
    end_barcode: BarcodeAssignment = Field(default_factory=BarcodeAssignment)
    pairs: list[DemuxPairEntry] = Field(default_factory=list)

    @field_validator("mode")
    @classmethod
    def validate_mode(cls, v: str) -> str:
        allowed = {"start_only", "end_only", "dual"}
        if v not in allowed:
            raise ValueError(f"mode must be one of {allowed}")
        return v


# ---------------------------------------------------------------------------
# Truncation config (SMA-seq overlay)
# ---------------------------------------------------------------------------

VALID_ASSIGNMENT_RULES = {
    "assign_to_target",
    "assign_by_start_barcode",
    "assign_by_end_barcode",
    "unclassified",
    "flag_for_review",
}


class TruncationRules(BaseModel):
    """Assignment rules for each truncation class."""
    full: str = "assign_to_target"
    trunc_3prime: str = "assign_to_target"
    trunc_target: str = "assign_by_start_barcode"
    trunc_barcode: str = "assign_by_start_barcode"
    adapter_only: str = "unclassified"
    chimeric: str = "flag_for_review"


class AutoReferences(BaseModel):
    """Auto-generated truncated reference config."""
    enabled: bool = True


class TruncationConfig(BaseModel):
    """Truncation detection and handling configuration."""
    min_barcode_confidence: float = 0.5
    min_target_fraction: float = 0.1
    adapter_search_window: int = 50
    rules: TruncationRules = Field(default_factory=TruncationRules)
    auto_references: AutoReferences = Field(default_factory=AutoReferences)


# ---------------------------------------------------------------------------
# Construct config (SMA-seq overlay)
# ---------------------------------------------------------------------------


class ConstructConfig(BaseModel):
    """Describes the expected molecular structure of a read."""
    adapter_5prime: str = ""
    adapter_3prime: str = ""
    insert_type: str = "amplicon"  # amplicon | genomic | synthetic
    notes: str = ""

    @field_validator("insert_type")
    @classmethod
    def validate_insert_type(cls, v: str) -> str:
        allowed = {"amplicon", "genomic", "synthetic"}
        if v not in allowed:
            raise ValueError(f"insert_type must be one of {allowed}")
        return v


# ---------------------------------------------------------------------------
# Assumptions
# ---------------------------------------------------------------------------


class Assumption(BaseModel):
    """A documented pipeline assumption."""
    key: str
    text: str
    why: str


# ---------------------------------------------------------------------------
# Quality metric documentation
# ---------------------------------------------------------------------------


class QualityMetrics(BaseModel):
    """Quality metric formula documentation."""
    q_bc: str = "Probability-averaged Phred: -10*log10(mean(10^(-Qi/10)))"
    q_ld: str = "-10*log10(min(max(1/L^2, ed/L), 1.0))"


# ---------------------------------------------------------------------------
# Classification config
# ---------------------------------------------------------------------------


class ClassificationConfig(BaseModel):
    """SMA-seq barcode classification parameters."""
    barcode_search_window: int = 100
    confidence_formula: str = "1.0 - (ed / barcode_length)"
    ambiguity_triggers_full_construct: bool = True


# ---------------------------------------------------------------------------
# Top-level experiment config (sma_experiment.toml)
# ---------------------------------------------------------------------------


class ExperimentConfig(BaseModel):
    """Top-level SMA-seq experiment configuration.

    This model maps to sma_experiment.toml and captures all
    SMA-seq-specific settings not covered by ONT-native formats.
    """
    description: str = ""
    construct: ConstructConfig = Field(default_factory=ConstructConfig)
    demultiplexing: DemultiplexingConfig = Field(
        default_factory=DemultiplexingConfig
    )
    classification: ClassificationConfig = Field(
        default_factory=ClassificationConfig
    )
    quality: QualityMetrics = Field(default_factory=QualityMetrics)
    truncation: TruncationConfig = Field(default_factory=TruncationConfig)
    assumptions: list[Assumption] = Field(default_factory=lambda: [
        Assumption(
            key="all_reads_ingested",
            text="All reads captured; filtering is post-hoc via SQL",
            why="Avoids premature filtering bias",
        ),
        Assumption(
            key="untrimmed_bams",
            text="BAMs must be untrimmed for barcode classification",
            why="Barcodes are in the read sequence",
        ),
        Assumption(
            key="native_barcodes_only",
            text="Only ONT native 24bp barcodes supported",
            why="Hardcoded in barcodes.py",
        ),
        Assumption(
            key="duplexed_required",
            text="Duplexed barcode pairs required",
            why="Single barcodes can't uniquely identify fwd/rev targets",
        ),
    ])
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_viz_models.py -v`
Expected: All tests PASS

**Step 5: Commit**

```bash
git add bin/viz/models.py tests/test_viz_models.py
git commit -m "feat: add Pydantic models for experiment config, arrangement, sample sheet, truncation"
```

---

### Task 3: Config Store (Read/Write ONT Files)

**Files:**
- Create: `bin/viz/config_store.py`
- Create: `tests/test_config_store.py`

**Step 1: Write tests for config store**

Create `tests/test_config_store.py`:
```python
"""Tests for config file read/write operations."""

from __future__ import annotations

from pathlib import Path

import pytest

from viz.config_store import ConfigStore


@pytest.fixture()
def experiment_dir(tmp_path: Path) -> Path:
    """Create a minimal experiment directory with config files."""
    exp = tmp_path / "experiment"
    exp.mkdir()
    (exp / "references").mkdir()

    # MinKNOW sample sheet
    (exp / "sample_sheet.csv").write_text(
        "flow_cell_id,kit,barcode,alias,type\n"
        "FAL12345,SQK-NBD114-96,barcode05--barcode10,CYP2D6_fwd,test_sample\n"
        "FAL12345,SQK-NBD114-96,barcode10--barcode05,CYP2D6_rev,test_sample\n"
    )

    # Reference FASTAs
    (exp / "references" / "CYP2D6_fwd.fasta").write_text(
        ">CYP2D6_fwd\nACGTACGTACGT\n"
    )
    (exp / "references" / "CYP2D6_rev.fasta").write_text(
        ">CYP2D6_rev\nTGCATGCATGCA\n"
    )

    return exp


class TestConfigStoreInit:
    def test_creates_default_toml_if_missing(self, experiment_dir: Path):
        store = ConfigStore(experiment_dir)
        assert (experiment_dir / "sma_experiment.toml").exists()

    def test_loads_existing_toml(self, experiment_dir: Path):
        (experiment_dir / "sma_experiment.toml").write_text(
            '[experiment]\ndescription = "My experiment"\n'
        )
        store = ConfigStore(experiment_dir)
        assert store.experiment_config.description == "My experiment"


class TestSampleSheet:
    def test_read_sample_sheet(self, experiment_dir: Path):
        store = ConfigStore(experiment_dir)
        entries = store.read_sample_sheet()
        assert len(entries) == 2
        assert entries[0].alias == "CYP2D6_fwd"
        assert entries[0].upstream_barcode == "nb05"

    def test_write_sample_sheet(self, experiment_dir: Path):
        store = ConfigStore(experiment_dir)
        entries = store.read_sample_sheet()
        entries.append(entries[0].model_copy(update={
            "barcode": "barcode01--barcode02",
            "alias": "GeneX_fwd",
        }))
        store.write_sample_sheet(entries)
        reloaded = store.read_sample_sheet()
        assert len(reloaded) == 3


class TestTargets:
    def test_read_targets(self, experiment_dir: Path):
        store = ConfigStore(experiment_dir)
        targets = store.read_targets()
        assert len(targets) == 2
        assert "CYP2D6_fwd" in {t.tgt_id for t in targets}

    def test_write_target(self, experiment_dir: Path):
        store = ConfigStore(experiment_dir)
        store.write_target("NewTarget", "AAACCCGGGTTT")
        targets = store.read_targets()
        ids = {t.tgt_id for t in targets}
        assert "NewTarget" in ids

    def test_delete_target(self, experiment_dir: Path):
        store = ConfigStore(experiment_dir)
        store.delete_target("CYP2D6_fwd")
        targets = store.read_targets()
        ids = {t.tgt_id for t in targets}
        assert "CYP2D6_fwd" not in ids


class TestArrangement:
    def test_read_arrangement_default(self, experiment_dir: Path):
        store = ConfigStore(experiment_dir)
        arr = store.read_arrangement()
        assert arr is None  # No arrangement.toml yet

    def test_write_and_read_arrangement(self, experiment_dir: Path):
        store = ConfigStore(experiment_dir)
        from viz.models import ArrangementConfig
        arr = ArrangementConfig(name="custom_kit", first_index=1, last_index=4)
        store.write_arrangement(arr)
        loaded = store.read_arrangement()
        assert loaded is not None
        assert loaded.name == "custom_kit"


class TestExperimentConfig:
    def test_save_and_reload(self, experiment_dir: Path):
        store = ConfigStore(experiment_dir)
        store.experiment_config.description = "Updated description"
        store.save_experiment_config()
        store2 = ConfigStore(experiment_dir)
        assert store2.experiment_config.description == "Updated description"


class TestValidation:
    def test_validate_consistent(self, experiment_dir: Path):
        store = ConfigStore(experiment_dir)
        errors = store.validate()
        assert isinstance(errors, list)
        # Should be valid - aliases match FASTA headers
        alias_errors = [e for e in errors if "alias" in e.lower() and "missing" in e.lower()]
        assert len(alias_errors) == 0

    def test_validate_missing_reference(self, experiment_dir: Path):
        # Add a sample sheet entry with no matching FASTA
        (experiment_dir / "sample_sheet.csv").write_text(
            "flow_cell_id,kit,barcode,alias,type\n"
            "FAL12345,SQK-NBD114-96,barcode05--barcode10,MISSING_TARGET,test_sample\n"
        )
        store = ConfigStore(experiment_dir)
        errors = store.validate()
        assert any("MISSING_TARGET" in e for e in errors)


class TestVersioning:
    def test_backup_created_on_write(self, experiment_dir: Path):
        store = ConfigStore(experiment_dir)
        store.write_sample_sheet(store.read_sample_sheet())
        backups = list((experiment_dir / ".backups").glob("sample_sheet_*.csv"))
        assert len(backups) >= 1
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_config_store.py -v`
Expected: ImportError

**Step 3: Implement config store**

Create `bin/viz/config_store.py`:
```python
"""Read/write ONT-native config files and SMA-seq overlay TOML.

Handles:
- sample_sheet.csv (MinKNOW CSV format)
- arrangement.toml (Dorado custom barcode arrangement)
- barcodes.fasta (Dorado barcode sequences)
- references/*.fasta (target reference sequences)
- sma_experiment.toml (SMA-seq overlay)
"""

from __future__ import annotations

import csv
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib  # type: ignore[no-redef]

import tomli_w

from viz.models import (
    ArrangementConfig,
    ExperimentConfig,
    SampleSheetEntry,
    ScoringParams,
    TargetRef,
)


class ConfigStore:
    """Manages all config files for an SMA-seq experiment directory."""

    MAX_BACKUPS = 10

    def __init__(self, experiment_dir: Path) -> None:
        self.dir = Path(experiment_dir).resolve()
        self._backups_dir = self.dir / ".backups"
        self.experiment_config = self._load_or_create_experiment_config()

    # ------------------------------------------------------------------
    # sma_experiment.toml
    # ------------------------------------------------------------------

    def _load_or_create_experiment_config(self) -> ExperimentConfig:
        toml_path = self.dir / "sma_experiment.toml"
        if toml_path.exists():
            with open(toml_path, "rb") as f:
                data = tomllib.load(f)
            # Flatten: top-level keys map to ExperimentConfig fields
            flat: dict[str, Any] = {}
            if "experiment" in data:
                flat.update(data["experiment"])
            for key in (
                "construct", "demultiplexing", "classification",
                "quality", "truncation", "assumptions",
            ):
                if key in data:
                    flat[key] = data[key]
            return ExperimentConfig(**flat)
        else:
            cfg = ExperimentConfig()
            self._write_experiment_toml(cfg)
            return cfg

    def save_experiment_config(self) -> None:
        self._backup("sma_experiment.toml")
        self._write_experiment_toml(self.experiment_config)

    def _write_experiment_toml(self, cfg: ExperimentConfig) -> None:
        data: dict[str, Any] = {}
        data["experiment"] = {"description": cfg.description}
        data["construct"] = cfg.construct.model_dump()
        data["demultiplexing"] = cfg.demultiplexing.model_dump()
        data["classification"] = cfg.classification.model_dump()
        data["quality"] = cfg.quality.model_dump()
        data["truncation"] = cfg.truncation.model_dump()
        data["assumptions"] = {"entries": [
            a.model_dump() for a in cfg.assumptions
        ]}
        toml_path = self.dir / "sma_experiment.toml"
        with open(toml_path, "wb") as f:
            tomli_w.dump(data, f)

    # ------------------------------------------------------------------
    # sample_sheet.csv (MinKNOW CSV)
    # ------------------------------------------------------------------

    def read_sample_sheet(self) -> list[SampleSheetEntry]:
        csv_path = self.dir / "sample_sheet.csv"
        if not csv_path.exists():
            return []
        entries: list[SampleSheetEntry] = []
        with open(csv_path, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                entries.append(SampleSheetEntry(
                    flow_cell_id=row.get("flow_cell_id", ""),
                    kit=row.get("kit", ""),
                    barcode=row["barcode"],
                    alias=row["alias"],
                    experiment_id=row.get("experiment_id", ""),
                    sample_id=row.get("sample_id", ""),
                    type=row.get("type", "test_sample"),
                ))
        return entries

    def write_sample_sheet(self, entries: list[SampleSheetEntry]) -> None:
        self._backup("sample_sheet.csv")
        csv_path = self.dir / "sample_sheet.csv"
        fieldnames = [
            "flow_cell_id", "kit", "barcode", "alias",
            "experiment_id", "sample_id", "type",
        ]
        with open(csv_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            for entry in entries:
                writer.writerow({
                    "flow_cell_id": entry.flow_cell_id,
                    "kit": entry.kit,
                    "barcode": entry.barcode,
                    "alias": entry.alias,
                    "experiment_id": entry.experiment_id,
                    "sample_id": entry.sample_id,
                    "type": entry.type,
                })

    # ------------------------------------------------------------------
    # arrangement.toml (Dorado custom barcodes)
    # ------------------------------------------------------------------

    def read_arrangement(self) -> ArrangementConfig | None:
        toml_path = self.dir / "arrangement.toml"
        if not toml_path.exists():
            return None
        with open(toml_path, "rb") as f:
            data = tomllib.load(f)
        arr_data = data.get("arrangement", {})
        scoring_data = data.get("scoring", {})
        if scoring_data:
            arr_data["scoring"] = scoring_data
        return ArrangementConfig(**arr_data)

    def write_arrangement(self, arr: ArrangementConfig) -> None:
        self._backup("arrangement.toml")
        toml_path = self.dir / "arrangement.toml"
        data: dict[str, Any] = {"arrangement": {}, "scoring": {}}
        arr_dict = arr.model_dump()
        scoring = arr_dict.pop("scoring", {})
        data["arrangement"] = arr_dict
        data["scoring"] = scoring
        with open(toml_path, "wb") as f:
            tomli_w.dump(data, f)

    # ------------------------------------------------------------------
    # barcodes.fasta (Dorado barcode sequences)
    # ------------------------------------------------------------------

    def read_barcode_sequences(self) -> dict[str, str]:
        fasta_path = self.dir / "barcodes.fasta"
        if not fasta_path.exists():
            return {}
        return self._read_fasta(fasta_path)

    def write_barcode_sequences(self, sequences: dict[str, str]) -> None:
        self._backup("barcodes.fasta")
        fasta_path = self.dir / "barcodes.fasta"
        self._write_fasta(fasta_path, sequences)

    # ------------------------------------------------------------------
    # references/*.fasta
    # ------------------------------------------------------------------

    def read_targets(self) -> list[TargetRef]:
        ref_dir = self.dir / "references"
        if not ref_dir.is_dir():
            return []
        targets: list[TargetRef] = []
        for pattern in ("*.fasta", "*.fa"):
            for fasta_path in ref_dir.glob(pattern):
                seqs = self._read_fasta(fasta_path)
                for tgt_id, seq in seqs.items():
                    targets.append(TargetRef(
                        tgt_id=tgt_id, sequence=seq, length=len(seq),
                    ))
        return sorted(targets, key=lambda t: t.tgt_id)

    def write_target(self, tgt_id: str, sequence: str) -> None:
        ref_dir = self.dir / "references"
        ref_dir.mkdir(exist_ok=True)
        fasta_path = ref_dir / f"{tgt_id}.fasta"
        self._write_fasta(fasta_path, {tgt_id: sequence})

    def delete_target(self, tgt_id: str) -> None:
        ref_dir = self.dir / "references"
        for pattern in ("*.fasta", "*.fa"):
            for fasta_path in ref_dir.glob(pattern):
                seqs = self._read_fasta(fasta_path)
                if tgt_id in seqs:
                    self._backup(f"references/{fasta_path.name}")
                    fasta_path.unlink()
                    return

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate(self) -> list[str]:
        errors: list[str] = []

        # Check sample sheet aliases have matching references
        entries = self.read_sample_sheet()
        targets = self.read_targets()
        target_ids = {t.tgt_id for t in targets}
        aliases = {e.alias for e in entries}

        for alias in aliases:
            if alias not in target_ids:
                errors.append(
                    f"Sample sheet alias '{alias}' has no matching reference FASTA"
                )

        for tgt_id in target_ids:
            if tgt_id not in aliases and entries:
                errors.append(
                    f"Reference '{tgt_id}' not referenced by any sample sheet entry"
                )

        return errors

    # ------------------------------------------------------------------
    # File versioning
    # ------------------------------------------------------------------

    def _backup(self, relative_path: str) -> None:
        src = self.dir / relative_path
        if not src.exists():
            return
        self._backups_dir.mkdir(exist_ok=True)
        stem = Path(relative_path).stem
        suffix = Path(relative_path).suffix
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        dst = self._backups_dir / f"{stem}_{ts}{suffix}"
        shutil.copy2(src, dst)

        # Prune old backups
        pattern = f"{stem}_*{suffix}"
        backups = sorted(self._backups_dir.glob(pattern))
        while len(backups) > self.MAX_BACKUPS:
            backups.pop(0).unlink()

    # ------------------------------------------------------------------
    # FASTA helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _read_fasta(path: Path) -> dict[str, str]:
        sequences: dict[str, str] = {}
        current_id: str | None = None
        parts: list[str] = []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_id is not None:
                        sequences[current_id] = "".join(parts)
                    current_id = line[1:].split()[0]
                    parts = []
                else:
                    parts.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(parts)
        return sequences

    @staticmethod
    def _write_fasta(path: Path, sequences: dict[str, str]) -> None:
        with open(path, "w") as f:
            for seq_id, seq in sequences.items():
                f.write(f">{seq_id}\n")
                # Wrap at 80 chars
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i + 80] + "\n")
```

**Step 4: Run tests to verify they pass**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_config_store.py -v`
Expected: All tests PASS

**Step 5: Commit**

```bash
git add bin/viz/config_store.py tests/test_config_store.py
git commit -m "feat: add config store for reading/writing ONT-native config files"
```

---

### Task 4: Base Template & Navigation

**Files:**
- Create: `bin/viz/templates/base.html`
- Create: `bin/viz/static/style.css`
- Modify: `bin/viz/templates/dashboard.html`
- Modify: `bin/viz/app.py`

**Step 1: Create base template with sidebar navigation**

Create `bin/viz/templates/base.html`:
```html
<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{% block title %}SMA-seq Config Visualizer{% endblock %}</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@picocss/pico@2/css/pico.min.css">
    <link rel="stylesheet" href="/static/style.css">
    <script src="https://unpkg.com/htmx.org@2.0.4"></script>
    {% block head %}{% endblock %}
</head>
<body>
  <div class="viz-layout">
    <nav class="viz-sidebar">
      <header>
        <h3>SMA-seq Viz</h3>
        <small id="validation-badge"
               hx-get="/api/validate/badge"
               hx-trigger="load, configChanged from:body"
               hx-swap="innerHTML">
        </small>
      </header>
      <ul>
        <li><a href="/" class="{% if active_page == 'dashboard' %}active{% endif %}">Dashboard</a></li>
        <li><a href="/sample-sheet" class="{% if active_page == 'sample_sheet' %}active{% endif %}">Sample Sheet</a></li>
        <li><a href="/barcodes" class="{% if active_page == 'barcodes' %}active{% endif %}">Barcodes</a></li>
        <li><a href="/construct" class="{% if active_page == 'construct' %}active{% endif %}">Construct</a></li>
        <li><a href="/targets" class="{% if active_page == 'targets' %}active{% endif %}">Targets</a></li>
        <li><a href="/assumptions" class="{% if active_page == 'assumptions' %}active{% endif %}">Assumptions</a></li>
      </ul>
      <footer>
        <small>{{ experiment_dir_name }}</small>
      </footer>
    </nav>
    <main class="viz-main container">
      {% block content %}{% endblock %}
    </main>
  </div>
</body>
</html>
```

**Step 2: Create CSS**

Create `bin/viz/static/style.css`:
```css
/* SMA-seq Config Visualizer */

.viz-layout {
  display: flex;
  min-height: 100vh;
}

.viz-sidebar {
  width: 220px;
  min-width: 220px;
  padding: 1rem;
  border-right: 1px solid var(--pico-muted-border-color);
  background: var(--pico-card-background-color);
}

.viz-sidebar header h3 {
  margin-bottom: 0.25rem;
}

.viz-sidebar ul {
  list-style: none;
  padding: 0;
  margin: 1rem 0;
}

.viz-sidebar ul li {
  margin: 0;
  padding: 0;
}

.viz-sidebar ul li a {
  display: block;
  padding: 0.5rem 0.75rem;
  border-radius: 4px;
  text-decoration: none;
}

.viz-sidebar ul li a.active {
  background: var(--pico-primary-background);
  color: var(--pico-primary-inverse);
}

.viz-sidebar footer {
  margin-top: auto;
  padding-top: 1rem;
  border-top: 1px solid var(--pico-muted-border-color);
}

.viz-main {
  flex: 1;
  padding: 2rem;
}

/* Validation badge */
.badge-ok { color: var(--pico-ins-color); }
.badge-error { color: var(--pico-del-color); }

/* Construct diagram */
.construct-region {
  cursor: pointer;
}

.construct-region:hover {
  opacity: 0.8;
}

/* Table helpers */
.table-actions {
  white-space: nowrap;
}

.inline-form {
  display: inline;
}

/* Truncation class colors */
.trunc-full { background: #4caf50; color: white; }
.trunc-3prime { background: #8bc34a; color: white; }
.trunc-target { background: #ff9800; color: white; }
.trunc-barcode { background: #ff5722; color: white; }
.trunc-adapter { background: #9e9e9e; color: white; }
.trunc-chimeric { background: #f44336; color: white; }

.trunc-label {
  display: inline-block;
  padding: 0.15rem 0.5rem;
  border-radius: 3px;
  font-size: 0.85rem;
}
```

**Step 3: Update dashboard template to extend base**

Rewrite `bin/viz/templates/dashboard.html`:
```html
{% extends "base.html" %}
{% block title %}Dashboard - SMA-seq Config Visualizer{% endblock %}
{% block content %}
<h1>Dashboard</h1>
<article>
  <header>Experiment</header>
  <p><strong>Directory:</strong> {{ experiment_dir }}</p>
  <p><strong>Description:</strong> {{ config.description or "Not set" }}</p>
</article>

<h2>File Status</h2>
<div class="grid">
  {% for name, exists in file_status.items() %}
  <article>
    <header>{{ name }}</header>
    {% if exists %}
      <p class="badge-ok">Present</p>
    {% else %}
      <p class="badge-error">Missing</p>
    {% endif %}
  </article>
  {% endfor %}
</div>

<h2>Quick Actions</h2>
<div class="grid">
  <button hx-post="/api/export" hx-swap="outerHTML">Export Static HTML</button>
</div>
{% endblock %}
```

**Step 4: Update app.py with config store integration and all page routes**

Update `bin/viz/app.py` to wire ConfigStore to routes, add page routes for sample-sheet, barcodes, construct, targets, assumptions. Each route renders its template with `active_page` and `experiment_dir_name` context.

**Step 5: Verify all pages render**

Run: `cd /tmp/ont-sma-seq && python -m bin.viz.app /tmp/test-experiment &`
Visit each URL manually. Kill: `kill %1`

**Step 6: Commit**

```bash
git add bin/viz/templates/ bin/viz/static/ bin/viz/app.py
git commit -m "feat: add base template with sidebar nav, CSS, and page routes"
```

---

### Task 5: API Endpoints - Sample Sheet CRUD

**Files:**
- Create: `bin/viz/api.py`
- Create: `bin/viz/templates/sample_sheet.html`
- Create: `tests/test_api.py`
- Modify: `bin/viz/app.py` (include API router)

**Step 1: Write API tests**

Create `tests/test_api.py`:
```python
"""Tests for the viz REST API."""

from __future__ import annotations

from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from viz.app import app, set_experiment_dir


@pytest.fixture()
def experiment_dir(tmp_path: Path) -> Path:
    exp = tmp_path / "experiment"
    exp.mkdir()
    (exp / "references").mkdir()
    (exp / "sample_sheet.csv").write_text(
        "flow_cell_id,kit,barcode,alias,type\n"
        "FAL12345,SQK-NBD114-96,barcode05--barcode10,CYP2D6_fwd,test_sample\n"
    )
    (exp / "references" / "CYP2D6_fwd.fasta").write_text(">CYP2D6_fwd\nACGT\n")
    return exp


@pytest.fixture()
def client(experiment_dir: Path) -> TestClient:
    set_experiment_dir(experiment_dir)
    return TestClient(app)


class TestSampleSheetAPI:
    def test_get_sample_sheet(self, client: TestClient):
        resp = client.get("/api/sample-sheet")
        assert resp.status_code == 200
        data = resp.json()
        assert len(data) == 1
        assert data[0]["alias"] == "CYP2D6_fwd"

    def test_put_sample_sheet(self, client: TestClient):
        entries = [
            {
                "flow_cell_id": "FAL12345",
                "kit": "SQK-NBD114-96",
                "barcode": "barcode05--barcode10",
                "alias": "CYP2D6_fwd",
                "type": "test_sample",
            },
            {
                "flow_cell_id": "FAL12345",
                "kit": "SQK-NBD114-96",
                "barcode": "barcode01--barcode02",
                "alias": "GeneX_fwd",
                "type": "test_sample",
            },
        ]
        resp = client.put("/api/sample-sheet", json=entries)
        assert resp.status_code == 200
        reloaded = client.get("/api/sample-sheet").json()
        assert len(reloaded) == 2


class TestTargetsAPI:
    def test_get_targets(self, client: TestClient):
        resp = client.get("/api/targets")
        assert resp.status_code == 200
        data = resp.json()
        assert len(data) == 1
        assert data[0]["tgt_id"] == "CYP2D6_fwd"

    def test_post_target(self, client: TestClient):
        resp = client.post("/api/targets", json={
            "tgt_id": "NewTarget",
            "sequence": "AAACCCGGGTTT",
        })
        assert resp.status_code == 200
        targets = client.get("/api/targets").json()
        assert len(targets) == 2

    def test_delete_target(self, client: TestClient):
        resp = client.delete("/api/targets/CYP2D6_fwd")
        assert resp.status_code == 200
        targets = client.get("/api/targets").json()
        assert len(targets) == 0


class TestValidationAPI:
    def test_validate(self, client: TestClient):
        resp = client.get("/api/validate")
        assert resp.status_code == 200
        data = resp.json()
        assert isinstance(data["errors"], list)


class TestExperimentAPI:
    def test_get_experiment(self, client: TestClient):
        resp = client.get("/api/experiment")
        assert resp.status_code == 200

    def test_put_experiment(self, client: TestClient):
        resp = client.put("/api/experiment", json={
            "description": "Updated",
        })
        assert resp.status_code == 200
        data = client.get("/api/experiment").json()
        assert data["description"] == "Updated"
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_api.py -v`
Expected: ImportError

**Step 3: Implement API router**

Create `bin/viz/api.py`:
```python
"""REST API endpoints for SMA-seq config CRUD."""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, HTTPException

from viz.config_store import ConfigStore
from viz.models import ArrangementConfig, SampleSheetEntry

router = APIRouter(prefix="/api")

# Set by app.py at startup
_store: ConfigStore | None = None


def set_store(store: ConfigStore) -> None:
    global _store
    _store = store


def get_store() -> ConfigStore:
    if _store is None:
        raise RuntimeError("ConfigStore not initialized")
    return _store


# ------------------------------------------------------------------
# Experiment
# ------------------------------------------------------------------


@router.get("/experiment")
async def get_experiment() -> dict[str, Any]:
    store = get_store()
    return store.experiment_config.model_dump()


@router.put("/experiment")
async def put_experiment(data: dict[str, Any]) -> dict[str, str]:
    store = get_store()
    for key, value in data.items():
        if hasattr(store.experiment_config, key):
            setattr(store.experiment_config, key, value)
    store.save_experiment_config()
    return {"status": "ok"}


# ------------------------------------------------------------------
# Sample Sheet
# ------------------------------------------------------------------


@router.get("/sample-sheet")
async def get_sample_sheet() -> list[dict[str, Any]]:
    store = get_store()
    entries = store.read_sample_sheet()
    return [e.model_dump() for e in entries]


@router.put("/sample-sheet")
async def put_sample_sheet(data: list[dict[str, Any]]) -> dict[str, str]:
    store = get_store()
    entries = [SampleSheetEntry(**d) for d in data]
    store.write_sample_sheet(entries)
    return {"status": "ok"}


@router.get("/sample-sheet/export")
async def export_sample_sheet():
    from fastapi.responses import FileResponse
    store = get_store()
    csv_path = store.dir / "sample_sheet.csv"
    if not csv_path.exists():
        raise HTTPException(404, "No sample sheet found")
    return FileResponse(csv_path, media_type="text/csv", filename="sample_sheet.csv")


# ------------------------------------------------------------------
# Barcodes
# ------------------------------------------------------------------


@router.get("/barcodes")
async def get_barcodes() -> dict[str, Any]:
    store = get_store()
    arr = store.read_arrangement()
    seqs = store.read_barcode_sequences()
    # Include built-in ONT barcodes from barcodes.py
    from barcodes import BARCODES
    return {
        "arrangement": arr.model_dump() if arr else None,
        "custom_sequences": seqs,
        "ont_barcodes": BARCODES,
    }


@router.put("/barcodes/arrangement")
async def put_arrangement(data: dict[str, Any]) -> dict[str, str]:
    store = get_store()
    arr = ArrangementConfig(**data)
    store.write_arrangement(arr)
    return {"status": "ok"}


@router.put("/barcodes/sequences")
async def put_barcode_sequences(data: dict[str, str]) -> dict[str, str]:
    store = get_store()
    store.write_barcode_sequences(data)
    return {"status": "ok"}


@router.get("/barcodes/distances")
async def get_barcode_distances() -> dict[str, Any]:
    """Compute pairwise edit distances for used barcodes."""
    import edlib
    store = get_store()
    entries = store.read_sample_sheet()
    used_ids: set[str] = set()
    for e in entries:
        used_ids.add(e.upstream_barcode)
        used_ids.add(e.downstream_barcode)

    from barcodes import BARCODES
    used_barcodes = {k: BARCODES[k] for k in sorted(used_ids) if k in BARCODES}

    ids = list(used_barcodes.keys())
    matrix: list[list[int]] = []
    for i, id_a in enumerate(ids):
        row: list[int] = []
        for j, id_b in enumerate(ids):
            if i == j:
                row.append(0)
            else:
                res = edlib.align(used_barcodes[id_a], used_barcodes[id_b], mode="NW", task="distance")
                row.append(res["editDistance"])
        matrix.append(row)

    return {"ids": ids, "matrix": matrix}


# ------------------------------------------------------------------
# Construct
# ------------------------------------------------------------------


@router.get("/construct")
async def get_construct() -> dict[str, Any]:
    store = get_store()
    return {
        "construct": store.experiment_config.construct.model_dump(),
        "demultiplexing": store.experiment_config.demultiplexing.model_dump(),
        "truncation": store.experiment_config.truncation.model_dump(),
        "classification": store.experiment_config.classification.model_dump(),
    }


@router.put("/construct")
async def put_construct(data: dict[str, Any]) -> dict[str, str]:
    store = get_store()
    if "construct" in data:
        for k, v in data["construct"].items():
            setattr(store.experiment_config.construct, k, v)
    if "demultiplexing" in data:
        from viz.models import DemultiplexingConfig
        store.experiment_config.demultiplexing = DemultiplexingConfig(**data["demultiplexing"])
    if "truncation" in data:
        from viz.models import TruncationConfig
        store.experiment_config.truncation = TruncationConfig(**data["truncation"])
    store.save_experiment_config()
    return {"status": "ok"}


@router.post("/construct/auto-refs")
async def generate_auto_refs() -> dict[str, Any]:
    """Generate truncated reference FASTAs from construct config."""
    store = get_store()
    from barcodes import BARCODES, reverse_complement

    entries = store.read_sample_sheet()
    targets = {t.tgt_id: t for t in store.read_targets()}
    construct = store.experiment_config.construct

    auto_dir = store.dir / "references" / "auto"
    auto_dir.mkdir(parents=True, exist_ok=True)

    generated: list[str] = []

    for entry in entries:
        alias = entry.alias
        if alias not in targets:
            continue

        target_seq = targets[alias].sequence
        bc1_seq = BARCODES.get(entry.upstream_barcode, "")
        bc2_rc_seq = reverse_complement(BARCODES.get(entry.downstream_barcode, ""))
        adapter_5 = construct.adapter_5prime
        adapter_3 = construct.adapter_3prime

        variants = {
            f"{alias}_full": adapter_5 + bc1_seq + target_seq + bc2_rc_seq + adapter_3,
            f"{alias}_no_adapter": bc1_seq + target_seq + bc2_rc_seq,
            f"{alias}_no_end_bc": adapter_5 + bc1_seq + target_seq,
            f"{alias}_bc_only": adapter_5 + bc1_seq,
        }

        for var_id, var_seq in variants.items():
            if var_seq:  # skip if all components empty
                fasta_path = auto_dir / f"{var_id}.fasta"
                store._write_fasta(fasta_path, {var_id: var_seq})
                generated.append(var_id)

    return {"generated": generated}


# ------------------------------------------------------------------
# Targets
# ------------------------------------------------------------------


@router.get("/targets")
async def get_targets() -> list[dict[str, Any]]:
    store = get_store()
    targets = store.read_targets()
    return [t.model_dump() for t in targets]


@router.post("/targets")
async def post_target(data: dict[str, str]) -> dict[str, str]:
    store = get_store()
    store.write_target(data["tgt_id"], data["sequence"])
    return {"status": "ok"}


@router.put("/targets/{alias}")
async def put_target(alias: str, data: dict[str, str]) -> dict[str, str]:
    store = get_store()
    store.delete_target(alias)
    store.write_target(data.get("tgt_id", alias), data["sequence"])
    return {"status": "ok"}


@router.delete("/targets/{alias}")
async def delete_target(alias: str) -> dict[str, str]:
    store = get_store()
    store.delete_target(alias)
    return {"status": "ok"}


# ------------------------------------------------------------------
# Assumptions
# ------------------------------------------------------------------


@router.get("/assumptions")
async def get_assumptions() -> list[dict[str, str]]:
    store = get_store()
    return [a.model_dump() for a in store.experiment_config.assumptions]


@router.put("/assumptions")
async def put_assumptions(data: list[dict[str, str]]) -> dict[str, str]:
    from viz.models import Assumption
    store = get_store()
    store.experiment_config.assumptions = [Assumption(**d) for d in data]
    store.save_experiment_config()
    return {"status": "ok"}


# ------------------------------------------------------------------
# Validation
# ------------------------------------------------------------------


@router.get("/validate")
async def validate() -> dict[str, Any]:
    store = get_store()
    errors = store.validate()
    return {"errors": errors, "valid": len(errors) == 0}


@router.get("/validate/badge")
async def validate_badge():
    from fastapi.responses import HTMLResponse
    store = get_store()
    errors = store.validate()
    if errors:
        return HTMLResponse(
            f'<span class="badge-error">{len(errors)} issue{"s" if len(errors) != 1 else ""}</span>'
        )
    return HTMLResponse('<span class="badge-ok">Valid</span>')


# ------------------------------------------------------------------
# Export
# ------------------------------------------------------------------


@router.post("/export")
async def export_html() -> dict[str, str]:
    # Placeholder - implemented in Task 9
    return {"status": "not_implemented"}
```

**Step 4: Update app.py to include API router and expose `set_experiment_dir`**

Add to `bin/viz/app.py`:
```python
from viz.api import router as api_router, set_store
from viz.config_store import ConfigStore

app.include_router(api_router)

def set_experiment_dir(path: Path) -> None:
    global experiment_dir
    experiment_dir = path.resolve()
    store = ConfigStore(experiment_dir)
    set_store(store)
```

Update `main()` to call `set_experiment_dir()`.

**Step 5: Run API tests**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_api.py -v`
Expected: All tests PASS

**Step 6: Commit**

```bash
git add bin/viz/api.py tests/test_api.py bin/viz/app.py
git commit -m "feat: add REST API endpoints for sample sheet, targets, barcodes, validation"
```

---

### Task 6: Page Templates (Sample Sheet, Barcodes, Targets, Assumptions)

**Files:**
- Create: `bin/viz/templates/sample_sheet.html`
- Create: `bin/viz/templates/barcodes.html`
- Create: `bin/viz/templates/targets.html`
- Create: `bin/viz/templates/assumptions.html`

**Step 1: Create sample sheet template with HTMX CRUD**

Create `bin/viz/templates/sample_sheet.html`:
```html
{% extends "base.html" %}
{% block title %}Sample Sheet - SMA-seq{% endblock %}
{% block content %}
<h1>Sample Sheet</h1>

<div id="validation-inline"
     hx-get="/api/validate/badge"
     hx-trigger="load, configChanged from:body"
     hx-swap="innerHTML"></div>

<table>
  <thead>
    <tr>
      <th>Barcode Pair</th>
      <th>Upstream</th>
      <th>Downstream</th>
      <th>Alias</th>
      <th>Kit</th>
      <th>Type</th>
      <th>Actions</th>
    </tr>
  </thead>
  <tbody id="sample-sheet-body"
         hx-get="/htmx/sample-sheet/rows"
         hx-trigger="load"
         hx-swap="innerHTML">
  </tbody>
</table>

<h3>Add Entry</h3>
<form hx-post="/htmx/sample-sheet/add"
      hx-target="#sample-sheet-body"
      hx-swap="innerHTML"
      hx-on::after-request="htmx.trigger(document.body, 'configChanged')">
  <div class="grid">
    <label>Upstream barcode
      <select name="upstream" required>
        {% for i in range(1, 97) %}
        <option value="barcode{{ '%02d'|format(i) }}">nb{{ '%02d'|format(i) }}</option>
        {% endfor %}
      </select>
    </label>
    <label>Downstream barcode
      <select name="downstream" required>
        {% for i in range(1, 97) %}
        <option value="barcode{{ '%02d'|format(i) }}">nb{{ '%02d'|format(i) }}</option>
        {% endfor %}
      </select>
    </label>
    <label>Alias
      <input type="text" name="alias" required pattern="[A-Za-z0-9_\.\-]+" placeholder="target_name">
    </label>
  </div>
  <button type="submit">Add</button>
</form>

<h3>Import / Export</h3>
<div class="grid">
  <a href="/api/sample-sheet/export" role="button" class="outline">Export CSV</a>
</div>
{% endblock %}
```

**Step 2: Create barcodes template**

Create `bin/viz/templates/barcodes.html` with a table of barcode sequences (ONT built-in + custom), arrangement editor form, and distance heatmap placeholder div for D3.

**Step 3: Create targets template**

Create `bin/viz/templates/targets.html` with CRUD table for reference FASTAs, upload form, and auto-reference preview section.

**Step 4: Create assumptions template**

Create `bin/viz/templates/assumptions.html` with editable list of assumptions rendered from `sma_experiment.toml`.

**Step 5: Add HTMX fragment routes to app.py**

Add routes like `/htmx/sample-sheet/rows`, `/htmx/sample-sheet/add`, `/htmx/sample-sheet/delete/{idx}` that return HTML fragments for HTMX to swap in.

**Step 6: Verify all pages render with data**

Run app with a test experiment directory, verify each page shows data and CRUD works.

**Step 7: Commit**

```bash
git add bin/viz/templates/ bin/viz/app.py
git commit -m "feat: add HTMX-powered page templates for sample sheet, barcodes, targets, assumptions"
```

---

### Task 7: Construct Diagram (D3.js)

**Files:**
- Create: `bin/viz/static/construct.js`
- Create: `bin/viz/templates/construct.html`

**Step 1: Create construct page template**

Create `bin/viz/templates/construct.html`:
```html
{% extends "base.html" %}
{% block title %}Construct - SMA-seq{% endblock %}
{% block head %}
<script src="https://d3js.org/d3.v7.min.js"></script>
<script src="/static/construct.js" defer></script>
{% endblock %}
{% block content %}
<h1>Read Construct</h1>

<h2>Full Library Product</h2>
<div id="construct-diagram"></div>
<div id="region-detail"></div>

<h2>Truncation Classes</h2>
<div id="truncation-diagram"></div>

<h2>Barcode Assignment</h2>
<div id="assignment-tables"
     hx-get="/htmx/construct/assignments"
     hx-trigger="load"
     hx-swap="innerHTML">
</div>

<h2>Truncation Rules</h2>
<div id="truncation-rules"
     hx-get="/htmx/construct/rules"
     hx-trigger="load"
     hx-swap="innerHTML">
</div>

<h2>Auto-Generated References</h2>
<button hx-post="/api/construct/auto-refs"
        hx-target="#auto-ref-result"
        hx-swap="innerHTML">
  Generate Truncated References
</button>
<div id="auto-ref-result"></div>

<script>
  // Pass construct data to D3
  fetch('/api/construct')
    .then(r => r.json())
    .then(data => {
      renderConstruct(data, '#construct-diagram', '#region-detail');
      renderTruncationLadder(data, '#truncation-diagram');
    });
</script>
{% endblock %}
```

**Step 2: Create D3 construct diagram**

Create `bin/viz/static/construct.js`:
```javascript
/**
 * SMA-seq Construct Diagram
 *
 * Renders an annotated linear molecule diagram showing:
 * - 5' adapter, upstream barcode, target insert, RC downstream barcode, 3' RC adapter
 * - Color-coded regions with labels
 * - Click to inspect region details
 *
 * Also renders the truncation class ladder showing 6 truncation classes.
 */

const COLORS = {
  adapter: '#607d8b',
  barcode: '#2196f3',
  target: '#4caf50',
  barcode_rc: '#ff9800',
  adapter_rc: '#9e9e9e',
};

const TRUNC_COLORS = {
  full: '#4caf50',
  trunc_3prime: '#8bc34a',
  trunc_target: '#ff9800',
  trunc_barcode: '#ff5722',
  adapter_only: '#9e9e9e',
  chimeric: '#f44336',
};

function renderConstruct(data, containerSel, detailSel) {
  const container = d3.select(containerSel);
  container.selectAll('*').remove();

  const width = 800;
  const height = 100;
  const margin = { top: 20, right: 20, bottom: 40, left: 20 };

  const svg = container.append('svg')
    .attr('viewBox', `0 0 ${width} ${height + margin.top + margin.bottom}`)
    .attr('width', '100%');

  const g = svg.append('g')
    .attr('transform', `translate(${margin.left}, ${margin.top})`);

  const regions = [
    { id: 'adapter_5', label: "5' Adapter", color: COLORS.adapter, width: 0.1 },
    { id: 'barcode1', label: 'Barcode 1', color: COLORS.barcode, width: 0.1 },
    { id: 'target', label: 'Target Insert', color: COLORS.target, width: 0.5 },
    { id: 'rc_barcode2', label: 'RC Barcode 2', color: COLORS.barcode_rc, width: 0.1 },
    { id: 'rc_adapter_3', label: "RC 3' Adapter", color: COLORS.adapter_rc, width: 0.1 },
  ];

  // Normalize widths
  const totalW = width - margin.left - margin.right;
  let x = 0;

  regions.forEach(r => {
    const w = r.width * totalW;

    g.append('rect')
      .attr('class', 'construct-region')
      .attr('x', x)
      .attr('y', 20)
      .attr('width', w)
      .attr('height', 40)
      .attr('fill', r.color)
      .attr('rx', 3)
      .on('click', () => showDetail(data, r, detailSel));

    g.append('text')
      .attr('x', x + w / 2)
      .attr('y', 45)
      .attr('text-anchor', 'middle')
      .attr('fill', 'white')
      .attr('font-size', '11px')
      .text(r.label);

    x += w + 4;
  });

  // Direction arrow
  svg.append('text')
    .attr('x', width / 2)
    .attr('y', height + margin.top + 15)
    .attr('text-anchor', 'middle')
    .attr('font-size', '12px')
    .text("5' ──────────────────────────────────→ 3'");
}

function renderTruncationLadder(data, containerSel) {
  const container = d3.select(containerSel);
  container.selectAll('*').remove();

  const width = 800;
  const rowH = 35;
  const margin = { top: 10, left: 20, right: 200 };

  const classes = [
    { id: 'full', label: 'Full', regions: [1, 1, 1, 1, 1] },
    { id: 'trunc_3prime', label: 'Trunc-3prime', regions: [0, 1, 1, 1, 0] },
    { id: 'trunc_target', label: 'Trunc-target', regions: [1, 1, 0.6, 0, 0] },
    { id: 'trunc_barcode', label: 'Trunc-barcode', regions: [1, 1, 0, 0, 0] },
    { id: 'adapter_only', label: 'Adapter-only', regions: [1, 0, 0, 0, 0] },
    { id: 'chimeric', label: 'Chimeric', regions: [0.5, 0.5, 0.3, 0.5, 0] },
  ];

  const regionColors = [COLORS.adapter, COLORS.barcode, COLORS.target, COLORS.barcode_rc, COLORS.adapter_rc];
  const totalW = width - margin.left - margin.right;
  const segW = totalW / 5;

  const svg = container.append('svg')
    .attr('viewBox', `0 0 ${width} ${classes.length * rowH + margin.top + 10}`)
    .attr('width', '100%');

  classes.forEach((cls, i) => {
    const y = margin.top + i * rowH;
    let x = margin.left;

    // Class label
    svg.append('text')
      .attr('x', width - margin.right + 10)
      .attr('y', y + 20)
      .attr('font-size', '12px')
      .attr('fill', TRUNC_COLORS[cls.id])
      .text(cls.label);

    // Draw regions
    cls.regions.forEach((frac, j) => {
      if (frac > 0) {
        const w = segW * frac;
        svg.append('rect')
          .attr('x', x)
          .attr('y', y + 5)
          .attr('width', w - 2)
          .attr('height', 22)
          .attr('fill', regionColors[j])
          .attr('opacity', frac < 1 ? 0.5 : 1)
          .attr('rx', 2);
      }
      x += segW;
    });

    // Assignment rule
    const rules = data.truncation?.rules || {};
    const rule = rules[cls.id] || '';
    svg.append('text')
      .attr('x', width - margin.right + 10)
      .attr('y', y + 32)
      .attr('font-size', '9px')
      .attr('fill', '#666')
      .text(rule);
  });
}

function showDetail(data, region, detailSel) {
  const detail = d3.select(detailSel);
  detail.selectAll('*').remove();

  const construct = data.construct || {};
  let info = '';

  switch (region.id) {
    case 'adapter_5':
      info = `<strong>5' Adapter:</strong> ${construct.adapter_5prime || 'Not set'}`;
      break;
    case 'barcode1':
      info = `<strong>Upstream Barcode:</strong> Forward orientation, classified in first ${data.classification?.barcode_search_window || 100}bp`;
      break;
    case 'target':
      info = `<strong>Target Insert:</strong> Type: ${construct.insert_type || 'amplicon'}. ${construct.notes || ''}`;
      break;
    case 'rc_barcode2':
      info = `<strong>Downstream Barcode (RC):</strong> Reverse complemented, classified in last ${data.classification?.barcode_search_window || 100}bp`;
      break;
    case 'rc_adapter_3':
      info = `<strong>3' Adapter (RC):</strong> ${construct.adapter_3prime || 'Not set'}`;
      break;
  }

  detail.html(`<article>${info}</article>`);
}
```

**Step 3: Add HTMX fragment routes for construct page**

Add `/htmx/construct/assignments` and `/htmx/construct/rules` routes to `app.py` that render HTML fragments showing the barcode assignment tables and truncation rule editor forms.

**Step 4: Verify construct page renders**

Run app, navigate to `/construct`, verify the D3 diagram renders with clickable regions and truncation ladder.

**Step 5: Commit**

```bash
git add bin/viz/static/construct.js bin/viz/templates/construct.html bin/viz/app.py
git commit -m "feat: add D3.js construct diagram with truncation ladder and region inspection"
```

---

### Task 8: Auto-Reference Generation

**Files:**
- Create: `tests/test_auto_refs.py`
- Modify: `bin/viz/api.py` (already has endpoint, just needs test coverage)

**Step 1: Write tests for auto-reference generation**

Create `tests/test_auto_refs.py`:
```python
"""Tests for auto-generated truncated reference sequences."""

from __future__ import annotations

from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from viz.app import app, set_experiment_dir


@pytest.fixture()
def experiment_with_adapter(tmp_path: Path) -> Path:
    exp = tmp_path / "experiment"
    exp.mkdir()
    (exp / "references").mkdir()

    (exp / "sample_sheet.csv").write_text(
        "flow_cell_id,kit,barcode,alias,type\n"
        "FAL12345,SQK-NBD114-96,barcode05--barcode10,target_A,test_sample\n"
    )
    (exp / "references" / "target_A.fasta").write_text(
        ">target_A\nACGTACGTACGTACGT\n"
    )
    # Set adapter sequences in config
    (exp / "sma_experiment.toml").write_text(
        '[experiment]\ndescription = "test"\n\n'
        '[construct]\nadapter_5prime = "AATTCCGG"\nadapter_3prime = "GGCCTTAA"\n'
        'insert_type = "amplicon"\nnotes = ""\n\n'
        '[demultiplexing]\nmode = "dual"\n\n'
        '[classification]\nbarcode_search_window = 100\n'
        'confidence_formula = "1.0 - (ed / barcode_length)"\n'
        'ambiguity_triggers_full_construct = true\n\n'
        '[quality]\nq_bc = "test"\nq_ld = "test"\n\n'
        '[truncation]\nmin_barcode_confidence = 0.5\n'
        'min_target_fraction = 0.1\nadapter_search_window = 50\n\n'
        '[truncation.rules]\nfull = "assign_to_target"\n'
        'trunc_3prime = "assign_to_target"\n'
        'trunc_target = "assign_by_start_barcode"\n'
        'trunc_barcode = "assign_by_start_barcode"\n'
        'adapter_only = "unclassified"\nchimeric = "flag_for_review"\n\n'
        '[truncation.auto_references]\nenabled = true\n'
    )
    return exp


@pytest.fixture()
def client(experiment_with_adapter: Path) -> TestClient:
    set_experiment_dir(experiment_with_adapter)
    return TestClient(app)


class TestAutoRefs:
    def test_generates_four_variants(self, client: TestClient, experiment_with_adapter: Path):
        resp = client.post("/api/construct/auto-refs")
        assert resp.status_code == 200
        data = resp.json()
        assert len(data["generated"]) == 4
        assert "target_A_full" in data["generated"]
        assert "target_A_no_adapter" in data["generated"]
        assert "target_A_no_end_bc" in data["generated"]
        assert "target_A_bc_only" in data["generated"]

    def test_auto_ref_files_exist(self, client: TestClient, experiment_with_adapter: Path):
        client.post("/api/construct/auto-refs")
        auto_dir = experiment_with_adapter / "references" / "auto"
        assert auto_dir.is_dir()
        fasta_files = list(auto_dir.glob("*.fasta"))
        assert len(fasta_files) == 4

    def test_full_variant_contains_all_components(self, client: TestClient, experiment_with_adapter: Path):
        client.post("/api/construct/auto-refs")
        full_fasta = experiment_with_adapter / "references" / "auto" / "target_A_full.fasta"
        content = full_fasta.read_text()
        # Should contain adapter + barcode + target + RC barcode + RC adapter
        assert ">target_A_full" in content
        assert "AATTCCGG" in content  # 5' adapter
        assert "ACGTACGTACGTACGT" in content  # target
```

**Step 2: Run tests**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_auto_refs.py -v`
Expected: All PASS

**Step 3: Commit**

```bash
git add tests/test_auto_refs.py
git commit -m "test: add tests for auto-generated truncated reference sequences"
```

---

### Task 9: Static HTML Export

**Files:**
- Create: `bin/viz/export.py`
- Create: `tests/test_export.py`

**Step 1: Write export tests**

Create `tests/test_export.py`:
```python
"""Tests for static HTML export."""

from __future__ import annotations

from pathlib import Path

import pytest

from viz.config_store import ConfigStore
from viz.export import export_experiment


@pytest.fixture()
def experiment_dir(tmp_path: Path) -> Path:
    exp = tmp_path / "experiment"
    exp.mkdir()
    (exp / "references").mkdir()
    (exp / "sample_sheet.csv").write_text(
        "flow_cell_id,kit,barcode,alias,type\n"
        "FAL12345,SQK-NBD114-96,barcode05--barcode10,CYP2D6_fwd,test_sample\n"
    )
    (exp / "references" / "CYP2D6_fwd.fasta").write_text(">CYP2D6_fwd\nACGT\n")
    return exp


class TestExport:
    def test_creates_output_directory(self, experiment_dir: Path, tmp_path: Path):
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        assert output_dir.is_dir()

    def test_creates_all_html_files(self, experiment_dir: Path, tmp_path: Path):
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        expected = [
            "index.html", "sample_sheet.html", "barcodes.html",
            "construct.html", "targets.html", "assumptions.html",
            "validation_report.html",
        ]
        for name in expected:
            assert (output_dir / name).exists(), f"Missing {name}"

    def test_html_is_self_contained(self, experiment_dir: Path, tmp_path: Path):
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        index = (output_dir / "index.html").read_text()
        # Should have inline styles, not external CSS links
        assert "<style>" in index
        # Should have data embedded
        assert "CYP2D6_fwd" in index

    def test_validation_report(self, experiment_dir: Path, tmp_path: Path):
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        report = (output_dir / "validation_report.html").read_text()
        assert "Validation" in report
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_export.py -v`
Expected: ImportError

**Step 3: Implement export**

Create `bin/viz/export.py`:
```python
"""Static HTML export for SMA-seq experiment configurations.

Generates self-contained HTML files that can be viewed without a server.
All CSS, JS, and data are inlined.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from viz.config_store import ConfigStore
from barcodes import BARCODES, reverse_complement


# The D3 construct.js code will be read and inlined
_JS_DIR = Path(__file__).resolve().parent / "static"
_CSS_DIR = Path(__file__).resolve().parent / "static"


def _inline_css() -> str:
    css_path = _CSS_DIR / "style.css"
    if css_path.exists():
        return css_path.read_text()
    return ""


def _inline_js() -> str:
    js_path = _JS_DIR / "construct.js"
    if js_path.exists():
        return js_path.read_text()
    return ""


def _base_html(title: str, body: str, data_json: str = "", nav_links: str = "") -> str:
    """Wrap body content in a self-contained HTML document."""
    css = _inline_css()
    js = _inline_js()
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
body {{ font-family: system-ui, -apple-system, sans-serif; margin: 0; padding: 2rem; color: #333; }}
table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; }}
th, td {{ border: 1px solid #ddd; padding: 0.5rem 0.75rem; text-align: left; }}
th {{ background: #f5f5f5; }}
h1, h2, h3 {{ color: #1a1a2e; }}
article {{ background: #f8f9fa; padding: 1rem; border-radius: 6px; margin: 1rem 0; }}
nav {{ margin-bottom: 2rem; padding-bottom: 1rem; border-bottom: 2px solid #eee; }}
nav a {{ margin-right: 1rem; text-decoration: none; color: #1976d2; }}
nav a:hover {{ text-decoration: underline; }}
.badge-ok {{ color: #4caf50; font-weight: bold; }}
.badge-error {{ color: #f44336; font-weight: bold; }}
{css}
@media print {{ body {{ padding: 1rem; }} nav {{ display: none; }} }}
</style>
</head>
<body>
<nav>{nav_links}</nav>
{body}
{f'<script>{js}</script>' if js else ''}
{f'<script>var DATA = {data_json};</script>' if data_json else ''}
</body>
</html>"""


def _nav_links() -> str:
    return (
        '<a href="index.html">Overview</a>'
        '<a href="sample_sheet.html">Sample Sheet</a>'
        '<a href="barcodes.html">Barcodes</a>'
        '<a href="construct.html">Construct</a>'
        '<a href="targets.html">Targets</a>'
        '<a href="assumptions.html">Assumptions</a>'
        '<a href="validation_report.html">Validation</a>'
    )


def export_experiment(
    experiment_dir: Path,
    output_dir: Path | None = None,
) -> Path:
    """Export an experiment configuration as self-contained static HTML.

    Parameters
    ----------
    experiment_dir : Path
        Path to experiment directory containing config files.
    output_dir : Path, optional
        Output directory. Defaults to exports/{exp_id}_{timestamp}/.

    Returns
    -------
    Path
        Path to the output directory.
    """
    store = ConfigStore(experiment_dir)
    cfg = store.experiment_config

    if output_dir is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = experiment_dir / "exports" / f"export_{ts}"

    output_dir.mkdir(parents=True, exist_ok=True)

    entries = store.read_sample_sheet()
    targets = store.read_targets()
    arrangement = store.read_arrangement()
    errors = store.validate()

    nav = _nav_links()

    # --- index.html ---
    body = f"<h1>SMA-seq Experiment Configuration</h1>"
    body += f"<article><strong>Description:</strong> {cfg.description or 'Not set'}</article>"
    body += f"<article><strong>Directory:</strong> {experiment_dir}</article>"
    body += f"<article><strong>Exported:</strong> {datetime.now().isoformat()}</article>"
    body += f"<h2>Summary</h2>"
    body += f"<p>Sample sheet entries: {len(entries)}</p>"
    body += f"<p>Target references: {len(targets)}</p>"
    body += f"<p>Demux mode: {cfg.demultiplexing.mode}</p>"
    body += f"<p>Validation: {'<span class=\"badge-ok\">Valid</span>' if not errors else f'<span class=\"badge-error\">{len(errors)} issues</span>'}</p>"

    # Construct diagram placeholder
    body += '<h2>Construct</h2>'
    body += '<div id="construct-diagram"></div>'
    body += '<div id="region-detail"></div>'
    body += '<h2>Truncation Classes</h2>'
    body += '<div id="truncation-diagram"></div>'

    construct_data = {
        "construct": cfg.construct.model_dump(),
        "truncation": cfg.truncation.model_dump(),
        "classification": cfg.classification.model_dump(),
    }
    body += f"""<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
document.addEventListener('DOMContentLoaded', function() {{
  var data = {json.dumps(construct_data)};
  renderConstruct(data, '#construct-diagram', '#region-detail');
  renderTruncationLadder(data, '#truncation-diagram');
}});
</script>"""

    (output_dir / "index.html").write_text(
        _base_html("SMA-seq Configuration Overview", body, nav_links=nav)
    )

    # --- sample_sheet.html ---
    body = "<h1>Sample Sheet</h1>"
    if entries:
        body += "<table><thead><tr><th>Barcode Pair</th><th>Upstream</th><th>Downstream</th><th>Alias</th><th>Kit</th><th>Type</th></tr></thead><tbody>"
        for e in entries:
            body += f"<tr><td>{e.barcode}</td><td>{e.upstream_barcode}</td><td>{e.downstream_barcode}</td><td>{e.alias}</td><td>{e.kit}</td><td>{e.type}</td></tr>"
        body += "</tbody></table>"
    else:
        body += "<p>No sample sheet entries.</p>"

    (output_dir / "sample_sheet.html").write_text(
        _base_html("Sample Sheet", body, nav_links=nav)
    )

    # --- barcodes.html ---
    body = "<h1>Barcodes</h1>"
    used_ids: set[str] = set()
    for e in entries:
        used_ids.add(e.upstream_barcode)
        used_ids.add(e.downstream_barcode)

    body += "<h2>Used Barcodes</h2>"
    body += "<table><thead><tr><th>ID</th><th>Sequence</th><th>RC</th><th>GC%</th></tr></thead><tbody>"
    for bc_id in sorted(used_ids):
        if bc_id in BARCODES:
            seq = BARCODES[bc_id]
            rc = reverse_complement(seq)
            gc = sum(1 for c in seq if c in "GC") / len(seq) * 100
            body += f"<tr><td>{bc_id}</td><td><code>{seq}</code></td><td><code>{rc}</code></td><td>{gc:.0f}%</td></tr>"
    body += "</tbody></table>"

    if arrangement:
        body += "<h2>Arrangement</h2>"
        body += f"<article><pre>{json.dumps(arrangement.model_dump(), indent=2)}</pre></article>"

    (output_dir / "barcodes.html").write_text(
        _base_html("Barcodes", body, nav_links=nav)
    )

    # --- construct.html ---
    body = "<h1>Read Construct</h1>"
    body += '<div id="construct-diagram"></div><div id="region-detail"></div>'
    body += "<h2>Truncation Classes</h2>"
    body += '<div id="truncation-diagram"></div>'

    body += "<h2>Barcode Assignment Mode</h2>"
    body += f"<article><strong>Mode:</strong> {cfg.demultiplexing.mode}</article>"

    if cfg.demultiplexing.pairs:
        body += "<h3>Pair Assignments</h3>"
        body += "<table><thead><tr><th>Start</th><th>End</th><th>Alias</th></tr></thead><tbody>"
        for p in cfg.demultiplexing.pairs:
            body += f"<tr><td>{p.start}</td><td>{p.end}</td><td>{p.alias}</td></tr>"
        body += "</tbody></table>"

    body += "<h2>Truncation Rules</h2>"
    body += "<table><thead><tr><th>Class</th><th>Rule</th></tr></thead><tbody>"
    rules = cfg.truncation.rules
    for field in ["full", "trunc_3prime", "trunc_target", "trunc_barcode", "adapter_only", "chimeric"]:
        body += f"<tr><td>{field}</td><td>{getattr(rules, field)}</td></tr>"
    body += "</tbody></table>"

    body += f"""<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
document.addEventListener('DOMContentLoaded', function() {{
  var data = {json.dumps(construct_data)};
  renderConstruct(data, '#construct-diagram', '#region-detail');
  renderTruncationLadder(data, '#truncation-diagram');
}});
</script>"""

    (output_dir / "construct.html").write_text(
        _base_html("Construct", body, nav_links=nav)
    )

    # --- targets.html ---
    body = "<h1>Target References</h1>"
    if targets:
        body += "<table><thead><tr><th>Target ID</th><th>Length</th><th>GC%</th></tr></thead><tbody>"
        for t in targets:
            body += f"<tr><td>{t.tgt_id}</td><td>{t.length} bp</td><td>{t.gc_content}%</td></tr>"
        body += "</tbody></table>"

        body += "<h2>Sequences</h2>"
        for t in targets:
            trunc_seq = t.sequence[:200] + ("..." if len(t.sequence) > 200 else "")
            body += f"<article><strong>{t.tgt_id}</strong> ({t.length} bp)<br><code style='word-break:break-all'>{trunc_seq}</code></article>"
    else:
        body += "<p>No target references found.</p>"

    (output_dir / "targets.html").write_text(
        _base_html("Targets", body, nav_links=nav)
    )

    # --- assumptions.html ---
    body = "<h1>Pipeline Assumptions</h1>"
    if cfg.assumptions:
        body += "<table><thead><tr><th>Key</th><th>Assumption</th><th>Rationale</th></tr></thead><tbody>"
        for a in cfg.assumptions:
            body += f"<tr><td><code>{a.key}</code></td><td>{a.text}</td><td>{a.why}</td></tr>"
        body += "</tbody></table>"
    else:
        body += "<p>No assumptions documented.</p>"

    body += "<h2>Quality Metrics</h2>"
    body += f"<article><strong>Q_bc:</strong> {cfg.quality.q_bc}</article>"
    body += f"<article><strong>Q_ld:</strong> {cfg.quality.q_ld}</article>"

    body += "<h2>Classification</h2>"
    body += f"<article><strong>Search window:</strong> {cfg.classification.barcode_search_window} bp</article>"
    body += f"<article><strong>Confidence formula:</strong> <code>{cfg.classification.confidence_formula}</code></article>"

    (output_dir / "assumptions.html").write_text(
        _base_html("Assumptions", body, nav_links=nav)
    )

    # --- validation_report.html ---
    body = "<h1>Validation Report</h1>"
    if errors:
        body += f'<p class="badge-error">{len(errors)} issue(s) found:</p><ul>'
        for e in errors:
            body += f"<li>{e}</li>"
        body += "</ul>"
    else:
        body += '<p class="badge-ok">All configurations are consistent.</p>'

    (output_dir / "validation_report.html").write_text(
        _base_html("Validation Report", body, nav_links=nav)
    )

    return output_dir


def main():
    """CLI entry point for export."""
    import sys

    if len(sys.argv) < 2:
        print("Usage: python -m bin.viz.export <experiment_dir> [--output <dir>]")
        sys.exit(1)

    experiment_dir = Path(sys.argv[1]).resolve()
    output_dir = None

    if "--output" in sys.argv:
        idx = sys.argv.index("--output")
        output_dir = Path(sys.argv[idx + 1]).resolve()

    result = export_experiment(experiment_dir, output_dir)
    print(f"[export] Static HTML exported to: {result}")


if __name__ == "__main__":
    main()
```

**Step 4: Wire export endpoint in api.py**

Update `/api/export` in `bin/viz/api.py`:
```python
@router.post("/export")
async def export_html() -> dict[str, str]:
    from viz.export import export_experiment
    store = get_store()
    output = export_experiment(store.dir)
    return {"status": "ok", "path": str(output)}
```

**Step 5: Run export tests**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_export.py -v`
Expected: All tests PASS

**Step 6: Commit**

```bash
git add bin/viz/export.py tests/test_export.py bin/viz/api.py
git commit -m "feat: add static HTML export with self-contained pages"
```

---

### Task 10: Integration Test & Polish

**Files:**
- Create: `tests/test_viz_integration.py`
- Modify: `bin/viz/app.py` (final wiring)

**Step 1: Write integration test**

Create `tests/test_viz_integration.py`:
```python
"""Integration test: full workflow from config creation to export."""

from __future__ import annotations

from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from viz.app import app, set_experiment_dir


@pytest.fixture()
def empty_experiment(tmp_path: Path) -> Path:
    exp = tmp_path / "new_experiment"
    exp.mkdir()
    return exp


@pytest.fixture()
def client(empty_experiment: Path) -> TestClient:
    set_experiment_dir(empty_experiment)
    return TestClient(app)


class TestFullWorkflow:
    """Test creating an experiment config from scratch and exporting."""

    def test_create_and_export(self, client: TestClient, empty_experiment: Path):
        # 1. Set experiment description
        resp = client.put("/api/experiment", json={
            "description": "Integration test experiment",
        })
        assert resp.status_code == 200

        # 2. Add sample sheet entries
        resp = client.put("/api/sample-sheet", json=[
            {
                "flow_cell_id": "FAL99999",
                "kit": "SQK-NBD114-96",
                "barcode": "barcode05--barcode10",
                "alias": "target_A",
                "type": "test_sample",
            },
        ])
        assert resp.status_code == 200

        # 3. Add target reference
        resp = client.post("/api/targets", json={
            "tgt_id": "target_A",
            "sequence": "ACGTACGTACGTACGTACGTACGT",
        })
        assert resp.status_code == 200

        # 4. Validate - should be clean
        resp = client.get("/api/validate")
        assert resp.status_code == 200
        assert resp.json()["valid"] is True

        # 5. Export
        resp = client.post("/api/export")
        assert resp.status_code == 200
        export_path = Path(resp.json()["path"])
        assert export_path.is_dir()
        assert (export_path / "index.html").exists()
        assert (export_path / "sample_sheet.html").exists()
        assert (export_path / "construct.html").exists()

        # 6. Verify exported HTML contains our data
        index = (export_path / "index.html").read_text()
        assert "Integration test experiment" in index
        assert "target_A" in index
```

**Step 2: Run integration tests**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_viz_integration.py -v`
Expected: All tests PASS

**Step 3: Run full test suite**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/ -v`
Expected: All tests PASS (existing + new)

**Step 4: Commit**

```bash
git add tests/test_viz_integration.py
git commit -m "test: add integration test for full create-validate-export workflow"
```

---

### Task 11: MinKNOW Output Directory Import

**Files:**
- Modify: `bin/viz/config_store.py`
- Create: `tests/test_minknow_import.py`

**Step 1: Write import test**

Create `tests/test_minknow_import.py`:
```python
"""Tests for importing config from MinKNOW output directory."""

from __future__ import annotations

from pathlib import Path

import pytest

from viz.config_store import ConfigStore


@pytest.fixture()
def minknow_output(tmp_path: Path) -> Path:
    """Create a fake MinKNOW output directory structure."""
    run_dir = tmp_path / "data" / "experiment_group" / "sample1" / "20260214_1200_MN12345_FAL12345_abcdef01"
    run_dir.mkdir(parents=True)

    (run_dir / "sample_sheet_FAL12345_20260214_1200_abcdef01.csv").write_text(
        "protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,flow_cell_product_code,kit,barcode,alias,type\n"
        "abc-def,MN12345,FAL12345,sample1,exp_group,FLO-PRO114,SQK-NBD114-96,barcode05--barcode10,target_A,test_sample\n"
    )

    (run_dir / "final_summary_FAL12345_abcdef01_12345678.txt").write_text(
        "protocol_run_id=abc-def\nflow_cell_id=FAL12345\nsample_id=sample1\n"
    )

    return run_dir


class TestMinKNOWImport:
    def test_import_from_minknow_dir(self, minknow_output: Path, tmp_path: Path):
        exp_dir = tmp_path / "imported_experiment"
        exp_dir.mkdir()
        store = ConfigStore(exp_dir)
        store.import_from_minknow(minknow_output)

        entries = store.read_sample_sheet()
        assert len(entries) == 1
        assert entries[0].alias == "target_A"
```

**Step 2: Run tests to verify they fail**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_minknow_import.py -v`
Expected: AttributeError - `import_from_minknow` not found

**Step 3: Implement MinKNOW import**

Add to `bin/viz/config_store.py`:
```python
def import_from_minknow(self, minknow_dir: Path) -> None:
    """Import config files from a MinKNOW output directory.

    Discovers and copies:
    - sample_sheet_*.csv -> sample_sheet.csv
    - final_summary_*.txt -> parses for metadata
    """
    minknow_dir = Path(minknow_dir).resolve()

    # Find sample sheet
    ss_files = list(minknow_dir.glob("sample_sheet_*.csv"))
    if ss_files:
        import shutil
        dest = self.dir / "sample_sheet.csv"
        shutil.copy2(ss_files[0], dest)

    # Find final summary for metadata
    fs_files = list(minknow_dir.glob("final_summary_*.txt"))
    if fs_files:
        metadata: dict[str, str] = {}
        with open(fs_files[0]) as f:
            for line in f:
                if "=" in line:
                    key, _, val = line.strip().partition("=")
                    metadata[key.strip()] = val.strip()
        if "sample_id" in metadata:
            self.experiment_config.description = (
                f"Imported from MinKNOW: {metadata.get('sample_id', '')}"
            )
            self.save_experiment_config()
```

**Step 4: Run tests**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/test_minknow_import.py -v`
Expected: All PASS

**Step 5: Add import endpoint to API**

Add to `bin/viz/api.py`:
```python
@router.post("/import/minknow")
async def import_minknow(data: dict[str, str]) -> dict[str, str]:
    from pathlib import Path
    store = get_store()
    minknow_path = Path(data["path"])
    if not minknow_path.is_dir():
        raise HTTPException(400, f"Not a directory: {minknow_path}")
    store.import_from_minknow(minknow_path)
    return {"status": "ok"}
```

**Step 6: Commit**

```bash
git add bin/viz/config_store.py tests/test_minknow_import.py bin/viz/api.py
git commit -m "feat: add MinKNOW output directory import"
```

---

### Task 12: Final Review & Documentation

**Step 1: Run full test suite**

Run: `cd /tmp/ont-sma-seq && PYTHONPATH=bin python -m pytest tests/ -v --tb=short`
Expected: All tests PASS

**Step 2: Verify app launches and all pages render**

Run: `cd /tmp/ont-sma-seq && python -m bin.viz.app /tmp/test-experiment`
Manually check: dashboard, sample sheet, barcodes, construct, targets, assumptions.

**Step 3: Verify static export works end-to-end**

Run: `cd /tmp/ont-sma-seq && python -m bin.viz.export /tmp/test-experiment --output /tmp/export-test`
Open `/tmp/export-test/index.html` in browser, verify all links work.

**Step 4: Commit any final fixes**

```bash
git add -A
git commit -m "chore: final polish and fixes for config visualizer"
```
