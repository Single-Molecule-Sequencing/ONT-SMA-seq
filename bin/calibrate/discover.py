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
    """Parse a final_summary_*.txt key=value file into RunInfo."""
    kv: dict[str, str] = {}
    for line in path.read_text().strip().splitlines():
        if "=" in line:
            key, _, value = line.partition("=")
            kv[key.strip()] = value.strip()

    run_dir = path.parent
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
    Returns all discovered runs sorted by start time.
    """
    runs: list[RunInfo] = []
    for fs_path in sorted(root.rglob("final_summary_*.txt")):
        runs.append(parse_final_summary(fs_path))
    runs.sort(key=lambda r: r.started)
    return runs
