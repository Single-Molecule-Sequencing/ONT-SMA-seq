"""Discover MinKNOW run directories and extract metadata from POD5 files."""
from __future__ import annotations

from pathlib import Path

import pod5 as p5

from sma_merge.models import RunInfo


def discover_runs(experiment_path: Path) -> list[RunInfo]:
    """Scan experiment_path for MinKNOW run directories.

    Looks for directories containing pod5_pass/ with at least one .pod5 file.
    Reads run metadata from the first POD5 file found in each run.
    """
    runs: list[RunInfo] = []

    for pod5_pass in sorted(experiment_path.rglob("pod5_pass")):
        if not pod5_pass.is_dir():
            continue

        run_dir = pod5_pass.parent
        pod5_files = sorted(pod5_pass.rglob("*.pod5"))
        if not pod5_files:
            continue

        with p5.Reader(pod5_files[0]) as reader:
            for read in reader.reads():
                ri = read.run_info
                ctx = ri.context_tags
                trk = ri.tracking_id
                runs.append(RunInfo(
                    run_dir=run_dir,
                    flow_cell_id=trk.get("flow_cell_id", ""),
                    device_id=trk.get("device_id", ""),
                    protocol_group_id=trk.get("protocol_group_id", ""),
                    basecall_model=ctx.get("basecall_model_simplex", ""),
                    sample_id=trk.get("sample_id", ""),
                    run_id=trk.get("run_id", ""),
                    sample_rate=ri.sample_rate,
                    pod5_dir=pod5_pass,
                    pod5_count=len(pod5_files),
                    mod_base_models=ctx.get("basecall_models_modified", ""),
                ))
                break

    return runs


def format_discovery_table(runs: list[RunInfo]) -> str:
    """Format discovered runs as a human-readable table."""
    if not runs:
        return "No runs found."

    lines = [
        f"{'FlowCell':<12} {'Device':<12} {'Protocol Group':<35} {'Model':<20} {'POD5s':>6}",
        "-" * 90,
    ]
    for r in runs:
        model_short = r.basecall_model.split("@")[-1] if "@" in r.basecall_model else r.basecall_model
        lines.append(
            f"{r.flow_cell_id:<12} {r.device_id:<12} {r.protocol_group_id:<35} {model_short:<20} {r.pod5_count:>6}"
        )
    return "\n".join(lines)
