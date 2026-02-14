"""Run grouping and BAM merge logic for barcode calibration."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pysam

from calibrate.discover import RunInfo


class ValidationError(Exception):
    """Raised when basecalling validation fails."""


@dataclass
class MergeGroup:
    """A group of runs sharing the same flow cell and sample."""

    flow_cell_id: str
    sample_id: str
    runs: list[RunInfo] = field(default_factory=list)


def _parse_started(ts: str) -> datetime:
    """Parse an ISO timestamp, handling trailing Z."""
    return datetime.fromisoformat(ts.replace("Z", "+00:00"))


def group_runs(runs: list[RunInfo], max_gap_hours: float = 24.0) -> list[MergeGroup]:
    """Group RunInfo objects by (flow_cell_id, sample_id).

    Within each key, runs are sorted by start time and split into
    separate groups when the gap between consecutive runs exceeds
    *max_gap_hours*.
    """
    buckets: dict[tuple[str, str], list[RunInfo]] = defaultdict(list)
    for run in runs:
        buckets[(run.flow_cell_id, run.sample_id)].append(run)

    groups: list[MergeGroup] = []
    for (fc, sid), bucket_runs in buckets.items():
        bucket_runs.sort(key=lambda r: _parse_started(r.started))
        current_group = MergeGroup(flow_cell_id=fc, sample_id=sid, runs=[bucket_runs[0]])
        for prev, cur in zip(bucket_runs, bucket_runs[1:]):
            gap = (_parse_started(cur.started) - _parse_started(prev.started)).total_seconds() / 3600
            if gap <= max_gap_hours:
                current_group.runs.append(cur)
            else:
                groups.append(current_group)
                current_group = MergeGroup(flow_cell_id=fc, sample_id=sid, runs=[cur])
        groups.append(current_group)

    return groups


def find_bam_files(run_dir: Path) -> list[Path]:
    """Find all .bam files recursively under run_dir/bam_pass/.

    Returns a sorted list of Path objects (empty if bam_pass/ does
    not exist or contains no BAMs).
    """
    bam_pass = run_dir / "bam_pass"
    if not bam_pass.is_dir():
        return []
    return sorted(bam_pass.rglob("*.bam"))


def validate_basecalling(bam_paths: list[Path]) -> str:
    """Validate that all BAMs used the same basecalling model.

    Opens each BAM, reads @PG header lines, and extracts the
    basecalling command line.

    Raises:
        ValidationError: if different models are detected across BAMs
            or if ``--trim`` is used without ``--no-trim``.

    Returns:
        The common basecalling model string.
    """
    models: set[str] = set()

    for bam_path in bam_paths:
        with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as af:
            header = af.header.to_dict()
            for pg in header.get("PG", []):
                cl = pg.get("CL", "")
                if not cl:
                    continue
                # Check for --trim without --no-trim
                parts = cl.split()
                if "--trim" in parts and "--no-trim" not in parts:
                    raise ValidationError(
                        f"BAM {bam_path.name} uses --trim without --no-trim"
                    )
                # Extract model: look for --model or -m flag, or
                # basecaller model token (e.g. dna_r10.4.1_...)
                model = _extract_model(parts)
                if model:
                    models.add(model)

    if len(models) > 1:
        raise ValidationError(
            f"Multiple basecalling models detected: {sorted(models)}"
        )

    return models.pop() if models else ""


def _extract_model(parts: list[str]) -> str | None:
    """Extract the basecalling model from a command-line token list."""
    for i, token in enumerate(parts):
        if token in ("--model", "-m") and i + 1 < len(parts):
            return parts[i + 1]
    # Fallback: look for a token matching common ONT model patterns
    for token in parts:
        if token.startswith("dna_r") or token.startswith("rna_"):
            return token
    return None


def merge_bams(group: MergeGroup, output_path: Path) -> dict[str, Any]:
    """Merge all BAMs from runs in *group* into *output_path*.

    Returns a provenance dict with merge metadata.
    """
    all_bams: list[Path] = []
    for run in group.runs:
        all_bams.extend(find_bam_files(run.run_dir))

    output_path.parent.mkdir(parents=True, exist_ok=True)

    pysam.merge("-f", str(output_path), *[str(b) for b in all_bams])

    return {
        "flow_cell_id": group.flow_cell_id,
        "sample_id": group.sample_id,
        "num_runs": len(group.runs),
        "num_bams": len(all_bams),
        "output": str(output_path),
        "source_bams": [str(b) for b in all_bams],
    }
