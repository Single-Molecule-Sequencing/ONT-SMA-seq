"""Data models for sma-merge."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class RunInfo:
    """Metadata for a single MinKNOW run directory."""

    run_dir: Path
    flow_cell_id: str
    device_id: str
    protocol_group_id: str
    basecall_model: str
    sample_id: str
    run_id: str
    sample_rate: int
    pod5_dir: Path
    pod5_count: int
    mod_base_models: str


@dataclass
class RunGroup:
    """A group of runs that should be merged (same flowcell/sample)."""

    flow_cell_id: str
    runs: list[RunInfo]
    basecall_model: str
    is_consistent: bool
    issues: list[str] = field(default_factory=list)


@dataclass
class MergeResult:
    """Result of a merge or subsample operation."""

    merged_pod5: Path
    output_bam: Path
    total_reads: int
    reads_tagged: int
