"""Validate and group MinKNOW runs for merging."""
from __future__ import annotations

from collections import defaultdict

from sma_merge.models import RunInfo, RunGroup


def validate_runs(runs: list[RunInfo]) -> list[RunGroup]:
    """Group runs by flow_cell_id and validate consistency within each group."""
    by_flowcell: dict[str, list[RunInfo]] = defaultdict(list)
    for run in runs:
        by_flowcell[run.flow_cell_id].append(run)

    groups: list[RunGroup] = []
    for fc_id, fc_runs in sorted(by_flowcell.items()):
        issues: list[str] = []

        models = {r.basecall_model for r in fc_runs}
        if len(models) > 1:
            issues.append(f"Inconsistent basecall models: {models}")

        protocols = {r.protocol_group_id for r in fc_runs}
        if len(protocols) > 1:
            issues.append(f"Multiple protocol groups on same flowcell: {protocols}")

        if all(r.sample_id == "" for r in fc_runs):
            issues.append("Warning: sample_id is empty for all runs (common with MinKNOW)")

        groups.append(RunGroup(
            flow_cell_id=fc_id,
            runs=fc_runs,
            basecall_model=fc_runs[0].basecall_model,
            is_consistent=not any(
                issue for issue in issues if not issue.startswith("Warning:")
            ),
            issues=issues,
        ))

    return groups


def format_validation(groups: list[RunGroup]) -> str:
    """Format validation results for display."""
    lines: list[str] = []
    for g in groups:
        status = "OK" if g.is_consistent else "FAIL"
        lines.append(f"[{status}] FlowCell {g.flow_cell_id}: {len(g.runs)} run(s), model={g.basecall_model}")
        for issue in g.issues:
            prefix = "  WARN:" if issue.startswith("Warning:") else "  ERROR:"
            lines.append(f"{prefix} {issue}")
    return "\n".join(lines)
