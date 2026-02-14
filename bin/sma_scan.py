#!/usr/bin/env python3
"""sma_scan.py — Scan experiment directory, discover MinKNOW runs, output manifest."""

import csv
import json as _json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


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


# ---------------------------------------------------------------------------
# Run discovery
# ---------------------------------------------------------------------------


def discover_runs(experiment_dir: Path) -> list[dict]:
    """Recursively discover MinKNOW run folders within an experiment directory.

    A run folder is identified by containing pod5*/ or pod5_pass/ directories,
    or final_summary*.txt files. Returns list of run info dicts.
    """
    experiment_dir = Path(experiment_dir)
    runs = []
    seen: set[str] = set()

    for candidate in sorted(experiment_dir.rglob("*")):
        if not candidate.is_dir():
            continue
        # Skip directories we already identified as children of a run
        if any(candidate.is_relative_to(Path(s)) for s in seen):
            continue

        try:
            children = list(candidate.iterdir())
        except PermissionError:
            continue

        has_pod5 = any(
            d.is_dir() and d.name in ("pod5", "pod5_pass")
            for d in children
        )
        has_final_summary = any(
            f.name.startswith("final_summary") and f.name.endswith(".txt")
            for f in children if f.is_file()
        )
        if not (has_pod5 or has_final_summary):
            continue

        seen.add(str(candidate))
        run_info = _extract_run_info(candidate, experiment_dir)
        runs.append(run_info)

    return runs


def _extract_run_info(run_dir: Path, experiment_dir: Path) -> dict:
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
        if f.name.startswith("final_summary") and f.name.endswith(".txt") and f.is_file():
            try:
                info["final_summary"] = str(f.relative_to(experiment_dir))
            except ValueError:
                info["final_summary"] = f.name
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
        if f.name.startswith("sample_sheet") and f.name.endswith(".csv") and f.is_file():
            try:
                info["sample_sheet"] = str(f.relative_to(experiment_dir))
            except ValueError:
                info["sample_sheet"] = f.name
            ss_data = parse_minknow_sample_sheet(f)
            info["experiment_id"] = ss_data.get("experiment_id")
            info["kit"] = ss_data.get("kit")
            info["flow_cell_product_code"] = ss_data.get("flow_cell_product_code")
            info["barcode_aliases"] = ss_data.get("barcode_aliases", {})
            break

    # --- report.json ---
    for f in sorted(run_dir.iterdir()):
        if f.name.startswith("report_") and f.name.endswith(".json") and f.is_file():
            try:
                info["report_json"] = str(f.relative_to(experiment_dir))
            except ValueError:
                info["report_json"] = f.name
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
            info["bam_count"] = sum(1 for f in bam_dir.rglob("*.bam"))

    # --- sequencing_summary ---
    for f in sorted(run_dir.iterdir()):
        if f.name.startswith("sequencing_summary") and f.name.endswith(".txt") and f.is_file():
            try:
                info["sequencing_summary"] = str(f.relative_to(experiment_dir))
            except ValueError:
                info["sequencing_summary"] = f.name
            break

    # --- Extract flow_cell_id from dir name if not in final_summary ---
    if info["flow_cell_id"] is None and len(parts) >= 4:
        info["flow_cell_id"] = parts[3]
    if info["instrument"] is None and len(parts) >= 3:
        info["instrument"] = parts[2]

    return info


# ---------------------------------------------------------------------------
# Merge decisions
# ---------------------------------------------------------------------------


def _parse_timestamp(ts: str | None) -> datetime | None:
    """Parse an ISO timestamp string, handling timezone offsets."""
    if not ts:
        return None
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
                # Handle timezone-aware vs naive comparison
                if prev_stop.tzinfo and not curr_start.tzinfo:
                    prev_stop = prev_stop.replace(tzinfo=None)
                elif curr_start.tzinfo and not prev_stop.tzinfo:
                    curr_start = curr_start.replace(tzinfo=None)
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


# ---------------------------------------------------------------------------
# Manifest builder
# ---------------------------------------------------------------------------


def build_manifest(experiment_dir: Path, runs: list[dict],
                   max_gap_hours: float = 24.0) -> dict:
    """Build the complete manifest JSON from discovered runs."""
    groups = decide_merges(runs, max_gap_hours)

    # Aggregate experiment-level info from first run
    first_run = runs[0] if runs else {}
    barcode_aliases: dict[str, str] = {}
    seq_summary_paths: list[str] = []
    pod5_sources: list[str] = []
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
    import sys

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
    main()
