"""CLI entry point for the SMA-seq barcode calibration pipeline."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from calibrate.discover import discover_runs
from calibrate.merge import (
    ValidationError,
    find_bam_files,
    group_runs,
    merge_bams,
    validate_basecalling,
)


def _cmd_discover(args: argparse.Namespace) -> int:
    """List discovered MinKNOW runs under root_dir."""
    root = Path(args.root_dir)
    if not root.is_dir():
        print(f"Error: {root} is not a directory", file=sys.stderr)
        return 1

    runs = discover_runs(root)
    if not runs:
        print("No runs found.")
        return 0

    print(f"Discovered {len(runs)} run(s):\n")
    for run in runs:
        print(f"  Flow cell: {run.flow_cell_id}")
        print(f"  Device:    {run.device_id}")
        print(f"  Sample:    {run.sample_id}")
        print(f"  Started:   {run.started}")
        print(f"  Directory: {run.run_dir}")
        bams = find_bam_files(run.run_dir)
        print(f"  BAM files: {len(bams)}")
        print()

    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    """Validate basecalling consistency across grouped runs."""
    root = Path(args.root_dir)
    if not root.is_dir():
        print(f"Error: {root} is not a directory", file=sys.stderr)
        return 1

    runs = discover_runs(root)
    if not runs:
        print("No runs found.")
        return 0

    groups = group_runs(runs)
    print(f"Found {len(groups)} merge group(s):\n")

    all_ok = True
    for i, group in enumerate(groups, 1):
        print(f"  Group {i}: FC={group.flow_cell_id}  sample={group.sample_id}  runs={len(group.runs)}")
        all_bams: list[Path] = []
        for run in group.runs:
            all_bams.extend(find_bam_files(run.run_dir))

        if not all_bams:
            print("    No BAM files found in this group.")
            continue

        try:
            model = validate_basecalling(all_bams)
            print(f"    Basecalling model: {model or '(unknown)'}")
            print(f"    BAM files: {len(all_bams)}")
            print("    Status: OK")
        except ValidationError as exc:
            print(f"    VALIDATION ERROR: {exc}")
            all_ok = False
        print()

    if all_ok:
        print("All groups validated successfully.")
    else:
        print("Some groups have validation errors.", file=sys.stderr)
    return 0 if all_ok else 1


def _cmd_ingest(args: argparse.Namespace) -> int:
    """Full pipeline: discover -> validate -> merge."""
    root = Path(args.root_dir)
    output_dir = Path(args.output_dir)

    if not root.is_dir():
        print(f"Error: {root} is not a directory", file=sys.stderr)
        return 1

    runs = discover_runs(root)
    if not runs:
        print("No runs found.")
        return 0

    groups = group_runs(runs)
    print(f"Processing {len(groups)} merge group(s)...")

    output_dir.mkdir(parents=True, exist_ok=True)

    for i, group in enumerate(groups, 1):
        label = f"{group.flow_cell_id}_{group.sample_id}"
        print(f"\n--- Group {i}: {label} ({len(group.runs)} run(s)) ---")

        # Collect BAMs
        all_bams: list[Path] = []
        for run in group.runs:
            all_bams.extend(find_bam_files(run.run_dir))

        if not all_bams:
            print("  Skipping: no BAM files found.")
            continue

        # Validate
        try:
            model = validate_basecalling(all_bams)
            print(f"  Basecalling model: {model or '(unknown)'}")
        except ValidationError as exc:
            print(f"  VALIDATION ERROR: {exc}", file=sys.stderr)
            return 1

        # Merge
        merged_path = output_dir / f"{label}_merged.bam"
        print(f"  Merging {len(all_bams)} BAM(s) -> {merged_path}")
        prov = merge_bams(group, merged_path)
        print(f"  Merge complete: {prov['num_bams']} BAMs from {prov['num_runs']} runs")

    print("\nPipeline complete.")
    return 0


def _cmd_viz(args: argparse.Namespace) -> int:
    """Launch visualization app (placeholder)."""
    print("Launching visualization...")
    return 0


def _cmd_export(args: argparse.Namespace) -> int:
    """Export results (placeholder)."""
    print("Exporting results...")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser with all subcommands."""
    parser = argparse.ArgumentParser(
        prog="calibrate",
        description="SMA-seq barcode calibration pipeline",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # discover
    p_discover = subparsers.add_parser(
        "discover",
        help="Discover MinKNOW output directories",
    )
    p_discover.add_argument("root_dir", help="Root directory to search for runs")
    p_discover.set_defaults(func=_cmd_discover)

    # validate
    p_validate = subparsers.add_parser(
        "validate",
        help="Validate basecalling consistency across runs",
    )
    p_validate.add_argument("root_dir", help="Root directory to search for runs")
    p_validate.set_defaults(func=_cmd_validate)

    # ingest
    p_ingest = subparsers.add_parser(
        "ingest",
        help="Full pipeline: discover, validate, merge",
    )
    p_ingest.add_argument("root_dir", help="Root directory to search for runs")
    p_ingest.add_argument(
        "-c", "--construct", dest="construct_toml",
        help="Path to construct TOML config",
    )
    p_ingest.add_argument(
        "-o", "--output", dest="output_dir", required=True,
        help="Output directory for merged BAMs",
    )
    p_ingest.set_defaults(func=_cmd_ingest)

    # viz
    p_viz = subparsers.add_parser(
        "viz",
        help="Launch visualization app",
    )
    p_viz.add_argument("db_path", nargs="+", help="Path(s) to database file(s)")
    p_viz.set_defaults(func=_cmd_viz)

    # export
    p_export = subparsers.add_parser(
        "export",
        help="Export results",
    )
    p_export.add_argument("db_path", help="Path to database file")
    p_export.add_argument("--output", dest="output_dir", help="Output directory")
    p_export.add_argument("--compare", action="store_true", help="Enable comparison mode")
    p_export.set_defaults(func=_cmd_export)

    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point for the calibrate CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
