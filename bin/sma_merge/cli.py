"""CLI entry point for sma-merge."""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

from sma_merge.discover import discover_runs, format_discovery_table
from sma_merge.validate import validate_runs, format_validation


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sma-merge",
        description="SMA-seq experiment merge tool: discover, validate, merge, or subsample MinKNOW runs.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # discover
    disc = sub.add_parser("discover", help="Discover MinKNOW runs and show metadata")
    disc.add_argument("path", type=Path, help="Path to experiment directory")

    # merge
    merge = sub.add_parser("merge", help="Full merge: combine POD5s + re-basecall")
    merge.add_argument("path", type=Path, help="Path to experiment directory")
    merge.add_argument("-o", "--output", type=Path, required=True, help="Output directory")
    merge.add_argument("--device", default="auto", help="Dorado device (default: auto)")
    merge.add_argument("--dorado", type=Path, default=None, help="Path to dorado binary")

    # subsample
    sample = sub.add_parser("subsample", help="Subsample N reads + basecall")
    sample.add_argument("path", type=Path, help="Path to experiment directory")
    sample.add_argument("-n", "--num-reads", type=int, required=True, help="Number of reads to subsample")
    sample.add_argument("-o", "--output", type=Path, required=True, help="Output directory")
    sample.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    sample.add_argument("--device", default="auto", help="Dorado device (default: auto)")
    sample.add_argument("--dorado", type=Path, default=None, help="Path to dorado binary")

    return parser


def cmd_discover(experiment_path: Path) -> None:
    """Run the discover subcommand."""
    print(f"Scanning {experiment_path} for MinKNOW runs...\n")
    runs = discover_runs(experiment_path)

    if not runs:
        print("No runs found.")
        return

    print(format_discovery_table(runs))
    print()

    groups = validate_runs(runs)
    print(format_validation(groups))


def cmd_merge(experiment_path: Path, output_dir: Path, device: str, dorado_path: Path | None) -> None:
    """Run the full merge pipeline."""
    from sma_merge.basecall import merge_pod5s, basecall
    from sma_merge.tag import build_pod5_lookup, tag_bam

    print(f"Scanning {experiment_path}...")
    runs = discover_runs(experiment_path)
    if not runs:
        print("No runs found.")
        sys.exit(1)

    groups = validate_runs(runs)
    print(format_validation(groups))

    for g in groups:
        if not g.is_consistent:
            print(f"\nERROR: Inconsistent runs in flowcell {g.flow_cell_id}. Aborting.")
            sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    for g in groups:
        name = g.flow_cell_id
        pod5_dirs = [r.pod5_dir for r in g.runs]
        model = g.basecall_model

        merged_pod5 = output_dir / f"{name}_merged.pod5"
        print(f"\nMerging POD5 files -> {merged_pod5}")
        merge_pod5s(pod5_dirs, merged_pod5)

        raw_bam = output_dir / f"{name}_raw.bam"
        print(f"Basecalling with {model} -> {raw_bam}")
        basecall(
            pod5_path=merged_pod5,
            output_bam=raw_bam,
            model=model,
            device=device,
            dorado_path=str(dorado_path) if dorado_path else None,
        )

        print("Building end_reason lookup from POD5...")
        lookup = build_pod5_lookup(merged_pod5)
        output_bam = output_dir / f"{name}_merged.bam"
        print(f"Tagging BAM -> {output_bam}")
        tagged = tag_bam(raw_bam, output_bam, lookup)
        print(f"Tagged {tagged}/{len(lookup)} reads with end_reason")

        raw_bam.unlink()

    print("\nDone!")


def cmd_subsample(
    experiment_path: Path, output_dir: Path, n_reads: int,
    seed: int, device: str, dorado_path: Path | None,
) -> None:
    """Run the subsample pipeline."""
    from sma_merge.basecall import basecall
    from sma_merge.subsample import subsample_pod5_multi
    from sma_merge.tag import build_pod5_lookup, tag_bam

    print(f"Scanning {experiment_path}...")
    runs = discover_runs(experiment_path)
    if not runs:
        print("No runs found.")
        sys.exit(1)

    groups = validate_runs(runs)
    print(format_validation(groups))

    output_dir.mkdir(parents=True, exist_ok=True)

    for g in groups:
        if not g.is_consistent:
            print(f"\nWARNING: Inconsistent runs in flowcell {g.flow_cell_id}, skipping.")
            continue

        name = g.flow_cell_id
        model = g.basecall_model
        pod5_dirs = [r.pod5_dir for r in g.runs]

        sub_pod5 = output_dir / f"{name}_sub{n_reads}.pod5"
        print(f"\nSubsampling {n_reads} reads from {len(pod5_dirs)} source dir(s) -> {sub_pod5}")
        selected = subsample_pod5_multi(pod5_dirs, sub_pod5, n_reads=n_reads, seed=seed)
        print(f"Selected {len(selected)} reads")

        raw_bam = output_dir / f"{name}_sub{n_reads}_raw.bam"
        print(f"Basecalling with {model} -> {raw_bam}")
        basecall(
            pod5_path=sub_pod5,
            output_bam=raw_bam,
            model=model,
            device=device,
            dorado_path=str(dorado_path) if dorado_path else None,
        )

        print("Building end_reason lookup from POD5...")
        lookup = build_pod5_lookup(sub_pod5)
        output_bam = output_dir / f"{name}_sub{n_reads}.bam"
        print(f"Tagging BAM -> {output_bam}")
        tagged = tag_bam(raw_bam, output_bam, lookup)
        print(f"Tagged {tagged}/{len(selected)} reads with end_reason")

        raw_bam.unlink()

    print("\nDone!")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "discover":
        cmd_discover(args.path)
    elif args.command == "merge":
        cmd_merge(args.path, args.output, args.device, args.dorado)
    elif args.command == "subsample":
        cmd_subsample(
            args.path, args.output, args.num_reads,
            args.seed, args.device, args.dorado,
        )


if __name__ == "__main__":
    main()
