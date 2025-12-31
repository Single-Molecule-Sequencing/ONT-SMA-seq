#!/usr/bin/env python3
"""Orchestrate the unified nanopore pipeline.

Watches for:
1. New experiments -> ready for basecalling
2. Completed basecall runs -> ready for SMA-seq analysis
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from lib import CentralDB


def show_status(db_path: str):
    """Show current pipeline status."""
    with CentralDB(db_path) as db:
        experiments = db.list_experiments()
        print(f"\n{'='*60}")
        print("PIPELINE STATUS")
        print(f"{'='*60}")
        print(f"\nExperiments: {len(experiments)}")

        # Count by status
        all_runs = db.list_basecall_runs()
        pending = len([r for r in all_runs if r['status'] == 'pending'])
        running = len([r for r in all_runs if r['status'] == 'running'])
        complete = len([r for r in all_runs if r['status'] == 'complete'])
        analyzed = len([r for r in all_runs if r['status'] == 'analyzed'])

        print(f"\nBasecall Runs:")
        print(f"  Pending:  {pending}")
        print(f"  Running:  {running}")
        print(f"  Complete: {complete}")
        print(f"  Analyzed: {analyzed}")

        # Show pending SMA analysis
        pending_sma = db.get_pending_sma_analysis()
        if pending_sma:
            print(f"\nReady for SMA Analysis ({len(pending_sma)}):")
            for run in pending_sma[:5]:
                print(f"  - {run['run_id']}: {run['exp_id']} ({run.get('model_tier', 'unknown')})")
            if len(pending_sma) > 5:
                print(f"  ... and {len(pending_sma) - 5} more")


def list_pending_sma(db_path: str):
    """List basecall runs ready for SMA-seq analysis."""
    with CentralDB(db_path) as db:
        pending = db.get_pending_sma_analysis()

        if not pending:
            print("No basecall runs pending SMA analysis.")
            return

        print(f"{'Run ID':<40} {'Exp ID':<20} {'Tier':<6} {'BAM Path'}")
        print("-" * 100)
        for run in pending:
            print(f"{run['run_id']:<40} {run['exp_id']:<20} {run.get('model_tier', 'N/A'):<6} {run.get('bam_path') or 'N/A'}")


def main():
    parser = argparse.ArgumentParser(
        description="Orchestrate the unified nanopore pipeline"
    )
    parser.add_argument(
        "-d", "--database",
        default="nanopore_unified.db",
        help="Path to central database [%(default)s]"
    )

    subparsers = parser.add_subparsers(dest='command', help='Commands')
    subparsers.add_parser('status', help='Show pipeline status')
    subparsers.add_parser('pending-sma', help='List runs ready for SMA analysis')

    args = parser.parse_args()

    if not Path(args.database).exists():
        print(f"[orchestrate] Error: Database not found: {args.database}", file=sys.stderr)
        sys.exit(1)

    if args.command == 'status':
        show_status(args.database)
    elif args.command == 'pending-sma':
        list_pending_sma(args.database)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
