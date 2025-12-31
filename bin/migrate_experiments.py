#!/usr/bin/env python3
"""Migrate experiments from existing nanopore_experiments.db to unified database."""

import argparse
import sqlite3
import sys
from pathlib import Path

# Add lib to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from lib.db_schema import create_central_db
from lib.db_lookups import populate_lookups
from lib.db_ops import CentralDB


DEFAULT_SOURCE = "/nfs/turbo/umms-atheylab/nanopore_experiments.db"
DEFAULT_TARGET = "nanopore_unified.db"


def extract_exp_id(unique_id: str, experiment_path: str) -> str:
    """Generate exp_id from unique_id or experiment_path.

    Args:
        unique_id: UUID from source database
        experiment_path: Full path to experiment directory

    Returns:
        Short experiment identifier derived from unique_id or path
    """
    # Prefer unique_id (UUID) - use first 8 chars
    if unique_id:
        return unique_id[:8]

    # Fall back to extracting from experiment_path
    # Path format: /data1/project/sample/YYYYMMDD_HHMM_Instrument_FlowCell_RunID
    if experiment_path:
        path_parts = Path(experiment_path).name.split("_")
        if len(path_parts) >= 5:
            # Use RunID (last part, typically UUID prefix)
            return path_parts[-1][:8]

    # Last resort: hash the path
    return f"exp_{hash(experiment_path) % 100000:05d}"


def extract_pod5_dir(experiment_path: str) -> str:
    """Extract POD5 directory path from experiment path.

    Args:
        experiment_path: Full path to experiment directory

    Returns:
        Path to pod5_pass directory (standard ONT output location)
    """
    if not experiment_path:
        return None

    # Standard ONT directory structure: experiment_path/pod5_pass/
    pod5_dir = Path(experiment_path) / "pod5_pass"
    return str(pod5_dir)


def read_source_experiments(source_path: str, limit: int = None) -> list:
    """Read experiments from source database.

    Args:
        source_path: Path to source nanopore_experiments.db
        limit: Optional limit on number of rows to read

    Returns:
        List of experiment dictionaries
    """
    conn = sqlite3.connect(source_path)
    conn.row_factory = sqlite3.Row

    sql = """
        SELECT
            id,
            experiment_path,
            data_root,
            unique_id,
            instrument,
            flow_cell_id,
            sample_id,
            protocol_group_id,
            protocol,
            protocol_run_id,
            acquisition_run_id,
            started,
            acquisition_stopped,
            processing_stopped,
            basecalling_enabled,
            pod5_files,
            fastq_files,
            bam_files,
            created_at
        FROM experiments
        ORDER BY id
    """

    if limit:
        sql += f" LIMIT {limit}"

    cursor = conn.execute(sql)
    rows = cursor.fetchall()
    conn.close()

    return [dict(row) for row in rows]


def migrate_experiment(db: CentralDB, exp: dict, dry_run: bool = False) -> dict:
    """Migrate a single experiment to the unified database.

    Args:
        db: CentralDB instance for target database
        exp: Source experiment dictionary
        dry_run: If True, only return what would be inserted

    Returns:
        Dictionary of mapped fields for the experiment
    """
    # Generate exp_id from source data
    exp_id = extract_exp_id(exp.get("unique_id"), exp.get("experiment_path"))

    # Map source fields to target schema
    mapped = {
        "exp_id": exp_id,
        "experiment_path": exp.get("experiment_path"),
        "instrument": exp.get("instrument"),
        "flow_cell_id": exp.get("flow_cell_id"),
        "sample_id": exp.get("sample_id"),
        "protocol": exp.get("protocol"),
        "started": exp.get("started"),
        "pod5_count": exp.get("pod5_files"),
        "pod5_dir": extract_pod5_dir(exp.get("experiment_path")),
    }

    # Remove None values for cleaner output
    mapped_clean = {k: v for k, v in mapped.items() if v is not None}

    if not dry_run:
        # Insert into target database
        kwargs = {k: v for k, v in mapped_clean.items() if k != "exp_id"}
        db.insert_experiment(exp_id, **kwargs)

    return mapped_clean


def main():
    parser = argparse.ArgumentParser(
        description="Migrate experiments from existing nanopore_experiments.db to unified database"
    )
    parser.add_argument(
        "-s", "--source",
        default=DEFAULT_SOURCE,
        help=f"Source database path [%(default)s]"
    )
    parser.add_argument(
        "-t", "--target",
        default=DEFAULT_TARGET,
        help=f"Target database path [%(default)s]"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be migrated without inserting"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit number of experiments to migrate (default: all)"
    )
    args = parser.parse_args()

    # Validate source database exists
    source_path = Path(args.source)
    if not source_path.exists():
        print(f"[migrate_experiments] Error: Source database not found: {source_path}", file=sys.stderr)
        sys.exit(1)

    target_path = Path(args.target)

    # Read source experiments
    print(f"[migrate_experiments] Reading from source: {source_path}")

    # For dry-run, limit to 5 by default
    limit = args.limit
    if args.dry_run and limit is None:
        limit = 5

    experiments = read_source_experiments(str(source_path), limit=limit)
    print(f"[migrate_experiments] Found {len(experiments)} experiments to migrate")

    if args.dry_run:
        print(f"\n[migrate_experiments] DRY RUN - showing first {len(experiments)} experiments:\n")
        print("-" * 80)

        for exp in experiments:
            # Show source data
            print(f"Source ID: {exp.get('id')}")
            print(f"  unique_id:       {exp.get('unique_id')}")
            print(f"  experiment_path: {exp.get('experiment_path')}")
            print(f"  instrument:      {exp.get('instrument')}")
            print(f"  flow_cell_id:    {exp.get('flow_cell_id')}")
            print(f"  sample_id:       {exp.get('sample_id')}")
            print(f"  pod5_files:      {exp.get('pod5_files')}")

            # Show mapped data
            mapped = migrate_experiment(None, exp, dry_run=True)
            print(f"\n  --> Would insert as:")
            print(f"      exp_id:    {mapped.get('exp_id')}")
            print(f"      pod5_dir:  {mapped.get('pod5_dir')}")
            print("-" * 80)

        print(f"\n[migrate_experiments] DRY RUN complete. No data inserted.")
        print(f"[migrate_experiments] Run without --dry-run to perform migration.")
        return

    # Initialize target database if needed
    if not target_path.exists():
        print(f"[migrate_experiments] Creating target database: {target_path}")
        conn = create_central_db(str(target_path))
        populate_lookups(conn)
        conn.close()

    # Perform migration
    print(f"[migrate_experiments] Migrating to target: {target_path}")

    success_count = 0
    skip_count = 0
    error_count = 0

    with CentralDB(str(target_path)) as db:
        for exp in experiments:
            exp_id = extract_exp_id(exp.get("unique_id"), exp.get("experiment_path"))

            # Check if already exists
            existing = db.get_experiment(exp_id)
            if existing:
                print(f"  [SKIP] {exp_id} already exists")
                skip_count += 1
                continue

            try:
                migrate_experiment(db, exp, dry_run=False)
                print(f"  [OK] {exp_id}")
                success_count += 1
            except Exception as e:
                print(f"  [ERROR] {exp_id}: {e}", file=sys.stderr)
                error_count += 1

    # Summary
    print(f"\n[migrate_experiments] Migration complete:")
    print(f"  - Migrated:  {success_count}")
    print(f"  - Skipped:   {skip_count}")
    print(f"  - Errors:    {error_count}")
    print(f"  - Target DB: {target_path}")


if __name__ == "__main__":
    main()
