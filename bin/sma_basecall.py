#!/usr/bin/env python3
"""sma_basecall.py — Basecall from pod5 or merge existing BAMs."""

import json as _json
import random
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


DORADO_BIN = Path.home() / "dorado" / "bin" / "dorado"


# ---------------------------------------------------------------------------
# BAM collection
# ---------------------------------------------------------------------------


def collect_bam_files(bam_dirs: list[str]) -> list[Path]:
    """Collect all .bam files from directories, including subdirectories."""
    bams = []
    for d in bam_dirs:
        dp = Path(d)
        if not dp.is_dir():
            continue
        for f in sorted(dp.rglob("*.bam")):
            if f.is_file() and f.stat().st_size > 0:
                bams.append(f)
    return bams


def merge_bams(bam_files: list[Path], output: Path) -> None:
    """Merge multiple BAM files into one sorted, indexed BAM."""
    if not bam_files:
        raise ValueError("No BAM files to merge")

    output.parent.mkdir(parents=True, exist_ok=True)

    if len(bam_files) == 1:
        subprocess.run(
            ["samtools", "sort", "-o", str(output), str(bam_files[0])],
            check=True,
        )
    else:
        merged_tmp = output.with_suffix(".unsorted.bam")
        subprocess.run(
            ["samtools", "merge", "-f", str(merged_tmp)]
            + [str(b) for b in bam_files],
            check=True,
        )
        subprocess.run(
            ["samtools", "sort", "-o", str(output), str(merged_tmp)],
            check=True,
        )
        merged_tmp.unlink()

    # Index
    subprocess.run(["samtools", "index", str(output)], check=True)


# ---------------------------------------------------------------------------
# Dorado basecalling
# ---------------------------------------------------------------------------


def build_dorado_command(
    model: str,
    pod5_dir: str,
    dorado_bin: str = str(DORADO_BIN),
    read_ids_file: str | None = None,
) -> list[str]:
    """Build Dorado basecaller command line."""
    cmd = [
        dorado_bin, "basecaller", model, pod5_dir,
        "--no-trim", "--emit-moves",
    ]
    if read_ids_file:
        cmd.extend(["--read-ids", read_ids_file])
    return cmd


def run_dorado(
    model: str,
    pod5_dir: str,
    output_bam: Path,
    dorado_bin: str = str(DORADO_BIN),
    read_ids_file: str | None = None,
) -> None:
    """Run Dorado basecaller and write output BAM."""
    cmd = build_dorado_command(model, pod5_dir, dorado_bin, read_ids_file)
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    with open(output_bam, "wb") as out_fh:
        subprocess.run(cmd, stdout=out_fh, check=True)


# ---------------------------------------------------------------------------
# Pod5 subsampling
# ---------------------------------------------------------------------------


def subsample_read_ids(all_ids: list[str], n: int, seed: int = 42) -> list[str]:
    """Randomly subsample n read IDs from a list."""
    if n >= len(all_ids):
        return list(all_ids)
    rng = random.Random(seed)
    return rng.sample(all_ids, n)


def get_pod5_read_ids(pod5_dir: str) -> list[str]:
    """Get all read IDs from pod5 files in a directory."""
    import pod5

    ids = []
    pod5_path = Path(pod5_dir)
    for pod5_file in sorted(pod5_path.glob("*.pod5")):
        with pod5.Reader(pod5_file) as reader:
            for read in reader.reads():
                ids.append(str(read.read_id))
    return ids


# ---------------------------------------------------------------------------
# Sequencing summary merge
# ---------------------------------------------------------------------------


def merge_sequencing_summaries(summary_paths: list[str], output: str) -> int:
    """Merge multiple sequencing_summary.txt files into one TSV.

    Returns the number of data lines written.
    """
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    header_written = False
    count = 0
    with open(out_path, "w") as out_fh:
        for sp in summary_paths:
            with open(sp) as in_fh:
                for i, line in enumerate(in_fh):
                    if i == 0:
                        if not header_written:
                            out_fh.write(line)
                            header_written = True
                        continue
                    out_fh.write(line)
                    count += 1
    return count


def parse_sequencing_summary(summary_path: str) -> dict[str, dict]:
    """Parse a sequencing summary TSV into a dict keyed by read_id.

    Returns {read_id: {end_reason, duration, mean_qscore_template, ...}}.
    """
    result: dict[str, dict] = {}
    with open(summary_path) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) != len(header):
                continue
            row = dict(zip(header, fields))
            read_id = row.get("read_id", "")
            if not read_id:
                continue
            entry: dict = {}
            if "end_reason" in row:
                entry["end_reason"] = row["end_reason"]
            if "duration" in row:
                try:
                    entry["duration"] = float(row["duration"])
                except ValueError:
                    pass
            if "mean_qscore_template" in row:
                try:
                    entry["mean_qscore"] = float(row["mean_qscore_template"])
                except ValueError:
                    pass
            result[read_id] = entry
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Basecall from pod5 or merge existing BAMs"
    )
    parser.add_argument("manifest", type=Path,
                        help="Path to manifest JSON from sma_scan.py")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--from-pod5", action="store_true",
                      help="Re-basecall from pod5 files using Dorado")
    mode.add_argument("--from-bam", action="store_true",
                      help="Merge existing BAM files")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Output directory")
    parser.add_argument("--model", default="sup",
                        help="Dorado model [%(default)s]")
    parser.add_argument("--subsample", type=int, default=None,
                        help="Subsample N reads from pod5 (--from-pod5 only)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for subsampling [%(default)s]")
    parser.add_argument("--dorado-bin", default=str(DORADO_BIN),
                        help="Path to dorado binary")
    args = parser.parse_args()

    # Load manifest
    manifest = _json.loads(args.manifest.read_text())
    exp_id = manifest.get("experiment_id", "experiment")
    args.output.mkdir(parents=True, exist_ok=True)

    if args.from_bam:
        _run_from_bam(manifest, exp_id, args.output)
    elif args.from_pod5:
        _run_from_pod5(manifest, exp_id, args.output, args.model,
                       args.subsample, args.seed, args.dorado_bin)


def _run_from_bam(manifest: dict, exp_id: str, output_dir: Path) -> None:
    """Merge existing BAM files from all runs."""
    bam_dirs = []
    for fc in manifest.get("flow_cells", []):
        for run in fc.get("runs", []):
            run_dir = run.get("run_dir", "")
            bam_dir = run.get("bam_dir")
            if run_dir and bam_dir:
                bam_dirs.append(str(Path(run_dir) / bam_dir))

    print(f"[sma_basecall] Collecting BAMs from {len(bam_dirs)} directories")
    bam_files = collect_bam_files(bam_dirs)
    print(f"[sma_basecall] Found {len(bam_files)} BAM files")

    if not bam_files:
        sys.exit("[sma_basecall] Error: no BAM files found")

    merged_bam = output_dir / f"{exp_id}_merged.bam"
    print(f"[sma_basecall] Merging to {merged_bam}")
    merge_bams(bam_files, merged_bam)

    # Count reads in output
    import pysam
    with pysam.AlignmentFile(str(merged_bam), check_sq=False) as f:
        read_count = sum(1 for _ in f)
    print(f"[sma_basecall] Merged BAM: {read_count} reads")

    # Merge sequencing summaries
    seq_paths = manifest.get("action", {}).get("sequencing_summary_paths", [])
    exp_dir = manifest.get("experiment_dir", "")
    full_seq_paths = []
    for sp in seq_paths:
        full = Path(exp_dir) / sp
        if full.exists():
            full_seq_paths.append(str(full))

    summary_out = output_dir / f"{exp_id}_summary.tsv"
    if full_seq_paths:
        n = merge_sequencing_summaries(full_seq_paths, str(summary_out))
        print(f"[sma_basecall] Merged {n} summary records to {summary_out}")

    # Write provenance
    provenance = {
        "mode": "from_bam",
        "bam_files": len(bam_files),
        "total_reads": read_count,
        "output_bam": str(merged_bam),
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    prov_path = output_dir / f"{exp_id}_basecall.json"
    prov_path.write_text(_json.dumps(provenance, indent=2))
    print(f"[sma_basecall] Provenance written to {prov_path}")


def _run_from_pod5(manifest: dict, exp_id: str, output_dir: Path,
                   model: str, subsample_n: int | None, seed: int,
                   dorado_bin: str) -> None:
    """Re-basecall from pod5 files using Dorado."""
    pod5_sources = manifest.get("action", {}).get("pod5_sources", [])
    if not pod5_sources:
        sys.exit("[sma_basecall] Error: no pod5 sources in manifest")

    suffix = f"_subsample_{subsample_n}" if subsample_n else "_merged"
    output_bam = output_dir / f"{exp_id}{suffix}.bam"

    read_ids_file = None
    if subsample_n:
        print(f"[sma_basecall] Subsampling {subsample_n} reads from pod5")
        all_ids = []
        for pod5_dir in pod5_sources:
            all_ids.extend(get_pod5_read_ids(pod5_dir))
        print(f"[sma_basecall] Total reads in pod5: {len(all_ids)}")
        sampled = subsample_read_ids(all_ids, subsample_n, seed)
        ids_path = output_dir / f"{exp_id}_sampled_ids.txt"
        ids_path.write_text("\n".join(sampled) + "\n")
        read_ids_file = str(ids_path)
        print(f"[sma_basecall] Sampled {len(sampled)} read IDs")

    # Basecall each pod5 source
    temp_bams: list[Path] = []
    for i, pod5_dir in enumerate(pod5_sources):
        print(f"[sma_basecall] Basecalling {pod5_dir}")
        if len(pod5_sources) > 1:
            temp_bam = output_dir / f"_temp_bc_{i}.bam"
        else:
            temp_bam = output_bam.with_suffix(".unsorted.bam")
        run_dorado(model, pod5_dir, temp_bam, dorado_bin, read_ids_file)
        temp_bams.append(temp_bam)

    # Merge if multiple sources
    if len(temp_bams) > 1:
        merge_bams(temp_bams, output_bam)
        for t in temp_bams:
            t.unlink(missing_ok=True)
    else:
        # Sort and index single output
        subprocess.run(
            ["samtools", "sort", "-o", str(output_bam), str(temp_bams[0])],
            check=True,
        )
        temp_bams[0].unlink(missing_ok=True)
        subprocess.run(["samtools", "index", str(output_bam)], check=True)

    # Count reads
    import pysam
    with pysam.AlignmentFile(str(output_bam)) as f:
        read_count = sum(1 for _ in f)
    print(f"[sma_basecall] Output BAM: {read_count} reads")

    # Merge sequencing summaries
    seq_paths = manifest.get("action", {}).get("sequencing_summary_paths", [])
    exp_dir = manifest.get("experiment_dir", "")
    full_seq_paths = []
    for sp in seq_paths:
        full = Path(exp_dir) / sp
        if full.exists():
            full_seq_paths.append(str(full))

    summary_out = output_dir / f"{exp_id}_summary.tsv"
    if full_seq_paths:
        n = merge_sequencing_summaries(full_seq_paths, str(summary_out))
        print(f"[sma_basecall] Merged {n} summary records")

    # Write provenance
    provenance = {
        "mode": "from_pod5",
        "model": model,
        "subsample": subsample_n,
        "seed": seed,
        "pod5_sources": pod5_sources,
        "total_reads": read_count,
        "output_bam": str(output_bam),
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    prov_path = output_dir / f"{exp_id}_basecall.json"
    prov_path.write_text(_json.dumps(provenance, indent=2))
    print(f"[sma_basecall] Provenance written to {prov_path}")


if __name__ == "__main__":
    main()
