"""Subsample reads from POD5 files."""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pod5 as p5


def collect_read_ids(pod5_dir: Path) -> list[str]:
    """Collect all read IDs from POD5 files under a directory."""
    read_ids: list[str] = []
    pod5_files = sorted(pod5_dir.rglob("*.pod5"))
    for pod5_file in pod5_files:
        with p5.Reader(pod5_file) as reader:
            for read in reader.reads():
                read_ids.append(str(read.read_id))
    return read_ids


def collect_read_ids_multi(pod5_dirs: list[Path]) -> list[str]:
    """Collect all read IDs from multiple POD5 directories."""
    read_ids: list[str] = []
    for d in pod5_dirs:
        read_ids.extend(collect_read_ids(d))
    return read_ids


def subsample_pod5(
    pod5_dir: Path,
    output_pod5: Path,
    n_reads: int,
    seed: int = 42,
) -> list[str]:
    """Subsample N reads from POD5 files in a single directory.

    Returns the list of selected read IDs.
    """
    return subsample_pod5_multi([pod5_dir], output_pod5, n_reads, seed)


def subsample_pod5_multi(
    pod5_dirs: list[Path],
    output_pod5: Path,
    n_reads: int,
    seed: int = 42,
) -> list[str]:
    """Subsample N reads from POD5 files across multiple directories.

    Collects read IDs from all directories, randomly selects N, then
    calls pod5 filter on the original source files (no intermediate merge).

    Returns the list of selected read IDs.
    """
    all_ids = collect_read_ids_multi(pod5_dirs)
    rng = np.random.default_rng(seed)
    n = min(n_reads, len(all_ids))
    indices = rng.choice(len(all_ids), size=n, replace=False)
    selected = sorted([all_ids[i] for i in indices])

    # Gather all source pod5 files
    all_pod5_files: list[Path] = []
    for d in pod5_dirs:
        all_pod5_files.extend(sorted(d.rglob("*.pod5")))

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".txt", delete=False
    ) as f:
        for rid in selected:
            f.write(f"{rid}\n")
        id_file = Path(f.name)

    try:
        output_pod5.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                "pod5", "filter",
                *[str(f) for f in all_pod5_files],
                "--output", str(output_pod5),
                "--ids", str(id_file),
            ],
            check=True,
        )
    finally:
        id_file.unlink(missing_ok=True)

    return selected
