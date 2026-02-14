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


def subsample_pod5(
    pod5_dir: Path,
    output_pod5: Path,
    n_reads: int,
    seed: int = 42,
) -> list[str]:
    """Subsample N reads from POD5 files and write to a new POD5.

    Returns the list of selected read IDs.
    """
    all_ids = collect_read_ids(pod5_dir)
    rng = np.random.default_rng(seed)
    n = min(n_reads, len(all_ids))
    indices = rng.choice(len(all_ids), size=n, replace=False)
    selected = sorted([all_ids[i] for i in indices])

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".txt", delete=False
    ) as f:
        for rid in selected:
            f.write(f"{rid}\n")
        id_file = Path(f.name)

    try:
        pod5_files = sorted(pod5_dir.rglob("*.pod5"))
        output_pod5.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                "pod5", "filter",
                *[str(f) for f in pod5_files],
                "--output", str(output_pod5),
                "--ids", str(id_file),
            ],
            check=True,
        )
    finally:
        id_file.unlink(missing_ok=True)

    return selected
