"""POD5 merging and dorado basecalling wrappers."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def find_dorado() -> str:
    """Find the dorado binary, checking PATH then common locations."""
    dorado_in_path = shutil.which("dorado")
    if dorado_in_path:
        return dorado_in_path

    home = Path.home()
    candidates = [
        home / "dorado" / "bin" / "dorado",
        home / ".local" / "bin" / "dorado",
        Path("/usr/local/bin/dorado"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)

    return "dorado"


def merge_pod5s(pod5_dirs: list[Path], output_pod5: Path) -> None:
    """Merge all POD5 files from multiple directories into one file."""
    all_files: list[Path] = []
    for d in pod5_dirs:
        all_files.extend(sorted(d.rglob("*.pod5")))

    if not all_files:
        raise ValueError("No POD5 files found in provided directories")

    output_pod5.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "pod5", "merge",
            *[str(f) for f in all_files],
            "--output", str(output_pod5),
        ],
        check=True,
    )


def basecall(
    pod5_path: Path,
    output_bam: Path,
    model: str,
    device: str = "auto",
    emit_moves: bool = True,
    dorado_path: str | None = None,
) -> None:
    """Run dorado basecaller on a POD5 file.

    Produces an unaligned BAM (no trimming, no demultiplexing).
    Dorado outputs unaligned BAM to stdout by default.
    """
    dorado = dorado_path or find_dorado()
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        dorado, "basecaller",
        model,
        str(pod5_path),
        "--no-trim",
        "--device", device,
    ]
    if emit_moves:
        cmd.append("--emit-moves")

    with open(output_bam, "wb") as f:
        subprocess.run(cmd, stdout=f, check=True)
