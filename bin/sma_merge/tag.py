"""Extract end_reason from POD5 and tag BAM reads."""
from __future__ import annotations

from pathlib import Path

import pod5 as p5
import pysam


def build_pod5_lookup(pod5_path: Path) -> dict[str, tuple[str, int]]:
    """Build a {read_id: (end_reason, num_samples)} lookup from POD5 file(s).

    Parameters
    ----------
    pod5_path : Path
        A single POD5 file or a directory containing POD5 files.
    """
    lookup: dict[str, tuple[str, int]] = {}
    if pod5_path.is_dir():
        pod5_files = sorted(pod5_path.rglob("*.pod5"))
    else:
        pod5_files = [pod5_path]

    for f in pod5_files:
        with p5.Reader(f) as reader:
            for read in reader.reads():
                er = read.end_reason.reason.name.lower()
                lookup[str(read.read_id)] = (er, read.num_samples)

    return lookup


def tag_bam(
    input_bam: Path,
    output_bam: Path,
    pod5_lookup: dict[str, tuple[str, int]],
) -> int:
    """Add end_reason (er) and signal_length (sl) tags to BAM reads.

    Returns the number of reads successfully tagged.
    """
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    tagged = 0
    with pysam.AlignmentFile(str(input_bam), check_sq=False) as in_af:
        with pysam.AlignmentFile(str(output_bam), "wb", header=in_af.header) as out_af:
            for read in in_af:
                rid = read.query_name
                if rid in pod5_lookup:
                    er, sl = pod5_lookup[rid]
                    read.set_tag("er", er, "Z")
                    read.set_tag("sl", sl, "i")
                    tagged += 1
                out_af.write(read)
    return tagged
