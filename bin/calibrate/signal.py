"""Load per-read metadata from MinKNOW sequencing summary files."""
from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_sequencing_summary(path: Path) -> dict[str, dict]:
    """Load sequencing summary TSV into a dict keyed by read_id.

    Reads the MinKNOW sequencing_summary_*.txt file and extracts
    duration, end_reason, mean_qscore, and sequence_length per read.

    Parameters
    ----------
    path : Path
        Path to a sequencing_summary_*.txt file.

    Returns
    -------
    dict[str, dict]
        Mapping of read_id to dict with keys: duration, end_reason,
        mean_qscore, sequence_length.
    """
    usecols = ["read_id", "duration", "end_reason",
               "sequence_length_template", "mean_qscore_template"]
    df = pd.read_csv(path, sep="\t", usecols=usecols)
    result: dict[str, dict] = {}
    for _, row in df.iterrows():
        result[row["read_id"]] = {
            "duration": float(row["duration"]),
            "end_reason": str(row["end_reason"]),
            "mean_qscore": float(row["mean_qscore_template"]),
            "sequence_length": int(row["sequence_length_template"]),
        }
    return result
