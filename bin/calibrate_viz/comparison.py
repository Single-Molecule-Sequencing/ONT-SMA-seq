"""Cross-experiment comparison computation."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np

from calibrate_viz.distributions import compute_kde, find_kde_peaks, load_distribution_data


def compute_experiment_summary(db_path: Path) -> dict:
    """Compute summary metrics for a single experiment database.

    Returns
    -------
    dict
        Keys: total_reads, mean_readlen, mean_signal_duration, mean_qscore,
        trunc_proportions (dict), end_reason_proportions (dict),
        barcode_count, target_count.
    """
    conn = sqlite3.connect(db_path)

    total = conn.execute("SELECT COUNT(*) FROM Reads").fetchone()[0]

    # Mean values -- handle NULL/missing columns gracefully
    mean_readlen = 0.0
    try:
        row = conn.execute("SELECT AVG(readlen) FROM Reads WHERE readlen IS NOT NULL").fetchone()
        mean_readlen = float(row[0]) if row[0] else 0.0
    except Exception:
        pass

    mean_signal = 0.0
    try:
        row = conn.execute("SELECT AVG(signal_duration_s) FROM Reads WHERE signal_duration_s IS NOT NULL").fetchone()
        mean_signal = float(row[0]) if row[0] else 0.0
    except Exception:
        pass

    mean_qscore = 0.0
    try:
        row = conn.execute("SELECT AVG(mean_qscore) FROM Reads WHERE mean_qscore IS NOT NULL").fetchone()
        mean_qscore = float(row[0]) if row[0] else 0.0
    except Exception:
        pass

    # Truncation proportions
    trunc_proportions: dict[str, float] = {}
    try:
        rows = conn.execute(
            "SELECT trunc_level, COUNT(*) FROM Reads WHERE trunc_level IS NOT NULL GROUP BY trunc_level"
        ).fetchall()
        trunc_total = sum(r[1] for r in rows) or 1
        trunc_proportions = {r[0]: round(r[1] / trunc_total, 4) for r in rows}
    except Exception:
        pass

    # End reason proportions
    er_proportions: dict[str, float] = {}
    try:
        rows = conn.execute(
            "SELECT ER, COUNT(*) FROM Reads WHERE ER IS NOT NULL GROUP BY ER"
        ).fetchall()
        er_total = sum(r[1] for r in rows) or 1
        er_proportions = {r[0]: round(r[1] / er_total, 4) for r in rows}
    except Exception:
        pass

    # Unique counts
    barcode_count = 0
    target_count = 0
    try:
        barcode_count = conn.execute(
            "SELECT COUNT(DISTINCT bc_start_id) FROM Reads WHERE bc_start_id IS NOT NULL"
        ).fetchone()[0]
        target_count = conn.execute(
            "SELECT COUNT(DISTINCT tgt_id) FROM Reads WHERE tgt_id IS NOT NULL"
        ).fetchone()[0]
    except Exception:
        pass

    conn.close()

    return {
        "total_reads": total,
        "mean_readlen": round(mean_readlen, 1),
        "mean_signal_duration": round(mean_signal, 3),
        "mean_qscore": round(mean_qscore, 2),
        "trunc_proportions": trunc_proportions,
        "end_reason_proportions": er_proportions,
        "barcode_count": barcode_count,
        "target_count": target_count,
    }


def compute_overlay_distributions(
    db_paths: list[Path],
    column: str = "readlen",
    n_points: int = 512,
) -> list[dict]:
    """Compute KDE distributions for the same column across multiple databases.

    Returns one KDE per database for overlay plotting on shared axes.

    Returns
    -------
    list[dict]
        Each dict has: db_name, x (list), y (list), count, peaks.
    """
    results = []
    for db_path in db_paths:
        data = load_distribution_data(db_path, column)
        if isinstance(data, dict):
            # Shouldn't happen for ungrouped, but handle
            continue
        if len(data) < 2:
            results.append({
                "db_name": db_path.stem,
                "x": [],
                "y": [],
                "count": len(data),
                "peaks": [],
            })
            continue
        x, y = compute_kde(data, n_points=n_points)
        peaks = find_kde_peaks(x, y)
        results.append({
            "db_name": db_path.stem,
            "x": x.tolist(),
            "y": y.tolist(),
            "count": len(data),
            "peaks": peaks,
        })
    return results
