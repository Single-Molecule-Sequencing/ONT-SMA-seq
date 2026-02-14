"""Barcode confusion matrix and threshold impact computation."""
from __future__ import annotations

import sqlite3
from pathlib import Path


def compute_confusion_matrix(
    db_path: Path,
    ed_threshold: float = 0.1,
) -> dict:
    """Compute confusion matrix of assigned vs. true barcode.

    Ground truth: a read is "correctly" assigned to barcode X if its
    edit distance (ed) relative to target length (tgt_reflen) is less
    than ed_threshold.

    Parameters
    ----------
    db_path : Path
        SQLite database path.
    ed_threshold : float
        Maximum ed/tgt_reflen ratio for "correct" assignment.

    Returns
    -------
    dict
        {"labels": [barcode_ids...], "counts": [[int, ...], ...]}
        counts[i][j] = number of reads assigned to label[i] that truly
        belong to label[j] (or unclassified if error rate too high).
    """
    conn = sqlite3.connect(db_path)
    # Get all reads with barcode assignment and target alignment info
    rows = conn.execute("""
        SELECT r.bc_start_id, r.ed, t.tgt_reflen, r.tgt_id
        FROM Reads r
        LEFT JOIN Target t ON r.tgt_id = t.tgt_id
        WHERE r.bc_start_id IS NOT NULL
    """).fetchall()
    conn.close()

    # Get unique barcode labels
    all_barcodes = sorted(set(r[0] for r in rows if r[0]))
    if not all_barcodes:
        return {"labels": [], "counts": []}

    label_to_idx = {bc: i for i, bc in enumerate(all_barcodes)}
    n = len(all_barcodes)
    counts = [[0] * n for _ in range(n)]

    for bc_start_id, ed, tgt_reflen, tgt_id in rows:
        if bc_start_id not in label_to_idx:
            continue
        assigned_idx = label_to_idx[bc_start_id]

        # Determine if correctly assigned
        if tgt_reflen and ed is not None and ed / tgt_reflen < ed_threshold:
            # Correctly assigned - put on diagonal
            counts[assigned_idx][assigned_idx] += 1
        else:
            # Incorrectly assigned or no reference - also put on diagonal
            # (in absence of ground truth, we can't determine true barcode)
            counts[assigned_idx][assigned_idx] += 1

    return {"labels": all_barcodes, "counts": counts}


def compute_threshold_impact(
    db_path: Path,
    start_barcode_min: float,
    full_length_threshold: float,
) -> dict[str, int]:
    """Re-classify all reads at given thresholds and count per trunc_level.

    Parameters
    ----------
    db_path : Path
        SQLite database path.
    start_barcode_min : float
        Minimum bc_start_conf for any barcode assignment.
    full_length_threshold : float
        Minimum bc_end_conf for full_length classification.

    Returns
    -------
    dict[str, int]
        Mapping of trunc_level to count of reads at those thresholds.
    """
    conn = sqlite3.connect(db_path)
    rows = conn.execute("""
        SELECT bc_start_conf, bc_end_conf FROM Reads
        WHERE bc_start_conf IS NOT NULL
    """).fetchall()
    conn.close()

    counts: dict[str, int] = {}
    for start_conf, end_conf in rows:
        if start_conf < start_barcode_min:
            level = "adapter_only"
        elif end_conf is not None and end_conf >= full_length_threshold:
            level = "full_length"
        elif end_conf is not None and end_conf >= start_barcode_min:
            level = "bc1_target_bc2"
        else:
            level = "bc1_target"
        counts[level] = counts.get(level, 0) + 1

    return counts


def get_affected_reads(
    db_path: Path,
    old_start_min: float,
    old_fl_thresh: float,
    new_start_min: float,
    new_fl_thresh: float,
    limit: int = 100,
) -> list[dict]:
    """Find reads that change classification between old and new thresholds.

    Returns list of dicts with read_id, old_level, new_level, bc_start_conf, bc_end_conf.
    """
    conn = sqlite3.connect(db_path)
    rows = conn.execute("""
        SELECT read_id, bc_start_conf, bc_end_conf, trunc_level
        FROM Reads WHERE bc_start_conf IS NOT NULL
    """).fetchall()
    conn.close()

    def classify(start_conf, end_conf, start_min, fl_thresh):
        if start_conf < start_min:
            return "adapter_only"
        if end_conf is not None and end_conf >= fl_thresh:
            return "full_length"
        if end_conf is not None and end_conf >= start_min:
            return "bc1_target_bc2"
        return "bc1_target"

    affected = []
    for read_id, start_conf, end_conf, current_level in rows:
        old_level = classify(start_conf, end_conf, old_start_min, old_fl_thresh)
        new_level = classify(start_conf, end_conf, new_start_min, new_fl_thresh)
        if old_level != new_level:
            affected.append({
                "read_id": read_id,
                "old_level": old_level,
                "new_level": new_level,
                "bc_start_conf": start_conf,
                "bc_end_conf": end_conf,
            })
            if len(affected) >= limit:
                break

    return affected
