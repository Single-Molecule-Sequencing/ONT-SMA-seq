"""Barcode separation analysis: pairwise edit distances and separation metrics."""
from __future__ import annotations

import edlib


def compute_pairwise_distances(barcodes: dict[str, str]) -> dict[str, dict[str, int]]:
    """Compute pairwise edit distance matrix for a set of barcodes.

    Uses edlib NW (global alignment) mode for comparing barcode sequences
    of the same length.

    Parameters
    ----------
    barcodes : dict[str, str]
        Mapping of barcode_id to sequence.

    Returns
    -------
    dict[str, dict[str, int]]
        Nested dict: matrix[id_a][id_b] = edit_distance.
        Symmetric with zero diagonal.
    """
    ids = sorted(barcodes.keys())
    matrix: dict[str, dict[str, int]] = {bc: {} for bc in ids}

    for i, id_a in enumerate(ids):
        matrix[id_a][id_a] = 0
        for j in range(i + 1, len(ids)):
            id_b = ids[j]
            result = edlib.align(
                barcodes[id_a], barcodes[id_b], mode="NW", task="distance"
            )
            dist = result["editDistance"]
            matrix[id_a][id_b] = dist
            matrix[id_b][id_a] = dist

    return matrix


def compute_separation_metrics(
    read_eds: dict[str, dict[str, list[float]]],
) -> dict[str, dict]:
    """Compute per-barcode separation metrics from empirical edit distances.

    Parameters
    ----------
    read_eds : dict
        Structure: {"bc_id": {"true_match": [ed1, ed2, ...], "next_best": [ed1, ...]}}
        true_match = edit distances of reads to their assigned barcode.
        next_best = edit distances of reads to the next-best-matching barcode.

    Returns
    -------
    dict[str, dict]
        Per-barcode metrics: mean_true_ed, mean_next_best_ed, separation_gap,
        estimated_error_rate (fraction of true_match > min(next_best)).
    """
    metrics: dict[str, dict] = {}
    for bc_id, eds in read_eds.items():
        true_eds = eds["true_match"]
        next_eds = eds["next_best"]
        if not true_eds or not next_eds:
            continue

        mean_true = sum(true_eds) / len(true_eds)
        mean_next = sum(next_eds) / len(next_eds)
        gap = mean_next - mean_true

        # Estimated error rate: what fraction of true_match eds exceed
        # the minimum next_best ed?
        min_next = min(next_eds)
        error_count = sum(1 for ed in true_eds if ed >= min_next)
        error_rate = error_count / len(true_eds)

        metrics[bc_id] = {
            "mean_true_ed": round(mean_true, 2),
            "mean_next_best_ed": round(mean_next, 2),
            "separation_gap": round(gap, 2),
            "estimated_error_rate": round(error_rate, 4),
            "n_reads": len(true_eds),
        }

    return metrics


def load_barcode_separation_data(
    db_path,
    barcodes: dict[str, str],
) -> dict[str, dict[str, list[float]]]:
    """Load empirical per-read edit distances from database.

    For each read assigned to barcode X, compute:
    - true_match: edit distance to barcode X
    - next_best: minimum edit distance to any other barcode in the set

    Parameters
    ----------
    db_path : Path
        SQLite database path.
    barcodes : dict[str, str]
        The barcode set used in this experiment.

    Returns
    -------
    dict[str, dict[str, list[float]]]
        Per-barcode read edit distance data.
    """
    import sqlite3

    conn = sqlite3.connect(db_path)
    rows = conn.execute("""
        SELECT bc_start_id, bc_start_ed FROM Reads
        WHERE bc_start_id IS NOT NULL AND bc_start_ed IS NOT NULL
    """).fetchall()
    conn.close()

    result: dict[str, dict[str, list[float]]] = {}
    for bc_id, ed in rows:
        if bc_id not in barcodes:
            continue
        if bc_id not in result:
            result[bc_id] = {"true_match": [], "next_best": []}
        result[bc_id]["true_match"].append(float(ed))

        # Compute next-best edit distance
        best_alt = float("inf")
        for alt_id, alt_seq in barcodes.items():
            if alt_id == bc_id:
                continue
            alt_result = edlib.align(
                barcodes[bc_id], alt_seq, mode="NW", task="distance"
            )
            best_alt = min(best_alt, alt_result["editDistance"])
        if best_alt < float("inf"):
            result[bc_id]["next_best"].append(float(best_alt))

    return result
