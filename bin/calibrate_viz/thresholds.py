"""Threshold optimization using ROC-style analysis."""
from __future__ import annotations

from pathlib import Path

import numpy as np


def compute_roc(
    true_confs: np.ndarray,
    false_confs: np.ndarray,
    n_points: int = 100,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute ROC curve for a confidence threshold.

    Parameters
    ----------
    true_confs : np.ndarray
        Confidence scores of correctly classified reads (true positives when above threshold).
    false_confs : np.ndarray
        Confidence scores of incorrectly classified reads (false positives when above threshold).
    n_points : int
        Number of threshold points to evaluate.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        (fpr, tpr, thresholds) — sorted by ascending FPR.
    """
    all_scores = np.concatenate([true_confs, false_confs])
    thresholds = np.linspace(float(all_scores.min()), float(all_scores.max()), n_points)

    tpr_list = []
    fpr_list = []
    for t in thresholds:
        tp = np.sum(true_confs >= t)
        fn = np.sum(true_confs < t)
        fp = np.sum(false_confs >= t)
        tn = np.sum(false_confs < t)

        tpr_val = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        fpr_val = fp / (fp + tn) if (fp + tn) > 0 else 0.0
        tpr_list.append(tpr_val)
        fpr_list.append(fpr_val)

    fpr = np.array(fpr_list)
    tpr = np.array(tpr_list)

    # Sort by FPR ascending
    order = np.argsort(fpr)
    return fpr[order], tpr[order], thresholds[order]


def recommend_threshold(
    true_confs: np.ndarray,
    false_confs: np.ndarray,
    target_fpr: float = 0.05,
) -> dict:
    """Find optimal threshold for a target false positive rate.

    Returns dict with: threshold, sensitivity, specificity, f1, fpr, reads_affected.
    """
    fpr, tpr, thresholds = compute_roc(true_confs, false_confs)

    # Find threshold where FPR is closest to target without exceeding it
    valid = fpr <= target_fpr
    if not np.any(valid):
        # If we can't achieve the target FPR, use the lowest FPR available
        idx = np.argmin(fpr)
    else:
        # Among valid thresholds, pick the one with highest TPR
        valid_indices = np.where(valid)[0]
        idx = valid_indices[np.argmax(tpr[valid_indices])]

    t = float(thresholds[idx])
    tp = int(np.sum(true_confs >= t))
    fn = int(np.sum(true_confs < t))
    fp = int(np.sum(false_confs >= t))
    tn = int(np.sum(false_confs < t))

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * precision * sensitivity / (precision + sensitivity) if (precision + sensitivity) > 0 else 0.0

    return {
        "threshold": round(t, 4),
        "sensitivity": round(sensitivity, 4),
        "specificity": round(specificity, 4),
        "f1": round(f1, 4),
        "fpr": round(float(fpr[idx]), 4),
        "true_positives": tp,
        "false_positives": fp,
        "reads_affected": fn + fp,
    }


def export_thresholds_to_toml(
    output_path: Path,
    start_barcode_min: float,
    full_length_threshold: float,
    flank_max_error_rate: float,
) -> None:
    """Write selected thresholds to a construct TOML confidence section."""
    # Try tomli_w first, fall back to manual
    try:
        import tomli_w
        data = {
            "sma": {
                "confidence": {
                    "start_barcode_min": start_barcode_min,
                    "full_length_threshold": full_length_threshold,
                    "flank_max_error_rate": flank_max_error_rate,
                }
            }
        }
        with open(output_path, "wb") as f:
            tomli_w.dump(data, f)
    except ImportError:
        with open(output_path, "w") as f:
            f.write("[sma.confidence]\n")
            f.write(f"start_barcode_min = {start_barcode_min}\n")
            f.write(f"full_length_threshold = {full_length_threshold}\n")
            f.write(f"flank_max_error_rate = {flank_max_error_rate}\n")
