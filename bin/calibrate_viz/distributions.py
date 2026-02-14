"""KDE distribution computation with peak detection for barcode calibration."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks, savgol_filter
from scipy.stats import gaussian_kde


def compute_kde(
    data: np.ndarray,
    n_points: int = 512,
    bw_method: str | float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute a Gaussian KDE over the data range.

    Parameters
    ----------
    data : np.ndarray
        1-D array of numeric values.
    n_points : int
        Number of evaluation points.
    bw_method : str, float, or None
        Bandwidth method passed to scipy gaussian_kde.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (x, y) arrays suitable for plotting.
    """
    if len(data) < 2:
        return np.array([]), np.array([])
    kde = gaussian_kde(data, bw_method=bw_method)
    x_min, x_max = float(data.min()), float(data.max())
    margin = (x_max - x_min) * 0.1
    x = np.linspace(x_min - margin, x_max + margin, n_points)
    y = kde(x)
    return x, y


def find_kde_peaks(
    x: np.ndarray,
    y: np.ndarray,
    window_length: int = 31,
    polyorder: int = 3,
    prominence: float = 0.01,
) -> list[float]:
    """Find peaks in a KDE curve using Savitzky-Golay smoothing.

    Parameters
    ----------
    x : np.ndarray
        X-axis values from compute_kde.
    y : np.ndarray
        Y-axis (density) values from compute_kde.
    window_length : int
        Savitzky-Golay window length (must be odd).
    polyorder : int
        Polynomial order for smoothing.
    prominence : float
        Minimum peak prominence for scipy find_peaks.

    Returns
    -------
    list[float]
        X-values of detected peaks.
    """
    if len(y) < window_length:
        return []
    smoothed = savgol_filter(y, window_length, polyorder)
    peak_indices, _ = find_peaks(smoothed, prominence=prominence)
    return [float(x[i]) for i in peak_indices]


def load_distribution_data(
    db_path: Path,
    column: str,
    group_by: str | None = None,
    where: str | None = None,
) -> dict[str, np.ndarray] | np.ndarray:
    """Load a numeric column from the Reads table for distribution analysis.

    Parameters
    ----------
    db_path : Path
        Path to SQLite database.
    column : str
        Column name to load (e.g., 'readlen', 'signal_duration_s').
    group_by : str or None
        If provided, group results by this column and return a dict.
    where : str or None
        Additional SQL WHERE clause.

    Returns
    -------
    dict[str, np.ndarray] or np.ndarray
        If group_by is provided, dict mapping group values to arrays.
        Otherwise a single array.
    """
    conn = sqlite3.connect(db_path)
    try:
        if group_by:
            query = f"SELECT [{group_by}], [{column}] FROM Reads WHERE [{column}] IS NOT NULL"
            if where:
                query += f" AND ({where})"
            rows = conn.execute(query).fetchall()
            groups: dict[str, list[float]] = {}
            for grp, val in rows:
                key = str(grp) if grp is not None else "unknown"
                groups.setdefault(key, []).append(float(val))
            return {k: np.array(v) for k, v in groups.items()}
        else:
            query = f"SELECT [{column}] FROM Reads WHERE [{column}] IS NOT NULL"
            if where:
                query += f" AND ({where})"
            rows = conn.execute(query).fetchall()
            return np.array([float(r[0]) for r in rows])
    finally:
        conn.close()


def compute_grouped_kde(
    db_path: Path,
    column: str,
    group_by: str | None = None,
    n_points: int = 512,
) -> dict:
    """Load data and compute KDE for each group.

    Returns
    -------
    dict
        Structure: {"groups": {"label": {"x": [...], "y": [...]}, ...},
                     "peaks": {"label": [peak_x, ...], ...}}
    """
    if group_by:
        grouped = load_distribution_data(db_path, column, group_by=group_by)
        result: dict = {"groups": {}, "peaks": {}}
        for label, data in grouped.items():
            if len(data) >= 2:
                x, y = compute_kde(data, n_points=n_points)
                result["groups"][label] = {
                    "x": x.tolist(),
                    "y": y.tolist(),
                    "count": len(data),
                }
                result["peaks"][label] = find_kde_peaks(x, y)
        return result
    else:
        data = load_distribution_data(db_path, column)
        if len(data) < 2:
            return {"groups": {}, "peaks": {}}
        x, y = compute_kde(data, n_points=n_points)
        return {
            "groups": {"all": {"x": x.tolist(), "y": y.tolist(), "count": len(data)}},
            "peaks": {"all": find_kde_peaks(x, y)},
        }
