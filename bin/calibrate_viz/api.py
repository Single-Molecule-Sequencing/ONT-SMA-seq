"""REST API router for SMA-seq barcode calibration visualizer."""

from __future__ import annotations

import sqlite3
from pathlib import Path

from fastapi import APIRouter, Query
from fastapi.responses import HTMLResponse, JSONResponse

router = APIRouter(prefix="/api")

# ---------------------------------------------------------------------------
# Module-level state
# ---------------------------------------------------------------------------

_db_paths: list[Path] = []


def set_db_paths(paths: list[Path]) -> None:
    """Set the list of calibration database paths."""
    global _db_paths
    _db_paths = paths


def _get_db(index: int = 0) -> Path:
    """Return a DB path by index, raising on invalid index."""
    if not _db_paths:
        raise ValueError("No databases loaded. Call set_db_paths() first.")
    if index < 0 or index >= len(_db_paths):
        raise ValueError(f"DB index {index} out of range (0..{len(_db_paths) - 1}).")
    return _db_paths[index]


def _get_columns(db_path: Path) -> list[str]:
    """Return column names of the Reads table."""
    conn = sqlite3.connect(db_path)
    try:
        cursor = conn.execute("PRAGMA table_info(Reads)")
        return [row[1] for row in cursor.fetchall()]
    finally:
        conn.close()


def _get_read_count(db_path: Path) -> int:
    """Return total row count in the Reads table."""
    conn = sqlite3.connect(db_path)
    try:
        result = conn.execute("SELECT COUNT(*) FROM Reads").fetchone()
        return result[0] if result else 0
    finally:
        conn.close()


def _get_tables(db_path: Path) -> list[str]:
    """Return list of table names in the database."""
    conn = sqlite3.connect(db_path)
    try:
        rows = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
        ).fetchall()
        return [r[0] for r in rows]
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------


@router.get("/overview")
async def api_overview():
    """Summary stats per loaded DB, returned as HTML fragment."""
    if not _db_paths:
        return HTMLResponse("<p>No databases loaded.</p>")

    html = ""
    for i, db_path in enumerate(_db_paths):
        try:
            count = _get_read_count(db_path)
            columns = _get_columns(db_path)
            tables = _get_tables(db_path)
        except Exception as exc:
            html += f"""
            <article>
              <header><strong>{db_path.name}</strong> (DB {i})</header>
              <p class="badge-error">Error reading database: {exc}</p>
            </article>
            """
            continue

        col_list = ", ".join(f"<code>{c}</code>" for c in columns)
        table_list = ", ".join(f"<code>{t}</code>" for t in tables)
        html += f"""
        <article>
          <header><strong>{db_path.name}</strong> (DB {i})</header>
          <table>
            <tr><td><strong>Path</strong></td><td><code>{db_path}</code></td></tr>
            <tr><td><strong>Total reads</strong></td><td>{count:,}</td></tr>
            <tr><td><strong>Tables</strong></td><td>{table_list}</td></tr>
            <tr><td><strong>Reads columns</strong></td><td>{col_list}</td></tr>
          </table>
        </article>
        """

    return HTMLResponse(html)


@router.get("/distributions")
async def api_distributions(
    column: str = Query("readlen"),
    group_by: str = Query(None),
    db: int = Query(0),
):
    """Return grouped KDE data as JSON."""
    from calibrate_viz.distributions import compute_grouped_kde

    db_path = _get_db(db)
    result = compute_grouped_kde(db_path, column, group_by=group_by)
    return JSONResponse(result)


@router.get("/distributions/peaks")
async def api_distribution_peaks(
    column: str = Query("readlen"),
    db: int = Query(0),
):
    """Return peak positions for the overall distribution."""
    from calibrate_viz.distributions import compute_kde, find_kde_peaks, load_distribution_data

    db_path = _get_db(db)
    data = load_distribution_data(db_path, column)
    if hasattr(data, "__len__") and len(data) < 2:
        return JSONResponse({"peaks": []})
    x, y = compute_kde(data)
    peaks = find_kde_peaks(x, y)
    return JSONResponse({"peaks": peaks})


@router.get("/threshold-impact", response_class=HTMLResponse)
async def api_threshold_impact(
    column: str = Query("readlen"),
    threshold: float = Query(...),
    group_by: str = Query(None),
    db: int = Query(0),
):
    """Return HTML fragment with classification counts at threshold.

    Shows count of reads above/below the threshold, optionally per group.
    """
    db_path = _get_db(db)
    conn = sqlite3.connect(db_path)
    try:
        if group_by:
            query = f"""
                SELECT [{group_by}],
                       SUM(CASE WHEN [{column}] < ? THEN 1 ELSE 0 END) as below,
                       SUM(CASE WHEN [{column}] >= ? THEN 1 ELSE 0 END) as above,
                       COUNT(*) as total
                FROM Reads
                WHERE [{column}] IS NOT NULL
                GROUP BY [{group_by}]
                ORDER BY [{group_by}]
            """
            rows = conn.execute(query, (threshold, threshold)).fetchall()
        else:
            query = f"""
                SELECT 'all' as grp,
                       SUM(CASE WHEN [{column}] < ? THEN 1 ELSE 0 END) as below,
                       SUM(CASE WHEN [{column}] >= ? THEN 1 ELSE 0 END) as above,
                       COUNT(*) as total
                FROM Reads
                WHERE [{column}] IS NOT NULL
            """
            rows = conn.execute(query, (threshold, threshold)).fetchall()
    finally:
        conn.close()

    html = f"""
    <table>
      <caption>Threshold: {column} = {threshold:,.1f}</caption>
      <thead>
        <tr>
          <th>Group</th>
          <th>Below</th>
          <th>Above</th>
          <th>Total</th>
          <th>% Above</th>
        </tr>
      </thead>
      <tbody>
    """
    for row in rows:
        grp, below, above, total = row
        pct = (above / total * 100) if total > 0 else 0
        grp_label = str(grp) if grp is not None else "unknown"
        html += f"""
        <tr>
          <td>{grp_label}</td>
          <td>{below:,}</td>
          <td>{above:,}</td>
          <td>{total:,}</td>
          <td>{pct:.1f}%</td>
        </tr>
        """

    html += "</tbody></table>"
    return HTMLResponse(html)
