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


# ---------------------------------------------------------------------------
# Separation endpoints
# ---------------------------------------------------------------------------


@router.get("/separation/pairwise")
async def api_pairwise_distances(db: int = Query(0)):
    """Return pairwise barcode edit distance matrix.

    Discovers barcodes used in the experiment from the Reads table,
    then computes all-vs-all NW edit distances.
    """
    from barcodes import BARCODES
    from calibrate_viz.separation import compute_pairwise_distances

    db_path = _get_db(db)
    conn = sqlite3.connect(db_path)
    try:
        rows = conn.execute(
            "SELECT DISTINCT bc_start_id FROM Reads "
            "WHERE bc_start_id IS NOT NULL ORDER BY bc_start_id"
        ).fetchall()
    finally:
        conn.close()

    # Build a subset of barcodes actually used in this experiment
    used_ids = [r[0] for r in rows if r[0] in BARCODES]
    if not used_ids:
        return JSONResponse({"ids": [], "matrix": {}})

    subset = {bc_id: BARCODES[bc_id] for bc_id in used_ids}
    matrix = compute_pairwise_distances(subset)
    return JSONResponse({"ids": sorted(subset.keys()), "matrix": matrix})


@router.get("/separation/metrics")
async def api_separation_metrics(db: int = Query(0)):
    """Return per-barcode separation metrics.

    Loads empirical edit distances from the database and computes
    separation gap and estimated error rate per barcode.
    """
    from barcodes import BARCODES
    from calibrate_viz.separation import (
        compute_separation_metrics,
        load_barcode_separation_data,
    )

    db_path = _get_db(db)

    # Discover which barcodes are used
    conn = sqlite3.connect(db_path)
    try:
        rows = conn.execute(
            "SELECT DISTINCT bc_start_id FROM Reads "
            "WHERE bc_start_id IS NOT NULL ORDER BY bc_start_id"
        ).fetchall()
    finally:
        conn.close()

    used_ids = [r[0] for r in rows if r[0] in BARCODES]
    if not used_ids:
        return JSONResponse({})

    subset = {bc_id: BARCODES[bc_id] for bc_id in used_ids}
    read_eds = load_barcode_separation_data(db_path, subset)
    metrics = compute_separation_metrics(read_eds)
    return JSONResponse(metrics)


# ---------------------------------------------------------------------------
# Confidence / confusion endpoints
# ---------------------------------------------------------------------------


@router.get("/confidence")
async def api_confidence_data(db: int = Query(0)):
    """Return confidence data for scatter and KDE plots.

    Returns JSON with:
      - reads: list of {bc_start_conf, bc_end_conf, trunc_level} for scatter
      - start_kde: grouped KDE of bc_start_conf by trunc_level
      - end_kde: grouped KDE of bc_end_conf by trunc_level
    """
    from calibrate_viz.distributions import compute_grouped_kde

    db_path = _get_db(db)

    # Scatter data
    conn = sqlite3.connect(db_path)
    try:
        rows = conn.execute("""
            SELECT bc_start_conf, bc_end_conf, trunc_level
            FROM Reads
            WHERE bc_start_conf IS NOT NULL
        """).fetchall()
    finally:
        conn.close()

    reads = [
        {
            "bc_start_conf": r[0],
            "bc_end_conf": r[1],
            "trunc_level": r[2] or "unknown",
        }
        for r in rows
    ]

    # KDE for bc_start_conf grouped by trunc_level
    start_kde = compute_grouped_kde(db_path, "bc_start_conf", group_by="trunc_level")

    # KDE for bc_end_conf grouped by trunc_level
    end_kde = compute_grouped_kde(db_path, "bc_end_conf", group_by="trunc_level")

    return JSONResponse({
        "reads": reads,
        "start_kde": start_kde,
        "end_kde": end_kde,
    })


@router.get("/confusion-matrix")
async def api_confusion_matrix(
    ed_threshold: float = Query(0.1),
    db: int = Query(0),
):
    """Return confusion matrix data as HTML table."""
    from calibrate_viz.confusion import compute_confusion_matrix

    db_path = _get_db(db)
    matrix = compute_confusion_matrix(db_path, ed_threshold=ed_threshold)

    labels = matrix["labels"]
    counts = matrix["counts"]

    if not labels:
        return HTMLResponse("<p>No barcode assignments found.</p>")

    html = """
    <table>
      <caption>Confusion Matrix (assigned vs. true barcode)</caption>
      <thead>
        <tr>
          <th>Assigned \\ True</th>
    """
    for label in labels:
        html += f"<th>{label}</th>"
    html += "<th>Total</th></tr></thead><tbody>"

    for i, label in enumerate(labels):
        row_total = sum(counts[i])
        html += f"<tr><td><strong>{label}</strong></td>"
        for j in range(len(labels)):
            val = counts[i][j]
            cls = ' class="badge-ok"' if i == j and val > 0 else ""
            html += f"<td{cls}>{val:,}</td>"
        html += f"<td>{row_total:,}</td></tr>"

    html += "</tbody></table>"
    return HTMLResponse(html)


@router.get("/threshold-classify", response_class=HTMLResponse)
async def api_threshold_classify(
    start_barcode_min: float = Query(0.6),
    full_length_threshold: float = Query(0.75),
    db: int = Query(0),
):
    """Return classification counts at given thresholds as HTML."""
    from calibrate_viz.confusion import compute_threshold_impact

    db_path = _get_db(db)
    counts = compute_threshold_impact(
        db_path,
        start_barcode_min=start_barcode_min,
        full_length_threshold=full_length_threshold,
    )

    total = sum(counts.values())
    html = f"""
    <table>
      <caption>Classification at start_min={start_barcode_min:.2f},
               fl_thresh={full_length_threshold:.2f}</caption>
      <thead>
        <tr>
          <th>Classification</th>
          <th>Count</th>
          <th>% of Total</th>
        </tr>
      </thead>
      <tbody>
    """

    level_order = ["full_length", "bc1_target_bc2", "bc1_target", "adapter_only"]
    for level in level_order:
        count = counts.get(level, 0)
        pct = (count / total * 100) if total > 0 else 0
        html += f"""
        <tr>
          <td>{level}</td>
          <td>{count:,}</td>
          <td>{pct:.1f}%</td>
        </tr>
        """

    html += f"""
      </tbody>
      <tfoot>
        <tr>
          <td><strong>Total</strong></td>
          <td><strong>{total:,}</strong></td>
          <td><strong>100%</strong></td>
        </tr>
      </tfoot>
    </table>
    """
    return HTMLResponse(html)


@router.get("/affected-reads", response_class=HTMLResponse)
async def api_affected_reads(
    old_start_min: float = Query(0.6),
    old_fl_thresh: float = Query(0.75),
    new_start_min: float = Query(0.6),
    new_fl_thresh: float = Query(0.75),
    db: int = Query(0),
):
    """Return affected reads table as HTML."""
    from calibrate_viz.confusion import get_affected_reads

    db_path = _get_db(db)
    affected = get_affected_reads(
        db_path,
        old_start_min=old_start_min,
        old_fl_thresh=old_fl_thresh,
        new_start_min=new_start_min,
        new_fl_thresh=new_fl_thresh,
    )

    if not affected:
        return HTMLResponse("<p>No reads changed classification.</p>")

    html = f"""
    <table>
      <caption>{len(affected)} read(s) changed classification</caption>
      <thead>
        <tr>
          <th>Read ID</th>
          <th>Old Level</th>
          <th>New Level</th>
          <th>bc_start_conf</th>
          <th>bc_end_conf</th>
        </tr>
      </thead>
      <tbody>
    """
    for r in affected:
        end_conf = f"{r['bc_end_conf']:.3f}" if r["bc_end_conf"] is not None else "N/A"
        html += f"""
        <tr>
          <td><code>{r['read_id']}</code></td>
          <td>{r['old_level']}</td>
          <td>{r['new_level']}</td>
          <td>{r['bc_start_conf']:.3f}</td>
          <td>{end_conf}</td>
        </tr>
        """

    html += "</tbody></table>"
    return HTMLResponse(html)
