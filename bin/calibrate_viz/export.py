"""Static HTML and SVG export for calibration reports."""
from __future__ import annotations

import html as _html
import json
import sqlite3
from pathlib import Path

import numpy as np

from calibrate_viz.distributions import compute_grouped_kde, load_distribution_data, compute_kde, find_kde_peaks
from calibrate_viz.comparison import compute_experiment_summary


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_PAGES = [
    ("index.html", "Overview"),
    ("distributions.html", "Distributions"),
    ("confidence.html", "Confidence"),
    ("separation.html", "Separation"),
    ("thresholds.html", "Thresholds"),
]

# Colour palette for KDE curves (10-colour categorical)
_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------


def _esc(text: str) -> str:
    """HTML-escape a string."""
    return _html.escape(str(text))


def _nav_links() -> str:
    """Return HTML navigation links for all pages."""
    links = [f'<a href="{fn}">{_esc(label)}</a>' for fn, label in _PAGES]
    return " | ".join(links)


def _base_html(
    title: str,
    body: str,
    data_json: str = "",
    include_d3: bool = False,
) -> str:
    """Wrap body content in a self-contained HTML document."""
    d3_tag = ""
    if include_d3:
        d3_tag = '<script src="https://d3js.org/d3.v7.min.js"></script>'

    data_block = ""
    if data_json:
        data_block = (
            '<script id="page-data" type="application/json">\n'
            f"{data_json}\n"
            "</script>"
        )

    nav = _nav_links()

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{_esc(title)}</title>
<style>
body {{
    font-family: system-ui, -apple-system, sans-serif;
    max-width: 1100px;
    margin: 0 auto;
    padding: 2rem;
    color: #333;
    line-height: 1.6;
}}
nav.export-nav {{
    background: #f5f5f5;
    padding: 0.75rem 1rem;
    border-radius: 0.5rem;
    margin-bottom: 2rem;
    font-size: 0.9rem;
}}
nav.export-nav a {{
    text-decoration: none;
    color: #1976d2;
    margin: 0 0.25rem;
}}
nav.export-nav a:hover {{
    text-decoration: underline;
}}
table {{
    width: 100%;
    border-collapse: collapse;
    margin: 1rem 0;
}}
th, td {{
    padding: 0.5rem 0.75rem;
    border: 1px solid #ddd;
    text-align: left;
}}
th {{
    background: #f5f5f5;
    font-weight: 600;
}}
tr:nth-child(even) {{
    background: #fafafa;
}}
code {{
    background: #f0f0f0;
    padding: 0.1rem 0.3rem;
    border-radius: 3px;
    font-size: 0.9em;
}}
.summary-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 1rem;
    margin: 1rem 0;
}}
.summary-card {{
    background: #f5f5f5;
    padding: 1rem;
    border-radius: 0.5rem;
    text-align: center;
}}
.summary-card .value {{
    font-size: 1.5rem;
    font-weight: 700;
    color: #1976d2;
}}
.summary-card .label {{
    font-size: 0.85rem;
    color: #666;
}}
h1 {{ border-bottom: 2px solid #1976d2; padding-bottom: 0.5rem; }}
h2 {{ border-bottom: 1px solid #ddd; padding-bottom: 0.3rem; margin-top: 2rem; }}

@media print {{
    nav.export-nav {{ display: none; }}
    body {{ padding: 0; }}
}}
</style>
{d3_tag}
</head>
<body>
<nav class="export-nav">{nav}</nav>
{body}
{data_block}
</body>
</html>"""


# ---------------------------------------------------------------------------
# SVG generation (pure Python, no matplotlib)
# ---------------------------------------------------------------------------


def _build_svg(
    kde_data: dict,
    x_label: str,
    width: int = 800,
    height: int = 400,
    margin_top: int = 20,
    margin_right: int = 20,
    margin_bottom: int = 50,
    margin_left: int = 60,
) -> str:
    """Build an SVG string from KDE data.

    Parameters
    ----------
    kde_data : dict
        Output of compute_grouped_kde with "groups" and "peaks" keys.
    x_label : str
        Label for the x-axis.
    width, height : int
        Total SVG dimensions.
    """
    groups = kde_data.get("groups", {})
    if not groups:
        return f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"><text x="50%" y="50%" text-anchor="middle">No data</text></svg>'

    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom

    # Compute global x/y ranges
    all_x: list[float] = []
    all_y: list[float] = []
    for g in groups.values():
        all_x.extend(g["x"])
        all_y.extend(g["y"])

    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = 0.0, max(all_y)

    x_range = x_max - x_min if x_max != x_min else 1.0
    y_range = y_max - y_min if y_max != y_min else 1.0

    def sx(v: float) -> float:
        return margin_left + (v - x_min) / x_range * plot_w

    def sy(v: float) -> float:
        return margin_top + plot_h - (v - y_min) / y_range * plot_h

    parts: list[str] = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
    )

    # Background
    parts.append(f'<rect width="{width}" height="{height}" fill="white"/>')

    # Axes
    # X-axis line
    parts.append(
        f'<line x1="{margin_left}" y1="{margin_top + plot_h}" '
        f'x2="{margin_left + plot_w}" y2="{margin_top + plot_h}" '
        f'stroke="#333" stroke-width="1"/>'
    )
    # Y-axis line
    parts.append(
        f'<line x1="{margin_left}" y1="{margin_top}" '
        f'x2="{margin_left}" y2="{margin_top + plot_h}" '
        f'stroke="#333" stroke-width="1"/>'
    )

    # X-axis ticks (5 ticks)
    for i in range(6):
        val = x_min + i * x_range / 5
        px = sx(val)
        py = margin_top + plot_h
        parts.append(
            f'<line x1="{px:.1f}" y1="{py}" x2="{px:.1f}" y2="{py + 5}" stroke="#333" stroke-width="1"/>'
        )
        label = f"{val:.1f}" if abs(val) < 1000 else f"{val:.0f}"
        parts.append(
            f'<text x="{px:.1f}" y="{py + 18}" text-anchor="middle" '
            f'font-size="11" fill="#333">{_esc(label)}</text>'
        )

    # Y-axis ticks (4 ticks)
    for i in range(5):
        val = y_min + i * y_range / 4
        px = margin_left
        py = sy(val)
        parts.append(
            f'<line x1="{px - 5}" y1="{py:.1f}" x2="{px}" y2="{py:.1f}" stroke="#333" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{px - 8}" y="{py + 4:.1f}" text-anchor="end" '
            f'font-size="11" fill="#333">{val:.4f}</text>'
        )

    # Axis labels
    parts.append(
        f'<text x="{margin_left + plot_w / 2}" y="{height - 5}" '
        f'text-anchor="middle" font-size="13" fill="#333">{_esc(x_label)}</text>'
    )
    parts.append(
        f'<text x="15" y="{margin_top + plot_h / 2}" '
        f'text-anchor="middle" font-size="13" fill="#333" '
        f'transform="rotate(-90, 15, {margin_top + plot_h / 2})">Density</text>'
    )

    # KDE curves
    legend_y = margin_top + 15
    for idx, (label, g) in enumerate(sorted(groups.items())):
        color = _COLORS[idx % len(_COLORS)]
        xs = g["x"]
        ys = g["y"]
        if not xs:
            continue

        points = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in zip(xs, ys))
        parts.append(
            f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2"/>'
        )

        # Legend entry
        lx = margin_left + plot_w - 150
        parts.append(
            f'<rect x="{lx}" y="{legend_y - 8}" width="12" height="12" '
            f'fill="{color}" rx="2"/>'
        )
        count = g.get("count", "")
        legend_text = f"{_esc(label)} (n={count})" if count else _esc(label)
        parts.append(
            f'<text x="{lx + 16}" y="{legend_y + 2}" font-size="11" fill="#333">'
            f'{legend_text}</text>'
        )
        legend_y += 18

    parts.append("</svg>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Page generators
# ---------------------------------------------------------------------------


def _write_index(
    path: Path,
    experiment_id: str,
    summary: dict,
    signal_kde: dict,
    readlen_kde: dict,
) -> None:
    """Generate index.html -- overview page."""
    trunc_rows = ""
    for level, prop in summary.get("trunc_proportions", {}).items():
        trunc_rows += f"<tr><td>{_esc(level)}</td><td>{prop:.2%}</td></tr>\n"

    er_rows = ""
    for reason, prop in summary.get("end_reason_proportions", {}).items():
        er_rows += f"<tr><td>{_esc(reason)}</td><td>{prop:.2%}</td></tr>\n"

    # Inline SVG for a quick visual
    signal_svg = _build_svg(signal_kde, "Signal Duration (s)", width=500, height=250)
    readlen_svg = _build_svg(readlen_kde, "Read Length (bp)", width=500, height=250)

    body = f"""
<h1>Calibration Report: {_esc(experiment_id)}</h1>

<h2>Summary</h2>
<div class="summary-grid">
  <div class="summary-card">
    <div class="value">{summary.get('total_reads', 0):,}</div>
    <div class="label">Total Reads</div>
  </div>
  <div class="summary-card">
    <div class="value">{summary.get('mean_readlen', 0):.0f}</div>
    <div class="label">Mean Read Length</div>
  </div>
  <div class="summary-card">
    <div class="value">{summary.get('mean_signal_duration', 0):.3f} s</div>
    <div class="label">Mean Signal Duration</div>
  </div>
  <div class="summary-card">
    <div class="value">{summary.get('mean_qscore', 0):.2f}</div>
    <div class="label">Mean Q-Score</div>
  </div>
  <div class="summary-card">
    <div class="value">{summary.get('barcode_count', 0)}</div>
    <div class="label">Unique Barcodes</div>
  </div>
  <div class="summary-card">
    <div class="value">{summary.get('target_count', 0)}</div>
    <div class="label">Unique Targets</div>
  </div>
</div>

<h2>Signal Duration Distribution</h2>
{signal_svg}

<h2>Read Length Distribution</h2>
{readlen_svg}

<h2>Truncation Proportions</h2>
<table>
<thead><tr><th>Level</th><th>Proportion</th></tr></thead>
<tbody>{trunc_rows if trunc_rows else '<tr><td colspan="2">No data</td></tr>'}</tbody>
</table>

<h2>End Reason Proportions</h2>
<table>
<thead><tr><th>End Reason</th><th>Proportion</th></tr></thead>
<tbody>{er_rows if er_rows else '<tr><td colspan="2">No data</td></tr>'}</tbody>
</table>
"""

    data = {
        "summary": summary,
        "signal_kde": _serialize_kde(signal_kde),
        "readlen_kde": _serialize_kde(readlen_kde),
    }
    html = _base_html(
        title=f"Calibration Report - {experiment_id}",
        body=body,
        data_json=json.dumps(data, indent=2),
    )
    path.write_text(html, encoding="utf-8")


def _write_distributions(
    path: Path,
    experiment_id: str,
    signal_kde: dict,
    readlen_kde: dict,
) -> None:
    """Generate distributions.html -- KDE distribution page."""
    signal_svg = _build_svg(signal_kde, "Signal Duration (s)")
    readlen_svg = _build_svg(readlen_kde, "Read Length (bp)")

    body = f"""
<h1>Distributions: {_esc(experiment_id)}</h1>

<h2>Signal Duration KDE</h2>
{signal_svg}

<h2>Read Length KDE</h2>
{readlen_svg}
"""
    data = {
        "signal_kde": _serialize_kde(signal_kde),
        "readlen_kde": _serialize_kde(readlen_kde),
    }
    html = _base_html(
        title=f"Distributions - {experiment_id}",
        body=body,
        data_json=json.dumps(data, indent=2),
        include_d3=True,
    )
    path.write_text(html, encoding="utf-8")


def _write_confidence(
    path: Path,
    experiment_id: str,
    db_path: Path,
) -> None:
    """Generate confidence.html -- barcode confidence analysis."""
    # Load confidence data
    conn = sqlite3.connect(db_path)
    rows = conn.execute("""
        SELECT bc_start_id, bc_start_conf, bc_end_conf
        FROM Reads
        WHERE bc_start_conf IS NOT NULL
    """).fetchall()
    conn.close()

    # Group confidence by barcode
    bc_confs: dict[str, list[float]] = {}
    end_confs: list[float] = []
    for bc_id, start_conf, end_conf in rows:
        bc_confs.setdefault(str(bc_id), []).append(float(start_conf))
        if end_conf is not None:
            end_confs.append(float(end_conf))

    # Build confidence table
    table_rows = ""
    for bc_id in sorted(bc_confs.keys()):
        confs = bc_confs[bc_id]
        mean_c = sum(confs) / len(confs) if confs else 0
        min_c = min(confs) if confs else 0
        max_c = max(confs) if confs else 0
        table_rows += (
            f"<tr><td>{_esc(bc_id)}</td>"
            f"<td>{len(confs)}</td>"
            f"<td>{mean_c:.4f}</td>"
            f"<td>{min_c:.4f}</td>"
            f"<td>{max_c:.4f}</td></tr>\n"
        )

    body = f"""
<h1>Confidence Analysis: {_esc(experiment_id)}</h1>

<h2>Start Barcode Confidence by Barcode</h2>
<table>
<thead><tr><th>Barcode</th><th>Count</th><th>Mean Conf</th><th>Min</th><th>Max</th></tr></thead>
<tbody>{table_rows if table_rows else '<tr><td colspan="5">No data</td></tr>'}</tbody>
</table>

<h2>End Barcode Confidence</h2>
<p>Total reads with end confidence: {len(end_confs)}</p>
<p>Mean end confidence: {sum(end_confs) / len(end_confs):.4f}</p>
"""

    data = {
        "bc_confs": {k: {"mean": sum(v) / len(v), "count": len(v)} for k, v in bc_confs.items()},
        "end_conf_count": len(end_confs),
        "end_conf_mean": sum(end_confs) / len(end_confs) if end_confs else 0,
    }
    html = _base_html(
        title=f"Confidence - {experiment_id}",
        body=body,
        data_json=json.dumps(data, indent=2),
        include_d3=True,
    )
    path.write_text(html, encoding="utf-8")


def _write_separation(
    path: Path,
    experiment_id: str,
    db_path: Path,
) -> None:
    """Generate separation.html -- barcode separation metrics."""
    # Load per-barcode edit distance stats from DB
    conn = sqlite3.connect(db_path)
    rows = conn.execute("""
        SELECT bc_start_id, AVG(ed), MIN(ed), MAX(ed), COUNT(*)
        FROM Reads
        WHERE bc_start_id IS NOT NULL AND ed IS NOT NULL
        GROUP BY bc_start_id
    """).fetchall()
    conn.close()

    table_rows = ""
    sep_data: dict[str, dict] = {}
    for bc_id, avg_ed, min_ed, max_ed, count in rows:
        table_rows += (
            f"<tr><td>{_esc(str(bc_id))}</td>"
            f"<td>{count}</td>"
            f"<td>{avg_ed:.2f}</td>"
            f"<td>{min_ed}</td>"
            f"<td>{max_ed}</td></tr>\n"
        )
        sep_data[str(bc_id)] = {
            "avg_ed": round(float(avg_ed), 2),
            "min_ed": int(min_ed),
            "max_ed": int(max_ed),
            "count": int(count),
        }

    body = f"""
<h1>Barcode Separation: {_esc(experiment_id)}</h1>

<h2>Edit Distance by Barcode</h2>
<table>
<thead><tr><th>Barcode</th><th>Count</th><th>Mean ED</th><th>Min ED</th><th>Max ED</th></tr></thead>
<tbody>{table_rows if table_rows else '<tr><td colspan="5">No data</td></tr>'}</tbody>
</table>
"""

    html = _base_html(
        title=f"Separation - {experiment_id}",
        body=body,
        data_json=json.dumps({"separation": sep_data}, indent=2),
        include_d3=True,
    )
    path.write_text(html, encoding="utf-8")


def _write_thresholds(
    path: Path,
    experiment_id: str,
    db_path: Path,
) -> None:
    """Generate thresholds.html -- threshold impact analysis."""
    from calibrate_viz.confusion import compute_threshold_impact

    # Compute threshold impact at a few reference thresholds
    thresholds_data: list[dict] = []
    for start_min in [0.5, 0.6, 0.7, 0.8, 0.9]:
        for fl_thresh in [0.5, 0.6, 0.7, 0.8, 0.9]:
            impact = compute_threshold_impact(db_path, start_min, fl_thresh)
            if impact:
                thresholds_data.append({
                    "start_min": start_min,
                    "fl_thresh": fl_thresh,
                    "counts": impact,
                })

    # Build a summary table of a few key thresholds
    table_rows = ""
    for entry in thresholds_data:
        counts = entry["counts"]
        total = sum(counts.values())
        fl_pct = counts.get("full_length", 0) / total * 100 if total else 0
        table_rows += (
            f"<tr><td>{entry['start_min']:.1f}</td>"
            f"<td>{entry['fl_thresh']:.1f}</td>"
            f"<td>{total}</td>"
            f"<td>{fl_pct:.1f}%</td></tr>\n"
        )

    body = f"""
<h1>Threshold Analysis: {_esc(experiment_id)}</h1>

<h2>Threshold Impact Grid</h2>
<table>
<thead><tr><th>Start BC Min</th><th>FL Threshold</th><th>Total Reads</th><th>Full-Length %</th></tr></thead>
<tbody>{table_rows if table_rows else '<tr><td colspan="4">No data</td></tr>'}</tbody>
</table>
"""

    html = _base_html(
        title=f"Thresholds - {experiment_id}",
        body=body,
        data_json=json.dumps({"thresholds": thresholds_data}, indent=2),
        include_d3=True,
    )
    path.write_text(html, encoding="utf-8")


def _write_svg(path: Path, kde_data: dict, x_label: str) -> None:
    """Write an SVG file from KDE data."""
    svg = _build_svg(kde_data, x_label)
    path.write_text(svg, encoding="utf-8")


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------


def _serialize_kde(kde_data: dict) -> dict:
    """Convert KDE data to JSON-safe format (lists instead of numpy arrays)."""
    result: dict = {"groups": {}, "peaks": {}}
    for label, g in kde_data.get("groups", {}).items():
        result["groups"][label] = {
            "x": [float(v) for v in g.get("x", [])],
            "y": [float(v) for v in g.get("y", [])],
            "count": g.get("count", 0),
        }
    for label, peaks in kde_data.get("peaks", {}).items():
        result["peaks"][label] = [float(p) for p in peaks]
    return result


# ---------------------------------------------------------------------------
# Main export function
# ---------------------------------------------------------------------------


def export_calibration_report(
    db_path: Path,
    output_dir: Path,
    experiment_id: str | None = None,
) -> Path:
    """Export a self-contained calibration report as static HTML + SVG.

    Parameters
    ----------
    db_path : Path
        SQLite database to export from.
    output_dir : Path
        Directory to write output files into. Created if needed.
    experiment_id : str or None
        Label for the report. Defaults to db filename stem.

    Returns
    -------
    Path
        The output directory.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    if experiment_id is None:
        experiment_id = db_path.stem

    # Compute data for all pages
    summary = compute_experiment_summary(db_path)

    # Compute KDE distributions for signal duration and read length
    signal_kde = compute_grouped_kde(db_path, "signal_duration_s")
    readlen_kde = compute_grouped_kde(db_path, "readlen")

    # Generate index.html
    _write_index(output_dir / "index.html", experiment_id, summary, signal_kde, readlen_kde)

    # Generate distributions.html
    _write_distributions(output_dir / "distributions.html", experiment_id, signal_kde, readlen_kde)

    # Generate confidence.html
    _write_confidence(output_dir / "confidence.html", experiment_id, db_path)

    # Generate separation.html
    _write_separation(output_dir / "separation.html", experiment_id, db_path)

    # Generate thresholds.html
    _write_thresholds(output_dir / "thresholds.html", experiment_id, db_path)

    # Generate SVG figures
    _write_svg(figures_dir / "signal_length_kde.svg", signal_kde, "Signal Duration (s)")
    _write_svg(figures_dir / "read_length_kde.svg", readlen_kde, "Read Length (bp)")

    return output_dir
