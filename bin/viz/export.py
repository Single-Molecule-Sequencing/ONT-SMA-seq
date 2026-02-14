"""Static HTML export engine for SMA-seq config visualizer.

Generates 7 self-contained HTML files from an experiment directory:
  1. index.html          - Overview with construct diagram
  2. sample_sheet.html   - Sample sheet table
  3. barcodes.html       - Barcode sequences and arrangement
  4. construct.html      - Construct diagram + truncation classes
  5. targets.html        - Target reference sequences
  6. assumptions.html    - Documented assumptions and quality metrics
  7. validation_report.html - Cross-validation results
"""

from __future__ import annotations

import html
import json
from datetime import datetime, timezone
from pathlib import Path

from viz.config_store import ConfigStore

# Try importing barcodes module for GC/RC calculations
try:
    import barcodes as _barcodes_mod
except ImportError:
    _barcodes_mod = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_PAGES = [
    ("index.html", "Overview"),
    ("sample_sheet.html", "Sample Sheet"),
    ("barcodes.html", "Barcodes"),
    ("construct.html", "Construct"),
    ("targets.html", "Targets"),
    ("assumptions.html", "Assumptions"),
    ("validation_report.html", "Validation Report"),
]


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _inline_css() -> str:
    """Read and return the contents of static/style.css."""
    css_path = Path(__file__).resolve().parent / "static" / "style.css"
    if css_path.exists():
        return css_path.read_text()
    return ""


def _inline_js() -> str:
    """Read and return the contents of static/construct.js."""
    js_path = Path(__file__).resolve().parent / "static" / "construct.js"
    if js_path.exists():
        return js_path.read_text()
    return ""


def _nav_links() -> str:
    """Return an HTML string with navigation links to all 7 pages."""
    links: list[str] = []
    for filename, label in _PAGES:
        links.append(f'<a href="{filename}">{_esc(label)}</a>')
    return " | ".join(links)


def _esc(text: str) -> str:
    """HTML-escape a string."""
    return html.escape(str(text))


def _gc_content(seq: str) -> float:
    """Compute GC content as a percentage."""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in ("G", "C"))
    return round((gc / len(seq)) * 100.0, 2)


def _base_html(
    title: str,
    body: str,
    data_json: str = "",
    nav_links: str = "",
    include_d3: bool = False,
    include_construct_js: bool = False,
) -> str:
    """Wrap body content in a complete self-contained HTML document.

    Parameters
    ----------
    title : str
        Page title.
    body : str
        HTML body content.
    data_json : str
        Optional JSON data block to embed inline.
    nav_links : str
        HTML string for navigation links.
    include_d3 : bool
        Whether to include the D3.js CDN script tag.
    include_construct_js : bool
        Whether to inline construct.js.
    """
    css = _inline_css()
    d3_tag = ""
    if include_d3:
        d3_tag = '<script src="https://d3js.org/d3.v7.min.js"></script>'

    js_block = ""
    if include_construct_js:
        construct_js = _inline_js()
        js_block = f"<script>\n{construct_js}\n</script>"

    data_block = ""
    if data_json:
        data_block = (
            f'<script id="page-data" type="application/json">\n'
            f"{data_json}\n"
            f"</script>"
        )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{_esc(title)}</title>
<style>
{css}

/* Export-specific styles */
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
.valid {{ color: #2e7d32; font-weight: bold; }}
.invalid {{ color: #c62828; font-weight: bold; }}
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
<nav class="export-nav">{nav_links}</nav>
{body}
{data_block}
{js_block}
</body>
</html>"""


# ---------------------------------------------------------------------------
# Page generators
# ---------------------------------------------------------------------------


def _gen_index_entries_summary(entries: list) -> str:
    """Generate a brief list of sample sheet aliases for the index page."""
    if not entries:
        return "<p>No sample sheet entries.</p>"
    items = "\n".join(
        f"<li>{_esc(e.alias)} ({_esc(e.barcode)})</li>" for e in entries
    )
    return f"<ul>\n{items}\n</ul>"


def _gen_index_targets_summary(targets: list) -> str:
    """Generate a brief list of targets for the index page."""
    if not targets:
        return "<p>No targets.</p>"
    items = "\n".join(
        f"<li>{_esc(t.tgt_id)} ({t.length} bp, GC {t.gc_content:.1f}%)</li>"
        for t in targets
    )
    return f"<ul>\n{items}\n</ul>"


def _gen_index(
    store: ConfigStore,
    entries: list,
    targets: list,
    errors: list[str],
    timestamp: str,
    construct_data: dict,
) -> str:
    """Generate index.html - Overview page."""
    cfg = store.experiment_config
    description = cfg.description or "(no description)"
    demux_mode = cfg.demultiplexing.mode
    valid_status = (
        '<span class="valid">All valid</span>'
        if not errors
        else f'<span class="invalid">{len(errors)} error(s)</span>'
    )

    body = f"""
<h1>SMA-seq Experiment Export</h1>
<p><strong>Description:</strong> {_esc(description)}</p>
<p><strong>Directory:</strong> <code>{_esc(str(store.dir))}</code></p>
<p><strong>Exported:</strong> {_esc(timestamp)}</p>

<h2>Summary</h2>
<div class="summary-grid">
  <div class="summary-card">
    <div class="value">{len(entries)}</div>
    <div class="label">Sample Sheet Entries</div>
  </div>
  <div class="summary-card">
    <div class="value">{len(targets)}</div>
    <div class="label">Targets</div>
  </div>
  <div class="summary-card">
    <div class="value">{_esc(demux_mode)}</div>
    <div class="label">Demux Mode</div>
  </div>
  <div class="summary-card">
    <div class="value">{valid_status}</div>
    <div class="label">Validation Status</div>
  </div>
</div>

<h2>Sample Sheet Entries</h2>
{_gen_index_entries_summary(entries)}

<h2>Targets</h2>
{_gen_index_targets_summary(targets)}

<h2>Construct Diagram</h2>
<div id="construct-diagram"></div>
<div id="construct-detail"></div>

<h2>Truncation Classes</h2>
<div id="truncation-ladder"></div>

<script>
document.addEventListener("DOMContentLoaded", function() {{
  var dataEl = document.getElementById("page-data");
  if (dataEl && typeof renderConstruct === "function") {{
    var data = JSON.parse(dataEl.textContent);
    renderConstruct(data, "#construct-diagram", "#construct-detail");
    renderTruncationLadder(data, "#truncation-ladder");
  }}
}});
</script>
"""
    data_json = json.dumps(construct_data, indent=2)
    return _base_html(
        title="SMA-seq Export - Overview",
        body=body,
        data_json=data_json,
        nav_links=_nav_links(),
        include_d3=True,
        include_construct_js=True,
    )


def _gen_sample_sheet(entries: list) -> str:
    """Generate sample_sheet.html - Sample sheet table."""
    rows = ""
    for e in entries:
        d = e.model_dump()
        rows += (
            f"<tr>"
            f"<td>{_esc(d.get('barcode', ''))}</td>"
            f"<td>{_esc(d.get('upstream_barcode', ''))}</td>"
            f"<td>{_esc(d.get('downstream_barcode', ''))}</td>"
            f"<td>{_esc(d.get('alias', ''))}</td>"
            f"<td>{_esc(d.get('kit', ''))}</td>"
            f"<td>{_esc(d.get('type', ''))}</td>"
            f"</tr>\n"
        )

    body = f"""
<h1>Sample Sheet</h1>
<p>{len(entries)} entries</p>
<table>
<thead>
<tr>
  <th>Barcode Pair</th>
  <th>Upstream</th>
  <th>Downstream</th>
  <th>Alias</th>
  <th>Kit</th>
  <th>Type</th>
</tr>
</thead>
<tbody>
{rows}
</tbody>
</table>
"""
    return _base_html(
        title="SMA-seq Export - Sample Sheet",
        body=body,
        nav_links=_nav_links(),
    )


def _gen_barcodes(store: ConfigStore, entries: list) -> str:
    """Generate barcodes.html - Barcode table and arrangement config."""
    # Collect used barcode IDs
    used_ids: set[str] = set()
    for entry in entries:
        if entry.upstream_barcode:
            used_ids.add(entry.upstream_barcode)
        if entry.downstream_barcode:
            used_ids.add(entry.downstream_barcode)

    # Resolve sequences
    custom_seqs = store.read_barcode_sequences()
    all_barcodes: dict[str, str] = {}
    if _barcodes_mod is not None:
        all_barcodes.update(_barcodes_mod.BARCODES)
    all_barcodes.update(custom_seqs)

    rows = ""
    for bc_id in sorted(used_ids):
        seq = all_barcodes.get(bc_id, "")
        rc = ""
        if seq and _barcodes_mod is not None:
            rc = _barcodes_mod.reverse_complement(seq)
        gc = _gc_content(seq) if seq else 0.0
        rows += (
            f"<tr>"
            f"<td>{_esc(bc_id)}</td>"
            f"<td><code>{_esc(seq)}</code></td>"
            f"<td><code>{_esc(rc)}</code></td>"
            f"<td>{gc:.1f}%</td>"
            f"</tr>\n"
        )

    # Arrangement config
    arrangement = store.read_arrangement()
    arrangement_block = ""
    if arrangement:
        arr_json = json.dumps(arrangement.model_dump(), indent=2)
        arrangement_block = f"""
<h2>Arrangement Configuration</h2>
<pre><code>{_esc(arr_json)}</code></pre>
"""

    body = f"""
<h1>Barcodes</h1>
<h2>Used Barcodes</h2>
<table>
<thead>
<tr>
  <th>ID</th>
  <th>Sequence</th>
  <th>Reverse Complement</th>
  <th>GC%</th>
</tr>
</thead>
<tbody>
{rows}
</tbody>
</table>
{arrangement_block}
"""
    return _base_html(
        title="SMA-seq Export - Barcodes",
        body=body,
        nav_links=_nav_links(),
    )


def _gen_construct(store: ConfigStore, construct_data: dict) -> str:
    """Generate construct.html - Construct diagram, barcode assignments, truncation rules."""
    cfg = store.experiment_config
    demux = cfg.demultiplexing
    truncation = cfg.truncation

    # Barcode assignment mode and pairs table
    pairs_rows = ""
    for pair in demux.pairs:
        d = pair.model_dump()
        pairs_rows += (
            f"<tr>"
            f"<td>{_esc(d.get('start', ''))}</td>"
            f"<td>{_esc(d.get('end', ''))}</td>"
            f"<td>{_esc(d.get('alias', ''))}</td>"
            f"</tr>\n"
        )

    pairs_table = ""
    if pairs_rows:
        pairs_table = f"""
<h2>Barcode Pairs</h2>
<table>
<thead><tr><th>Start</th><th>End</th><th>Alias</th></tr></thead>
<tbody>{pairs_rows}</tbody>
</table>
"""

    # Truncation rules table
    rules = truncation.rules.model_dump()
    rules_rows = ""
    for key, action in rules.items():
        rules_rows += f"<tr><td>{_esc(key)}</td><td>{_esc(action)}</td></tr>\n"

    body = f"""
<h1>Construct</h1>

<h2>Construct Diagram</h2>
<div id="construct-diagram"></div>
<div id="construct-detail"></div>

<h2>Truncation Classes</h2>
<div id="truncation-ladder"></div>

<h2>Barcode Assignment</h2>
<p><strong>Mode:</strong> {_esc(demux.mode)}</p>
{pairs_table}

<h2>Truncation Rules</h2>
<table>
<thead><tr><th>Class</th><th>Action</th></tr></thead>
<tbody>{rules_rows}</tbody>
</table>

<script>
document.addEventListener("DOMContentLoaded", function() {{
  var dataEl = document.getElementById("page-data");
  if (dataEl && typeof renderConstruct === "function") {{
    var data = JSON.parse(dataEl.textContent);
    renderConstruct(data, "#construct-diagram", "#construct-detail");
    renderTruncationLadder(data, "#truncation-ladder");
  }}
}});
</script>
"""
    data_json = json.dumps(construct_data, indent=2)
    return _base_html(
        title="SMA-seq Export - Construct",
        body=body,
        data_json=data_json,
        nav_links=_nav_links(),
        include_d3=True,
        include_construct_js=True,
    )


def _gen_targets(targets: list) -> str:
    """Generate targets.html - Target reference table."""
    rows = ""
    for t in targets:
        d = t.model_dump()
        seq = d.get("sequence", "")
        display_seq = seq[:200] + "..." if len(seq) > 200 else seq
        gc = d.get("gc_content", 0.0)
        rows += (
            f"<tr>"
            f"<td>{_esc(d.get('tgt_id', ''))}</td>"
            f"<td>{d.get('length', 0)}</td>"
            f"<td>{gc:.2f}%</td>"
            f"<td><code>{_esc(display_seq)}</code></td>"
            f"</tr>\n"
        )

    body = f"""
<h1>Targets</h1>
<p>{len(targets)} target reference(s)</p>
<table>
<thead>
<tr>
  <th>ID</th>
  <th>Length</th>
  <th>GC%</th>
  <th>Sequence</th>
</tr>
</thead>
<tbody>
{rows}
</tbody>
</table>
"""
    return _base_html(
        title="SMA-seq Export - Targets",
        body=body,
        nav_links=_nav_links(),
    )


def _gen_assumptions(store: ConfigStore) -> str:
    """Generate assumptions.html - Assumptions table and quality metrics."""
    cfg = store.experiment_config
    assumptions = cfg.assumptions
    quality = cfg.quality
    classification = cfg.classification

    assumption_rows = ""
    for a in assumptions:
        d = a.model_dump()
        assumption_rows += (
            f"<tr>"
            f"<td><code>{_esc(d.get('key', ''))}</code></td>"
            f"<td>{_esc(d.get('text', ''))}</td>"
            f"<td>{_esc(d.get('why', ''))}</td>"
            f"</tr>\n"
        )

    body = f"""
<h1>Assumptions</h1>

<h2>Documented Assumptions</h2>
<table>
<thead>
<tr><th>Key</th><th>Assumption</th><th>Rationale</th></tr>
</thead>
<tbody>
{assumption_rows}
</tbody>
</table>

<h2>Quality Metrics</h2>
<table>
<tr><th>Metric</th><th>Formula</th></tr>
<tr><td>Q_bc (barcode quality)</td><td><code>{_esc(quality.q_bc)}</code></td></tr>
<tr><td>Q_ld (alignment quality)</td><td><code>{_esc(quality.q_ld)}</code></td></tr>
</table>

<h2>Classification Parameters</h2>
<table>
<tr><th>Parameter</th><th>Value</th></tr>
<tr><td>Barcode search window</td><td>{classification.barcode_search_window}</td></tr>
<tr><td>Confidence formula</td><td><code>{_esc(classification.confidence_formula)}</code></td></tr>
<tr><td>Ambiguity triggers full construct</td><td>{"Yes" if classification.ambiguity_triggers_full_construct else "No"}</td></tr>
</table>
"""
    return _base_html(
        title="SMA-seq Export - Assumptions",
        body=body,
        nav_links=_nav_links(),
    )


def _gen_validation_report(errors: list[str]) -> str:
    """Generate validation_report.html - Validation results."""
    if errors:
        error_items = "\n".join(
            f"<li>{_esc(e)}</li>" for e in errors
        )
        content = f"""
<p class="invalid">Found {len(errors)} validation error(s):</p>
<ul>
{error_items}
</ul>
"""
    else:
        content = '<p class="valid">All valid - no errors found.</p>'

    body = f"""
<h1>Validation Report</h1>
{content}
"""
    return _base_html(
        title="SMA-seq Export - Validation Report",
        body=body,
        nav_links=_nav_links(),
    )


# ---------------------------------------------------------------------------
# Main export function
# ---------------------------------------------------------------------------


def export_experiment(
    experiment_dir: Path | str,
    output_dir: Path | str | None = None,
) -> Path:
    """Export an experiment as 7 self-contained HTML files.

    Parameters
    ----------
    experiment_dir : Path | str
        Root directory of the experiment (must contain config files).
    output_dir : Path | str | None
        Output directory for the HTML files. If *None*, creates
        ``exports/export_<timestamp>/`` inside experiment_dir.

    Returns
    -------
    Path
        The output directory path.
    """
    experiment_dir = Path(experiment_dir).resolve()
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

    if output_dir is None:
        output_dir = experiment_dir / "exports" / f"export_{timestamp}"
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build the ConfigStore and read all data
    store = ConfigStore(experiment_dir)
    entries = store.read_sample_sheet()
    targets = store.read_targets()
    errors = store.validate()
    cfg = store.experiment_config

    # Construct data for D3 diagrams (matches API /api/construct shape)
    construct_data = {
        "construct": cfg.construct.model_dump(),
        "demultiplexing": cfg.demultiplexing.model_dump(),
        "truncation": cfg.truncation.model_dump(),
        "classification": cfg.classification.model_dump(),
    }

    # Generate and write all pages
    pages = {
        "index.html": _gen_index(
            store, entries, targets, errors, timestamp, construct_data
        ),
        "sample_sheet.html": _gen_sample_sheet(entries),
        "barcodes.html": _gen_barcodes(store, entries),
        "construct.html": _gen_construct(store, construct_data),
        "targets.html": _gen_targets(targets),
        "assumptions.html": _gen_assumptions(store),
        "validation_report.html": _gen_validation_report(errors),
    }

    for filename, html_content in pages.items():
        (output_dir / filename).write_text(html_content, encoding="utf-8")

    return output_dir


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main() -> None:
    import sys

    if len(sys.argv) < 2:
        print("Usage: python -m bin.viz.export <experiment_dir> [--output <dir>]")
        sys.exit(1)

    experiment_dir = Path(sys.argv[1]).resolve()
    output_dir = None
    if "--output" in sys.argv:
        idx = sys.argv.index("--output")
        output_dir = Path(sys.argv[idx + 1]).resolve()

    result = export_experiment(experiment_dir, output_dir)
    print(f"[export] Static HTML exported to: {result}")


if __name__ == "__main__":
    main()
