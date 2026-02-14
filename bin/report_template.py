"""HTML report generator for barcode classification results.

Generates an interactive dark-themed HTML report from analysis data produced
by ``report_analysis.analyze_classification``.  The HTML is fully self-contained
with no external dependencies.

The single public entry-point is :func:`generate_html`.
"""

from __future__ import annotations

import json
from typing import Any


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_html(analysis: dict[str, Any], exp_metadata: dict[str, Any]) -> str:
    """Build a standalone HTML report from classification analysis data.

    Parameters
    ----------
    analysis : dict
        Output of ``report_analysis.analyze_classification``.  Expected keys:

        * ``summary`` -- {total, full_length, truncated, matched, unmatched}
        * ``per_target_stats`` -- tgt_id -> stat dict
        * ``confidence_distributions`` -- {start_confs, end_confs}
        * ``pair_matrix`` -- {counts (Counter), used_barcodes (sorted list)}
        * ``reads_table`` -- list of per-read summary dicts
        * ``detailed_reads`` -- list of per-read detail dicts
        * ``barcode_info`` -- bc_id -> {fwd, rc}
        * ``pairing_table`` -- list of {upstream, downstream, alias}
        * ``references`` -- optional tgt_id -> {seq, length}
        * ``flank_front``, ``flank_rear`` -- optional flanking sequences

    exp_metadata : dict
        Experiment metadata with keys ``exp_id``, ``flow_cell_id``,
        ``sample_id`` (any may be None).

    Returns
    -------
    str
        Complete standalone HTML document.
    """
    # Unpack analysis sections
    summary = analysis["summary"]
    per_target_stats = analysis["per_target_stats"]
    conf_dists = analysis["confidence_distributions"]
    pair_matrix = analysis["pair_matrix"]
    reads_table = analysis["reads_table"]
    detailed_reads = analysis["detailed_reads"]
    barcode_info = analysis["barcode_info"]
    pairing_table = analysis["pairing_table"]
    references = analysis.get("references") or {}
    flank_front = analysis.get("flank_front")
    flank_rear = analysis.get("flank_rear")

    # Metadata
    exp_id = exp_metadata.get("exp_id") or "unknown"
    flow_cell_id = exp_metadata.get("flow_cell_id") or "unknown"
    sample_id = exp_metadata.get("sample_id")

    # Build set of expected pairs from pairing table for matrix highlighting
    expected_pairs: set[tuple[str, str]] = set()
    pair_alias_map: dict[tuple[str, str], str] = {}
    for entry in pairing_table:
        key = (entry["upstream"], entry["downstream"])
        expected_pairs.add(key)
        pair_alias_map[key] = entry["alias"]

    # Subtitle parts
    subtitle_parts = [f"Experiment: {exp_id}", f"Flow Cell: {flow_cell_id}"]
    if sample_id:
        subtitle_parts.append(f"Sample: {sample_id}")
    subtitle_parts.append(f"{summary['total']} reads")
    subtitle = " &middot; ".join(subtitle_parts)

    # Percentage helper
    def _pct(num: int | float, denom: int | float) -> str:
        if denom == 0:
            return "0"
        return f"{num / denom * 100:.0f}"

    # Confidence badge helper
    def _conf_cls(val: float | None) -> str:
        if val is None:
            return "conf-low"
        if val >= 0.9:
            return "conf-high"
        if val >= 0.75:
            return "conf-med"
        return "conf-low"

    # ------------------------------------------------------------------
    # Start building HTML
    # ------------------------------------------------------------------
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>SMA-seq Barcode Classification Report</title>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{
    font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
    background: #0f1117;
    color: #e0e0e0;
    line-height: 1.6;
    padding: 20px;
}}
.container {{ max-width: 1400px; margin: 0 auto; }}
h1 {{
    font-size: 28px;
    font-weight: 700;
    margin-bottom: 8px;
    background: linear-gradient(90deg, #60a5fa, #a78bfa);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
}}
h2 {{
    font-size: 20px;
    font-weight: 600;
    color: #93c5fd;
    margin: 30px 0 15px;
    border-bottom: 1px solid #1e293b;
    padding-bottom: 8px;
}}
h3 {{
    font-size: 16px;
    font-weight: 600;
    color: #c4b5fd;
    margin: 20px 0 10px;
}}
.subtitle {{ color: #94a3b8; font-size: 14px; margin-bottom: 25px; }}

/* Cards */
.card-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 15px;
    margin-bottom: 25px;
}}
.card {{
    background: #1e1e2e;
    border-radius: 12px;
    padding: 20px;
    border: 1px solid #2d2d44;
}}
.card-label {{ font-size: 12px; color: #94a3b8; text-transform: uppercase; letter-spacing: 1px; }}
.card-value {{ font-size: 32px; font-weight: 700; margin: 5px 0; }}
.card-detail {{ font-size: 13px; color: #64748b; }}
.card-value.green {{ color: #4ade80; }}
.card-value.amber {{ color: #fbbf24; }}
.card-value.red {{ color: #f87171; }}
.card-value.blue {{ color: #60a5fa; }}

/* Construct diagram */
.construct-diagram {{
    background: #1e1e2e;
    border-radius: 12px;
    padding: 20px;
    border: 1px solid #2d2d44;
    margin-bottom: 25px;
    overflow-x: auto;
}}
.construct-row {{
    display: flex;
    align-items: center;
    margin: 8px 0;
    font-family: 'Cascadia Code', 'Fira Code', monospace;
    font-size: 12px;
    white-space: nowrap;
}}
.construct-block {{
    padding: 6px 10px;
    border-radius: 6px;
    margin: 0 2px;
    text-align: center;
    min-width: 60px;
    position: relative;
}}
.construct-block.adapter {{ background: #374151; color: #9ca3af; }}
.construct-block.flank {{ background: #1e3a5f; color: #60a5fa; }}
.construct-block.barcode {{ background: #3b2f63; color: #c4b5fd; }}
.construct-block.target {{ background: #1a3f2e; color: #4ade80; flex: 1; min-width: 200px; }}
.construct-block small {{ display: block; font-size: 10px; opacity: 0.7; }}
.arrow {{ color: #475569; margin: 0 4px; font-size: 18px; }}
.search-window {{
    border: 2px dashed #f59e0b;
    border-radius: 8px;
    padding: 4px;
    display: flex;
    align-items: center;
}}
.search-label {{
    font-size: 10px;
    color: #f59e0b;
    position: absolute;
    top: -18px;
    left: 0;
    white-space: nowrap;
}}

/* Tables */
table {{
    width: 100%;
    border-collapse: collapse;
    margin: 10px 0;
    font-size: 13px;
}}
th {{
    background: #1a1a2e;
    color: #93c5fd;
    padding: 10px 12px;
    text-align: left;
    font-weight: 600;
    position: sticky;
    top: 0;
    z-index: 10;
    cursor: pointer;
    user-select: none;
}}
th:hover {{ background: #222244; }}
td {{
    padding: 8px 12px;
    border-bottom: 1px solid #1e293b;
}}
tr:hover td {{ background: #1a1a2e; }}
.table-container {{
    background: #1e1e2e;
    border-radius: 12px;
    border: 1px solid #2d2d44;
    overflow: hidden;
    margin-bottom: 25px;
}}
.table-scroll {{
    max-height: 500px;
    overflow-y: auto;
}}

/* Confidence badges */
.conf-badge {{
    display: inline-block;
    padding: 2px 8px;
    border-radius: 10px;
    font-size: 12px;
    font-weight: 600;
}}
.conf-high {{ background: #065f46; color: #34d399; }}
.conf-med {{ background: #78350f; color: #fbbf24; }}
.conf-low {{ background: #7f1d1d; color: #fca5a5; }}

/* Sequence viewer */
.seq-viewer {{
    background: #111827;
    border-radius: 8px;
    padding: 15px;
    margin: 10px 0;
    overflow-x: auto;
    font-family: 'Cascadia Code', 'Fira Code', monospace;
    font-size: 12px;
    line-height: 1.8;
    position: relative;
}}
.seq-ruler {{
    color: #475569;
    margin-bottom: 5px;
    user-select: none;
}}
.seq-text {{
    letter-spacing: 0.5px;
    word-break: break-all;
}}
.seq-text .bc-start {{ background: #5b21b6; color: #e9d5ff; border-radius: 2px; }}
.seq-text .bc-end {{ background: #9d174d; color: #fce7f3; border-radius: 2px; }}
.seq-text .flank {{ background: #1e3a5f; color: #93c5fd; border-radius: 2px; }}
.seq-text .target-region {{ background: #064e3b; color: #6ee7b7; border-radius: 2px; }}
.seq-text .search-zone {{ text-decoration: underline; text-decoration-color: #f59e0b; text-decoration-style: dashed; }}

/* Read detail panel */
.read-detail {{
    background: #1e1e2e;
    border-radius: 12px;
    border: 1px solid #2d2d44;
    padding: 20px;
    margin-bottom: 20px;
}}
.read-header {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    margin-bottom: 15px;
    flex-wrap: wrap;
    gap: 10px;
}}
.read-meta {{
    display: flex;
    gap: 15px;
    flex-wrap: wrap;
}}
.meta-chip {{
    display: inline-flex;
    align-items: center;
    gap: 5px;
    padding: 4px 10px;
    border-radius: 20px;
    font-size: 12px;
    font-weight: 500;
}}
.meta-chip.full {{ background: #065f46; color: #34d399; }}
.meta-chip.trunc {{ background: #78350f; color: #fbbf24; }}
.meta-chip.matched {{ background: #1e3a5f; color: #60a5fa; }}
.meta-chip.unmatched {{ background: #7f1d1d; color: #fca5a5; }}
.meta-chip.bin {{ background: #1a1a2e; color: #94a3b8; border: 1px solid #334155; }}

/* Competitor tables */
.competitor-grid {{
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 15px;
    margin: 15px 0;
}}
@media (max-width: 900px) {{
    .competitor-grid {{ grid-template-columns: 1fr; }}
}}
.competitor-table {{ font-size: 12px; }}
.competitor-table th {{ padding: 6px 10px; }}
.competitor-table td {{ padding: 5px 10px; }}
.competitor-table tr:first-child td {{ font-weight: 700; }}

/* Read structure diagram (per read) */
.read-structure {{
    display: flex;
    align-items: stretch;
    margin: 10px 0;
    height: 40px;
    border-radius: 6px;
    overflow: hidden;
    font-size: 10px;
    font-weight: 600;
}}
.read-structure .segment {{
    display: flex;
    align-items: center;
    justify-content: center;
    padding: 0 4px;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    min-width: 2px;
}}
.read-structure .seg-start-bc {{ background: #5b21b6; color: #e9d5ff; }}
.read-structure .seg-flank {{ background: #1e3a5f; color: #60a5fa; }}
.read-structure .seg-target {{ background: #064e3b; color: #6ee7b7; flex: 1; }}
.read-structure .seg-end-bc {{ background: #9d174d; color: #fce7f3; }}
.read-structure .seg-unknown {{ background: #374151; color: #9ca3af; flex: 1; }}

/* Tabs */
.tab-bar {{
    display: flex;
    gap: 5px;
    margin-bottom: 20px;
    border-bottom: 1px solid #2d2d44;
    padding-bottom: 5px;
    flex-wrap: wrap;
}}
.tab-btn {{
    padding: 8px 16px;
    border-radius: 8px 8px 0 0;
    border: 1px solid transparent;
    background: transparent;
    color: #94a3b8;
    cursor: pointer;
    font-size: 14px;
    font-weight: 500;
    transition: all 0.2s;
}}
.tab-btn:hover {{ background: #1e1e2e; color: #e0e0e0; }}
.tab-btn.active {{ background: #1e1e2e; color: #60a5fa; border-color: #2d2d44; border-bottom-color: #1e1e2e; }}
.tab-panel {{ display: none; }}
.tab-panel.active {{ display: block; }}

/* Pair matrix */
.pair-matrix {{
    display: grid;
    gap: 2px;
    margin: 15px 0;
}}
.pair-cell {{
    padding: 8px;
    text-align: center;
    border-radius: 4px;
    font-size: 12px;
    font-weight: 600;
}}
.pair-cell.header {{ background: #1a1a2e; color: #93c5fd; }}
.pair-cell.high {{ background: #065f46; color: #34d399; }}
.pair-cell.med {{ background: #78350f; color: #fbbf24; }}
.pair-cell.low {{ background: #1a1a2e; color: #475569; }}
.pair-cell.expected {{ border: 2px solid #60a5fa; }}

/* Histograms */
.histogram {{
    display: flex;
    align-items: flex-end;
    gap: 2px;
    height: 120px;
    padding: 10px 0;
}}
.hist-bar {{
    flex: 1;
    min-width: 8px;
    border-radius: 3px 3px 0 0;
    position: relative;
    cursor: pointer;
    transition: opacity 0.2s;
}}
.hist-bar:hover {{ opacity: 0.8; }}
.hist-bar .tooltip {{
    display: none;
    position: absolute;
    bottom: 100%;
    left: 50%;
    transform: translateX(-50%);
    background: #1a1a2e;
    border: 1px solid #334155;
    padding: 4px 8px;
    border-radius: 4px;
    font-size: 11px;
    white-space: nowrap;
    z-index: 20;
}}
.hist-bar:hover .tooltip {{ display: block; }}
.hist-labels {{
    display: flex;
    gap: 2px;
    font-size: 10px;
    color: #475569;
}}
.hist-labels span {{ flex: 1; text-align: center; min-width: 8px; }}

/* Filter controls */
.filter-bar {{
    display: flex;
    gap: 10px;
    margin-bottom: 15px;
    align-items: center;
    flex-wrap: wrap;
}}
.filter-bar select, .filter-bar input {{
    background: #111827;
    border: 1px solid #334155;
    color: #e0e0e0;
    padding: 6px 12px;
    border-radius: 6px;
    font-size: 13px;
}}
.filter-bar label {{ font-size: 13px; color: #94a3b8; }}

/* Legend */
.legend {{
    display: flex;
    gap: 15px;
    flex-wrap: wrap;
    margin: 10px 0;
    font-size: 12px;
}}
.legend-item {{
    display: flex;
    align-items: center;
    gap: 5px;
}}
.legend-swatch {{
    width: 14px;
    height: 14px;
    border-radius: 3px;
}}

.note {{ font-size: 12px; color: #64748b; font-style: italic; margin: 10px 0; }}
</style>
</head>
<body>
<div class="container">

<h1>SMA-seq Barcode Classification Report</h1>
<p class="subtitle">{subtitle}</p>

<!-- Tabs -->
<div class="tab-bar">
    <button class="tab-btn active" onclick="showTab('construct')">Construct Structure</button>
    <button class="tab-btn" onclick="showTab('summary')">Classification Summary</button>
    <button class="tab-btn" onclick="showTab('reads')">All Reads Table</button>
    <button class="tab-btn" onclick="showTab('details')">Read Details</button>
    <button class="tab-btn" onclick="showTab('barcodes')">Barcode Sequences</button>
</div>
"""

    # ==================================================================
    # Tab 1: Construct Structure
    # ==================================================================
    html += """
<!-- Tab 1: Construct Structure -->
<div id="tab-construct" class="tab-panel active">
<h2>Expected Construct Structure</h2>
<p class="note">The classifier searches the first and last 100bp of each read for barcode sequences using semi-global alignment (edlib HW mode).</p>

<div class="construct-diagram">
    <h3>Full-Length Read (both barcodes present)</h3>
    <div class="construct-row" style="margin-top: 15px;">
        <div style="position: relative; display: inline-flex;">
            <span class="search-label">Search window: first 100bp</span>
            <div class="search-window">
"""

    if flank_front is not None:
        html += f'                <div class="construct-block flank"><small>flank_front</small>{flank_front}</div>\n'

    html += '                <div class="construct-block barcode"><small>upstream BC</small>24bp barcode</div>\n'

    if flank_rear is not None:
        html += f'                <div class="construct-block flank"><small>flank_rear</small>{flank_rear}</div>\n'

    html += """            </div>
        </div>
        <span class="arrow">&rarr;</span>
        <div class="construct-block target"><small>target sequence</small>variable length</div>
        <span class="arrow">&rarr;</span>
        <div style="position: relative; display: inline-flex;">
            <span class="search-label" style="left: auto; right: 0;">Search window: last 100bp</span>
            <div class="search-window">
"""

    if flank_rear is not None:
        html += f'                <div class="construct-block flank"><small>RC(flank_rear)</small>{_reverse_complement_simple(flank_rear)}</div>\n'

    html += '                <div class="construct-block barcode" style="background: #831843;"><small>RC(downstream)</small>24bp RC barcode</div>\n'

    if flank_front is not None:
        html += f'                <div class="construct-block flank"><small>RC(flank_front)</small>{_reverse_complement_simple(flank_front)}</div>\n'

    html += """            </div>
        </div>
    </div>

    <h3 style="margin-top: 25px;">Truncated Read (only start barcode)</h3>
    <div class="construct-row" style="margin-top: 15px;">
        <div class="search-window">
"""

    if flank_front is not None:
        html += f'            <div class="construct-block flank"><small>flank_front</small>{flank_front}</div>\n'

    html += '            <div class="construct-block barcode"><small>upstream BC</small>24bp barcode</div>\n'

    if flank_rear is not None:
        html += f'            <div class="construct-block flank"><small>flank_rear</small>{flank_rear}</div>\n'

    html += """        </div>
        <span class="arrow">&rarr;</span>
        <div class="construct-block target" style="min-width: 120px;"><small>partial target</small>truncated</div>
        <span class="arrow" style="color: #f87171;">&#10007; truncated</span>
    </div>
    <p class="note" style="margin-top: 10px;">Truncated reads lack the downstream barcode. The end barcode confidence score reliably detects these (typically &lt;0.75).</p>
</div>

<div class="construct-diagram">
    <h3>Barcode Pairing (Sample Sheet)</h3>
    <table>
        <tr><th>Upstream BC</th><th>Downstream BC</th><th>Alias</th><th>Target Length</th></tr>
"""

    for entry in pairing_table:
        tgt_len_str = ""
        ref_info = references.get(entry["alias"])
        if ref_info is not None:
            tgt_len_str = f"{ref_info['length']}bp"
        else:
            tgt_len_str = "N/A"
        html += f'        <tr><td>{entry["upstream"]}</td><td>{entry["downstream"]}</td><td>{entry["alias"]}</td><td>{tgt_len_str}</td></tr>\n'

    html += """    </table>
    <p class="note" style="margin-top: 10px;">Barcodes are searched in 100bp windows at the start and end of each read.</p>
</div>

<div class="construct-diagram">
    <h3>Classification Algorithm</h3>
    <ol style="font-size: 14px; line-height: 2; padding-left: 20px;">
        <li><strong>Extract segments:</strong> first 100bp (start) and last 100bp (end) of each read</li>
        <li><strong>Align barcodes:</strong> Semi-global alignment (edlib HW mode) of each expected barcode against the segment</li>
        <li><strong>Select best match:</strong> Lowest edit distance wins. Confidence = 1.0 - (ED / 24)</li>
        <li><strong>End barcode uses RC:</strong> The downstream barcode appears as its reverse complement in the read</li>
        <li><strong>Look up pair:</strong> (start_bc, end_bc) matched against sample sheet for target assignment</li>
        <li><strong>Align to target:</strong> Matched reads aligned (NW mode) to assigned reference for quality metrics</li>
    </ol>
</div>
</div>
"""

    # ==================================================================
    # Tab 2: Classification Summary
    # ==================================================================
    total = summary["total"]
    full_length = summary["full_length"]
    truncated = summary["truncated"]
    matched = summary["matched"]
    unmatched = summary["unmatched"]

    html += f"""
<!-- Tab 2: Classification Summary -->
<div id="tab-summary" class="tab-panel">
<h2>Classification Summary</h2>

<!-- Summary Cards -->
<div class="card-grid">
    <div class="card">
        <div class="card-label">Total Reads</div>
        <div class="card-value blue">{total}</div>
        <div class="card-detail">All classified reads</div>
    </div>
    <div class="card">
        <div class="card-label">Full-Length</div>
        <div class="card-value green">{full_length}</div>
        <div class="card-detail">{_pct(full_length, total)}% of reads</div>
    </div>
    <div class="card">
        <div class="card-label">Truncated</div>
        <div class="card-value amber">{truncated}</div>
        <div class="card-detail">{_pct(truncated, total)}% of reads</div>
    </div>
    <div class="card">
        <div class="card-label">Pair Matched</div>
        <div class="card-value green">{matched}</div>
        <div class="card-detail">{_pct(matched, total)}% assigned target</div>
    </div>
    <div class="card">
        <div class="card-label">Unmatched</div>
        <div class="card-value red">{unmatched}</div>
        <div class="card-detail">{_pct(unmatched, total)}% no target</div>
    </div>
</div>

<h3>Per-Target Statistics</h3>
<div class="table-container">
<table>
<tr><th>Target</th><th>Reads</th><th>Avg Length</th><th>Start Conf</th><th>End Conf</th><th>Full-Length %</th><th>Avg ED</th><th>Avg Q_LD</th></tr>
"""

    # Sort targets: named targets first alphabetically, then "unmatched"
    sorted_targets = sorted(
        per_target_stats.keys(),
        key=lambda t: (t == "unmatched", t),
    )

    for tgt_id in sorted_targets:
        stats = per_target_stats[tgt_id]
        sc_cls = _conf_cls(stats["avg_start_conf"])
        ec_cls = _conf_cls(stats["avg_end_conf"])
        html += (
            f'<tr><td><strong>{tgt_id}</strong></td>'
            f'<td>{stats["count"]}</td>'
            f'<td>{stats["avg_length"]:.0f}bp</td>'
            f'<td><span class="conf-badge {sc_cls}">{stats["avg_start_conf"]:.3f}</span></td>'
            f'<td><span class="conf-badge {ec_cls}">{stats["avg_end_conf"]:.3f}</span></td>'
            f'<td>{stats["full_length_pct"]:.0f}%</td>'
            f'<td>{stats["avg_ed"]:.1f}</td>'
            f'<td>{stats["avg_q_ld"]:.1f} dB</td></tr>\n'
        )

    html += """</table>
</div>

<h3>Barcode Pair Matrix</h3>
<p class="note">Rows = start barcode, Columns = end barcode. Expected pairs are highlighted with a blue border.</p>
<div class="construct-diagram">
"""

    # Pair matrix grid
    used_barcodes = pair_matrix["used_barcodes"]
    pair_counts = pair_matrix["counts"]
    n_bc = len(used_barcodes)

    html += f'<div class="pair-matrix" style="grid-template-columns: 80px repeat({n_bc}, 1fr);">\n'
    html += '<div class="pair-cell header"></div>'
    for bc in used_barcodes:
        html += f'<div class="pair-cell header">end:{bc}</div>'

    for bc_start in used_barcodes:
        html += f'\n<div class="pair-cell header">start:{bc_start}</div>'
        for bc_end in used_barcodes:
            count = pair_counts.get((bc_start, bc_end), 0)
            is_exp = (bc_start, bc_end) in expected_pairs
            cls = "expected " if is_exp else ""
            if count >= 20:
                cls += "high"
            elif count >= 5:
                cls += "med"
            else:
                cls += "low"
            alias_label = pair_alias_map.get((bc_start, bc_end), "")
            tooltip = f"{alias_label} " if alias_label else ""
            html += f'<div class="pair-cell {cls}" title="{tooltip}{bc_start}--{bc_end}: {count} reads">{count}</div>'

    html += "\n</div></div>\n"

    # Confidence histograms
    html += """
<h3>Confidence Score Distributions</h3>
<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;">
<div class="construct-diagram">
    <h3>Start Barcode Confidence</h3>
    <div class="histogram" id="start-hist"></div>
    <p class="note">High confidence (&gt;0.9) = barcode clearly present. Low (~0.58) = no real match.</p>
</div>
<div class="construct-diagram">
    <h3>End Barcode Confidence</h3>
    <div class="histogram" id="end-hist"></div>
    <p class="note">Bimodal: high = full-length read, low = truncated (barcode absent).</p>
</div>
</div>
</div>
"""

    # ==================================================================
    # Tab 3: All Reads Table
    # ==================================================================

    # Build filter options for target dropdown
    target_options = ""
    for tgt_id in sorted_targets:
        target_options += f'        <option value="{tgt_id}">{tgt_id}</option>\n'

    html += f"""
<!-- Tab 3: All Reads -->
<div id="tab-reads" class="tab-panel">
<h2>All Reads ({total})</h2>

<div class="filter-bar">
    <label>Filter by target:</label>
    <select id="filter-target" onchange="filterTable()">
        <option value="">All</option>
{target_options}    </select>
    <label>Filter by status:</label>
    <select id="filter-status" onchange="filterTable()">
        <option value="">All</option>
        <option value="full">Full-length</option>
        <option value="trunc">Truncated</option>
        <option value="matched">Matched</option>
        <option value="unmatched">Unmatched</option>
    </select>
</div>

<div class="table-container">
<div class="table-scroll">
<table id="reads-table">
<thead>
<tr>
    <th onclick="sortTable(0)">#</th>
    <th onclick="sortTable(1)">Read ID</th>
    <th onclick="sortTable(2)">Length</th>
    <th onclick="sortTable(3)">Target</th>
    <th onclick="sortTable(4)">Start BC</th>
    <th onclick="sortTable(5)">Start ED</th>
    <th onclick="sortTable(6)">Start Conf</th>
    <th onclick="sortTable(7)">End BC</th>
    <th onclick="sortTable(8)">End ED</th>
    <th onclick="sortTable(9)">End Conf</th>
    <th onclick="sortTable(10)">Pair</th>
    <th onclick="sortTable(11)">Status</th>
    <th onclick="sortTable(12)">Target ED</th>
    <th onclick="sortTable(13)">Q_LD</th>
</tr>
</thead>
<tbody>
"""

    for i, r in enumerate(reads_table):
        sc_cls = _conf_cls(r.get("bc_start_conf"))
        ec_cls = _conf_cls(r.get("bc_end_conf"))
        status = "full" if r["is_full_length"] else "trunc"
        tgt_id_val = r.get("tgt_id", "")
        is_matched = not tgt_id_val.startswith("unmatched")
        matched_str = "matched" if is_matched else "unmatched"

        sc_val = r["bc_start_conf"] if r["bc_start_conf"] is not None else ""
        ec_val = r["bc_end_conf"] if r["bc_end_conf"] is not None else ""
        ed_val = r["ed"] if r["ed"] is not None else ""
        qld_val = r["q_ld"] if r["q_ld"] is not None else ""

        fl_label = "Full" if r["is_full_length"] else "Trunc"
        fl_cls = "full" if r["is_full_length"] else "trunc"

        html += (
            f'<tr data-target="{tgt_id_val}" data-status="{status}" data-matched="{matched_str}">'
            f'<td>{i + 1}</td>'
            f'<td title="{r["read_id_full"]}">{r["read_id"]}</td>'
            f'<td>{r["length"]}</td>'
            f'<td>{tgt_id_val}</td>'
            f'<td>{r["bc_start"]}</td>'
            f'<td>{r.get("bc_start_ed", "")}</td>'
            f'<td><span class="conf-badge {sc_cls}">{sc_val}</span></td>'
            f'<td>{r["bc_end"]}</td>'
            f'<td>{r.get("bc_end_ed", "")}</td>'
            f'<td><span class="conf-badge {ec_cls}">{ec_val}</span></td>'
            f'<td>{r["pair"]}</td>'
            f'<td><span class="meta-chip {fl_cls}">{fl_label}</span></td>'
            f'<td>{ed_val}</td>'
            f'<td>{qld_val}</td>'
            f'</tr>\n'
        )

    html += """</tbody></table></div></div>
</div>
"""

    # ==================================================================
    # Tab 4: Read Details
    # ==================================================================
    html += """
<!-- Tab 4: Read Details -->
<div id="tab-details" class="tab-panel">
<h2>Individual Read Analysis</h2>
<p class="note">Showing representative reads per target (full-length and truncated where available).</p>

<div class="legend">
    <div class="legend-item"><div class="legend-swatch" style="background:#5b21b6;"></div>Start barcode</div>
    <div class="legend-item"><div class="legend-swatch" style="background:#9d174d;"></div>End barcode (RC)</div>
    <div class="legend-item"><div class="legend-swatch" style="background:#1e3a5f;"></div>Flanking sequence</div>
    <div class="legend-item"><div class="legend-swatch" style="background:#374151;"></div>Unknown/unassigned</div>
</div>
"""

    for ex in detailed_reads:
        is_fl = ex["is_full_length"]
        fl_cls = "full" if is_fl else "trunc"
        fl_label = "Full-length" if is_fl else "Truncated"
        tgt_id_val = ex.get("tgt_id", "")
        is_matched = tgt_id_val and not tgt_id_val.startswith("unmatched")
        m_cls = "matched" if is_matched else "unmatched"
        pair_val = ex.get("pair", "")

        # Truncate read_id for display
        rid_display = ex["read_id"]
        if len(rid_display) > 20:
            rid_display = rid_display[:20] + "..."

        html += f"""
<div class="read-detail">
    <div class="read-header">
        <div>
            <strong style="font-size: 14px;">{rid_display}</strong>
            <span style="color: #64748b; font-size: 12px; margin-left: 10px;">{ex['length']}bp</span>
        </div>
        <div class="read-meta">
            <span class="meta-chip {fl_cls}">{fl_label}</span>
            <span class="meta-chip {'matched' if is_matched else 'unmatched'}">Target: {tgt_id_val if tgt_id_val else 'none'}</span>
            <span class="meta-chip {m_cls}">{'Pair: ' + pair_val if is_matched else 'Unmatched: ' + pair_val}</span>
        </div>
    </div>
"""

        # Read structure bar
        seq_len = ex["length"]
        regions = sorted(ex.get("regions", []), key=lambda r: r["start"])

        html += '    <div class="read-structure">\n'
        pos = 0
        for reg in regions:
            if reg["start"] > pos:
                gap_pct = (reg["start"] - pos) / seq_len * 100
                html += f'        <div class="segment seg-unknown" style="width:{gap_pct:.1f}%"></div>\n'
            w_pct = max(0.5, (reg["end"] - reg["start"]) / seq_len * 100)
            reg_type = reg.get("type", "")
            if reg_type in ("start_barcode", "barcode_start"):
                seg_cls = "seg-start-bc"
            elif reg_type in ("end_barcode", "barcode_end"):
                seg_cls = "seg-end-bc"
            elif reg_type == "flank":
                seg_cls = "seg-flank"
            else:
                seg_cls = "seg-unknown"
            label = reg.get("label", "")[:15]
            html += f'        <div class="segment {seg_cls}" style="width:{w_pct:.1f}%" title="{reg.get("label", "")} pos:{reg["start"]}-{reg["end"]} ed:{reg.get("ed", "")}">{label}</div>\n'
            pos = reg["end"]
        if pos < seq_len:
            gap_pct = (seq_len - pos) / seq_len * 100
            html += f'        <div class="segment seg-unknown" style="width:{gap_pct:.1f}%"></div>\n'
        html += '    </div>\n'

        # Annotated sequence viewer
        seq = ex.get("seq", "")
        if seq:
            html += '    <div class="seq-viewer">\n'

            # Ruler
            ruler = ""
            for p in range(0, seq_len, 10):
                ruler += f"{p:<10d}"
            html += f'        <div class="seq-ruler">{ruler[:seq_len]}</div>\n'

            # Build annotation map
            ann_map: list[str | None] = [None] * seq_len
            for reg in regions:
                reg_type = reg.get("type", "")
                if reg_type in ("start_barcode", "barcode_start"):
                    cls = "bc-start"
                elif reg_type in ("end_barcode", "barcode_end"):
                    cls = "bc-end"
                elif reg_type == "flank":
                    cls = "flank"
                else:
                    cls = ""
                if cls:
                    for p in range(reg["start"], min(reg["end"], seq_len)):
                        ann_map[p] = cls

            # Build HTML spans
            seq_html = ""
            current_cls: str | None = None
            chunk = ""
            for p in range(seq_len):
                c = ann_map[p]
                if c != current_cls:
                    if chunk:
                        if current_cls:
                            seq_html += f'<span class="{current_cls}">{chunk}</span>'
                        else:
                            seq_html += chunk
                    chunk = seq[p]
                    current_cls = c
                else:
                    chunk += seq[p]
            if chunk:
                if current_cls:
                    seq_html += f'<span class="{current_cls}">{chunk}</span>'
                else:
                    seq_html += chunk

            html += f'        <div class="seq-text">{seq_html}</div>\n'
            html += '    </div>\n'

        # Competitor tables
        start_competitors = ex.get("start_competitors", [])
        end_competitors = ex.get("end_competitors", [])

        html += '    <div class="competitor-grid">\n'
        html += '        <div><h3>Start Barcode Alignment (first 100bp)</h3>\n'
        html += '        <table class="competitor-table"><tr><th>Barcode</th><th>Edit Dist</th><th>Confidence</th><th>Position</th></tr>\n'
        for idx, comp in enumerate(start_competitors):
            bold = ' style="font-weight:700;"' if idx == 0 else ""
            bc_id = comp.get("bc_id", comp.get("bc", ""))
            ed = comp.get("edit_distance", comp.get("ed", ""))
            conf = comp.get("confidence", comp.get("conf", ""))
            if isinstance(conf, float):
                conf = f"{conf:.3f}"
            start_pos = comp.get("start", "")
            end_pos = comp.get("end", "")
            pos_str = f"{start_pos}-{end_pos}"
            html += f'        <tr{bold}><td>{bc_id}</td><td>{ed}</td><td>{conf}</td><td>{pos_str}</td></tr>\n'
        html += '        </table></div>\n'

        html += '        <div><h3>End Barcode Alignment (last 100bp, RC)</h3>\n'
        html += '        <table class="competitor-table"><tr><th>Barcode</th><th>Edit Dist</th><th>Confidence</th><th>Position</th></tr>\n'
        for idx, comp in enumerate(end_competitors):
            bold = ' style="font-weight:700;"' if idx == 0 else ""
            bc_id = comp.get("bc_id", comp.get("bc", ""))
            ed = comp.get("edit_distance", comp.get("ed", ""))
            conf = comp.get("confidence", comp.get("conf", ""))
            if isinstance(conf, float):
                conf = f"{conf:.3f}"
            start_pos = comp.get("start", "")
            end_pos = comp.get("end", "")
            pos_str = f"{start_pos}-{end_pos}"
            html += f'        <tr{bold}><td>{bc_id}</td><td>{ed}</td><td>{conf}</td><td>{pos_str}</td></tr>\n'
        html += '        </table></div>\n'
        html += '    </div>\n'

        # Alignment metrics
        db_ed = ex.get("db_ed")
        db_q_ld = ex.get("db_q_ld")
        if db_ed is not None:
            q_ld_display = f"{db_q_ld:.1f}" if db_q_ld is not None else "N/A"
            html += f'    <div style="margin-top:10px; font-size:13px; color:#94a3b8;">Target alignment: ED={db_ed}, Q_LD={q_ld_display} dB</div>\n'

        html += '</div>\n'

    html += "</div>\n"

    # ==================================================================
    # Tab 5: Barcode Sequences
    # ==================================================================
    html += """
<!-- Tab 5: Barcode Sequences -->
<div id="tab-barcodes" class="tab-panel">
<h2>Barcode Sequences Used</h2>
<div class="table-container">
<table>
<tr><th>Barcode ID</th><th>Forward Sequence (5'&rarr;3')</th><th>Reverse Complement</th><th>Used As</th></tr>
"""

    for bc_id in sorted(barcode_info.keys()):
        info = barcode_info[bc_id]
        usage_parts: list[str] = []
        for entry in pairing_table:
            if entry["upstream"] == bc_id:
                usage_parts.append(f"upstream in {entry['alias']}")
            if entry["downstream"] == bc_id:
                usage_parts.append(f"downstream in {entry['alias']}")
        usage_str = ", ".join(usage_parts) if usage_parts else "unused"
        html += (
            f'<tr><td><strong>{bc_id}</strong></td>'
            f'<td style="font-family:monospace;font-size:12px;">{info["fwd"]}</td>'
            f'<td style="font-family:monospace;font-size:12px;">{info["rc"]}</td>'
            f'<td>{usage_str}</td></tr>\n'
        )

    html += """</table></div>

<h3>Classification Algorithm</h3>
<div class="construct-diagram">
<ol style="font-size: 14px; line-height: 2; padding-left: 20px;">
    <li><strong>Extract segments:</strong> first 100bp (start) and last 100bp (end) of each read</li>
    <li><strong>Align barcodes:</strong> Semi-global alignment (edlib HW mode) of each expected barcode against the segment</li>
    <li><strong>Select best match:</strong> Lowest edit distance wins. Confidence = 1.0 - (ED / 24)</li>
    <li><strong>End barcode uses RC:</strong> The downstream barcode appears as its reverse complement in the read</li>
    <li><strong>Look up pair:</strong> (start_bc, end_bc) matched against sample sheet for target assignment</li>
    <li><strong>Align to target:</strong> Matched reads aligned (NW mode) to assigned reference for quality metrics</li>
</ol>
</div>
</div>
"""

    # ==================================================================
    # JavaScript
    # ==================================================================
    start_confs_json = json.dumps(conf_dists["start_confs"])
    end_confs_json = json.dumps(conf_dists["end_confs"])

    html += f"""
<script>
const startConfs = {start_confs_json};
const endConfs = {end_confs_json};

function showTab(name) {{
    document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
    document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
    document.getElementById('tab-' + name).classList.add('active');
    event.target.classList.add('active');
}}

function filterTable() {{
    const targetVal = document.getElementById('filter-target').value;
    const statusVal = document.getElementById('filter-status').value;
    document.querySelectorAll('#reads-table tbody tr').forEach(tr => {{
        let show = true;
        if (targetVal && tr.dataset.target !== targetVal) show = false;
        if (statusVal === 'full' && tr.dataset.status !== 'full') show = false;
        if (statusVal === 'trunc' && tr.dataset.status !== 'trunc') show = false;
        if (statusVal === 'matched' && tr.dataset.matched !== 'matched') show = false;
        if (statusVal === 'unmatched' && tr.dataset.matched !== 'unmatched') show = false;
        tr.style.display = show ? '' : 'none';
    }});
}}

let sortDir = {{}};
function sortTable(col) {{
    const table = document.getElementById('reads-table');
    const tbody = table.querySelector('tbody');
    const rows = Array.from(tbody.querySelectorAll('tr'));
    sortDir[col] = !(sortDir[col] || false);

    rows.sort((a, b) => {{
        let va = a.cells[col].textContent.trim();
        let vb = b.cells[col].textContent.trim();
        let na = parseFloat(va), nb = parseFloat(vb);
        if (!isNaN(na) && !isNaN(nb)) {{
            return sortDir[col] ? na - nb : nb - na;
        }}
        return sortDir[col] ? va.localeCompare(vb) : vb.localeCompare(va);
    }});
    rows.forEach(r => tbody.appendChild(r));
}}

function drawHist(containerId, values, color) {{
    const container = document.getElementById(containerId);
    const bins = 20;
    const counts = new Array(bins).fill(0);
    values.forEach(v => {{
        const idx = Math.min(Math.floor(v * bins), bins - 1);
        counts[idx]++;
    }});
    const max = Math.max(...counts, 1);

    let html = '';
    for (let i = 0; i < bins; i++) {{
        const pct = counts[i] / max * 100;
        const lo = (i / bins).toFixed(2);
        const hi = ((i + 1) / bins).toFixed(2);
        html += '<div class="hist-bar" style="height:' + pct + '%;background:' + color + ';">' +
                '<div class="tooltip">' + lo + '-' + hi + ': ' + counts[i] + ' reads</div></div>';
    }}
    container.innerHTML = html;

    let labels = '<div class="hist-labels">';
    for (let i = 0; i < bins; i++) {{
        labels += '<span>' + (i % 5 === 0 ? (i / bins).toFixed(1) : '') + '</span>';
    }}
    labels += '</div>';
    container.insertAdjacentHTML('afterend', labels);
}}

drawHist('start-hist', startConfs, '#7c3aed');
drawHist('end-hist', endConfs, '#be185d');
</script>

</div>
</body>
</html>
"""

    return html


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _reverse_complement_simple(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    This is a local helper so the module has no dependency on ``barcodes.py``
    for the construct diagram display.
    """
    return seq.translate(_COMPLEMENT)[::-1]
