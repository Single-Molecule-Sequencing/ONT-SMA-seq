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

        Optional extended keys (for Statistics tab):

        * ``length_stats`` -- {min, max, mean, median, n50, lengths,
          lengths_full, lengths_trunc}
        * ``quality_stats`` -- {q_bc_values, q_ld_values, ed_values,
          mean_q_bc, mean_q_ld, mean_ed}
        * ``confidence_scatter`` -- [{x, y, tgt, fl, len}, ...]
        * ``end_reason_counts`` -- {reason: count, ...}
        * ``classification_by_confidence`` -- [{threshold, matched, pct}, ...]

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

    # Extended stats (optional, with safe defaults)
    length_stats = analysis.get("length_stats") or {}
    quality_stats = analysis.get("quality_stats") or {}
    confidence_scatter = analysis.get("confidence_scatter") or []
    end_reason_counts = analysis.get("end_reason_counts") or {}
    classification_by_confidence = analysis.get("classification_by_confidence") or []
    truncation_counts = analysis.get("truncation_counts") or {}

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

    # Summary values
    total = summary["total"]
    full_length = summary["full_length"]
    truncated = summary["truncated"]
    matched = summary["matched"]
    unmatched = summary["unmatched"]

    # Stats bar values
    len_min = length_stats.get("min", "?")
    len_max = length_stats.get("max", "?")
    len_mean = length_stats.get("mean", "?")
    len_n50 = length_stats.get("n50", "?")
    mean_q_bc = quality_stats.get("mean_q_bc", "?")
    mean_q_ld = quality_stats.get("mean_q_ld", "?")

    # Format stats bar values
    def _fmt(v: Any, suffix: str = "", fmt: str = ".0f") -> str:
        if v == "?" or v is None:
            return "?"
        try:
            return f"{float(v):{fmt}}{suffix}"
        except (ValueError, TypeError):
            return str(v)

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
.subtitle {{ color: #94a3b8; font-size: 14px; margin-bottom: 15px; }}

/* Cards */
.card-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 15px;
    margin-bottom: 15px;
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

/* Stats bar */
.stats-bar {{
    display: flex;
    gap: 20px;
    flex-wrap: wrap;
    padding: 10px 20px;
    background: #1a1a2e;
    border-radius: 8px;
    margin-bottom: 20px;
    font-size: 13px;
    color: #94a3b8;
}}
.stats-bar .stat-val {{ color: #e0e0e0; font-weight: 600; }}

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
tr.row-selected td {{ background: #1e3a5f !important; }}
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
    margin-bottom: 20px;
    overflow: hidden;
}}
.read-detail-header {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 15px 20px;
    flex-wrap: wrap;
    gap: 10px;
    cursor: pointer;
    user-select: none;
}}
.read-detail-header:hover {{
    background: #252540;
}}
.read-detail-body {{
    padding: 0 20px 20px;
    overflow: hidden;
    transition: max-height 0.3s ease, padding 0.3s ease;
}}
.read-detail.collapsed .read-detail-body {{
    max-height: 0 !important;
    padding-top: 0;
    padding-bottom: 0;
}}
.collapse-icon {{
    transition: transform 0.2s;
    display: inline-block;
    margin-right: 8px;
    color: #94a3b8;
}}
.read-detail.collapsed .collapse-icon {{
    transform: rotate(-90deg);
}}
.read-header {{
    display: flex;
    justify-content: space-between;
    align-items: center;
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
.pair-cell.clickable {{ cursor: pointer; transition: opacity 0.2s; }}
.pair-cell.clickable:hover {{ opacity: 0.75; box-shadow: 0 0 0 2px #60a5fa; }}
.pair-cell.high {{ background: #065f46; color: #34d399; }}
.pair-cell.med {{ background: #78350f; color: #fbbf24; }}
.pair-cell.low {{ background: #1a1a2e; color: #475569; }}
.pair-cell.expected {{ border: 2px solid #60a5fa; }}

/* Histograms (legacy div-based) */
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
.row-count {{ font-size: 13px; color: #94a3b8; margin-left: auto; }}

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

/* Canvas chart styles */
canvas {{ display: block; margin: 10px auto; }}
.chart-container {{ position: relative; }}
.chart-tooltip {{
    position: absolute;
    display: none;
    background: #1a1a2e;
    border: 1px solid #334155;
    padding: 6px 10px;
    border-radius: 4px;
    font-size: 11px;
    pointer-events: none;
    z-index: 30;
    color: #e0e0e0;
    white-space: nowrap;
}}
.chart-section {{
    background: #1e1e2e;
    border-radius: 12px;
    padding: 20px;
    border: 1px solid #2d2d44;
    margin-bottom: 20px;
}}
.chart-stats {{
    display: flex;
    gap: 20px;
    flex-wrap: wrap;
    margin-top: 10px;
    font-size: 13px;
    color: #94a3b8;
}}
.chart-stats .cs-val {{ color: #e0e0e0; font-weight: 600; }}
</style>
</head>
<body>
<div class="container">

<h1>SMA-seq Barcode Classification Report</h1>
<p class="subtitle">{subtitle}</p>

<!-- Always-visible summary cards (above tabs) -->
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

<!-- Compact stats bar -->
<div class="stats-bar">
    <span>Read lengths: <span class="stat-val">{_fmt(len_min)}</span>&ndash;<span class="stat-val">{_fmt(len_max)}</span> bp, mean <span class="stat-val">{_fmt(len_mean)}</span> bp, N50 <span class="stat-val">{_fmt(len_n50)}</span> bp</span>
    <span>|</span>
    <span>Quality: mean Q_BC <span class="stat-val">{_fmt(mean_q_bc, fmt=".1f")}</span>, mean Q_LD <span class="stat-val">{_fmt(mean_q_ld, " dB", ".1f")}</span></span>
</div>

<!-- Tabs -->
<div class="tab-bar">
    <button class="tab-btn active" onclick="showTab('construct', this)">Construct Structure</button>
    <button class="tab-btn" onclick="showTab('summary', this)">Classification Summary</button>
    <button class="tab-btn" onclick="showTab('statistics', this)">Statistics</button>
    <button class="tab-btn" onclick="showTab('reads', this)">All Reads Table</button>
    <button class="tab-btn" onclick="showTab('details', this)">Read Details</button>
    <button class="tab-btn" onclick="showTab('barcodes', this)">Barcode Sequences</button>
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
    html += f"""
<!-- Tab 2: Classification Summary -->
<div id="tab-summary" class="tab-panel">
<h2>Classification Summary</h2>

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
"""

    # Truncation breakdown table (if available)
    if truncation_counts:
        html += """
<h3>Truncation Level Breakdown</h3>
<div class="table-container">
<table>
<tr><th>Level</th><th>Count</th><th>Percentage</th></tr>
"""
        # Order levels: full_length, bc1_target_bc2, bc1_target, bc1_only, adapter_only
        level_order = ["full_length", "bc1_target_bc2", "bc1_target", "bc1_only", "adapter_only"]
        for level in level_order:
            if level in truncation_counts:
                count = truncation_counts[level]
                pct = (count / total * 100.0) if total else 0.0
                html += f'<tr><td><strong>{level}</strong></td><td>{count}</td><td>{pct:.1f}%</td></tr>\n'

        # Add any other levels not in the ordered list
        for level in sorted(truncation_counts.keys()):
            if level not in level_order:
                count = truncation_counts[level]
                pct = (count / total * 100.0) if total else 0.0
                html += f'<tr><td><strong>{level}</strong></td><td>{count}</td><td>{pct:.1f}%</td></tr>\n'

        html += """</table>
</div>
"""

    html += """
<h3>Barcode Pair Matrix</h3>
<p class="note">Rows = start barcode, Columns = end barcode. Expected pairs are highlighted with a blue border. Click a cell to filter the reads table to that pair.</p>
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
            # Make cells with reads clickable
            click_attr = ""
            if count > 0:
                pair_str = f"{bc_start}--{bc_end}"
                click_attr = f' onclick="filterByPair(\'{pair_str}\')"'
                cls += " clickable"
            html += f'<div class="pair-cell {cls}" title="{tooltip}{bc_start}--{bc_end}: {count} reads"{click_attr}>{count}</div>'

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
    # Tab 3: Statistics (NEW)
    # ==================================================================
    # Prepare JSON data for JS charts
    lengths_json = json.dumps(length_stats.get("lengths", []))
    lengths_full_json = json.dumps(length_stats.get("lengths_full", []))
    lengths_trunc_json = json.dumps(length_stats.get("lengths_trunc", []))
    q_bc_json = json.dumps(quality_stats.get("q_bc_values", []))
    q_ld_json = json.dumps(quality_stats.get("q_ld_values", []))
    ed_values_json = json.dumps(quality_stats.get("ed_values", []))
    scatter_json = json.dumps(confidence_scatter)
    end_reason_json = json.dumps(end_reason_counts)
    conf_curve_json = json.dumps(classification_by_confidence)

    # Per-target lengths for box plots
    per_tgt_lengths: dict[str, list] = {}
    for tgt_id in sorted_targets:
        tgt_stats = per_target_stats.get(tgt_id, {})
        if isinstance(tgt_stats, dict):
            per_tgt_lengths[tgt_id] = tgt_stats.get("lengths", [])
    per_tgt_lengths_json = json.dumps(per_tgt_lengths)

    # Stats values for display
    len_median = length_stats.get("median", "?")
    mean_ed = quality_stats.get("mean_ed", "?")

    html += f"""
<!-- Tab 3: Statistics -->
<div id="tab-statistics" class="tab-panel">
<h2>Statistics</h2>

<!-- (a) Read Length Distribution -->
<div class="chart-section">
    <h3>Read Length Distribution</h3>
    <div class="chart-container">
        <canvas id="canvas-length-hist"></canvas>
        <div class="chart-tooltip" id="tip-length-hist"></div>
    </div>
    <div class="legend" style="margin-top:8px;">
        <div class="legend-item"><div class="legend-swatch" style="background:#4ade80;"></div>Full-length</div>
        <div class="legend-item"><div class="legend-swatch" style="background:#fbbf24;"></div>Truncated</div>
        <div class="legend-item"><div class="legend-swatch" style="background:#60a5fa; width:2px; height:14px;"></div>Mean</div>
        <div class="legend-item"><div class="legend-swatch" style="background:#a78bfa; width:2px; height:14px;"></div>Median</div>
        <div class="legend-item"><div class="legend-swatch" style="background:#f87171; width:2px; height:14px;"></div>N50</div>
    </div>
    <div class="chart-stats">
        <span>Min: <span class="cs-val">{_fmt(len_min)} bp</span></span>
        <span>Max: <span class="cs-val">{_fmt(len_max)} bp</span></span>
        <span>Mean: <span class="cs-val">{_fmt(len_mean)} bp</span></span>
        <span>Median: <span class="cs-val">{_fmt(len_median)} bp</span></span>
        <span>N50: <span class="cs-val">{_fmt(len_n50)} bp</span></span>
    </div>
</div>

<!-- (b) Per-Target Length Comparison -->
<div class="chart-section">
    <h3>Per-Target Length Comparison</h3>
    <div class="chart-container">
        <canvas id="canvas-target-box"></canvas>
        <div class="chart-tooltip" id="tip-target-box"></div>
    </div>
    <p class="note">Box plot showing median (line), Q1-Q3 (box), and min-max (whiskers) for each target.</p>
</div>

<!-- (c) Quality Score Distributions -->
<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;">
<div class="chart-section">
    <h3>Q_BC Distribution</h3>
    <div class="chart-container">
        <canvas id="canvas-qbc-hist"></canvas>
        <div class="chart-tooltip" id="tip-qbc-hist"></div>
    </div>
    <div class="chart-stats">
        <span>Mean Q_BC: <span class="cs-val">{_fmt(mean_q_bc, fmt=".2f")}</span></span>
    </div>
</div>
<div class="chart-section">
    <h3>Q_LD Distribution (dB)</h3>
    <div class="chart-container">
        <canvas id="canvas-qld-hist"></canvas>
        <div class="chart-tooltip" id="tip-qld-hist"></div>
    </div>
    <div class="chart-stats">
        <span>Mean Q_LD: <span class="cs-val">{_fmt(mean_q_ld, " dB", ".1f")}</span></span>
    </div>
</div>
</div>

<!-- (d) Edit Distance Distribution -->
<div class="chart-section">
    <h3>Edit Distance Distribution</h3>
    <div class="chart-container">
        <canvas id="canvas-ed-hist"></canvas>
        <div class="chart-tooltip" id="tip-ed-hist"></div>
    </div>
    <div class="chart-stats">
        <span>Mean ED: <span class="cs-val">{_fmt(mean_ed, fmt=".1f")}</span></span>
    </div>
</div>

<!-- (e) Start vs End Confidence Scatter -->
<div class="chart-section">
    <h3>Start vs End Confidence Scatter Plot</h3>
    <div class="chart-container">
        <canvas id="canvas-conf-scatter"></canvas>
        <div class="chart-tooltip" id="tip-conf-scatter"></div>
    </div>
    <div class="legend" style="margin-top:8px;">
        <div class="legend-item"><span style="display:inline-block;width:10px;height:10px;border-radius:50%;background:#4ade80;margin-right:4px;"></span>Full-length (filled)</div>
        <div class="legend-item"><span style="display:inline-block;width:10px;height:10px;border-radius:50%;border:2px solid #fbbf24;margin-right:4px;"></span>Truncated (open)</div>
    </div>
    <p class="note">Points colored by target. Hover for read info.</p>
</div>

<!-- (f) Confidence Threshold Curve -->
<div class="chart-section">
    <h3>Confidence Threshold Curve</h3>
    <div class="chart-container">
        <canvas id="canvas-conf-curve"></canvas>
        <div class="chart-tooltip" id="tip-conf-curve"></div>
    </div>
    <div id="conf-curve-info" class="chart-stats">
        <span>Click on the curve to see matched reads at a given threshold.</span>
    </div>
</div>

<!-- (g) End Reason Breakdown -->
<div class="chart-section">
    <h3>End Reason Breakdown</h3>
    <div class="chart-container">
        <canvas id="canvas-end-reason"></canvas>
        <div class="chart-tooltip" id="tip-end-reason"></div>
    </div>
</div>

</div>
"""

    # ==================================================================
    # Tab 4: All Reads Table
    # ==================================================================

    # Build filter options for target dropdown
    target_options = ""
    for tgt_id in sorted_targets:
        target_options += f'        <option value="{tgt_id}">{tgt_id}</option>\n'

    html += f"""
<!-- Tab 4: All Reads -->
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
    <label>Filter by pair:</label>
    <select id="filter-pair" onchange="filterTable()">
        <option value="">All</option>
    </select>
    <label>Search Read ID:</label>
    <input type="text" id="search-read-id" placeholder="Type to search..." oninput="filterTable()" style="width:180px;">
    <span class="row-count" id="row-count">Showing {total} of {total} reads</span>
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

    # Collect unique pairs for the pair filter dropdown
    unique_pairs: set[str] = set()

    for i, r in enumerate(reads_table):
        sc_cls = _conf_cls(r.get("bc_start_conf"))
        ec_cls = _conf_cls(r.get("bc_end_conf"))
        status = "full" if r["is_full_length"] else "trunc"
        tgt_id_val = r.get("tgt_id", "")
        is_matched = not tgt_id_val.startswith("unmatched")
        matched_str = "matched" if is_matched else "unmatched"
        pair_val = r.get("pair", "")
        if pair_val:
            unique_pairs.add(pair_val)

        sc_val = r["bc_start_conf"] if r["bc_start_conf"] is not None else ""
        ec_val = r["bc_end_conf"] if r["bc_end_conf"] is not None else ""
        ed_val = r["ed"] if r["ed"] is not None else ""
        qld_val = r["q_ld"] if r["q_ld"] is not None else ""

        fl_label = "Full" if r["is_full_length"] else "Trunc"
        fl_cls = "full" if r["is_full_length"] else "trunc"

        html += (
            f'<tr data-target="{tgt_id_val}" data-status="{status}" data-matched="{matched_str}" data-pair="{pair_val}" data-readid="{r["read_id_full"]}" onclick="selectRow(this)">'
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
            f'<td>{pair_val}</td>'
            f'<td><span class="meta-chip {fl_cls}">{fl_label}</span></td>'
            f'<td>{ed_val}</td>'
            f'<td>{qld_val}</td>'
            f'</tr>\n'
        )

    html += """</tbody></table></div></div>
</div>
"""

    # ==================================================================
    # Tab 5: Read Details (collapsible panels)
    # ==================================================================
    html += """
<!-- Tab 5: Read Details -->
<div id="tab-details" class="tab-panel">
<h2>Individual Read Analysis</h2>
<p class="note">Showing representative reads per target (full-length and truncated where available). Click a header to expand/collapse.</p>

<div class="legend">
    <div class="legend-item"><div class="legend-swatch" style="background:#5b21b6;"></div>Start barcode</div>
    <div class="legend-item"><div class="legend-swatch" style="background:#9d174d;"></div>End barcode (RC)</div>
    <div class="legend-item"><div class="legend-swatch" style="background:#1e3a5f;"></div>Flanking sequence</div>
    <div class="legend-item"><div class="legend-swatch" style="background:#374151;"></div>Unknown/unassigned</div>
</div>
"""

    for idx_detail, ex in enumerate(detailed_reads):
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

        # First panel expanded, rest collapsed
        collapsed_cls = "" if idx_detail == 0 else " collapsed"

        html += f"""
<div class="read-detail{collapsed_cls}" id="read-detail-{idx_detail}">
    <div class="read-detail-header" onclick="toggleReadDetail('read-detail-{idx_detail}')">
        <div>
            <span class="collapse-icon">&#9660;</span>
            <strong style="font-size: 14px;">{rid_display}</strong>
            <span style="color: #64748b; font-size: 12px; margin-left: 10px;">{ex['length']}bp</span>
        </div>
        <div class="read-meta">
            <span class="meta-chip {fl_cls}">{fl_label}</span>
            <span class="meta-chip {'matched' if is_matched else 'unmatched'}">Target: {tgt_id_val if tgt_id_val else 'none'}</span>
            <span class="meta-chip {m_cls}">{'Pair: ' + pair_val if is_matched else 'Unmatched: ' + pair_val}</span>
        </div>
    </div>
    <div class="read-detail-body">
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

        html += '    </div>\n'  # close read-detail-body
        html += '</div>\n'  # close read-detail

    html += "</div>\n"

    # ==================================================================
    # Tab 6: Barcode Sequences
    # ==================================================================
    html += """
<!-- Tab 6: Barcode Sequences -->
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
    # Populate pair filter dropdown via JS
    # ==================================================================
    unique_pairs_json = json.dumps(sorted(unique_pairs))

    # ==================================================================
    # JavaScript
    # ==================================================================
    start_confs_json = json.dumps(conf_dists["start_confs"])
    end_confs_json = json.dumps(conf_dists["end_confs"])

    html += f"""
<script>
// ===================== DATA =====================
const startConfs = {start_confs_json};
const endConfs = {end_confs_json};
const lengthsAll = {lengths_json};
const lengthsFull = {lengths_full_json};
const lengthsTrunc = {lengths_trunc_json};
const qBcValues = {q_bc_json};
const qLdValues = {q_ld_json};
const edValues = {ed_values_json};
const scatterData = {scatter_json};
const endReasonCounts = {end_reason_json};
const confCurveData = {conf_curve_json};
const perTargetLengths = {per_tgt_lengths_json};
const uniquePairs = {unique_pairs_json};
const lengthMean = {json.dumps(length_stats.get("mean"))};
const lengthMedian = {json.dumps(length_stats.get("median"))};
const lengthN50 = {json.dumps(length_stats.get("n50"))};
const meanQBc = {json.dumps(quality_stats.get("mean_q_bc"))};
const meanQLd = {json.dumps(quality_stats.get("mean_q_ld"))};
const meanEd = {json.dumps(quality_stats.get("mean_ed"))};

// ===================== TAB NAVIGATION =====================
function showTab(name, btn) {{
    document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
    document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
    document.getElementById('tab-' + name).classList.add('active');
    if (btn) btn.classList.add('active');
    // Render charts when Statistics tab becomes visible
    if (name === 'statistics') {{
        requestAnimationFrame(renderAllCharts);
    }}
}}

// ===================== TABLE FILTERING =====================
function filterTable() {{
    const targetVal = document.getElementById('filter-target').value;
    const statusVal = document.getElementById('filter-status').value;
    const pairVal = document.getElementById('filter-pair').value;
    const searchVal = document.getElementById('search-read-id').value.toLowerCase();
    let shown = 0;
    let total = 0;
    document.querySelectorAll('#reads-table tbody tr').forEach(tr => {{
        total++;
        let show = true;
        if (targetVal && tr.dataset.target !== targetVal) show = false;
        if (statusVal === 'full' && tr.dataset.status !== 'full') show = false;
        if (statusVal === 'trunc' && tr.dataset.status !== 'trunc') show = false;
        if (statusVal === 'matched' && tr.dataset.matched !== 'matched') show = false;
        if (statusVal === 'unmatched' && tr.dataset.matched !== 'unmatched') show = false;
        if (pairVal && tr.dataset.pair !== pairVal) show = false;
        if (searchVal && !tr.dataset.readid.toLowerCase().includes(searchVal)) show = false;
        tr.style.display = show ? '' : 'none';
        if (show) shown++;
    }});
    document.getElementById('row-count').textContent = 'Showing ' + shown + ' of ' + total + ' reads';
}}

function selectRow(tr) {{
    document.querySelectorAll('#reads-table tbody tr').forEach(r => r.classList.remove('row-selected'));
    tr.classList.add('row-selected');
}}

function filterByPair(pairStr) {{
    // Set pair filter dropdown and switch to reads tab
    const pairSelect = document.getElementById('filter-pair');
    // Find or create the option
    let found = false;
    for (const opt of pairSelect.options) {{
        if (opt.value === pairStr) {{ found = true; break; }}
    }}
    if (!found) {{
        const opt = document.createElement('option');
        opt.value = pairStr;
        opt.textContent = pairStr;
        pairSelect.appendChild(opt);
    }}
    pairSelect.value = pairStr;
    filterTable();
    // Switch to reads tab
    document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
    document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
    document.getElementById('tab-reads').classList.add('active');
    document.querySelectorAll('.tab-btn').forEach(b => {{
        if (b.textContent.includes('All Reads')) b.classList.add('active');
    }});
}}

// ===================== TABLE SORTING =====================
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

// ===================== COLLAPSIBLE READ DETAILS =====================
function toggleReadDetail(id) {{
    const el = document.getElementById(id);
    el.classList.toggle('collapsed');
}}

// ===================== LEGACY DIV-BASED HISTOGRAMS =====================
function drawHist(containerId, values, color) {{
    const container = document.getElementById(containerId);
    if (!container || values.length === 0) return;
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

// ===================== CANVAS CHART UTILITIES =====================
const CHART_BG = '#0f1117';
const CHART_SURFACE = '#1e1e2e';
const CHART_GRID = '#2d2d44';
const CHART_TEXT = '#94a3b8';
const CHART_TEXT_LIGHT = '#e0e0e0';
const TARGET_PALETTE = ['#4ade80','#60a5fa','#f87171','#fbbf24','#a78bfa','#fb923c','#2dd4bf','#f472b6','#818cf8','#facc15'];

function getTargetColor(tgt, targets) {{
    if (!targets) targets = [];
    const idx = targets.indexOf(tgt);
    return TARGET_PALETTE[idx >= 0 ? idx % TARGET_PALETTE.length : 0];
}}

function setupCanvas(canvasId) {{
    const canvas = document.getElementById(canvasId);
    if (!canvas) return null;
    const container = canvas.parentElement;
    const rect = container.getBoundingClientRect();
    const w = Math.max(rect.width - 40, 300);
    const dpr = window.devicePixelRatio || 1;
    canvas.width = w * dpr;
    canvas.height = 220 * dpr;
    canvas.style.width = w + 'px';
    canvas.style.height = '220px';
    const ctx = canvas.getContext('2d');
    ctx.scale(dpr, dpr);
    ctx.clearRect(0, 0, w, 220);
    return {{ canvas, ctx, w: w, h: 220 }};
}}

function showTooltip(tipId, canvas, x, y, text) {{
    const tip = document.getElementById(tipId);
    if (!tip) return;
    const cr = canvas.getBoundingClientRect();
    tip.style.display = 'block';
    tip.innerHTML = text;
    // Position relative to chart-container
    let tx = x + 10;
    let ty = y - 10;
    if (tx + 150 > cr.width) tx = x - 140;
    if (ty < 0) ty = 10;
    tip.style.left = tx + 'px';
    tip.style.top = ty + 'px';
}}

function hideTooltip(tipId) {{
    const tip = document.getElementById(tipId);
    if (tip) tip.style.display = 'none';
}}

// Draw axes with labels
function drawAxes(ctx, margin, w, h, xLabel, yLabel) {{
    ctx.strokeStyle = CHART_GRID;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(margin.left, margin.top);
    ctx.lineTo(margin.left, h - margin.bottom);
    ctx.lineTo(w - margin.right, h - margin.bottom);
    ctx.stroke();

    ctx.fillStyle = CHART_TEXT;
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'center';
    if (xLabel) ctx.fillText(xLabel, margin.left + (w - margin.left - margin.right) / 2, h - 2);
    ctx.save();
    ctx.translate(10, margin.top + (h - margin.top - margin.bottom) / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center';
    if (yLabel) ctx.fillText(yLabel, 0, 0);
    ctx.restore();
}}

// Draw tick labels on x-axis
function drawXTicks(ctx, margin, w, h, ticks) {{
    ctx.fillStyle = CHART_TEXT;
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'center';
    ticks.forEach(t => {{
        ctx.fillText(t.label, t.x, h - margin.bottom + 14);
        ctx.strokeStyle = CHART_GRID;
        ctx.beginPath();
        ctx.moveTo(t.x, h - margin.bottom);
        ctx.lineTo(t.x, h - margin.bottom + 4);
        ctx.stroke();
    }});
}}

// Draw tick labels on y-axis
function drawYTicks(ctx, margin, h, ticks) {{
    ctx.fillStyle = CHART_TEXT;
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'right';
    ticks.forEach(t => {{
        ctx.fillText(t.label, margin.left - 5, t.y + 3);
        ctx.strokeStyle = CHART_GRID;
        ctx.lineWidth = 0.5;
        ctx.beginPath();
        ctx.moveTo(margin.left, t.y);
        ctx.lineTo(margin.left - 3, t.y);
        ctx.stroke();
    }});
}}

// Generate nice tick values
function niceTicksForRange(minVal, maxVal, desiredCount) {{
    if (maxVal <= minVal) return [{{ val: minVal, label: String(minVal) }}];
    const range = maxVal - minVal;
    const rough = range / desiredCount;
    const mag = Math.pow(10, Math.floor(Math.log10(rough)));
    let step = mag;
    if (rough / mag >= 5) step = mag * 5;
    else if (rough / mag >= 2) step = mag * 2;
    const ticks = [];
    let v = Math.ceil(minVal / step) * step;
    while (v <= maxVal) {{
        ticks.push({{ val: v, label: step >= 1 ? String(Math.round(v)) : v.toFixed(2) }});
        v += step;
    }}
    return ticks;
}}

// ===================== CANVAS HISTOGRAM =====================
function drawCanvasHistogram(canvasId, tipId, values, color, meanVal, meanColor, extraVerticals) {{
    const setup = setupCanvas(canvasId);
    if (!setup || values.length === 0) return;
    const {{ canvas, ctx, w, h }} = setup;
    const margin = {{ top: 15, right: 20, bottom: 30, left: 55 }};
    const plotW = w - margin.left - margin.right;
    const plotH = h - margin.top - margin.bottom;

    // Compute bins
    const binCount = Math.min(40, Math.max(10, Math.ceil(Math.sqrt(values.length))));
    const minV = Math.min(...values);
    const maxV = Math.max(...values);
    const range = maxV - minV || 1;
    const binWidth = range / binCount;
    const bins = new Array(binCount).fill(0);
    values.forEach(v => {{
        let idx = Math.floor((v - minV) / binWidth);
        if (idx >= binCount) idx = binCount - 1;
        bins[idx]++;
    }});
    const maxCount = Math.max(...bins, 1);

    // Draw axes
    drawAxes(ctx, margin, w, h, '', '');

    // Y ticks
    const yTicks = niceTicksForRange(0, maxCount, 4);
    drawYTicks(ctx, margin, h, yTicks.map(t => ({{
        label: t.label,
        y: margin.top + plotH - (t.val / maxCount) * plotH
    }})));

    // X ticks
    const xTickVals = niceTicksForRange(minV, maxV, 5);
    drawXTicks(ctx, margin, w, h, xTickVals.map(t => ({{
        label: t.label,
        x: margin.left + ((t.val - minV) / range) * plotW
    }})));

    // Draw bars
    const barW = plotW / binCount;
    const barRects = [];
    for (let i = 0; i < binCount; i++) {{
        const bx = margin.left + i * barW;
        const bh = (bins[i] / maxCount) * plotH;
        const by = margin.top + plotH - bh;
        ctx.fillStyle = color;
        ctx.fillRect(bx + 1, by, barW - 2, bh);
        barRects.push({{ x: bx, y: by, w: barW, h: bh, lo: minV + i * binWidth, hi: minV + (i + 1) * binWidth, count: bins[i] }});
    }}

    // Draw vertical lines (mean, median, n50)
    const verticals = [];
    if (meanVal != null) verticals.push({{ val: meanVal, color: meanColor || '#60a5fa', label: 'Mean' }});
    if (extraVerticals) verticals.push(...extraVerticals);
    verticals.forEach(vl => {{
        const vx = margin.left + ((vl.val - minV) / range) * plotW;
        if (vx >= margin.left && vx <= margin.left + plotW) {{
            ctx.strokeStyle = vl.color;
            ctx.lineWidth = 2;
            ctx.setLineDash([4, 3]);
            ctx.beginPath();
            ctx.moveTo(vx, margin.top);
            ctx.lineTo(vx, margin.top + plotH);
            ctx.stroke();
            ctx.setLineDash([]);
            ctx.fillStyle = vl.color;
            ctx.font = '10px sans-serif';
            ctx.textAlign = 'center';
            ctx.fillText(vl.label + ': ' + (Number.isInteger(vl.val) ? vl.val : vl.val.toFixed(1)), vx, margin.top - 3);
        }}
    }});

    // Hover
    canvas.addEventListener('mousemove', function(e) {{
        const rect = canvas.getBoundingClientRect();
        const mx = (e.clientX - rect.left);
        const my = (e.clientY - rect.top);
        let hit = null;
        for (const br of barRects) {{
            if (mx >= br.x && mx <= br.x + br.w && my >= br.y && my <= br.y + br.h) {{
                hit = br; break;
            }}
        }}
        if (hit) {{
            showTooltip(tipId, canvas, mx, my,
                hit.lo.toFixed(1) + ' - ' + hit.hi.toFixed(1) + ': <b>' + hit.count + '</b> reads');
        }} else {{
            hideTooltip(tipId);
        }}
    }});
    canvas.addEventListener('mouseleave', () => hideTooltip(tipId));
}}

// ===================== STACKED HISTOGRAM (full vs trunc) =====================
function drawStackedLengthHistogram(canvasId, tipId) {{
    const setup = setupCanvas(canvasId);
    if (!setup || lengthsAll.length === 0) return;
    const {{ canvas, ctx, w, h }} = setup;
    const margin = {{ top: 20, right: 20, bottom: 30, left: 55 }};
    const plotW = w - margin.left - margin.right;
    const plotH = h - margin.top - margin.bottom;

    const allVals = lengthsAll;
    const minV = Math.min(...allVals);
    const maxV = Math.max(...allVals);
    const range = maxV - minV || 1;
    const binCount = Math.min(50, Math.max(10, Math.ceil(Math.sqrt(allVals.length))));
    const binWidth = range / binCount;

    const binsFull = new Array(binCount).fill(0);
    const binsTrunc = new Array(binCount).fill(0);
    lengthsFull.forEach(v => {{
        let idx = Math.floor((v - minV) / binWidth);
        if (idx >= binCount) idx = binCount - 1;
        binsFull[idx]++;
    }});
    lengthsTrunc.forEach(v => {{
        let idx = Math.floor((v - minV) / binWidth);
        if (idx >= binCount) idx = binCount - 1;
        binsTrunc[idx]++;
    }});

    const maxCount = Math.max(...binsFull.map((f, i) => f + binsTrunc[i]), 1);

    drawAxes(ctx, margin, w, h, 'Read Length (bp)', 'Count');

    const yTicks = niceTicksForRange(0, maxCount, 4);
    drawYTicks(ctx, margin, h, yTicks.map(t => ({{
        label: t.label,
        y: margin.top + plotH - (t.val / maxCount) * plotH
    }})));

    const xTickVals = niceTicksForRange(minV, maxV, 6);
    drawXTicks(ctx, margin, w, h, xTickVals.map(t => ({{
        label: t.label,
        x: margin.left + ((t.val - minV) / range) * plotW
    }})));

    const barW = plotW / binCount;
    const barRects = [];
    for (let i = 0; i < binCount; i++) {{
        const bx = margin.left + i * barW;
        const totalH = ((binsFull[i] + binsTrunc[i]) / maxCount) * plotH;
        const fullH = (binsFull[i] / maxCount) * plotH;
        const truncH = (binsTrunc[i] / maxCount) * plotH;

        // Draw trunc on bottom, full on top (stacked)
        const baseY = margin.top + plotH;
        // Truncated portion
        ctx.fillStyle = '#fbbf24';
        ctx.fillRect(bx + 1, baseY - truncH, barW - 2, truncH);
        // Full-length portion
        ctx.fillStyle = '#4ade80';
        ctx.fillRect(bx + 1, baseY - truncH - fullH, barW - 2, fullH);

        barRects.push({{ x: bx, w: barW, y: baseY - totalH, h: totalH,
            lo: minV + i * binWidth, hi: minV + (i + 1) * binWidth,
            full: binsFull[i], trunc: binsTrunc[i] }});
    }}

    // Vertical lines
    const verticals = [];
    if (lengthMean != null) verticals.push({{ val: lengthMean, color: '#60a5fa', label: 'Mean' }});
    if (lengthMedian != null) verticals.push({{ val: lengthMedian, color: '#a78bfa', label: 'Median' }});
    if (lengthN50 != null) verticals.push({{ val: lengthN50, color: '#f87171', label: 'N50' }});
    verticals.forEach(vl => {{
        const vx = margin.left + ((vl.val - minV) / range) * plotW;
        if (vx >= margin.left && vx <= margin.left + plotW) {{
            ctx.strokeStyle = vl.color;
            ctx.lineWidth = 2;
            ctx.setLineDash([4, 3]);
            ctx.beginPath();
            ctx.moveTo(vx, margin.top);
            ctx.lineTo(vx, margin.top + plotH);
            ctx.stroke();
            ctx.setLineDash([]);
            ctx.fillStyle = vl.color;
            ctx.font = '10px sans-serif';
            ctx.textAlign = 'center';
            ctx.fillText(vl.label + ': ' + Math.round(vl.val), vx, margin.top - 5);
        }}
    }});

    canvas.addEventListener('mousemove', function(e) {{
        const rect = canvas.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;
        let hit = null;
        for (const br of barRects) {{
            if (mx >= br.x && mx <= br.x + br.w && my >= br.y && my <= br.y + br.h) {{
                hit = br; break;
            }}
        }}
        if (hit) {{
            showTooltip(tipId, canvas, mx, my,
                Math.round(hit.lo) + ' - ' + Math.round(hit.hi) + ' bp<br>' +
                'Full-length: <b>' + hit.full + '</b><br>Truncated: <b>' + hit.trunc + '</b>');
        }} else {{
            hideTooltip(tipId);
        }}
    }});
    canvas.addEventListener('mouseleave', () => hideTooltip(tipId));
}}

// ===================== BOX PLOT =====================
function drawBoxPlot(canvasId, tipId) {{
    const setup = setupCanvas(canvasId);
    if (!setup) return;
    const {{ canvas, ctx, w, h }} = setup;
    const margin = {{ top: 15, right: 20, bottom: 40, left: 55 }};
    const plotW = w - margin.left - margin.right;
    const plotH = h - margin.top - margin.bottom;

    const targets = Object.keys(perTargetLengths).filter(t => perTargetLengths[t].length > 0);
    if (targets.length === 0) return;

    // Compute stats for each target
    const stats = targets.map(t => {{
        const vals = perTargetLengths[t].slice().sort((a, b) => a - b);
        const n = vals.length;
        const q1 = vals[Math.floor(n * 0.25)];
        const median = vals[Math.floor(n * 0.5)];
        const q3 = vals[Math.floor(n * 0.75)];
        return {{ tgt: t, min: vals[0], max: vals[n - 1], q1, median, q3, n }};
    }});

    const globalMin = Math.min(...stats.map(s => s.min));
    const globalMax = Math.max(...stats.map(s => s.max));
    const range = globalMax - globalMin || 1;

    drawAxes(ctx, margin, w, h, 'Read Length (bp)', '');

    // X ticks
    const xTickVals = niceTicksForRange(globalMin, globalMax, 6);
    drawXTicks(ctx, margin, w, h, xTickVals.map(t => ({{
        label: t.label,
        x: margin.left + ((t.val - globalMin) / range) * plotW
    }})));

    const boxH = Math.min(24, (plotH - 10) / targets.length - 4);
    const boxRects = [];

    stats.forEach((s, i) => {{
        const cy = margin.top + (i + 0.5) * (plotH / targets.length);
        const toX = v => margin.left + ((v - globalMin) / range) * plotW;
        const color = TARGET_PALETTE[i % TARGET_PALETTE.length];

        // Target label
        ctx.fillStyle = CHART_TEXT_LIGHT;
        ctx.font = '11px sans-serif';
        ctx.textAlign = 'right';
        // Draw label on left if space, otherwise inside
        ctx.fillStyle = CHART_TEXT;

        // Whiskers (min to max)
        ctx.strokeStyle = color;
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(toX(s.min), cy);
        ctx.lineTo(toX(s.max), cy);
        ctx.stroke();
        // Min cap
        ctx.beginPath();
        ctx.moveTo(toX(s.min), cy - boxH / 4);
        ctx.lineTo(toX(s.min), cy + boxH / 4);
        ctx.stroke();
        // Max cap
        ctx.beginPath();
        ctx.moveTo(toX(s.max), cy - boxH / 4);
        ctx.lineTo(toX(s.max), cy + boxH / 4);
        ctx.stroke();

        // Box Q1-Q3
        const bx = toX(s.q1);
        const bw = toX(s.q3) - bx;
        ctx.fillStyle = color + '44';
        ctx.fillRect(bx, cy - boxH / 2, bw, boxH);
        ctx.strokeStyle = color;
        ctx.lineWidth = 1.5;
        ctx.strokeRect(bx, cy - boxH / 2, bw, boxH);

        // Median line
        ctx.strokeStyle = '#fff';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(toX(s.median), cy - boxH / 2);
        ctx.lineTo(toX(s.median), cy + boxH / 2);
        ctx.stroke();

        // Label
        ctx.fillStyle = CHART_TEXT_LIGHT;
        ctx.font = '11px sans-serif';
        ctx.textAlign = 'left';
        ctx.fillText(s.tgt + ' (n=' + s.n + ')', margin.left + 5, cy - boxH / 2 - 3);

        boxRects.push({{ y: cy - boxH / 2, h: boxH, s }});
    }});

    canvas.addEventListener('mousemove', function(e) {{
        const rect = canvas.getBoundingClientRect();
        const my = e.clientY - rect.top;
        const mx = e.clientX - rect.left;
        let hit = null;
        for (const br of boxRects) {{
            if (my >= br.y && my <= br.y + br.h) {{ hit = br.s; break; }}
        }}
        if (hit) {{
            showTooltip(tipId, canvas, mx, my,
                '<b>' + hit.tgt + '</b> (n=' + hit.n + ')<br>' +
                'Min: ' + hit.min + ' | Q1: ' + hit.q1 + '<br>' +
                'Median: ' + hit.median + ' | Q3: ' + hit.q3 + '<br>' +
                'Max: ' + hit.max);
        }} else {{
            hideTooltip(tipId);
        }}
    }});
    canvas.addEventListener('mouseleave', () => hideTooltip(tipId));
}}

// ===================== SCATTER PLOT =====================
function drawScatterPlot(canvasId, tipId) {{
    const setup = setupCanvas(canvasId);
    if (!setup || scatterData.length === 0) return;
    const {{ canvas, ctx, w, h }} = setup;
    const margin = {{ top: 15, right: 20, bottom: 30, left: 55 }};
    const plotW = w - margin.left - margin.right;
    const plotH = h - margin.top - margin.bottom;

    const targets = [...new Set(scatterData.map(d => d.tgt))];

    drawAxes(ctx, margin, w, h, 'Start Confidence', 'End Confidence');

    // Ticks (0..1 range)
    const ticks01 = [0, 0.2, 0.4, 0.6, 0.8, 1.0];
    drawXTicks(ctx, margin, w, h, ticks01.map(t => ({{
        label: t.toFixed(1),
        x: margin.left + t * plotW
    }})));
    drawYTicks(ctx, margin, h, ticks01.map(t => ({{
        label: t.toFixed(1),
        y: margin.top + plotH - t * plotH
    }})));

    // Grid lines
    ctx.strokeStyle = CHART_GRID;
    ctx.lineWidth = 0.5;
    ticks01.forEach(t => {{
        const gx = margin.left + t * plotW;
        const gy = margin.top + plotH - t * plotH;
        ctx.beginPath(); ctx.moveTo(gx, margin.top); ctx.lineTo(gx, margin.top + plotH); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(margin.left, gy); ctx.lineTo(margin.left + plotW, gy); ctx.stroke();
    }});

    // Points
    const points = [];
    scatterData.forEach(d => {{
        const px = margin.left + (d.x || 0) * plotW;
        const py = margin.top + plotH - (d.y || 0) * plotH;
        const color = getTargetColor(d.tgt, targets);
        const r = 3.5;

        ctx.beginPath();
        ctx.arc(px, py, r, 0, Math.PI * 2);
        if (d.fl) {{
            ctx.fillStyle = color;
            ctx.fill();
        }} else {{
            ctx.strokeStyle = color;
            ctx.lineWidth = 1.5;
            ctx.stroke();
        }}
        points.push({{ px, py, r, d }});
    }});

    canvas.addEventListener('mousemove', function(e) {{
        const rect = canvas.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;
        let hit = null;
        for (const p of points) {{
            const dx = mx - p.px, dy = my - p.py;
            if (dx * dx + dy * dy < (p.r + 4) * (p.r + 4)) {{ hit = p.d; break; }}
        }}
        if (hit) {{
            showTooltip(tipId, canvas, mx, my,
                'Target: <b>' + hit.tgt + '</b><br>' +
                'Start: ' + (hit.x || 0).toFixed(3) + ' | End: ' + (hit.y || 0).toFixed(3) + '<br>' +
                'Length: ' + (hit.len || '?') + 'bp | ' + (hit.fl ? 'Full-length' : 'Truncated'));
        }} else {{
            hideTooltip(tipId);
        }}
    }});
    canvas.addEventListener('mouseleave', () => hideTooltip(tipId));
}}

// ===================== CONFIDENCE THRESHOLD CURVE =====================
function drawConfCurve(canvasId, tipId) {{
    const setup = setupCanvas(canvasId);
    if (!setup || confCurveData.length === 0) return;
    const {{ canvas, ctx, w, h }} = setup;
    const margin = {{ top: 15, right: 20, bottom: 30, left: 55 }};
    const plotW = w - margin.left - margin.right;
    const plotH = h - margin.top - margin.bottom;

    const maxMatched = Math.max(...confCurveData.map(d => d.matched), 1);

    drawAxes(ctx, margin, w, h, 'Confidence Threshold', 'Matched Reads');

    // X ticks
    const xRange = confCurveData.length > 1 ?
        [confCurveData[0].threshold, confCurveData[confCurveData.length - 1].threshold] : [0, 1];
    const xTickVals = niceTicksForRange(xRange[0], xRange[1], 5);
    const xSpan = xRange[1] - xRange[0] || 1;
    drawXTicks(ctx, margin, w, h, xTickVals.map(t => ({{
        label: t.label,
        x: margin.left + ((t.val - xRange[0]) / xSpan) * plotW
    }})));

    // Y ticks
    const yTicks = niceTicksForRange(0, maxMatched, 4);
    drawYTicks(ctx, margin, h, yTicks.map(t => ({{
        label: t.label,
        y: margin.top + plotH - (t.val / maxMatched) * plotH
    }})));

    // Draw line
    ctx.strokeStyle = '#60a5fa';
    ctx.lineWidth = 2;
    ctx.beginPath();
    const linePoints = [];
    confCurveData.forEach((d, i) => {{
        const px = margin.left + ((d.threshold - xRange[0]) / xSpan) * plotW;
        const py = margin.top + plotH - (d.matched / maxMatched) * plotH;
        if (i === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
        linePoints.push({{ px, py, d }});
    }});
    ctx.stroke();

    // Fill area
    ctx.fillStyle = 'rgba(96,165,250,0.1)';
    ctx.beginPath();
    linePoints.forEach((lp, i) => {{
        if (i === 0) ctx.moveTo(lp.px, lp.py);
        else ctx.lineTo(lp.px, lp.py);
    }});
    ctx.lineTo(linePoints[linePoints.length - 1].px, margin.top + plotH);
    ctx.lineTo(linePoints[0].px, margin.top + plotH);
    ctx.closePath();
    ctx.fill();

    // Click to show threshold
    canvas.addEventListener('click', function(e) {{
        const rect = canvas.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        let closest = null;
        let minDist = Infinity;
        for (const lp of linePoints) {{
            const dist = Math.abs(mx - lp.px);
            if (dist < minDist) {{ minDist = dist; closest = lp; }}
        }}
        if (closest) {{
            const info = document.getElementById('conf-curve-info');
            info.innerHTML = 'At threshold <b>' + closest.d.threshold.toFixed(2) + '</b>: ' +
                '<b>' + closest.d.matched + '</b> matched reads (' + (closest.d.pct || 0).toFixed(1) + '%)';
        }}
    }});

    canvas.addEventListener('mousemove', function(e) {{
        const rect = canvas.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;
        let closest = null;
        let minDist = Infinity;
        for (const lp of linePoints) {{
            const dist = Math.abs(mx - lp.px);
            if (dist < minDist) {{ minDist = dist; closest = lp; }}
        }}
        if (closest && minDist < 20) {{
            showTooltip(tipId, canvas, mx, my,
                'Threshold: ' + closest.d.threshold.toFixed(2) + '<br>' +
                'Matched: <b>' + closest.d.matched + '</b> (' + (closest.d.pct || 0).toFixed(1) + '%)');
        }} else {{
            hideTooltip(tipId);
        }}
    }});
    canvas.addEventListener('mouseleave', () => hideTooltip(tipId));
}}

// ===================== END REASON BAR CHART =====================
function drawEndReasonChart(canvasId, tipId) {{
    const setup = setupCanvas(canvasId);
    if (!setup) return;
    const {{ canvas, ctx, w, h }} = setup;
    const margin = {{ top: 15, right: 20, bottom: 15, left: 200 }};
    const plotW = w - margin.left - margin.right;
    const plotH = h - margin.top - margin.bottom;

    const reasons = Object.keys(endReasonCounts);
    if (reasons.length === 0) return;
    reasons.sort((a, b) => endReasonCounts[b] - endReasonCounts[a]);

    const maxCount = Math.max(...reasons.map(r => endReasonCounts[r]), 1);
    const totalCount = reasons.reduce((s, r) => s + endReasonCounts[r], 0);
    const barH = Math.min(22, (plotH - 5) / reasons.length - 3);

    const barColors = ['#4ade80', '#60a5fa', '#fbbf24', '#f87171', '#a78bfa', '#fb923c', '#2dd4bf', '#f472b6'];
    const barRects = [];

    reasons.forEach((reason, i) => {{
        const count = endReasonCounts[reason];
        const cy = margin.top + (i + 0.5) * (plotH / reasons.length);
        const bw = (count / maxCount) * plotW;
        const color = barColors[i % barColors.length];

        // Label
        ctx.fillStyle = CHART_TEXT;
        ctx.font = '11px sans-serif';
        ctx.textAlign = 'right';
        const label = reason.length > 28 ? reason.substring(0, 26) + '..' : reason;
        ctx.fillText(label, margin.left - 8, cy + 4);

        // Bar
        ctx.fillStyle = color;
        ctx.fillRect(margin.left, cy - barH / 2, bw, barH);

        // Count label
        const pct = ((count / totalCount) * 100).toFixed(1);
        ctx.fillStyle = CHART_TEXT_LIGHT;
        ctx.font = '11px sans-serif';
        ctx.textAlign = 'left';
        ctx.fillText(count + ' (' + pct + '%)', margin.left + bw + 5, cy + 4);

        barRects.push({{ y: cy - barH / 2, h: barH, reason, count, pct }});
    }});

    canvas.addEventListener('mousemove', function(e) {{
        const rect = canvas.getBoundingClientRect();
        const my = e.clientY - rect.top;
        const mx = e.clientX - rect.left;
        let hit = null;
        for (const br of barRects) {{
            if (my >= br.y && my <= br.y + br.h) {{ hit = br; break; }}
        }}
        if (hit) {{
            showTooltip(tipId, canvas, mx, my,
                '<b>' + hit.reason + '</b>: ' + hit.count + ' (' + hit.pct + '%)');
        }} else {{
            hideTooltip(tipId);
        }}
    }});
    canvas.addEventListener('mouseleave', () => hideTooltip(tipId));
}}

// ===================== RENDER ALL CHARTS =====================
let chartsRendered = false;
function renderAllCharts() {{
    if (chartsRendered) return;
    chartsRendered = true;

    // (a) Read length distribution (stacked)
    drawStackedLengthHistogram('canvas-length-hist', 'tip-length-hist');

    // (b) Per-target box plot
    drawBoxPlot('canvas-target-box', 'tip-target-box');

    // (c) Quality distributions
    drawCanvasHistogram('canvas-qbc-hist', 'tip-qbc-hist', qBcValues, '#7c3aed', meanQBc, '#60a5fa', []);
    drawCanvasHistogram('canvas-qld-hist', 'tip-qld-hist', qLdValues, '#be185d', meanQLd, '#60a5fa', []);

    // (d) Edit distance distribution
    drawCanvasHistogram('canvas-ed-hist', 'tip-ed-hist', edValues, '#f59e0b', meanEd, '#60a5fa', []);

    // (e) Scatter plot
    drawScatterPlot('canvas-conf-scatter', 'tip-conf-scatter');

    // (f) Confidence threshold curve
    drawConfCurve('canvas-conf-curve', 'tip-conf-curve');

    // (g) End reason breakdown
    drawEndReasonChart('canvas-end-reason', 'tip-end-reason');
}}

// ===================== INIT =====================
document.addEventListener('DOMContentLoaded', function() {{
    // Legacy div-based histograms for Summary tab
    drawHist('start-hist', startConfs, '#7c3aed');
    drawHist('end-hist', endConfs, '#be185d');

    // Populate pair filter dropdown
    const pairSelect = document.getElementById('filter-pair');
    uniquePairs.forEach(p => {{
        const opt = document.createElement('option');
        opt.value = p;
        opt.textContent = p;
        pairSelect.appendChild(opt);
    }});

    // Initial row count
    filterTable();
}});
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
