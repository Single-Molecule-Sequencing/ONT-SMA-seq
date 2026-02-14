"""SMA-seq experiment configuration visualizer - FastAPI application."""

from __future__ import annotations

import sys
import webbrowser
from pathlib import Path

import edlib
import uvicorn
from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

import barcodes as _barcodes_mod
from viz.api import router as api_router, set_store as set_api_store
from viz.config_store import ConfigStore
from viz.models import Assumption, SampleSheetEntry

app = FastAPI(title="SMA-seq Config Visualizer")
app.include_router(api_router)

_HERE = Path(__file__).resolve().parent

templates = Jinja2Templates(directory=_HERE / "templates")
app.mount("/static", StaticFiles(directory=_HERE / "static"), name="static")

experiment_dir: Path | None = None
_store: ConfigStore | None = None


def set_experiment_dir(path: Path) -> None:
    """Set the global experiment directory and initialise the ConfigStore."""
    global experiment_dir, _store
    experiment_dir = path.resolve()
    _store = ConfigStore(experiment_dir)
    set_api_store(_store)


def _ctx(request: Request, active_page: str) -> dict:
    """Build common template context shared by all pages."""
    return {
        "request": request,
        "active_page": active_page,
        "experiment_dir": str(experiment_dir),
        "experiment_dir_name": experiment_dir.name if experiment_dir else "",
    }


# ---------------------------------------------------------------------------
# Page routes
# ---------------------------------------------------------------------------


@app.get("/")
async def dashboard(request: Request):
    """Render the main dashboard."""
    ctx = _ctx(request, "dashboard")
    ctx["config"] = _store.experiment_config if _store else None
    ctx["file_status"] = {
        "sample_sheet.csv": (experiment_dir / "sample_sheet.csv").exists(),
        "arrangement.toml": (experiment_dir / "arrangement.toml").exists(),
        "barcodes.fasta": (experiment_dir / "barcodes.fasta").exists(),
        "references/": (experiment_dir / "references").is_dir(),
        "sma_experiment.toml": (experiment_dir / "sma_experiment.toml").exists(),
    } if experiment_dir else {}
    return templates.TemplateResponse("dashboard.html", ctx)


@app.get("/sample-sheet")
async def sample_sheet_page(request: Request):
    """Render the sample sheet page."""
    ctx = _ctx(request, "sample_sheet")
    return templates.TemplateResponse("sample_sheet.html", ctx)


@app.get("/barcodes")
async def barcodes_page(request: Request):
    """Render the barcodes page."""
    ctx = _ctx(request, "barcodes")
    ctx["arrangement"] = _store.read_arrangement() if _store else None
    return templates.TemplateResponse("barcodes.html", ctx)


@app.get("/construct")
async def construct_page(request: Request):
    """Render the construct page."""
    ctx = _ctx(request, "construct")
    return templates.TemplateResponse("construct.html", ctx)


@app.get("/targets")
async def targets_page(request: Request):
    """Render the targets page."""
    ctx = _ctx(request, "targets")
    return templates.TemplateResponse("targets.html", ctx)


@app.get("/assumptions")
async def assumptions_page(request: Request):
    """Render the assumptions page."""
    ctx = _ctx(request, "assumptions")
    cfg = _store.experiment_config if _store else None
    ctx["classification"] = cfg.classification if cfg else None
    return templates.TemplateResponse("assumptions.html", ctx)


# ---------------------------------------------------------------------------
# HTMX fragment routes - Sample Sheet
# ---------------------------------------------------------------------------


@app.get("/htmx/sample-sheet/rows", response_class=HTMLResponse)
async def htmx_sample_sheet_rows():
    """Return sample sheet table rows as an HTML fragment."""
    entries = _store.read_sample_sheet()
    rows = ""
    for i, e in enumerate(entries):
        rows += f"""<tr>
            <td>{e.barcode}</td>
            <td>{e.upstream_barcode}</td>
            <td>{e.downstream_barcode}</td>
            <td>{e.alias}</td>
            <td>{e.kit}</td>
            <td>{e.type}</td>
            <td><button hx-delete="/htmx/sample-sheet/delete/{i}" hx-target="#sample-sheet-body" hx-swap="innerHTML" class="outline secondary">Delete</button></td>
        </tr>"""
    return rows or "<tr><td colspan='7'>No entries</td></tr>"


@app.post("/htmx/sample-sheet/add", response_class=HTMLResponse)
async def htmx_sample_sheet_add(request: Request):
    """Add a new sample sheet entry from form data and return updated rows."""
    form = await request.form()
    upstream = form["upstream"]
    downstream = form["downstream"]
    alias = form["alias"]
    kit = form.get("kit", "")
    entries = _store.read_sample_sheet()
    entries.append(SampleSheetEntry(
        barcode=f"{upstream}--{downstream}",
        alias=alias,
        kit=kit,
    ))
    _store.write_sample_sheet(entries)
    return await htmx_sample_sheet_rows()


@app.delete("/htmx/sample-sheet/delete/{idx}", response_class=HTMLResponse)
async def htmx_sample_sheet_delete(idx: int):
    """Delete a sample sheet entry by index and return updated rows."""
    entries = _store.read_sample_sheet()
    if 0 <= idx < len(entries):
        entries.pop(idx)
        _store.write_sample_sheet(entries)
    return await htmx_sample_sheet_rows()


# ---------------------------------------------------------------------------
# HTMX fragment routes - Barcodes
# ---------------------------------------------------------------------------


def _gc_percent(seq: str) -> float:
    """Compute GC percentage of a DNA sequence."""
    if not seq:
        return 0.0
    gc_count = sum(1 for base in seq.upper() if base in ("G", "C"))
    return round((gc_count / len(seq)) * 100.0, 2)


@app.get("/htmx/barcodes/used-table", response_class=HTMLResponse)
async def htmx_barcodes_used_table():
    """Return a table of ONT barcodes used in the sample sheet."""
    entries = _store.read_sample_sheet()

    # Collect unique barcode IDs referenced by sample sheet
    used_ids: set[str] = set()
    for entry in entries:
        if entry.upstream_barcode:
            used_ids.add(entry.upstream_barcode)
        if entry.downstream_barcode:
            used_ids.add(entry.downstream_barcode)

    all_barcodes = _barcodes_mod.BARCODES
    custom_seqs = _store.read_barcode_sequences()
    merged = {**all_barcodes, **custom_seqs}

    if not merged:
        return "<p>No barcodes available.</p>"

    rows = """<table>
        <thead>
            <tr>
                <th>ID</th>
                <th>Sequence</th>
                <th>Reverse Complement</th>
                <th>GC%</th>
                <th>Status</th>
            </tr>
        </thead>
        <tbody>"""

    for bc_id in sorted(merged.keys()):
        seq = merged[bc_id]
        rc = _barcodes_mod.reverse_complement(seq)
        gc = _gc_percent(seq)
        is_used = bc_id in used_ids
        badge = ('<span class="badge-ok">Used</span>'
                 if is_used
                 else '<span style="color: var(--pico-muted-color);">Unused</span>')
        rows += f"""<tr>
            <td>{bc_id}</td>
            <td><code>{seq}</code></td>
            <td><code>{rc}</code></td>
            <td>{gc:.1f}%</td>
            <td>{badge}</td>
        </tr>"""

    rows += "</tbody></table>"
    return rows


@app.put("/htmx/barcodes/arrangement", response_class=HTMLResponse)
async def htmx_barcodes_arrangement(request: Request):
    """Update arrangement config from form data."""
    from viz.models import ArrangementConfig, ScoringParams

    form = await request.form()

    scoring = ScoringParams(
        max_barcode_penalty=int(form.get("scoring_max_barcode_penalty", 11)),
        min_barcode_penalty_dist=int(form.get("scoring_min_barcode_penalty_dist", 3)),
        min_separation_only_dist=int(form.get("scoring_min_separation_only_dist", 6)),
        barcode_end_proximity=int(form.get("scoring_barcode_end_proximity", 75)),
        front_barcode_window=int(form.get("scoring_front_barcode_window", 175)),
        rear_barcode_window=int(form.get("scoring_rear_barcode_window", 175)),
    )

    arr = ArrangementConfig(
        name=form.get("name", ""),
        kit=form.get("kit", ""),
        mask1_front=form.get("mask1_front", ""),
        mask1_rear=form.get("mask1_rear", ""),
        mask2_front=form.get("mask2_front", ""),
        mask2_rear=form.get("mask2_rear", ""),
        barcode1_pattern=form.get("barcode1_pattern", "BC%02i"),
        barcode2_pattern=form.get("barcode2_pattern", "BC%02i"),
        first_index=int(form.get("first_index", 1)),
        last_index=int(form.get("last_index", 96)),
        scoring=scoring,
    )

    _store.write_arrangement(arr)
    return HTMLResponse(content="<p class='badge-ok'>Arrangement saved.</p>")


@app.get("/htmx/barcodes/distances", response_class=HTMLResponse)
async def htmx_barcodes_distances():
    """Render pairwise edit distance matrix as an HTML table."""
    entries = _store.read_sample_sheet()

    # Collect unique barcode IDs referenced by sample sheet
    used_ids: set[str] = set()
    for entry in entries:
        if entry.upstream_barcode:
            used_ids.add(entry.upstream_barcode)
        if entry.downstream_barcode:
            used_ids.add(entry.downstream_barcode)

    if len(used_ids) < 2:
        return "<p>Need at least 2 barcodes in the sample sheet to compute distances.</p>"

    # Resolve sequences (prefer custom, fall back to ONT built-in)
    custom_seqs = _store.read_barcode_sequences()
    all_barcodes = {**_barcodes_mod.BARCODES, **custom_seqs}

    bc_seqs: dict[str, str] = {}
    for bc_id in sorted(used_ids):
        if bc_id in all_barcodes:
            bc_seqs[bc_id] = all_barcodes[bc_id]

    ids = sorted(bc_seqs.keys())

    if len(ids) < 2:
        return "<p>Could not resolve sequences for enough barcodes.</p>"

    # Build distance matrix
    dist_map: dict[tuple[str, str], int] = {}
    for i, id_a in enumerate(ids):
        for j in range(i + 1, len(ids)):
            id_b = ids[j]
            result = edlib.align(
                bc_seqs[id_a], bc_seqs[id_b], mode="NW", task="distance"
            )
            dist_map[(id_a, id_b)] = result["editDistance"]
            dist_map[(id_b, id_a)] = result["editDistance"]

    # Render as HTML table
    html = "<table><thead><tr><th></th>"
    for bc_id in ids:
        html += f"<th>{bc_id}</th>"
    html += "</tr></thead><tbody>"

    for id_a in ids:
        html += f"<tr><th>{id_a}</th>"
        for id_b in ids:
            if id_a == id_b:
                html += "<td>-</td>"
            else:
                ed = dist_map.get((id_a, id_b), "?")
                html += f"<td>{ed}</td>"
        html += "</tr>"

    html += "</tbody></table>"
    return html


# ---------------------------------------------------------------------------
# HTMX fragment routes - Targets
# ---------------------------------------------------------------------------


@app.get("/htmx/targets/rows", response_class=HTMLResponse)
async def htmx_targets_rows():
    """Return target reference table rows as an HTML fragment."""
    targets = _store.read_targets()
    rows = ""
    for t in targets:
        rows += f"""<tr>
            <td>{t.tgt_id}</td>
            <td>{t.length}</td>
            <td>{t.gc_content:.1f}%</td>
            <td><button hx-delete="/htmx/targets/delete/{t.tgt_id}" hx-target="#targets-body" hx-swap="innerHTML" class="outline secondary">Delete</button></td>
        </tr>"""
    return rows or "<tr><td colspan='4'>No targets</td></tr>"


@app.post("/htmx/targets/add", response_class=HTMLResponse)
async def htmx_targets_add(request: Request):
    """Add a new target reference from form data and return updated rows."""
    form = await request.form()
    tgt_id = form["tgt_id"].strip()
    raw_sequence = form["sequence"].strip()

    # Parse sequence: strip FASTA header lines if present, join remaining
    lines = raw_sequence.splitlines()
    seq_lines = [
        line.strip()
        for line in lines
        if line.strip() and not line.strip().startswith(">")
    ]
    sequence = "".join(seq_lines).upper()

    _store.write_target(tgt_id, sequence)
    return await htmx_targets_rows()


@app.delete("/htmx/targets/delete/{alias}", response_class=HTMLResponse)
async def htmx_targets_delete(alias: str):
    """Delete a target reference and return updated rows."""
    _store.delete_target(alias)
    return await htmx_targets_rows()


@app.get("/htmx/targets/auto-refs", response_class=HTMLResponse)
async def htmx_targets_auto_refs():
    """Show contents of references/auto/ directory if it exists."""
    auto_dir = _store.dir / "references" / "auto"
    if not auto_dir.exists():
        return (
            "<p>No auto-generated references found. "
            "Use the Construct page to generate them.</p>"
        )

    fasta_files = sorted(auto_dir.glob("*.fasta"))
    if not fasta_files:
        return "<p>Auto-references directory exists but contains no FASTA files.</p>"

    html = "<table><thead><tr><th>File</th><th>Size</th></tr></thead><tbody>"
    for f in fasta_files:
        size = f.stat().st_size
        html += f"<tr><td><code>{f.name}</code></td><td>{size} bytes</td></tr>"
    html += "</tbody></table>"
    return html


# ---------------------------------------------------------------------------
# HTMX fragment routes - Assumptions
# ---------------------------------------------------------------------------


@app.get("/htmx/assumptions/rows", response_class=HTMLResponse)
async def htmx_assumptions_rows():
    """Return assumption table rows as an HTML fragment."""
    assumptions = _store.experiment_config.assumptions
    rows = ""
    for i, a in enumerate(assumptions):
        rows += f"""<tr>
            <td><code>{a.key}</code></td>
            <td>{a.text}</td>
            <td>{a.why}</td>
            <td><button hx-delete="/htmx/assumptions/delete/{i}" hx-target="#assumptions-body" hx-swap="innerHTML" class="outline secondary">Delete</button></td>
        </tr>"""
    return rows or "<tr><td colspan='4'>No assumptions</td></tr>"


@app.post("/htmx/assumptions/add", response_class=HTMLResponse)
async def htmx_assumptions_add(request: Request):
    """Add a new assumption from form data and return updated rows."""
    form = await request.form()
    key = form["key"].strip()
    text = form["text"].strip()
    why = form["why"].strip()
    _store.experiment_config.assumptions.append(
        Assumption(key=key, text=text, why=why)
    )
    _store.save_experiment_config()
    return await htmx_assumptions_rows()


@app.delete("/htmx/assumptions/delete/{idx}", response_class=HTMLResponse)
async def htmx_assumptions_delete(idx: int):
    """Delete an assumption by index and return updated rows."""
    assumptions = _store.experiment_config.assumptions
    if 0 <= idx < len(assumptions):
        assumptions.pop(idx)
        _store.save_experiment_config()
    return await htmx_assumptions_rows()


# ---------------------------------------------------------------------------
# HTMX fragment routes - Construct page
# ---------------------------------------------------------------------------


@app.get("/htmx/construct/assignments", response_class=HTMLResponse)
async def htmx_construct_assignments():
    """Return barcode assignment tables as an HTMX HTML fragment."""
    cfg = _store.experiment_config
    html = f"<p><strong>Mode:</strong> {cfg.demultiplexing.mode}</p>"
    if cfg.demultiplexing.start_barcode:
        html += f"<p><strong>Start barcode:</strong> {cfg.demultiplexing.start_barcode}</p>"
    if cfg.demultiplexing.end_barcode:
        html += f"<p><strong>End barcode:</strong> {cfg.demultiplexing.end_barcode}</p>"
    if cfg.demultiplexing.pairs:
        html += "<table><thead><tr><th>Start</th><th>End</th><th>Alias</th></tr></thead><tbody>"
        for p in cfg.demultiplexing.pairs:
            html += f"<tr><td>{p.start}</td><td>{p.end}</td><td>{p.alias}</td></tr>"
        html += "</tbody></table>"
    else:
        html += "<p><em>No barcode pairs configured.</em></p>"
    return html


@app.get("/htmx/construct/rules", response_class=HTMLResponse)
async def htmx_construct_rules():
    """Return truncation rules as an HTMX HTML fragment table."""
    rules = _store.experiment_config.truncation.rules
    html = "<table><thead><tr><th>Class</th><th>Rule</th></tr></thead><tbody>"
    for field in ["full", "trunc_3prime", "trunc_target", "trunc_barcode", "adapter_only", "chimeric"]:
        val = getattr(rules, field)
        html += f"<tr><td><span class='trunc-label trunc-{field}'>{field}</span></td><td>{val}</td></tr>"
    html += "</tbody></table>"
    return html


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Entry point: parse experiment dir from argv, launch server and browser."""
    if len(sys.argv) < 2:
        print("Usage: python -m bin.viz.app <experiment_dir>")
        sys.exit(1)

    set_experiment_dir(Path(sys.argv[1]))

    host = "127.0.0.1"
    port = 8050
    webbrowser.open(f"http://{host}:{port}")
    uvicorn.run(app, host=host, port=port)


if __name__ == "__main__":
    main()
