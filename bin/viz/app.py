"""SMA-seq experiment configuration visualizer - FastAPI application."""

from __future__ import annotations

import sys
import webbrowser
from pathlib import Path

import uvicorn
from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from viz.api import router as api_router, set_store as set_api_store
from viz.config_store import ConfigStore

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
    return templates.TemplateResponse("assumptions.html", ctx)


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
