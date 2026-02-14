"""SMA-seq barcode calibration visualizer - FastAPI application."""

from __future__ import annotations

import sys
import webbrowser
from pathlib import Path

import uvicorn
from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from calibrate_viz.api import router as api_router, set_db_paths

app = FastAPI(title="SMA-seq Barcode Calibration")
app.include_router(api_router)

_HERE = Path(__file__).resolve().parent

templates = Jinja2Templates(directory=_HERE / "templates")
app.mount("/static", StaticFiles(directory=_HERE / "static"), name="static")

_db_paths: list[Path] = []


def set_databases(paths: list[Path]) -> None:
    """Set the list of calibration database paths."""
    global _db_paths
    _db_paths = [p.resolve() for p in paths]
    set_db_paths(_db_paths)


def _ctx(request: Request, active_page: str) -> dict:
    """Build common template context shared by all pages."""
    return {
        "request": request,
        "active_page": active_page,
        "db_paths": _db_paths,
        "db_names": [p.stem for p in _db_paths],
    }


# ---------------------------------------------------------------------------
# Page routes
# ---------------------------------------------------------------------------


@app.get("/")
async def overview(request: Request):
    """Render the overview page."""
    ctx = _ctx(request, "overview")
    return templates.TemplateResponse("overview.html", ctx)


@app.get("/distributions")
async def distributions(request: Request):
    """Render the distributions page."""
    ctx = _ctx(request, "distributions")
    return templates.TemplateResponse("distributions.html", ctx)


@app.get("/confidence")
async def confidence(request: Request):
    """Render the confidence page."""
    ctx = _ctx(request, "confidence")
    return templates.TemplateResponse("confidence.html", ctx)


@app.get("/separation")
async def separation(request: Request):
    """Render the separation page."""
    ctx = _ctx(request, "separation")
    return templates.TemplateResponse("separation.html", ctx)


@app.get("/thresholds")
async def thresholds(request: Request):
    """Render the thresholds page."""
    ctx = _ctx(request, "thresholds")
    return templates.TemplateResponse("thresholds.html", ctx)


@app.get("/compare")
async def compare(request: Request):
    """Render the compare page."""
    ctx = _ctx(request, "compare")
    return templates.TemplateResponse("compare.html", ctx)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Entry point: parse DB paths from argv, launch server and browser."""
    if len(sys.argv) < 2:
        print("Usage: python -m calibrate_viz.app <db_path> [db_path ...]")
        sys.exit(1)

    db_paths = [Path(p) for p in sys.argv[1:]]
    for p in db_paths:
        if not p.exists():
            print(f"Error: database not found: {p}")
            sys.exit(1)

    set_databases(db_paths)

    host = "127.0.0.1"
    port = 8051
    webbrowser.open(f"http://{host}:{port}")
    uvicorn.run(app, host=host, port=port)


if __name__ == "__main__":
    main()
