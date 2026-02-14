"""SMA-seq experiment configuration visualizer - FastAPI application."""

from __future__ import annotations

import sys
import webbrowser
from pathlib import Path

import uvicorn
from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

app = FastAPI(title="SMA-seq Config Visualizer")

_HERE = Path(__file__).resolve().parent

templates = Jinja2Templates(directory=_HERE / "templates")
app.mount("/static", StaticFiles(directory=_HERE / "static"), name="static")

experiment_dir: Path | None = None


def set_experiment_dir(path: Path) -> None:
    """Set the global experiment directory."""
    global experiment_dir
    experiment_dir = Path(path).resolve()


@app.get("/")
async def dashboard(request: Request):
    """Render the main dashboard."""
    return templates.TemplateResponse(
        "dashboard.html",
        {"request": request, "experiment_dir": str(experiment_dir)},
    )


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
