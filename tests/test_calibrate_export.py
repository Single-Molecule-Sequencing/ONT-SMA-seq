"""Tests for static HTML export of calibration results."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from calibrate_viz.export import export_calibration_report


@pytest.fixture
def populated_db(tmp_path: Path) -> Path:
    """Create a DB with enough data for a full export."""
    db_path = tmp_path / "test.db"
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE Reads (
            read_id TEXT PRIMARY KEY, readlen INTEGER,
            signal_duration_s REAL, mean_qscore REAL, ER TEXT,
            bc_start_id TEXT, bc_start_conf REAL,
            bc_end_id TEXT, bc_end_conf REAL,
            tgt_id TEXT, trunc_level TEXT, ed INTEGER
        )
    """)
    conn.execute("""
        CREATE TABLE Target (tgt_id TEXT PRIMARY KEY, tgt_reflen INTEGER)
    """)
    conn.execute("INSERT INTO Target VALUES ('V04_2_fwd', 500)")
    for i in range(50):
        conn.execute(
            "INSERT INTO Reads VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (f"r{i}", 400 + i * 5, 1.0 + i * 0.05, 10.0 + i * 0.1,
             "signal_positive", "nb05", 0.7 + i * 0.005,
             "nb10", 0.6 + i * 0.005,
             "V04_2_fwd", "full_length", 5 + i),
        )
    conn.commit()
    conn.close()
    return db_path


class TestExportCalibrationReport:
    def test_creates_output_directory(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        assert out.is_dir()

    def test_generates_all_html_pages(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        expected = ["index.html", "distributions.html", "confidence.html",
                     "separation.html", "thresholds.html"]
        for name in expected:
            assert (out / name).exists(), f"Missing {name}"

    def test_html_is_self_contained(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        html = (out / "index.html").read_text()
        assert "<style>" in html
        assert "DATA" in html or "data" in html

    def test_generates_svg_figures(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        figures = out / "figures"
        assert figures.is_dir()
        svgs = list(figures.glob("*.svg"))
        assert len(svgs) >= 2  # at least signal + readlen KDEs

    def test_svg_is_valid(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        for svg_file in (out / "figures").glob("*.svg"):
            content = svg_file.read_text()
            assert content.startswith("<svg")
            assert "</svg>" in content

    def test_custom_experiment_id(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out, experiment_id="MY_EXP")
        html = (out / "index.html").read_text()
        assert "MY_EXP" in html

    def test_default_experiment_id_from_db_stem(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        html = (out / "index.html").read_text()
        assert "test" in html  # db stem is "test"

    def test_returns_output_dir(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        result = export_calibration_report(populated_db, out)
        assert result == out

    def test_distributions_page_has_data(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        html = (out / "distributions.html").read_text()
        assert "<style>" in html
        assert "d3" in html.lower() or "D3" in html or "polyline" in html or "path" in html

    def test_confidence_page_has_data(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        html = (out / "confidence.html").read_text()
        assert "<style>" in html

    def test_nav_links_present(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        for name in ["index.html", "distributions.html", "confidence.html",
                      "separation.html", "thresholds.html"]:
            html = (out / name).read_text()
            assert "index.html" in html, f"Nav link to index.html missing in {name}"

    def test_print_css_present(self, populated_db: Path, tmp_path: Path):
        out = tmp_path / "export"
        export_calibration_report(populated_db, out)
        html = (out / "index.html").read_text()
        assert "@media print" in html
