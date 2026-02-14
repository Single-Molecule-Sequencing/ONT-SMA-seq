"""Tests for the static HTML export engine."""

from __future__ import annotations

from pathlib import Path

import pytest

from viz.export import export_experiment


@pytest.fixture()
def experiment_dir(tmp_path: Path) -> Path:
    exp = tmp_path / "experiment"
    exp.mkdir()
    (exp / "references").mkdir()
    (exp / "sample_sheet.csv").write_text(
        "flow_cell_id,kit,barcode,alias,type\n"
        "FAL12345,SQK-NBD114-96,barcode05--barcode10,CYP2D6_fwd,test_sample\n"
    )
    (exp / "references" / "CYP2D6_fwd.fasta").write_text(">CYP2D6_fwd\nACGT\n")
    return exp


class TestExport:
    def test_creates_output_directory(self, experiment_dir: Path, tmp_path: Path) -> None:
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        assert output_dir.is_dir()

    def test_creates_all_html_files(self, experiment_dir: Path, tmp_path: Path) -> None:
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        for name in [
            "index.html",
            "sample_sheet.html",
            "barcodes.html",
            "construct.html",
            "targets.html",
            "assumptions.html",
            "validation_report.html",
        ]:
            assert (output_dir / name).exists(), f"Missing {name}"

    def test_html_is_self_contained(self, experiment_dir: Path, tmp_path: Path) -> None:
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        index = (output_dir / "index.html").read_text()
        assert "<style>" in index
        assert "CYP2D6_fwd" in index

    def test_validation_report(self, experiment_dir: Path, tmp_path: Path) -> None:
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        report = (output_dir / "validation_report.html").read_text()
        assert "Validation" in report

    def test_default_output_dir(self, experiment_dir: Path) -> None:
        """When no output_dir given, creates exports/export_<timestamp>/ in experiment_dir."""
        result = export_experiment(experiment_dir)
        result_path = Path(result)
        assert result_path.is_dir()
        assert result_path.parent.name == "exports"

    def test_sample_sheet_page_content(self, experiment_dir: Path, tmp_path: Path) -> None:
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        html = (output_dir / "sample_sheet.html").read_text()
        assert "CYP2D6_fwd" in html
        assert "barcode05--barcode10" in html
        assert "<style>" in html

    def test_targets_page_content(self, experiment_dir: Path, tmp_path: Path) -> None:
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        html = (output_dir / "targets.html").read_text()
        assert "CYP2D6_fwd" in html
        assert "ACGT" in html

    def test_nav_links_present(self, experiment_dir: Path, tmp_path: Path) -> None:
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        for name in [
            "index.html",
            "sample_sheet.html",
            "barcodes.html",
            "construct.html",
            "targets.html",
            "assumptions.html",
            "validation_report.html",
        ]:
            html = (output_dir / name).read_text()
            assert "index.html" in html, f"Nav link to index.html missing in {name}"

    def test_print_css_present(self, experiment_dir: Path, tmp_path: Path) -> None:
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        index = (output_dir / "index.html").read_text()
        assert "@media print" in index

    def test_barcodes_page(self, experiment_dir: Path, tmp_path: Path) -> None:
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        html = (output_dir / "barcodes.html").read_text()
        assert "Barcodes" in html

    def test_assumptions_page(self, experiment_dir: Path, tmp_path: Path) -> None:
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        html = (output_dir / "assumptions.html").read_text()
        assert "Assumptions" in html
        # Default assumptions should be present
        assert "strand_orientation" in html

    def test_construct_page(self, experiment_dir: Path, tmp_path: Path) -> None:
        output_dir = tmp_path / "export_output"
        export_experiment(experiment_dir, output_dir)
        html = (output_dir / "construct.html").read_text()
        assert "Construct" in html

    def test_returns_output_path(self, experiment_dir: Path, tmp_path: Path) -> None:
        output_dir = tmp_path / "export_output"
        result = export_experiment(experiment_dir, output_dir)
        assert str(output_dir) == str(result)
