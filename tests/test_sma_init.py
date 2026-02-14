"""Tests for sma_init.py database initialization and QC."""

from pathlib import Path
import json
import sqlite3
import subprocess
import sys

import pytest
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))

from sma_init import ingest_bam_to_db, compute_qc_metrics, generate_qc_html, create_database


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _make_bam(path: Path, reads: list[tuple[str, str]],
              pg_line: dict | None = None) -> None:
    """Create a BAM file with given (read_id, sequence) pairs."""
    pg = pg_line or {
        "ID": "basecaller", "PN": "dorado",
        "VN": "0.8.4", "CL": "dorado basecaller sup --no-trim",
    }
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "PG": [pg],
    })
    with pysam.AlignmentFile(str(path), "wb", header=header) as f:
        for read_id, seq in reads:
            seg = pysam.AlignedSegment(header)
            seg.query_name = read_id
            seg.query_sequence = seq
            seg.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            seg.flag = 4
            f.write(seg)


# ---------------------------------------------------------------------------
# ingest_bam_to_db tests
# ---------------------------------------------------------------------------


class TestIngestBamToDb:
    def test_creates_database_with_reads(self, tmp_path):
        bam = tmp_path / "merged.bam"
        _make_bam(bam, [
            ("read_1", "ACGTACGT" * 25),
            ("read_2", "TGCATGCA" * 30),
        ])
        db_path = tmp_path / "test.db"
        count = ingest_bam_to_db(bam, db_path, "test_exp", {})
        assert count == 2

        conn = sqlite3.connect(str(db_path))
        total = conn.execute("SELECT COUNT(*) FROM Reads").fetchone()[0]
        assert total == 2

        row = conn.execute(
            "SELECT read_id, readlen, q_bc FROM Reads WHERE read_id='read_1'"
        ).fetchone()
        assert row[0] == "read_1"
        assert row[1] == 200  # 8*25
        assert row[2] > 0  # some Q score
        conn.close()

    def test_attaches_end_reason_from_summary(self, tmp_path):
        bam = tmp_path / "merged.bam"
        _make_bam(bam, [("read_1", "ACGT" * 50)])
        db_path = tmp_path / "test.db"
        summary_map = {
            "read_1": {"end_reason": "signal_positive", "duration": 1.5}
        }
        ingest_bam_to_db(bam, db_path, "exp", summary_map)

        conn = sqlite3.connect(str(db_path))
        row = conn.execute(
            "SELECT ER, signal_duration_s FROM Reads WHERE read_id='read_1'"
        ).fetchone()
        assert row[0] == "signal_positive"
        assert row[1] == 1.5
        conn.close()

    def test_extracts_model_from_pg(self, tmp_path):
        bam = tmp_path / "merged.bam"
        _make_bam(bam, [("read_1", "ACGT" * 50)])
        db_path = tmp_path / "test.db"
        ingest_bam_to_db(bam, db_path, "exp", {})

        conn = sqlite3.connect(str(db_path))
        row = conn.execute(
            "SELECT model_tier, model_ver FROM Reads WHERE read_id='read_1'"
        ).fetchone()
        assert row[0] == "sup"
        assert row[1] == "0.8.4"
        conn.close()

    def test_creates_exp_record(self, tmp_path):
        bam = tmp_path / "merged.bam"
        _make_bam(bam, [("read_1", "ACGT")])
        db_path = tmp_path / "test.db"
        ingest_bam_to_db(bam, db_path, "my_exp", {})

        conn = sqlite3.connect(str(db_path))
        row = conn.execute("SELECT exp_id FROM Exp").fetchone()
        assert row[0] == "my_exp"
        conn.close()


# ---------------------------------------------------------------------------
# compute_qc_metrics tests
# ---------------------------------------------------------------------------


class TestComputeQcMetrics:
    def _make_db(self, tmp_path, reads):
        """Create a DB with given read data."""
        db_path = tmp_path / "test.db"
        bam = tmp_path / "test.bam"
        _make_bam(bam, [(r[0], r[1]) for r in reads])
        summary_map = {}
        for r in reads:
            if len(r) > 2:
                summary_map[r[0]] = r[2]
        ingest_bam_to_db(bam, db_path, "exp", summary_map)
        return db_path

    def test_computes_total_reads(self, tmp_path):
        db = self._make_db(tmp_path, [
            ("r1", "ACGT" * 50),
            ("r2", "TGCA" * 30),
        ])
        metrics = compute_qc_metrics(db)
        assert metrics["total_reads"] == 2

    def test_computes_length_stats(self, tmp_path):
        db = self._make_db(tmp_path, [
            ("r1", "ACGT" * 50),  # 200bp
            ("r2", "TGCA" * 100),  # 400bp
        ])
        metrics = compute_qc_metrics(db)
        assert len(metrics["read_lengths"]) == 2
        assert metrics["mean_length"] == 300.0
        assert metrics["total_bases"] == 600

    def test_computes_end_reason_counts(self, tmp_path):
        db = self._make_db(tmp_path, [
            ("r1", "ACGT" * 50, {"end_reason": "signal_positive"}),
            ("r2", "TGCA" * 30, {"end_reason": "signal_positive"}),
            ("r3", "AAAA" * 20, {"end_reason": "data_service_unblock_mux_change"}),
        ])
        metrics = compute_qc_metrics(db)
        assert metrics["end_reason_counts"]["signal_positive"] == 2
        assert metrics["end_reason_counts"]["data_service_unblock_mux_change"] == 1

    def test_computes_qscore_stats(self, tmp_path):
        db = self._make_db(tmp_path, [
            ("r1", "ACGT" * 50),
            ("r2", "TGCA" * 30),
        ])
        metrics = compute_qc_metrics(db)
        assert metrics["mean_qbc"] is not None
        assert metrics["mean_qbc"] > 0


# ---------------------------------------------------------------------------
# generate_qc_html tests
# ---------------------------------------------------------------------------


class TestGenerateQcHtml:
    def test_html_contains_tabs(self, tmp_path):
        db = TestComputeQcMetrics()._make_db(tmp_path, [
            ("r1", "ACGT" * 50, {"end_reason": "signal_positive", "duration": 1.0}),
            ("r2", "TGCA" * 30, {"end_reason": "unblock", "duration": 0.3}),
        ])
        metrics = compute_qc_metrics(db)
        manifest = {
            "experiment_id": "test_exp",
            "kit": "SQK-NBD114-24",
            "flow_cell_product": "FLO-MIN114",
            "software": {"basecaller_version": "0.8.4"},
            "flow_cells": [{"flow_cell_id": "FBD69411", "runs": [{}]}],
        }
        html = generate_qc_html(metrics, manifest)
        assert 'id="overview"' in html
        assert 'id="lengths"' in html
        assert 'id="quality"' in html
        assert 'id="endreasons"' in html
        assert 'id="duration"' in html
        assert "test_exp" in html
        assert len(html) > 1000


# ---------------------------------------------------------------------------
# CLI integration test
# ---------------------------------------------------------------------------


class TestInitCLI:
    def test_full_pipeline(self, tmp_path):
        """End-to-end: BAM → DB + QC report."""
        # Create test BAM
        bam = tmp_path / "merged.bam"
        _make_bam(bam, [
            ("read_1", "ACGTACGT" * 25),
            ("read_2", "TGCATGCA" * 30),
            ("read_3", "AAAACCCC" * 40),
        ])

        # Create summary
        summary = tmp_path / "summary.tsv"
        summary.write_text(
            "read_id\tend_reason\tduration\tmean_qscore_template\n"
            "read_1\tsignal_positive\t2.5\t15.0\n"
            "read_2\tdata_service_unblock_mux_change\t0.3\t10.0\n"
            "read_3\tsignal_positive\t3.0\t14.5\n"
        )

        # Create manifest
        manifest = {
            "experiment_id": "test_exp",
            "kit": "SQK-NBD114-24",
            "flow_cell_product": "FLO-MIN114",
            "software": {"basecaller_version": "0.8.4"},
            "flow_cells": [{"flow_cell_id": "FBD69411", "runs": [{}]}],
        }
        manifest_path = tmp_path / "manifest.json"
        manifest_path.write_text(json.dumps(manifest))

        out_dir = tmp_path / "output"
        result = subprocess.run(
            [sys.executable, "bin/sma_init.py",
             "--bam", str(bam),
             "--sequencing-summary", str(summary),
             "--manifest", str(manifest_path),
             "-o", str(out_dir)],
            capture_output=True, text=True,
            cwd="/tmp/ont-sma-seq",
        )
        assert result.returncode == 0, f"stderr: {result.stderr}\nstdout: {result.stdout}"

        # Check outputs
        db_path = out_dir / "test_exp.db"
        assert db_path.exists()
        report_path = out_dir / "test_exp_qc_report.html"
        assert report_path.exists()
        assert report_path.stat().st_size > 1000
        prov_path = out_dir / "test_exp_init.json"
        assert prov_path.exists()

        # Verify DB contents
        conn = sqlite3.connect(str(db_path))
        count = conn.execute("SELECT COUNT(*) FROM Reads").fetchone()[0]
        assert count == 3

        # Verify end reasons attached
        er = conn.execute(
            "SELECT ER FROM Reads WHERE read_id='read_1'"
        ).fetchone()[0]
        assert er == "signal_positive"
        conn.close()

        # Verify provenance
        prov = json.loads(prov_path.read_text())
        assert prov["total_reads"] == 3
