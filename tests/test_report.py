#!/usr/bin/env python3
"""Tests for the barcode classification report generator."""

import csv
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from report_analysis import (
    analyze_classification,
    classify_barcode_detailed,
    find_flank_position,
)

# ---------------------------------------------------------------------------
# Test data constants
# ---------------------------------------------------------------------------

NB05 = "AAGGTTACACAAACCCTGGACAAG"
NB10 = "GAGAGGACAAAGGTTTCAACGCTT"
FLANK_F = "AAGGTTAA"
FLANK_R = "CAGCACCT"
REV_FLANK_F = "AGGTGCTG"
REV_FLANK_R = "TTAACCTTAGCAAT"


def rc(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))


TARGET_SEQ = "ATCGATCGATCGATCGATCG" * 5  # 100bp


def make_full_length_read(bc_start=NB05, bc_end=NB10, target=TARGET_SEQ):
    """Build a full-length construct with both barcodes."""
    return FLANK_F + bc_start + FLANK_R + target + REV_FLANK_F + rc(bc_end) + REV_FLANK_R


def make_truncated_read(bc_start=NB05, target=TARGET_SEQ):
    """Build a truncated read with only the start barcode."""
    return FLANK_F + bc_start + FLANK_R + target


def make_db_row(read_id, seq, tgt_id="test_target",
                bc_start="nb05", bc_start_ed=0, bc_start_conf=1.0,
                bc_end="nb10", bc_end_ed=0, bc_end_conf=1.0,
                ed=5, q_ld=20.0):
    """Build a synthetic DB row dict."""
    return {
        "read_id": read_id,
        "readseq": seq,
        "readlen": len(seq),
        "tgt_id": tgt_id,
        "ed": ed,
        "q_bc": 15.0,
        "q_ld": q_ld,
        "ER": "signal_positive",
        "bc_start_id": bc_start,
        "bc_start_ed": bc_start_ed,
        "bc_start_conf": bc_start_conf,
        "bc_end_id": bc_end,
        "bc_end_ed": bc_end_ed,
        "bc_end_conf": bc_end_conf,
    }


# ---------------------------------------------------------------------------
# Unit tests: classify_barcode_detailed
# ---------------------------------------------------------------------------

class TestClassifyBarcodeDetailed:

    def test_returns_all_barcodes_sorted(self):
        segment = FLANK_F + NB05 + FLANK_R + "ATCGATCG" * 5
        expected = {"nb05": NB05, "nb10": NB10}
        results = classify_barcode_detailed(segment, expected)
        assert len(results) == 2
        assert results[0]["bc_id"] == "nb05"
        assert results[0]["edit_distance"] <= results[1]["edit_distance"]

    def test_perfect_match_has_zero_ed(self):
        segment = NB05 + "A" * 76
        expected = {"nb05": NB05}
        results = classify_barcode_detailed(segment, expected)
        assert results[0]["edit_distance"] == 0
        assert results[0]["confidence"] == 1.0

    def test_positions_returned(self):
        segment = "AAAAAA" + NB05 + "TTTTTT"
        expected = {"nb05": NB05}
        results = classify_barcode_detailed(segment, expected)
        assert results[0]["start"] == 6
        # end is inclusive (edlib convention) so barcode ends at position 29
        assert results[0]["end"] == 29


# ---------------------------------------------------------------------------
# Unit tests: find_flank_position
# ---------------------------------------------------------------------------

class TestFindFlankPosition:

    def test_finds_exact_flank(self):
        seq = "AAAA" + FLANK_F + "TTTT"
        result = find_flank_position(seq, FLANK_F, 0, len(seq))
        assert result is not None
        assert result["start"] == 4
        assert result["edit_distance"] == 0

    def test_returns_none_for_bad_match(self):
        seq = "GGGGGGGGGGGGGGGG"
        result = find_flank_position(seq, FLANK_F, 0, len(seq))
        # Should be None because no good match exists
        # (or very high ED depending on threshold)
        if result is not None:
            assert result["edit_distance"] > len(FLANK_F) // 2

    def test_respects_region_bounds(self):
        seq = "AAAA" + FLANK_F + "TTTT" + FLANK_F + "CCCC"
        # Only search second half
        mid = len(seq) // 2
        result = find_flank_position(seq, FLANK_F, mid, len(seq))
        assert result is not None
        assert result["start"] >= mid


# ---------------------------------------------------------------------------
# Unit tests: analyze_classification
# ---------------------------------------------------------------------------

class TestAnalyzeClassification:

    @pytest.fixture
    def sample_data(self):
        pair_to_alias = {("nb05", "nb10"): "test_target"}
        full_seq = make_full_length_read()
        trunc_seq = make_truncated_read()

        reads = [
            make_db_row("read_0001", full_seq, bc_end_conf=1.0),
            make_db_row("read_0002", full_seq, bc_end_conf=0.958),
            make_db_row("read_0003", trunc_seq, tgt_id="unmatched_nb05_nb09",
                        bc_end="nb09", bc_end_ed=10, bc_end_conf=0.583,
                        ed=None, q_ld=None),
            make_db_row("read_0004", trunc_seq, tgt_id="unmatched_nb05_nb09",
                        bc_end="nb09", bc_end_ed=10, bc_end_conf=0.583,
                        ed=None, q_ld=None),
            make_db_row("read_0005", trunc_seq, tgt_id="unmatched_nb05_nb09",
                        bc_end="nb09", bc_end_ed=10, bc_end_conf=0.583,
                        ed=None, q_ld=None),
        ]
        return reads, pair_to_alias

    def test_summary_counts(self, sample_data):
        reads, pairs = sample_data
        result = analyze_classification(reads, pairs)
        s = result["summary"]
        assert s["total"] == 5
        assert s["matched"] == 2
        assert s["unmatched"] == 3
        assert s["full_length"] == 2
        assert s["truncated"] == 3

    def test_per_target_stats(self, sample_data):
        reads, pairs = sample_data
        result = analyze_classification(reads, pairs)
        stats = result["per_target_stats"]
        assert "test_target" in stats
        assert stats["test_target"]["count"] == 2
        assert "unmatched" in stats
        assert stats["unmatched"]["count"] == 3

    def test_confidence_distributions(self, sample_data):
        reads, pairs = sample_data
        result = analyze_classification(reads, pairs)
        dists = result["confidence_distributions"]
        assert len(dists["start_confs"]) == 5
        assert len(dists["end_confs"]) == 5

    def test_pair_matrix(self, sample_data):
        reads, pairs = sample_data
        result = analyze_classification(reads, pairs)
        pm = result["pair_matrix"]
        assert ("nb05", "nb10") in pm["counts"]
        assert pm["counts"][("nb05", "nb10")] == 2
        assert ("nb05", "nb09") in pm["counts"]
        assert pm["counts"][("nb05", "nb09")] == 3
        assert "nb05" in pm["used_barcodes"]

    def test_reads_table_complete(self, sample_data):
        reads, pairs = sample_data
        result = analyze_classification(reads, pairs)
        table = result["reads_table"]
        assert len(table) == 5
        # Check structure of first row
        row = table[0]
        assert "read_id" in row
        assert "bc_start" in row
        assert "bc_end" in row
        assert "is_full_length" in row
        assert "pair" in row

    def test_detailed_reads_selected(self, sample_data):
        reads, pairs = sample_data
        result = analyze_classification(reads, pairs, max_detail_reads=4)
        detailed = result["detailed_reads"]
        assert len(detailed) >= 1
        assert len(detailed) <= 4
        # Should have both full-length and truncated representatives
        fl_count = sum(1 for d in detailed if d["is_full_length"])
        tr_count = sum(1 for d in detailed if not d["is_full_length"])
        assert fl_count >= 1
        assert tr_count >= 1

    def test_detailed_reads_have_competitors(self, sample_data):
        reads, pairs = sample_data
        result = analyze_classification(reads, pairs)
        for d in result["detailed_reads"]:
            assert "start_competitors" in d
            assert "end_competitors" in d
            assert len(d["start_competitors"]) > 0

    def test_barcode_info(self, sample_data):
        reads, pairs = sample_data
        result = analyze_classification(reads, pairs)
        info = result["barcode_info"]
        assert "nb05" in info
        assert "nb10" in info
        assert "fwd" in info["nb05"]
        assert "rc" in info["nb05"]

    def test_pairing_table(self, sample_data):
        reads, pairs = sample_data
        result = analyze_classification(reads, pairs)
        pt = result["pairing_table"]
        assert len(pt) == 1
        assert pt[0]["upstream"] == "nb05"
        assert pt[0]["downstream"] == "nb10"
        assert pt[0]["alias"] == "test_target"

    def test_flank_annotation_when_provided(self, sample_data):
        reads, pairs = sample_data
        result = analyze_classification(
            reads, pairs, flank_front=FLANK_F, flank_rear=FLANK_R
        )
        # Flanks should appear as region annotations in detailed reads
        detailed = result["detailed_reads"]
        assert len(detailed) >= 1
        # At least one detailed read should have flank regions
        all_labels = []
        for d in detailed:
            all_labels.extend(r["label"] for r in d.get("regions", []))
        assert any("flank" in label for label in all_labels)

    def test_references_passthrough(self, sample_data):
        reads, pairs = sample_data
        refs = {"test_target": (TARGET_SEQ, len(TARGET_SEQ))}
        result = analyze_classification(reads, pairs, references=refs)
        assert "references" in result


# ---------------------------------------------------------------------------
# Unit tests: generate_html
# ---------------------------------------------------------------------------

class TestGenerateHtml:

    def test_produces_valid_html(self, sample_data_for_html):
        from report_template import generate_html
        analysis, metadata = sample_data_for_html
        html = generate_html(analysis, metadata)
        assert html.startswith("<!DOCTYPE html>")
        assert "</html>" in html

    def test_contains_all_tabs(self, sample_data_for_html):
        from report_template import generate_html
        analysis, metadata = sample_data_for_html
        html = generate_html(analysis, metadata)
        assert 'id="tab-construct"' in html
        assert 'id="tab-summary"' in html
        assert 'id="tab-reads"' in html
        assert 'id="tab-details"' in html
        assert 'id="tab-barcodes"' in html

    def test_contains_experiment_metadata(self, sample_data_for_html):
        from report_template import generate_html
        analysis, metadata = sample_data_for_html
        html = generate_html(analysis, metadata)
        assert "TEST_EXP" in html

    def test_contains_javascript(self, sample_data_for_html):
        from report_template import generate_html
        analysis, metadata = sample_data_for_html
        html = generate_html(analysis, metadata)
        assert "function showTab" in html
        assert "function filterTable" in html
        assert "function sortTable" in html
        assert "function drawHist" in html

    @pytest.fixture
    def sample_data_for_html(self):
        pair_to_alias = {("nb05", "nb10"): "test_target"}
        full_seq = make_full_length_read()
        trunc_seq = make_truncated_read()
        reads = [
            make_db_row("read_0001", full_seq, bc_end_conf=1.0),
            make_db_row("read_0002", trunc_seq, tgt_id="unmatched_nb05_nb09",
                        bc_end="nb09", bc_end_ed=10, bc_end_conf=0.583,
                        ed=None, q_ld=None),
        ]
        analysis = analyze_classification(reads, pair_to_alias)
        metadata = {"exp_id": "TEST_EXP", "flow_cell_id": "FAL99999", "sample_id": "TEST"}
        return analysis, metadata


# ---------------------------------------------------------------------------
# Integration test: full pipeline → report
# ---------------------------------------------------------------------------

class TestReportIntegration:

    @pytest.fixture
    def full_pipeline_env(self, tmp_path):
        """Build DB via mkdb + ingest, then generate report."""
        target_seq = "ATCGATCGATCGATCGATCG" * 5

        # Build reads
        reads = []
        for i in range(3):
            seq = make_full_length_read()
            reads.append((f"read_{i:04d}", seq))
        for i in range(3, 5):
            seq = make_truncated_read()
            reads.append((f"read_{i:04d}", seq))

        # Write BAM
        bam_path = tmp_path / "FAL99999_20260214_TEST_sup_v5.2.0_trim0_0.bam"
        header = {"HD": {"VN": "1.6", "SO": "unsorted"}}
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            for read_id, seq in reads:
                a = pysam.AlignedSegment()
                a.query_name = read_id
                a.query_sequence = seq
                a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
                a.flag = 4
                bam.write(a)

        # Summary TSV
        summary = tmp_path / "summary.tsv"
        lines = ["read_id\tend_reason\n"]
        for read_id, _ in reads:
            lines.append(f"{read_id}\tsignal_positive\n")
        summary.write_text("".join(lines))

        # Sample sheet
        ss = tmp_path / "sample_sheet.csv"
        rows = [{"flow_cell_id": "FAL99999", "kit": "SQK-NBD114-96",
                 "sample_id": "20260214_TEST", "experiment_id": "exp_test",
                 "barcode": "barcode05--barcode10", "alias": "test_target"}]
        with open(ss, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

        # Reference FASTA
        ref_dir = tmp_path / "refs"
        ref_dir.mkdir()
        (ref_dir / "test_target.fasta").write_text(f">test_target\n{target_seq}\n")

        out_dir = tmp_path / "Output"

        # Step 1: mkdb
        r1 = subprocess.run(
            [sys.executable, "bin/mkdb.py",
             "-e", "FAL99999_20260214_TEST", "-o", str(out_dir)],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
        )
        assert r1.returncode == 0, r1.stderr

        db = out_dir / "SMA_FAL99999_20260214_TEST.db"

        # Step 2: ingest
        r2 = subprocess.run(
            [sys.executable, "bin/ingest.py",
             "-e", "FAL99999_20260214_TEST",
             "-b", str(bam_path), "-s", str(summary),
             "-d", str(db), "-o", str(tmp_path / "tagged.bam"),
             "-ss", str(ss), "-rd", str(ref_dir)],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
            env={**os.environ, "PYTHONPATH": "bin"},
        )
        assert r2.returncode == 0, f"{r2.stderr}\n{r2.stdout}"

        return {
            "db": db, "ss": ss, "out_dir": out_dir, "tmp": tmp_path,
        }

    def test_report_generates_html(self, full_pipeline_env):
        env = full_pipeline_env
        report_html = env["tmp"] / "report.html"

        r = subprocess.run(
            [sys.executable, "bin/report.py",
             "-d", str(env["db"]),
             "-ss", str(env["ss"]),
             "-o", str(report_html),
             "--flank-front", FLANK_F,
             "--flank-rear", FLANK_R],
            capture_output=True, text=True, cwd="/tmp/ont-sma-seq",
            env={**os.environ, "PYTHONPATH": "bin"},
        )
        assert r.returncode == 0, f"{r.stderr}\n{r.stdout}"
        assert report_html.exists()

        html = report_html.read_text()
        assert len(html) > 10_000  # at least 10KB

        # Verify key content
        assert "tab-construct" in html
        assert "tab-summary" in html
        assert "tab-reads" in html
        assert "tab-details" in html
        assert "tab-barcodes" in html
        assert "test_target" in html
        assert "nb05" in html
