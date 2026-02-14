"""Tests for sequencing summary parsing."""
from __future__ import annotations
from pathlib import Path
import pytest
from calibrate.signal import load_sequencing_summary


@pytest.fixture
def summary_tsv(tmp_path: Path) -> Path:
    p = tmp_path / "sequencing_summary.txt"
    p.write_text(
        "read_id\tduration\tend_reason\tsequence_length_template\t"
        "mean_qscore_template\tbarcode_arrangement\n"
        "read-001\t1.5\tsignal_positive\t500\t12.3\tbarcode05\n"
        "read-002\t0.8\tdata_service_unblock_mux_change\t200\t9.1\tunclassified\n"
        "read-003\t2.1\tsignal_positive\t750\t15.0\tbarcode10\n"
    )
    return p


class TestLoadSequencingSummary:
    def test_returns_dict_keyed_by_read_id(self, summary_tsv):
        data = load_sequencing_summary(summary_tsv)
        assert "read-001" in data
        assert "read-002" in data
        assert "read-003" in data

    def test_contains_duration(self, summary_tsv):
        data = load_sequencing_summary(summary_tsv)
        assert data["read-001"]["duration"] == pytest.approx(1.5)

    def test_contains_end_reason(self, summary_tsv):
        data = load_sequencing_summary(summary_tsv)
        assert data["read-001"]["end_reason"] == "signal_positive"

    def test_contains_mean_qscore(self, summary_tsv):
        data = load_sequencing_summary(summary_tsv)
        assert data["read-001"]["mean_qscore"] == pytest.approx(12.3)

    def test_contains_sequence_length(self, summary_tsv):
        data = load_sequencing_summary(summary_tsv)
        assert data["read-001"]["sequence_length"] == 500

    def test_reads_count(self, summary_tsv):
        data = load_sequencing_summary(summary_tsv)
        assert len(data) == 3
