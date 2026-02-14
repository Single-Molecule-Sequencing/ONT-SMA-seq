"""Tests for sma_scan.py experiment scanner."""

from pathlib import Path
import json
import subprocess
import sys

import pytest

# Module under test lives in bin/ — add to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))

from sma_scan import (
    parse_final_summary, parse_minknow_sample_sheet, parse_report_json,
    discover_runs, decide_merges,
)


class TestParseFinalSummary:
    def test_extracts_flow_cell_id(self, tmp_path):
        fs = tmp_path / "final_summary_FBD69411_abc_123.txt"
        fs.write_text(
            "instrument=MD-100098\n"
            "flow_cell_id=FBD69411\n"
            "protocol_group_id=my_experiment\n"
            "protocol=sequencing/sequencing_MIN114_DNA_e8_2_400K:FLO-MIN114:SQK-NBD114-24:400\n"
            "protocol_run_id=abc-def\n"
            "started=2025-12-29T10:56:12.341976-05:00\n"
            "acquisition_stopped=2026-01-01T10:56:12.413307-05:00\n"
            "pod5_files_in_final_dest=72\n"
            "bam_files_in_final_dest=2165\n"
        )
        result = parse_final_summary(fs)
        assert result["flow_cell_id"] == "FBD69411"
        assert result["instrument"] == "MD-100098"
        assert result["protocol_group_id"] == "my_experiment"
        assert result["protocol_run_id"] == "abc-def"
        assert result["started"] == "2025-12-29T10:56:12.341976-05:00"
        assert result["acquisition_stopped"] == "2026-01-01T10:56:12.413307-05:00"
        assert result["pod5_files_in_final_dest"] == 72
        assert result["bam_files_in_final_dest"] == 2165

    def test_missing_keys_return_none(self, tmp_path):
        fs = tmp_path / "final_summary_X_Y_Z.txt"
        fs.write_text("instrument=MD-100098\nflow_cell_id=FBD69411\n")
        result = parse_final_summary(fs)
        assert result["flow_cell_id"] == "FBD69411"
        assert result.get("started") is None

    def test_skips_blank_lines(self, tmp_path):
        fs = tmp_path / "final_summary.txt"
        fs.write_text("\n\ninstrument=MD-100098\n\nflow_cell_id=FBD69411\n\n")
        result = parse_final_summary(fs)
        assert result["flow_cell_id"] == "FBD69411"
        assert result["instrument"] == "MD-100098"


class TestParseMinkNowSampleSheet:
    def test_extracts_experiment_and_kit(self, tmp_path):
        ss = tmp_path / "sample_sheet_FBD69411_20251228_2101_abc.csv"
        ss.write_text(
            "protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,"
            "flow_cell_product_code,kit,barcode,alias,type,rowNumber\n"
            "abc-def,MD-100098,FBD69411,,my_experiment,FLO-MIN114,"
            "SQK-NBD114-24,barcode02,V04_2,test_sample,1\n"
            "abc-def,MD-100098,FBD69411,,my_experiment,FLO-MIN114,"
            "SQK-NBD114-24,barcode04,V04_4,test_sample,2\n"
        )
        result = parse_minknow_sample_sheet(ss)
        assert result["experiment_id"] == "my_experiment"
        assert result["kit"] == "SQK-NBD114-24"
        assert result["flow_cell_product_code"] == "FLO-MIN114"
        assert result["barcode_aliases"] == {
            "barcode02": "V04_2",
            "barcode04": "V04_4",
        }

    def test_empty_sample_sheet_returns_empty_aliases(self, tmp_path):
        ss = tmp_path / "sample_sheet.csv"
        ss.write_text(
            "protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,"
            "flow_cell_product_code,kit,barcode,alias,type,rowNumber\n"
        )
        result = parse_minknow_sample_sheet(ss)
        assert result["barcode_aliases"] == {}

    def test_handles_empty_barcode_field(self, tmp_path):
        ss = tmp_path / "sample_sheet.csv"
        ss.write_text(
            "protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,"
            "flow_cell_product_code,kit,barcode,alias,type,rowNumber\n"
            "abc-def,MD-100098,FBD69411,,my_experiment,FLO-MIN114,"
            "SQK-NBD114-24,,,test_sample,1\n"
        )
        result = parse_minknow_sample_sheet(ss)
        assert result["barcode_aliases"] == {}


class TestParseReportJson:
    def test_extracts_software_versions(self, tmp_path):
        report = tmp_path / "report_FBD69411_20251228_abc.json"
        report.write_text(json.dumps({
            "protocol_run_info": {
                "software_versions": {
                    "basecaller_build_version": "7.11.0+5d1db4a52",
                    "bream": "8.8.3",
                    "distribution_version": "25.09.16",
                    "minknow": {"full": "25.09.5"},
                },
                "flow_cell": {
                    "flow_cell_id": "FBD69411",
                    "channel_count": 512,
                    "product_code": "FLO-MIN114",
                },
                "device": {
                    "device_type": "MINION_MK1D",
                    "device_id": "MD-100098",
                },
            }
        }))
        result = parse_report_json(report)
        assert result["basecaller_version"] == "7.11.0+5d1db4a52"
        assert result["bream_version"] == "8.8.3"
        assert result["distribution_version"] == "25.09.16"
        assert result["channel_count"] == 512
        assert result["device_type"] == "MINION_MK1D"

    def test_missing_keys_return_none(self, tmp_path):
        report = tmp_path / "report.json"
        report.write_text(json.dumps({"protocol_run_info": {}}))
        result = parse_report_json(report)
        assert result["basecaller_version"] is None
        assert result["channel_count"] is None


# ---------------------------------------------------------------------------
# discover_runs tests
# ---------------------------------------------------------------------------


class TestDiscoverRuns:
    def _make_run_dir(self, base: Path, name: str, *,
                      flow_cell_id: str = "FBD69411",
                      instrument: str = "MD-100098",
                      run_id: str = "abc12345",
                      with_final_summary: bool = True,
                      with_pod5: bool = True,
                      with_bam: bool = True,
                      bam_demuxed: bool = False,
                      pod5_count: int = 3,
                      bam_count: int = 5) -> Path:
        """Helper to create a realistic MinKNOW run directory."""
        run_dir = base / name
        run_dir.mkdir(parents=True)
        if with_final_summary:
            fs = run_dir / f"final_summary_{flow_cell_id}_{run_id[:8]}_deadbeef.txt"
            fs.write_text(
                f"instrument={instrument}\n"
                f"flow_cell_id={flow_cell_id}\n"
                f"protocol_group_id=test_experiment\n"
                f"protocol=sequencing/seq:FLO-MIN114:SQK-NBD114-24:400\n"
                f"protocol_run_id={run_id}\n"
                f"started=2025-12-29T10:00:00-05:00\n"
                f"acquisition_stopped=2025-12-30T10:00:00-05:00\n"
                f"pod5_files_in_final_dest={pod5_count}\n"
                f"bam_files_in_final_dest={bam_count}\n"
            )
        if with_pod5:
            pod_dir = run_dir / "pod5_pass"
            pod_dir.mkdir()
            for i in range(pod5_count):
                (pod_dir / f"file_{i}.pod5").write_bytes(b"")
        if with_bam:
            bam_dir = run_dir / "bam_pass"
            bam_dir.mkdir()
            if bam_demuxed:
                for sub in ["unclassified", "V04_2", "V04_4"]:
                    (bam_dir / sub).mkdir()
                    (bam_dir / sub / "reads.bam").write_bytes(b"")
            else:
                for i in range(bam_count):
                    (bam_dir / f"reads_{i}.bam").write_bytes(b"")
        # Sample sheet
        ss = run_dir / f"sample_sheet_{flow_cell_id}_20251228_2101_{run_id[:8]}.csv"
        ss.write_text(
            "protocol_run_id,position_id,flow_cell_id,sample_id,experiment_id,"
            "flow_cell_product_code,kit,barcode,alias,type,rowNumber\n"
            f"{run_id},{instrument},{flow_cell_id},,test_experiment,FLO-MIN114,"
            "SQK-NBD114-24,barcode02,V04_2,test_sample,1\n"
        )
        return run_dir

    def test_finds_single_run(self, tmp_path):
        exp_dir = tmp_path / "experiment"
        exp_dir.mkdir()
        self._make_run_dir(
            exp_dir / "no_sample_id",
            "20251229_1055_MD-100098_FBD69411_abc12345",
        )
        runs = discover_runs(exp_dir)
        assert len(runs) == 1
        assert runs[0]["flow_cell_id"] == "FBD69411"

    def test_finds_multiple_runs(self, tmp_path):
        exp_dir = tmp_path / "experiment"
        exp_dir.mkdir()
        parent = exp_dir / "no_sample_id"
        self._make_run_dir(parent, "20251228_2219_MD-100098_FBD69411_run1xxxx",
                           run_id="run1xxxx")
        self._make_run_dir(parent, "20251229_1055_MD-100098_FBD69411_run2xxxx",
                           run_id="run2xxxx")
        runs = discover_runs(exp_dir)
        assert len(runs) == 2

    def test_detects_demuxed_bams(self, tmp_path):
        exp_dir = tmp_path / "experiment"
        exp_dir.mkdir()
        self._make_run_dir(
            exp_dir / "no_sample_id",
            "20251229_1055_MD-100098_FBD69411_abc12345",
            bam_demuxed=True,
        )
        runs = discover_runs(exp_dir)
        assert runs[0]["bam_demuxed"] is True

    def test_run_without_final_summary(self, tmp_path):
        exp_dir = tmp_path / "experiment"
        exp_dir.mkdir()
        self._make_run_dir(
            exp_dir / "no_sample_id",
            "20251228_2219_MD-100098_FBD69411_abc12345",
            with_final_summary=False,
        )
        runs = discover_runs(exp_dir)
        assert len(runs) == 1
        assert runs[0]["final_summary"] is None
        # Should still get flow_cell_id from dir name
        assert runs[0]["flow_cell_id"] == "FBD69411"

    def test_extracts_barcode_aliases(self, tmp_path):
        exp_dir = tmp_path / "experiment"
        exp_dir.mkdir()
        self._make_run_dir(
            exp_dir / "no_sample_id",
            "20251229_1055_MD-100098_FBD69411_abc12345",
        )
        runs = discover_runs(exp_dir)
        assert runs[0]["barcode_aliases"] == {"barcode02": "V04_2"}


# ---------------------------------------------------------------------------
# decide_merges tests
# ---------------------------------------------------------------------------


class TestDecideMerges:
    def test_auto_merge_same_flow_cell(self):
        runs = [
            {"flow_cell_id": "FBD69411", "instrument": "MD-100098",
             "experiment_id": "exp1", "kit": "SQK-NBD114-24",
             "protocol": "seq:FLO:SQK:400",
             "started": "2025-12-28T22:00:00-05:00",
             "stopped": "2025-12-29T08:00:00-05:00"},
            {"flow_cell_id": "FBD69411", "instrument": "MD-100098",
             "experiment_id": "exp1", "kit": "SQK-NBD114-24",
             "protocol": "seq:FLO:SQK:400",
             "started": "2025-12-29T10:00:00-05:00",
             "stopped": "2025-12-30T10:00:00-05:00"},
        ]
        groups = decide_merges(runs, max_gap_hours=24)
        assert len(groups) == 1
        assert groups[0]["merge_decision"] == "auto"
        assert len(groups[0]["runs"]) == 2
        assert groups[0]["time_gap_hours"] == pytest.approx(2.0, abs=0.1)

    def test_flag_when_gap_exceeds_threshold(self):
        runs = [
            {"flow_cell_id": "FBD69411", "instrument": "MD-100098",
             "experiment_id": "exp1", "kit": "SQK-NBD114-24",
             "protocol": "seq:FLO:SQK:400",
             "started": "2025-12-28T10:00:00-05:00",
             "stopped": "2025-12-28T20:00:00-05:00"},
            {"flow_cell_id": "FBD69411", "instrument": "MD-100098",
             "experiment_id": "exp1", "kit": "SQK-NBD114-24",
             "protocol": "seq:FLO:SQK:400",
             "started": "2025-12-31T10:00:00-05:00",
             "stopped": None},
        ]
        groups = decide_merges(runs, max_gap_hours=24)
        assert groups[0]["merge_decision"] == "review"
        assert "gap" in groups[0]["warnings"][0].lower()

    def test_separate_flow_cells(self):
        runs = [
            {"flow_cell_id": "FBD69411", "instrument": "MD-100098",
             "experiment_id": "exp1", "kit": "K1", "protocol": "P1",
             "started": "2025-12-28T10:00:00", "stopped": None},
            {"flow_cell_id": "FBD70599", "instrument": "MD-101527",
             "experiment_id": "exp1", "kit": "K1", "protocol": "P1",
             "started": "2025-12-28T10:00:00", "stopped": None},
        ]
        groups = decide_merges(runs, max_gap_hours=24)
        assert len(groups) == 2

    def test_single_run_no_merge(self):
        runs = [
            {"flow_cell_id": "FBD66244", "instrument": "MD-101527",
             "experiment_id": "exp1", "kit": "K1", "protocol": "P1",
             "started": "2025-12-30T17:10:00", "stopped": None},
        ]
        groups = decide_merges(runs, max_gap_hours=24)
        assert len(groups) == 1
        assert groups[0]["merge_decision"] == "single"
        assert len(groups[0]["runs"]) == 1

    def test_flag_different_experiment_id(self):
        runs = [
            {"flow_cell_id": "FBD69411", "instrument": "MD-100098",
             "experiment_id": "exp1", "kit": "K1", "protocol": "P1",
             "started": "2025-12-28T10:00:00", "stopped": "2025-12-28T20:00:00"},
            {"flow_cell_id": "FBD69411", "instrument": "MD-100098",
             "experiment_id": "exp2", "kit": "K1", "protocol": "P1",
             "started": "2025-12-28T21:00:00", "stopped": None},
        ]
        groups = decide_merges(runs, max_gap_hours=24)
        assert groups[0]["merge_decision"] == "review"
        assert any("experiment_id" in w for w in groups[0]["warnings"])
