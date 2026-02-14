"""Tests for sma_scan.py experiment scanner."""

from pathlib import Path
import json
import subprocess
import sys

import pytest

# Module under test lives in bin/ — add to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))

from sma_scan import parse_final_summary, parse_minknow_sample_sheet, parse_report_json


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
