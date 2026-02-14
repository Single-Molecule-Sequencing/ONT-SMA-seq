#!/usr/bin/env python3
"""sma_scan.py — Scan experiment directory, discover MinKNOW runs, output manifest."""

import csv
import json as _json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


def parse_final_summary(path: Path) -> dict:
    """Parse a MinKNOW final_summary*.txt file (key=value format).

    Integer-valued keys (file counts) are converted to int.
    All other values remain strings. Missing keys return None via dict.get().
    """
    INT_KEYS = {
        "pod5_files_in_final_dest", "pod5_files_in_fallback",
        "bam_files_in_final_dest", "bam_files_in_fallback",
        "fast5_files_in_final_dest", "fastq_files_in_final_dest",
        "bai_files_in_final_dest",
    }
    data: dict = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if "=" not in line:
                continue
            key, _, val = line.partition("=")
            if key in INT_KEYS:
                data[key] = int(val)
            else:
                data[key] = val
    return data


def parse_minknow_sample_sheet(path: Path) -> dict:
    """Parse MinKNOW-format sample_sheet*.csv.

    Returns dict with experiment_id, kit, flow_cell_product_code,
    and barcode_aliases mapping.
    """
    result: dict = {
        "experiment_id": None,
        "kit": None,
        "flow_cell_product_code": None,
        "barcode_aliases": {},
    }
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if result["experiment_id"] is None:
                result["experiment_id"] = row.get("experiment_id") or None
                result["kit"] = row.get("kit") or None
                result["flow_cell_product_code"] = row.get("flow_cell_product_code") or None
            bc = row.get("barcode", "")
            alias = row.get("alias", "")
            if bc and alias:
                result["barcode_aliases"][bc] = alias
    return result


def parse_report_json(path: Path) -> dict:
    """Parse MinKNOW report*.json for software versions and device info.

    Navigates nested JSON safely, returning None for missing keys.
    """
    with open(path) as fh:
        data = _json.load(fh)
    pri = data.get("protocol_run_info", {})
    sw = pri.get("software_versions", {})
    fc = pri.get("flow_cell", {})
    dev = pri.get("device", {})
    mk = sw.get("minknow", {})
    return {
        "basecaller_version": sw.get("basecaller_build_version"),
        "bream_version": sw.get("bream"),
        "distribution_version": sw.get("distribution_version"),
        "minknow_version": mk.get("full") if isinstance(mk, dict) else None,
        "channel_count": fc.get("channel_count"),
        "device_type": dev.get("device_type"),
    }
