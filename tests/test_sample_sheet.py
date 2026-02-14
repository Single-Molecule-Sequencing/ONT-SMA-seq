"""Tests for the sample_sheet module.

Tests cover:
- parse_sample_sheet: barcode pair extraction, alias mapping, expected count,
  error handling for missing columns and non-duplexed barcodes
- detect_barcode_ambiguity: no ambiguity case, ambiguity detected case
"""

from __future__ import annotations

from pathlib import Path

import pytest

from sample_sheet import detect_barcode_ambiguity, parse_sample_sheet


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def sample_sheet_csv(tmp_path: Path) -> Path:
    """Create a minimal MinKNOW-style sample sheet CSV with 3 rows.

    Rows:
      barcode05--barcode10  ->  CYP2D6_v04_fwd
      barcode10--barcode05  ->  CYP2D6_v04_rev
      barcode01--barcode02  ->  GeneX_fwd
    """
    csv_path = tmp_path / "sample_sheet.csv"
    csv_path.write_text(
        "flow_cell_id,kit,barcode,alias,type\n"
        "FC001,SQK-NBD114.96,barcode05--barcode10,CYP2D6_v04_fwd,test_sample\n"
        "FC001,SQK-NBD114.96,barcode10--barcode05,CYP2D6_v04_rev,test_sample\n"
        "FC001,SQK-NBD114.96,barcode01--barcode02,GeneX_fwd,test_sample\n"
    )
    return csv_path


@pytest.fixture()
def ambiguous_sheet(tmp_path: Path) -> Path:
    """Create a sample sheet where nb05 is upstream for two different targets.

    This is ambiguous because the same upstream barcode maps to different aliases.
    """
    csv_path = tmp_path / "ambiguous.csv"
    csv_path.write_text(
        "flow_cell_id,kit,barcode,alias,type\n"
        "FC001,SQK-NBD114.96,barcode05--barcode10,target_A,test_sample\n"
        "FC001,SQK-NBD114.96,barcode05--barcode12,target_B,test_sample\n"
    )
    return csv_path


# ---------------------------------------------------------------------------
# TestParseSampleSheet
# ---------------------------------------------------------------------------


class TestParseSampleSheet:
    """Validate parse_sample_sheet function."""

    def test_parses_barcode_pairs(self, sample_sheet_csv: Path):
        """All three expected barcode pairs must be present as keys."""
        mapping = parse_sample_sheet(sample_sheet_csv)
        assert ("nb05", "nb10") in mapping
        assert ("nb10", "nb05") in mapping
        assert ("nb01", "nb02") in mapping

    def test_alias_mapping(self, sample_sheet_csv: Path):
        """Each barcode pair must map to the correct alias string."""
        mapping = parse_sample_sheet(sample_sheet_csv)
        assert mapping[("nb05", "nb10")] == "CYP2D6_v04_fwd"
        assert mapping[("nb10", "nb05")] == "CYP2D6_v04_rev"
        assert mapping[("nb01", "nb02")] == "GeneX_fwd"

    def test_expected_barcodes(self, sample_sheet_csv: Path):
        """The fixture should produce exactly 3 entries."""
        mapping = parse_sample_sheet(sample_sheet_csv)
        assert len(mapping) == 3

    def test_missing_barcode_column(self, tmp_path: Path):
        """Missing 'barcode' column should raise ValueError mentioning 'barcode'."""
        bad_csv = tmp_path / "no_barcode.csv"
        bad_csv.write_text(
            "flow_cell_id,kit,alias,type\n"
            "FC001,SQK-NBD114.96,CYP2D6_v04_fwd,test_sample\n"
        )
        with pytest.raises(ValueError, match="barcode"):
            parse_sample_sheet(bad_csv)

    def test_invalid_barcode_format(self, tmp_path: Path):
        """Single barcode (not duplexed) should raise ValueError mentioning 'duplexed'."""
        bad_csv = tmp_path / "single_barcode.csv"
        bad_csv.write_text(
            "flow_cell_id,kit,barcode,alias,type\n"
            "FC001,SQK-NBD114.96,barcode05,CYP2D6_v04_fwd,test_sample\n"
        )
        with pytest.raises(ValueError, match="duplexed"):
            parse_sample_sheet(bad_csv)


# ---------------------------------------------------------------------------
# TestDetectBarcodeAmbiguity
# ---------------------------------------------------------------------------


class TestDetectBarcodeAmbiguity:
    """Validate detect_barcode_ambiguity function."""

    def test_no_ambiguity(self, sample_sheet_csv: Path):
        """Normal sample sheet with no shared upstream barcodes returns False."""
        mapping = parse_sample_sheet(sample_sheet_csv)
        assert detect_barcode_ambiguity(mapping) is False

    def test_ambiguity_detected(self, ambiguous_sheet: Path):
        """Sheet where nb05 is upstream for 2 different targets returns True."""
        mapping = parse_sample_sheet(ambiguous_sheet)
        assert detect_barcode_ambiguity(mapping) is True
