"""Tests for the construct_wizard module.

Tests cover:
- TestDiscoverReferences: FASTA file discovery and length reporting
- TestFormatConstructDiagram: ASCII diagram generation for dual and start-only modes
- TestBuildConstructToml: TOML file generation and round-trip validation
- TestKnownKits: Verifies the KNOWN_KITS dictionary contents
"""

from __future__ import annotations

from pathlib import Path

import pytest

from construct_wizard import (
    KNOWN_KITS,
    build_construct_toml,
    discover_references,
    format_construct_diagram,
)
from construct import parse_construct_toml


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def ref_dir(tmp_path: Path) -> Path:
    """Create a temporary directory with FASTA files and a non-FASTA file."""
    d = tmp_path / "references"
    d.mkdir()
    (d / "target_A.fasta").write_text(">target_A\nATCGATCG\n")
    (d / "target_B.fasta").write_text(">target_B\nGCTAGCTA\n")
    (d / "README.txt").write_text("Not a FASTA\n")
    return d


# ---------------------------------------------------------------------------
# TestDiscoverReferences
# ---------------------------------------------------------------------------


class TestDiscoverReferences:
    """Validate discover_references() finds FASTA files correctly."""

    def test_finds_fasta_files(self, ref_dir: Path):
        """Must find exactly 2 FASTA files in the fixture directory."""
        refs = discover_references(ref_dir)
        names = {r["name"] for r in refs}
        assert names == {"target_A", "target_B"}
        assert len(refs) == 2

    def test_reports_lengths(self, ref_dir: Path):
        """Each discovered reference must have a length field."""
        refs = discover_references(ref_dir)
        for ref in refs:
            assert "length" in ref
            assert isinstance(ref["length"], int)
            assert ref["length"] > 0

    def test_ignores_non_fasta(self, ref_dir: Path):
        """README.txt must not appear in the discovered references."""
        refs = discover_references(ref_dir)
        names = {r["name"] for r in refs}
        assert "README" not in names


# ---------------------------------------------------------------------------
# TestFormatConstructDiagram
# ---------------------------------------------------------------------------


class TestFormatConstructDiagram:
    """Validate format_construct_diagram() ASCII output."""

    def test_dual_mode(self):
        """Dual mode diagram must contain BC1, RC(BC2), and TARGET."""
        diagram = format_construct_diagram(
            mode="dual_independent",
            mask1_front="AAGGTTAA",
            mask1_rear="CAGCACCT",
            mask2_front="AGGTGCTG",
            mask2_rear="TTAACCTT",
        )
        assert "BC1" in diagram
        assert "RC(BC2)" in diagram
        assert "TARGET" in diagram

    def test_start_only_mode(self):
        """Start-only mode diagram must contain BC1 and TARGET but not BC2."""
        diagram = format_construct_diagram(
            mode="start_only",
            mask1_front="AAGGTTAA",
            mask1_rear="CAGCACCT",
        )
        assert "BC1" in diagram
        assert "TARGET" in diagram
        assert "BC2" not in diagram


# ---------------------------------------------------------------------------
# TestBuildConstructToml
# ---------------------------------------------------------------------------


class TestBuildConstructToml:
    """Validate build_construct_toml() produces valid TOML files."""

    def test_produces_valid_toml(self, tmp_path: Path):
        """Generated TOML must contain [arrangement] and [sma] sections."""
        output = tmp_path / "construct.toml"
        build_construct_toml(
            output_path=output,
            name="SMA_test_dual",
            kit="SQK-NBD114-96",
            mode="dual_independent",
            mask1_front="AAGGTTAA",
            mask1_rear="CAGCACCT",
            targets=[
                {
                    "barcode1": "NB01",
                    "barcode2": "NB02",
                    "alias": "test_target",
                    "reference": "references/test.fasta",
                }
            ],
            mask2_front="AGGTGCTG",
            mask2_rear="TTAACCTT",
        )
        content = output.read_text()
        assert "[arrangement]" in content
        assert "[sma]" in content

    def test_roundtrips_through_parser(self, tmp_path: Path):
        """Generated TOML must parse successfully through parse_construct_toml()."""
        output = tmp_path / "construct.toml"
        build_construct_toml(
            output_path=output,
            name="SMA_roundtrip",
            kit="SQK-NBD114-96",
            mode="dual_independent",
            mask1_front="AAGGTTAA",
            mask1_rear="CAGCACCT",
            targets=[
                {
                    "barcode1": "NB05",
                    "barcode2": "NB10",
                    "alias": "CYP2D6_fwd",
                    "reference": "references/CYP2D6_fwd.fasta",
                }
            ],
            mask2_front="AGGTGCTG",
            mask2_rear="TTAACCTT",
        )
        cfg = parse_construct_toml(output)
        assert cfg.mode == "dual_independent"
        assert cfg.arrangement.name == "SMA_roundtrip"
        assert cfg.arrangement.kit == "SQK-NBD114-96"
        assert len(cfg.targets) == 1
        assert cfg.targets[0].alias == "CYP2D6_fwd"


# ---------------------------------------------------------------------------
# TestKnownKits
# ---------------------------------------------------------------------------


class TestKnownKits:
    """Validate the KNOWN_KITS dictionary."""

    def test_has_common_kits(self):
        """SQK-NBD114-96 must be present with count=96."""
        assert "SQK-NBD114-96" in KNOWN_KITS
        assert KNOWN_KITS["SQK-NBD114-96"]["count"] == 96

    def test_all_kits_have_required_fields(self):
        """Every kit entry must have count, length, and desc."""
        for kit_name, kit_info in KNOWN_KITS.items():
            assert "count" in kit_info, f"{kit_name} missing 'count'"
            assert "length" in kit_info, f"{kit_name} missing 'length'"
            assert "desc" in kit_info, f"{kit_name} missing 'desc'"
