"""Tests for the mkrefs module (truncated reference generator).

Tests cover:
- TestGenerateTruncatedRefs: pure-logic generation of truncated reference sequences
  (full construct, bc1+target, bc1-only, three levels, start_only mode)
- TestWriteManifest: TSV output with header and entries
"""

from __future__ import annotations

from pathlib import Path

import pytest

from barcodes import reverse_complement
from mkrefs import generate_truncated_refs, write_manifest

# ---------------------------------------------------------------------------
# Test constants
# ---------------------------------------------------------------------------

NB05 = "CACAAAGACACCGACAACTTTCTT"
NB10 = "GAGAGGACAAAGGTTTCAACGCTT"
FLANK_F = "AAGGTTAA"
FLANK_R = "CAGCACCT"
REV_FLANK_F = "AGGTGCTG"
REV_FLANK_R = "TTAACCTT"
TARGET = "ATCGATCG" * 10  # 80bp


# ---------------------------------------------------------------------------
# TestGenerateTruncatedRefs
# ---------------------------------------------------------------------------


class TestGenerateTruncatedRefs:
    """Validate generate_truncated_refs() pure-logic output."""

    def test_generates_full_construct(self):
        """Full level must concatenate all flanks, barcodes, and target."""
        refs = generate_truncated_refs(
            alias="CYP2D6_fwd",
            target_seq=TARGET,
            bc1_seq=NB05,
            bc2_seq=NB10,
            mask1_front=FLANK_F,
            mask1_rear=FLANK_R,
            mask2_front=REV_FLANK_F,
            mask2_rear=REV_FLANK_R,
        )
        expected = (
            FLANK_F + NB05 + FLANK_R
            + TARGET
            + REV_FLANK_F + reverse_complement(NB10) + REV_FLANK_R
        )
        assert refs["full"] == expected

    def test_generates_bc1_target(self):
        """bc1_target level must include mask1_front + BC1 + mask1_rear + target."""
        refs = generate_truncated_refs(
            alias="CYP2D6_fwd",
            target_seq=TARGET,
            bc1_seq=NB05,
            bc2_seq=NB10,
            mask1_front=FLANK_F,
            mask1_rear=FLANK_R,
            mask2_front=REV_FLANK_F,
            mask2_rear=REV_FLANK_R,
        )
        expected = FLANK_F + NB05 + FLANK_R + TARGET
        assert refs["bc1_target"] == expected

    def test_generates_bc1_only(self):
        """bc1_only level must include mask1_front + BC1 + mask1_rear."""
        refs = generate_truncated_refs(
            alias="CYP2D6_fwd",
            target_seq=TARGET,
            bc1_seq=NB05,
            bc2_seq=NB10,
            mask1_front=FLANK_F,
            mask1_rear=FLANK_R,
            mask2_front=REV_FLANK_F,
            mask2_rear=REV_FLANK_R,
        )
        expected = FLANK_F + NB05 + FLANK_R
        assert refs["bc1_only"] == expected

    def test_returns_three_levels(self):
        """Dual mode must return exactly {full, bc1_target, bc1_only}."""
        refs = generate_truncated_refs(
            alias="CYP2D6_fwd",
            target_seq=TARGET,
            bc1_seq=NB05,
            bc2_seq=NB10,
            mask1_front=FLANK_F,
            mask1_rear=FLANK_R,
            mask2_front=REV_FLANK_F,
            mask2_rear=REV_FLANK_R,
        )
        assert set(refs.keys()) == {"full", "bc1_target", "bc1_only"}

    def test_start_only_mode(self):
        """No bc2 must return only {bc1_target, bc1_only} (no full)."""
        refs = generate_truncated_refs(
            alias="GeneX_fwd",
            target_seq=TARGET,
            bc1_seq=NB05,
            bc2_seq=None,
            mask1_front=FLANK_F,
            mask1_rear=FLANK_R,
            mask2_front=None,
            mask2_rear=None,
        )
        assert set(refs.keys()) == {"bc1_target", "bc1_only"}
        assert refs["bc1_target"] == FLANK_F + NB05 + FLANK_R + TARGET
        assert refs["bc1_only"] == FLANK_F + NB05 + FLANK_R


# ---------------------------------------------------------------------------
# TestWriteManifest
# ---------------------------------------------------------------------------


class TestWriteManifest:
    """Validate write_manifest() TSV output."""

    def test_writes_tsv(self, tmp_path: Path):
        """Manifest must have header + 2 data rows with correct columns."""
        manifest_path = tmp_path / "manifest.tsv"
        entries = [
            {
                "alias": "CYP2D6_fwd",
                "level": "full",
                "length": 200,
                "path": "truncated/CYP2D6_fwd_full.fasta",
            },
            {
                "alias": "CYP2D6_fwd",
                "level": "bc1_target",
                "length": 120,
                "path": "truncated/CYP2D6_fwd_bc1_target.fasta",
            },
        ]
        write_manifest(entries, manifest_path)

        lines = manifest_path.read_text().strip().splitlines()
        assert len(lines) == 3  # header + 2 entries

        # Verify header
        header = lines[0].split("\t")
        assert header == ["alias", "level", "length", "path"]

        # Verify first data row
        row1 = lines[1].split("\t")
        assert row1[0] == "CYP2D6_fwd"
        assert row1[1] == "full"
        assert row1[2] == "200"
        assert row1[3] == "truncated/CYP2D6_fwd_full.fasta"

        # Verify second data row
        row2 = lines[2].split("\t")
        assert row2[0] == "CYP2D6_fwd"
        assert row2[1] == "bc1_target"
        assert row2[2] == "120"
        assert row2[3] == "truncated/CYP2D6_fwd_bc1_target.fasta"
