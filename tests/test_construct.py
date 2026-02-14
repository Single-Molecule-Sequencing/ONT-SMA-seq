"""Tests for the construct module.

Tests cover:
- TestParseConstructToml: parsing dual mode, start_only mode, arrangement flanks,
  scoring defaults, confidence defaults, truncation defaults, barcode_pair_to_alias
  mapping, used_barcodes set, flank_sequences
- TestConstructValidation: rejects missing arrangement, missing sma, invalid mode,
  dual target without barcode2, invalid barcode ID, non-DNA flanks
"""

from __future__ import annotations

from pathlib import Path

import pytest

from construct import ConstructConfig, ValidationError, parse_construct_toml


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

DUAL_TOML = """\
[arrangement]
name = "SMA_CYP2D6_dual"
kit = "SQK-NBD114-96"
mask1_front = "AAGGTTAA"
mask1_rear  = "CAGCACCT"
mask2_front = "AGGTGCTG"
mask2_rear  = "TTAACCTTAGCAAT"
barcode1_pattern = "NB%02i"
barcode2_pattern = "NB%02i"
first_index = 1
last_index = 96

[scoring]
max_barcode_penalty = 11
min_barcode_penalty_dist = 3
front_barcode_window = 100
rear_barcode_window = 100
min_flank_score = 0.5

[sma]
mode = "dual_independent"

[sma.confidence]
full_length_threshold = 0.75
start_barcode_min = 0.6
flank_max_error_rate = 0.5

[sma.truncation]
auto_generate_refs = true
min_target_length = 20

[[sma.targets]]
barcode1 = "NB05"
barcode2 = "NB10"
alias = "CYP2D6_v04_fwd"
reference = "references/CYP2D6_v04_fwd.fasta"

[[sma.targets]]
barcode1 = "NB10"
barcode2 = "NB05"
alias = "CYP2D6_v04_rev"
reference = "references/CYP2D6_v04_rev.fasta"
"""

START_ONLY_TOML = """\
[arrangement]
name = "SMA_simple"
kit = "SQK-NBD114-96"
mask1_front = "AAGGTTAA"
mask1_rear  = "CAGCACCT"
barcode1_pattern = "NB%02i"
first_index = 1
last_index = 96

[sma]
mode = "start_only"

[[sma.targets]]
barcode1 = "NB01"
alias = "GeneX_fwd"
reference = "references/GeneX_fwd.fasta"
"""

MINIMAL_DUAL_TOML = """\
[arrangement]
name = "SMA_minimal"
kit = "SQK-NBD114-96"
mask1_front = "AAGGTTAA"
mask1_rear  = "CAGCACCT"
mask2_front = "AGGTGCTG"
mask2_rear  = "TTAACCTT"
barcode1_pattern = "NB%02i"
barcode2_pattern = "NB%02i"
first_index = 1
last_index = 96

[sma]
mode = "dual_independent"

[[sma.targets]]
barcode1 = "NB01"
barcode2 = "NB02"
alias = "target_A"
reference = "references/target_A.fasta"
"""


@pytest.fixture()
def dual_toml_path(tmp_path: Path) -> Path:
    """Write a full dual-mode construct TOML and return the path."""
    p = tmp_path / "dual_construct.toml"
    p.write_text(DUAL_TOML)
    return p


@pytest.fixture()
def start_only_toml_path(tmp_path: Path) -> Path:
    """Write a start_only-mode construct TOML and return the path."""
    p = tmp_path / "start_only_construct.toml"
    p.write_text(START_ONLY_TOML)
    return p


@pytest.fixture()
def minimal_dual_toml_path(tmp_path: Path) -> Path:
    """Write a minimal dual-mode construct TOML and return the path."""
    p = tmp_path / "minimal_dual.toml"
    p.write_text(MINIMAL_DUAL_TOML)
    return p


# ---------------------------------------------------------------------------
# TestParseConstructToml
# ---------------------------------------------------------------------------


class TestParseConstructToml:
    """Validate parse_construct_toml and ConstructConfig for valid inputs."""

    def test_parses_dual_mode(self, dual_toml_path: Path):
        """Dual-mode TOML must parse with mode='dual_independent'."""
        cfg = parse_construct_toml(dual_toml_path)
        assert cfg.mode == "dual_independent"

    def test_parses_start_only_mode(self, start_only_toml_path: Path):
        """Start-only TOML must parse with mode='start_only'."""
        cfg = parse_construct_toml(start_only_toml_path)
        assert cfg.mode == "start_only"

    def test_arrangement_flanks(self, dual_toml_path: Path):
        """Arrangement section flanking sequences must be parsed."""
        cfg = parse_construct_toml(dual_toml_path)
        assert cfg.arrangement.mask1_front == "AAGGTTAA"
        assert cfg.arrangement.mask1_rear == "CAGCACCT"
        assert cfg.arrangement.mask2_front == "AGGTGCTG"
        assert cfg.arrangement.mask2_rear == "TTAACCTTAGCAAT"

    def test_scoring_defaults(self, start_only_toml_path: Path):
        """When [scoring] is omitted, defaults must be populated."""
        cfg = parse_construct_toml(start_only_toml_path)
        assert cfg.scoring.max_barcode_penalty == 11
        assert cfg.scoring.min_barcode_penalty_dist == 3
        assert cfg.scoring.front_barcode_window == 100
        assert cfg.scoring.rear_barcode_window == 100
        assert cfg.scoring.min_flank_score == pytest.approx(0.5)

    def test_confidence_defaults(self, start_only_toml_path: Path):
        """When [sma.confidence] is omitted, defaults must be populated."""
        cfg = parse_construct_toml(start_only_toml_path)
        assert cfg.confidence.full_length_threshold == pytest.approx(0.75)
        assert cfg.confidence.start_barcode_min == pytest.approx(0.6)
        assert cfg.confidence.flank_max_error_rate == pytest.approx(0.5)

    def test_truncation_defaults(self, start_only_toml_path: Path):
        """When [sma.truncation] is omitted, defaults must be populated."""
        cfg = parse_construct_toml(start_only_toml_path)
        assert cfg.truncation.auto_generate_refs is True
        assert cfg.truncation.min_target_length == 20

    def test_barcode_pair_to_alias_dual(self, dual_toml_path: Path):
        """barcode_pair_to_alias() must return dict compatible with sample_sheet."""
        cfg = parse_construct_toml(dual_toml_path)
        mapping = cfg.barcode_pair_to_alias()
        assert mapping[("nb05", "nb10")] == "CYP2D6_v04_fwd"
        assert mapping[("nb10", "nb05")] == "CYP2D6_v04_rev"
        assert len(mapping) == 2

    def test_barcode_pair_to_alias_start_only(self, start_only_toml_path: Path):
        """start_only mode: barcode_pair_to_alias() uses barcode1 for both elements."""
        cfg = parse_construct_toml(start_only_toml_path)
        mapping = cfg.barcode_pair_to_alias()
        # In start_only mode, we use (barcode1, barcode1) as the key
        assert mapping[("nb01", "nb01")] == "GeneX_fwd"

    def test_used_barcode_ids(self, dual_toml_path: Path):
        """used_barcode_ids() must return all barcode IDs referenced in targets."""
        cfg = parse_construct_toml(dual_toml_path)
        used = cfg.used_barcode_ids()
        assert used == {"nb05", "nb10"}

    def test_flank_properties(self, dual_toml_path: Path):
        """flank_front and flank_rear must return mask1_rear and mask2_front."""
        cfg = parse_construct_toml(dual_toml_path)
        # flank_front = the trailing flank after barcode1 = mask1_rear
        assert cfg.flank_front == "CAGCACCT"
        # flank_rear = the leading flank before RC(barcode2) = mask2_front
        assert cfg.flank_rear == "AGGTGCTG"


# ---------------------------------------------------------------------------
# TestConstructValidation
# ---------------------------------------------------------------------------


class TestConstructValidation:
    """Validate that parse_construct_toml rejects invalid inputs."""

    def test_missing_arrangement(self, tmp_path: Path):
        """TOML without [arrangement] must raise ValidationError."""
        p = tmp_path / "no_arrangement.toml"
        p.write_text(
            '[sma]\nmode = "dual_independent"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\nbarcode2 = "NB02"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="arrangement"):
            parse_construct_toml(p)

    def test_missing_sma(self, tmp_path: Path):
        """TOML without [sma] must raise ValidationError."""
        p = tmp_path / "no_sma.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'barcode1_pattern = "NB%02i"\nfirst_index = 1\nlast_index = 96\n'
        )
        with pytest.raises(ValidationError, match="sma"):
            parse_construct_toml(p)

    def test_invalid_mode(self, tmp_path: Path):
        """mode must be 'dual_independent' or 'start_only'."""
        p = tmp_path / "bad_mode.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'barcode1_pattern = "NB%02i"\nfirst_index = 1\nlast_index = 96\n'
            '[sma]\nmode = "single_end"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="mode"):
            parse_construct_toml(p)

    def test_dual_target_without_barcode2(self, tmp_path: Path):
        """In dual_independent mode, targets must have barcode2."""
        p = tmp_path / "missing_bc2.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'mask2_front = "AGGTGCTG"\nmask2_rear = "TTAACCTT"\n'
            'barcode1_pattern = "NB%02i"\nbarcode2_pattern = "NB%02i"\n'
            "first_index = 1\nlast_index = 96\n"
            '[sma]\nmode = "dual_independent"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="barcode2"):
            parse_construct_toml(p)

    def test_invalid_barcode_id(self, tmp_path: Path):
        """Barcode ID outside NB01-NB96 must raise ValidationError."""
        p = tmp_path / "bad_bc.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'barcode1_pattern = "NB%02i"\nfirst_index = 1\nlast_index = 96\n'
            '[sma]\nmode = "start_only"\n'
            '[[sma.targets]]\nbarcode1 = "NB99"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="NB99"):
            parse_construct_toml(p)

    def test_non_dna_flanks(self, tmp_path: Path):
        """Flanking sequences with non-DNA characters must raise ValidationError."""
        p = tmp_path / "bad_flank.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAXGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'barcode1_pattern = "NB%02i"\nfirst_index = 1\nlast_index = 96\n'
            '[sma]\nmode = "start_only"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="DNA"):
            parse_construct_toml(p)

    def test_dual_mode_without_mask2_fails(self, tmp_path: Path):
        """dual_independent mode without mask2_front/mask2_rear must fail."""
        p = tmp_path / "dual_no_mask2.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'barcode1_pattern = "NB%02i"\nbarcode2_pattern = "NB%02i"\n'
            "first_index = 1\nlast_index = 96\n"
            '[sma]\nmode = "dual_independent"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\nbarcode2 = "NB02"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="mask2"):
            parse_construct_toml(p)

    def test_dual_mode_without_barcode2_pattern_fails(self, tmp_path: Path):
        """dual_independent mode without barcode2_pattern must fail."""
        p = tmp_path / "dual_no_bc2_pattern.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'mask2_front = "AGGTGCTG"\nmask2_rear = "TTAACCTT"\n'
            'barcode1_pattern = "NB%02i"\n'
            "first_index = 1\nlast_index = 96\n"
            '[sma]\nmode = "dual_independent"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\nbarcode2 = "NB02"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="barcode2_pattern"):
            parse_construct_toml(p)

    def test_start_only_with_mask2_fails(self, tmp_path: Path):
        """start_only mode with mask2_front or mask2_rear must fail."""
        p = tmp_path / "start_with_mask2.toml"
        p.write_text(
            "[arrangement]\n"
            'name = "test"\nkit = "SQK-NBD114-96"\n'
            'mask1_front = "AAGGTTAA"\nmask1_rear = "CAGCACCT"\n'
            'mask2_front = "AGGTGCTG"\nmask2_rear = "TTAACCTT"\n'
            'barcode1_pattern = "NB%02i"\nfirst_index = 1\nlast_index = 96\n'
            '[sma]\nmode = "start_only"\n'
            '[[sma.targets]]\nbarcode1 = "NB01"\n'
            'alias = "x"\nreference = "x.fasta"\n'
        )
        with pytest.raises(ValidationError, match="should not have mask2"):
            parse_construct_toml(p)
