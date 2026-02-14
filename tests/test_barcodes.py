"""Tests for the barcodes module.

Tests cover:
- BARCODES dictionary: 96 entries, all 24bp, correct naming, valid DNA
- reverse_complement: simple, palindrome, empty, roundtrip
- classify_barcode: perfect match, different barcode, imperfect match, poor match, subset
"""

import re

import pytest

from barcodes import BARCODES, BARCODE_LENGTH, classify_barcode, reverse_complement


# ---------------------------------------------------------------------------
# BARCODES dictionary tests
# ---------------------------------------------------------------------------


class TestBarcodesDict:
    """Validate the BARCODES constant."""

    def test_barcode_count(self):
        """There must be exactly 96 barcodes."""
        assert len(BARCODES) == 96

    def test_barcode_names(self):
        """Keys must be nb01 through nb96."""
        expected = {f"nb{i:02d}" for i in range(1, 97)}
        assert set(BARCODES.keys()) == expected

    def test_barcode_length_constant(self):
        """BARCODE_LENGTH must be 24."""
        assert BARCODE_LENGTH == 24

    def test_all_barcodes_are_24bp(self):
        """Every barcode sequence must be exactly 24 bases."""
        for name, seq in BARCODES.items():
            assert len(seq) == 24, f"{name} has length {len(seq)}, expected 24"

    def test_all_barcodes_uppercase_acgt(self):
        """Every barcode must contain only uppercase A, C, G, T."""
        pattern = re.compile(r"^[ACGT]+$")
        for name, seq in BARCODES.items():
            assert pattern.match(seq), f"{name} contains invalid characters: {seq}"

    def test_all_barcodes_unique(self):
        """All barcode sequences must be unique."""
        seqs = list(BARCODES.values())
        assert len(seqs) == len(set(seqs)), "Duplicate barcode sequences found"

    def test_spot_check_nb01(self):
        """Spot-check nb01 sequence."""
        assert BARCODES["nb01"] == "CACAAAGACACCGACAACTTTCTT"

    def test_spot_check_nb96(self):
        """Spot-check nb96 sequence."""
        assert BARCODES["nb96"] == "CTGAACGGTCATAGAGTCCACCAT"

    def test_spot_check_nb48(self):
        """Spot-check nb48 (middle of the range)."""
        assert BARCODES["nb48"] == "CATCTGGAACGTGGTACACCTGTA"


# ---------------------------------------------------------------------------
# reverse_complement tests
# ---------------------------------------------------------------------------


class TestReverseComplement:
    """Validate the reverse_complement function."""

    def test_simple_sequence(self):
        """Reverse complement of a simple sequence."""
        assert reverse_complement("ACGT") == "ACGT"  # palindrome

    def test_known_sequence(self):
        """Reverse complement of a non-palindromic sequence."""
        assert reverse_complement("AAACCC") == "GGGTTT"

    def test_single_base(self):
        """Single base reverse complement."""
        assert reverse_complement("A") == "T"
        assert reverse_complement("C") == "G"
        assert reverse_complement("G") == "C"
        assert reverse_complement("T") == "A"

    def test_palindrome(self):
        """A palindromic sequence should be its own reverse complement."""
        palindrome = "AATT"
        assert reverse_complement(palindrome) == palindrome

    def test_empty_string(self):
        """Empty string should return empty string."""
        assert reverse_complement("") == ""

    def test_roundtrip(self):
        """Applying reverse_complement twice returns the original."""
        seq = "CACAAAGACACCGACAACTTTCTT"
        assert reverse_complement(reverse_complement(seq)) == seq

    def test_all_same_base(self):
        """All-A should become all-T (reversed)."""
        assert reverse_complement("AAAA") == "TTTT"

    def test_longer_sequence(self):
        """Verify a longer known reverse complement."""
        # GATTACA -> reverse: ACATTAG -> complement: TGTAATC
        assert reverse_complement("GATTACA") == "TGTAATC"


# ---------------------------------------------------------------------------
# classify_barcode tests
# ---------------------------------------------------------------------------


class TestClassifyBarcode:
    """Validate the classify_barcode function."""

    def test_perfect_match_nb01(self):
        """A perfect match to nb01 should return ed=0, confidence=1.0."""
        segment = BARCODES["nb01"]
        result = classify_barcode(segment, BARCODES)
        assert result["barcode_id"] == "nb01"
        assert result["edit_distance"] == 0
        assert result["confidence"] == pytest.approx(1.0)

    def test_perfect_match_nb50(self):
        """A perfect match to nb50 should return ed=0, confidence=1.0."""
        segment = BARCODES["nb50"]
        result = classify_barcode(segment, BARCODES)
        assert result["barcode_id"] == "nb50"
        assert result["edit_distance"] == 0
        assert result["confidence"] == pytest.approx(1.0)

    def test_different_barcode(self):
        """nb02 sequence should classify as nb02, not nb01."""
        segment = BARCODES["nb02"]
        result = classify_barcode(segment, BARCODES)
        assert result["barcode_id"] == "nb02"
        assert result["edit_distance"] == 0

    def test_imperfect_match_2_mutations(self):
        """Introduce 2 mutations into nb01; should still classify as nb01 with ed=2."""
        # nb01: CACAAAGACACCGACAACTTTCTT
        # Mutate positions 0 and 5: C->T, A->G
        mutated = "TACAAGGACACCGACAACTTTCTT"
        result = classify_barcode(mutated, BARCODES)
        assert result["barcode_id"] == "nb01"
        assert result["edit_distance"] == 2
        assert result["confidence"] == pytest.approx(1.0 - 2.0 / 24.0)

    def test_poor_match_random_sequence(self):
        """A random-ish sequence should produce a high edit distance and low confidence."""
        random_seq = "AAAAAAAAAAAAAAAAAAAAAAAAA"  # 25 A's -- unlikely to match well
        result = classify_barcode(random_seq, BARCODES)
        # Should still return *some* barcode, but with poor confidence
        assert result["barcode_id"] in BARCODES
        assert result["edit_distance"] > 5
        assert result["confidence"] < 0.8

    def test_subset_of_barcodes(self):
        """Classify against a subset containing only nb01, nb02, nb03."""
        subset = {k: BARCODES[k] for k in ["nb01", "nb02", "nb03"]}
        segment = BARCODES["nb02"]
        result = classify_barcode(segment, subset)
        assert result["barcode_id"] == "nb02"
        assert result["edit_distance"] == 0

    def test_subset_excludes_true_barcode(self):
        """When the true barcode is not in the subset, should match the closest one."""
        subset = {k: BARCODES[k] for k in ["nb01", "nb03"]}
        segment = BARCODES["nb02"]  # not in subset
        result = classify_barcode(segment, subset)
        assert result["barcode_id"] in subset
        assert result["edit_distance"] > 0

    def test_result_keys(self):
        """Result dict must have exactly these keys."""
        result = classify_barcode(BARCODES["nb01"], BARCODES)
        assert set(result.keys()) == {"barcode_id", "edit_distance", "confidence"}

    def test_confidence_bounds(self):
        """Confidence should be between 0 and 1 inclusive."""
        for bc_id in ["nb01", "nb50", "nb96"]:
            result = classify_barcode(BARCODES[bc_id], BARCODES)
            assert 0.0 <= result["confidence"] <= 1.0

    def test_hw_mode_longer_segment(self):
        """HW mode should find the barcode within a longer segment (flanking DNA)."""
        # Embed nb10 in a longer context
        flanking = "GGGGGG" + BARCODES["nb10"] + "TTTTTT"
        result = classify_barcode(flanking, BARCODES)
        assert result["barcode_id"] == "nb10"
        assert result["edit_distance"] == 0
        assert result["confidence"] == pytest.approx(1.0)
