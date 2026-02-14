"""Tests for viz.models Pydantic models.

Test classes:
- TestExperimentConfig: minimal config, full roundtrip
- TestTruncationRules: defaults, valid rules
- TestSampleSheetEntry: duplexed barcode computed fields, single barcode rejection
- TestArrangementConfig: default scoring params
- TestTargetRef: GC content calculation
"""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from viz.models import (
    ArrangementConfig,
    Assumption,
    BarcodeAssignment,
    ClassificationConfig,
    ConstructConfig,
    DemultiplexingConfig,
    DemuxPairEntry,
    ExperimentConfig,
    QualityMetrics,
    SampleSheetEntry,
    ScoringParams,
    TargetRef,
    TruncationConfig,
    TruncationRules,
)


# ---------------------------------------------------------------------------
# TestExperimentConfig
# ---------------------------------------------------------------------------


class TestExperimentConfig:
    """Validate the top-level ExperimentConfig model."""

    def test_minimal_config(self):
        """A minimal ExperimentConfig should populate default assumptions."""
        cfg = ExperimentConfig(
            description="Minimal test experiment",
            construct=ConstructConfig(
                adapter_5prime="ATCG",
                adapter_3prime="GCAT",
                insert_type="amplicon",
            ),
            demultiplexing=DemultiplexingConfig(
                mode="dual",
                start_barcode="nb01",
                end_barcode="nb02",
                pairs=[],
            ),
            classification=ClassificationConfig(
                barcode_search_window=100,
                confidence_formula="1 - penalty / max_penalty",
                ambiguity_triggers_full_construct=True,
            ),
            quality=QualityMetrics(
                q_bc="Q_bc = -10 * log10(1 - confidence)",
                q_ld="Q_ld = -10 * log10(edit_distance / length)",
            ),
            truncation=TruncationConfig(
                min_barcode_confidence=0.6,
                min_target_fraction=0.8,
                adapter_search_window=150,
                rules=TruncationRules(),
            ),
        )
        # Should have 4 default assumptions
        assert len(cfg.assumptions) == 4
        assert cfg.description == "Minimal test experiment"

    def test_full_config_roundtrip(self):
        """An ExperimentConfig should survive JSON serialization round-trip."""
        cfg = ExperimentConfig(
            description="Roundtrip test",
            construct=ConstructConfig(
                adapter_5prime="AAAA",
                adapter_3prime="TTTT",
                insert_type="genomic",
                notes="Test notes",
            ),
            demultiplexing=DemultiplexingConfig(
                mode="start_only",
                start_barcode="nb05",
                end_barcode="nb10",
                pairs=[
                    DemuxPairEntry(start="nb05", end="nb10", alias="target_A"),
                ],
            ),
            classification=ClassificationConfig(
                barcode_search_window=200,
                confidence_formula="custom_formula",
                ambiguity_triggers_full_construct=False,
            ),
            quality=QualityMetrics(
                q_bc="formula_bc",
                q_ld="formula_ld",
            ),
            truncation=TruncationConfig(
                min_barcode_confidence=0.7,
                min_target_fraction=0.9,
                adapter_search_window=200,
                rules=TruncationRules(),
            ),
            assumptions=[
                Assumption(key="a1", text="Test assumption", why="because"),
            ],
        )
        json_str = cfg.model_dump_json()
        restored = ExperimentConfig.model_validate_json(json_str)
        assert restored.description == "Roundtrip test"
        assert restored.construct.insert_type == "genomic"
        assert len(restored.demultiplexing.pairs) == 1
        assert restored.demultiplexing.pairs[0].alias == "target_A"
        assert len(restored.assumptions) == 1


# ---------------------------------------------------------------------------
# TestTruncationRules
# ---------------------------------------------------------------------------


class TestTruncationRules:
    """Validate TruncationRules model defaults and custom values."""

    def test_defaults(self):
        """Default TruncationRules should have non-empty string values."""
        rules = TruncationRules()
        assert isinstance(rules.full, str)
        assert len(rules.full) > 0
        assert isinstance(rules.trunc_3prime, str)
        assert isinstance(rules.trunc_target, str)
        assert isinstance(rules.trunc_barcode, str)
        assert isinstance(rules.adapter_only, str)
        assert isinstance(rules.chimeric, str)

    def test_valid_rules(self):
        """Custom TruncationRules should accept all string fields."""
        rules = TruncationRules(
            full="both barcodes + target",
            trunc_3prime="missing 3' barcode",
            trunc_target="partial target",
            trunc_barcode="partial barcode",
            adapter_only="adapter only, no insert",
            chimeric="chimeric construct detected",
        )
        assert rules.full == "both barcodes + target"
        assert rules.chimeric == "chimeric construct detected"


# ---------------------------------------------------------------------------
# TestSampleSheetEntry
# ---------------------------------------------------------------------------


class TestSampleSheetEntry:
    """Validate SampleSheetEntry barcode parsing and computed fields."""

    def test_duplexed_barcode(self):
        """Duplexed barcode should parse into upstream and downstream nb-format."""
        entry = SampleSheetEntry(
            flow_cell_id="FC001",
            kit="SQK-NBD114.96",
            barcode="barcode05--barcode10",
            alias="CYP2D6_v04_fwd",
            experiment_id="EXP001",
            sample_id="SAMPLE001",
            type="test_sample",
        )
        assert entry.upstream_barcode == "nb05"
        assert entry.downstream_barcode == "nb10"

    def test_single_barcode_rejected(self):
        """A single (non-duplexed) barcode should be rejected by validation."""
        with pytest.raises(ValidationError, match="duplexed"):
            SampleSheetEntry(
                flow_cell_id="FC001",
                kit="SQK-NBD114.96",
                barcode="barcode05",
                alias="CYP2D6_v04_fwd",
                experiment_id="EXP001",
                sample_id="SAMPLE001",
                type="test_sample",
            )


# ---------------------------------------------------------------------------
# TestArrangementConfig
# ---------------------------------------------------------------------------


class TestArrangementConfig:
    """Validate ArrangementConfig and its nested ScoringParams defaults."""

    def test_defaults(self):
        """Default ArrangementConfig should have expected ScoringParams values."""
        arr = ArrangementConfig(name="test_arrangement")
        assert arr.scoring.max_barcode_penalty == 11
        assert arr.scoring.min_barcode_penalty_dist == 3
        assert arr.scoring.min_separation_only_dist == 6
        assert arr.scoring.barcode_end_proximity == 75
        assert arr.scoring.flank_left_pad == 5
        assert arr.scoring.flank_right_pad == 10
        assert arr.scoring.front_barcode_window == 175
        assert arr.scoring.rear_barcode_window == 175
        assert arr.scoring.midstrand_flank_score == pytest.approx(0.95)
        assert arr.first_index == 1
        assert arr.last_index == 96
        assert arr.barcode1_pattern == "BC%02i"
        assert arr.barcode2_pattern == "BC%02i"
        assert arr.rear_only_barcodes is False


# ---------------------------------------------------------------------------
# TestTargetRef
# ---------------------------------------------------------------------------


class TestTargetRef:
    """Validate TargetRef GC content computed field."""

    def test_gc_content(self):
        """AATTCCGG should give 50.0% GC content."""
        ref = TargetRef(tgt_id="test", sequence="AATTCCGG", length=8)
        assert ref.gc_content == pytest.approx(50.0)

    def test_gc_empty(self):
        """Empty sequence should give 0.0% GC content."""
        ref = TargetRef(tgt_id="empty", sequence="", length=0)
        assert ref.gc_content == pytest.approx(0.0)
