"""Pydantic v2 models for SMA-seq experiment configuration.

Defines the data structures used across the config visualizer: scoring
parameters, arrangement (custom barcode kits), sample sheets, target
references, demultiplexing, truncation handling, and the top-level
experiment configuration.
"""

from __future__ import annotations

import re
from typing import Literal

from pydantic import BaseModel, Field, computed_field, field_validator

# ---------------------------------------------------------------------------
# Barcode regex for duplexed format validation
# ---------------------------------------------------------------------------

_DUPLEXED_RE = re.compile(r"^(barcode\d+)--(barcode\d+)$")
_BARCODE_NUM_RE = re.compile(r"barcode(\d+)")


def _normalize_barcode(name: str) -> str:
    """Normalize ``barcodeNN`` to ``nbNN`` (zero-padded, range 1-96)."""
    m = _BARCODE_NUM_RE.fullmatch(name)
    if not m:
        raise ValueError(f"Invalid barcode name '{name}': expected 'barcodeNN'")
    num = int(m.group(1))
    if num < 1 or num > 96:
        raise ValueError(f"Barcode number {num} out of valid range 1-96")
    return f"nb{num:02d}"


# ---------------------------------------------------------------------------
# 1. ScoringParams - Dorado barcode scoring parameters
# ---------------------------------------------------------------------------


class ScoringParams(BaseModel):
    """Dorado barcode scoring parameters with sensible defaults."""

    max_barcode_penalty: int = 11
    min_barcode_penalty_dist: int = 3
    min_separation_only_dist: int = 6
    barcode_end_proximity: int = 75
    flank_left_pad: int = 5
    flank_right_pad: int = 10
    front_barcode_window: int = 175
    rear_barcode_window: int = 175
    midstrand_flank_score: float = 0.95


# ---------------------------------------------------------------------------
# 2. ArrangementConfig - Dorado custom barcode arrangement
# ---------------------------------------------------------------------------


class ArrangementConfig(BaseModel):
    """Dorado custom barcode arrangement file representation."""

    name: str
    kit: str = ""
    mask1_front: str = ""
    mask1_rear: str = ""
    mask2_front: str = ""
    mask2_rear: str = ""
    barcode1_pattern: str = "BC%02i"
    barcode2_pattern: str = "BC%02i"
    first_index: int = 1
    last_index: int = 96
    rear_only_barcodes: bool = False
    scoring: ScoringParams = Field(default_factory=ScoringParams)


# ---------------------------------------------------------------------------
# 3. SampleSheetEntry - MinKNOW CSV row
# ---------------------------------------------------------------------------


class SampleSheetEntry(BaseModel):
    """A single row from a MinKNOW sample sheet CSV.

    The barcode field must be in duplexed ``barcodeNN--barcodeNN`` format.
    Computed fields ``upstream_barcode`` and ``downstream_barcode`` normalize
    the barcode names to ``nbNN`` short format.
    """

    flow_cell_id: str = ""
    kit: str = ""
    barcode: str
    alias: str
    experiment_id: str = ""
    sample_id: str = ""
    type: str = "test_sample"

    @field_validator("barcode")
    @classmethod
    def _validate_duplexed_barcode(cls, v: str) -> str:
        if not _DUPLEXED_RE.match(v):
            raise ValueError(
                f"Barcode '{v}' is not in duplexed 'barcodeNN--barcodeNN' format. "
                f"Only duplexed barcodes are supported."
            )
        return v

    @computed_field  # type: ignore[prop-decorator]
    @property
    def upstream_barcode(self) -> str:
        """Upstream barcode normalized to ``nbNN`` format."""
        m = _DUPLEXED_RE.match(self.barcode)
        if not m:
            return ""
        return _normalize_barcode(m.group(1))

    @computed_field  # type: ignore[prop-decorator]
    @property
    def downstream_barcode(self) -> str:
        """Downstream barcode normalized to ``nbNN`` format."""
        m = _DUPLEXED_RE.match(self.barcode)
        if not m:
            return ""
        return _normalize_barcode(m.group(2))


# ---------------------------------------------------------------------------
# 4. TargetRef - Reference sequence
# ---------------------------------------------------------------------------


class TargetRef(BaseModel):
    """A target reference sequence with computed GC content."""

    tgt_id: str
    sequence: str
    length: int

    @computed_field  # type: ignore[prop-decorator]
    @property
    def gc_content(self) -> float:
        """GC content as a percentage (0.0-100.0)."""
        if not self.sequence:
            return 0.0
        gc_count = sum(1 for base in self.sequence.upper() if base in ("G", "C"))
        return round((gc_count / len(self.sequence)) * 100.0, 2)


# ---------------------------------------------------------------------------
# 5. BarcodeAssignment
# ---------------------------------------------------------------------------


class BarcodeAssignment(BaseModel):
    """Barcode-to-target assignment mapping."""

    assignments: dict[str, str]


# ---------------------------------------------------------------------------
# 6. DemuxPairEntry
# ---------------------------------------------------------------------------


class DemuxPairEntry(BaseModel):
    """A single demultiplexing pair entry."""

    start: str
    end: str
    alias: str


# ---------------------------------------------------------------------------
# 7. DemultiplexingConfig
# ---------------------------------------------------------------------------


class DemultiplexingConfig(BaseModel):
    """Demultiplexing configuration for barcode pair resolution."""

    mode: str = "dual"
    start_barcode: str = ""
    end_barcode: str = ""
    pairs: list[DemuxPairEntry] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# 8. TruncationRules
# ---------------------------------------------------------------------------


class TruncationRules(BaseModel):
    """Read classification rules for truncation analysis.

    Each field describes the action to take for a particular read category.
    """

    full: str = "assign_to_target"
    trunc_3prime: str = "assign_to_target"
    trunc_target: str = "assign_by_start_barcode"
    trunc_barcode: str = "assign_by_start_barcode"
    adapter_only: str = "unclassified"
    chimeric: str = "flag_for_review"


# ---------------------------------------------------------------------------
# 9. AutoReferences
# ---------------------------------------------------------------------------


class AutoReferences(BaseModel):
    """Configuration for automatic reference sequence detection."""

    enabled: bool = True


# ---------------------------------------------------------------------------
# 10. TruncationConfig
# ---------------------------------------------------------------------------


class TruncationConfig(BaseModel):
    """Truncation analysis configuration with thresholds and rules."""

    min_barcode_confidence: float = 0.5
    min_target_fraction: float = 0.1
    adapter_search_window: int = 50
    rules: TruncationRules = Field(default_factory=TruncationRules)
    auto_references: AutoReferences = Field(default_factory=AutoReferences)


# ---------------------------------------------------------------------------
# 11. ConstructConfig
# ---------------------------------------------------------------------------


class ConstructConfig(BaseModel):
    """Describes the physical construct structure (adapters + insert)."""

    adapter_5prime: str = ""
    adapter_3prime: str = ""
    insert_type: Literal["amplicon", "genomic", "synthetic"] = "amplicon"
    notes: str = ""


# ---------------------------------------------------------------------------
# 12. Assumption
# ---------------------------------------------------------------------------


class Assumption(BaseModel):
    """A documented assumption about the experiment configuration."""

    key: str
    text: str
    why: str


# ---------------------------------------------------------------------------
# 13. QualityMetrics
# ---------------------------------------------------------------------------


class QualityMetrics(BaseModel):
    """Quality metric formula documentation strings."""

    q_bc: str = "Probability-averaged Phred: -10*log10(mean(10^(-Qi/10)))"
    q_ld: str = "-10*log10(min(max(1/L^2, ed/L), 1.0))"


# ---------------------------------------------------------------------------
# 14. ClassificationConfig
# ---------------------------------------------------------------------------


class ClassificationConfig(BaseModel):
    """Read classification configuration."""

    barcode_search_window: int = 100
    confidence_formula: str = "1.0 - (ed / barcode_length)"
    ambiguity_triggers_full_construct: bool = True


# ---------------------------------------------------------------------------
# 15. ExperimentConfig - top-level
# ---------------------------------------------------------------------------


_DEFAULT_ASSUMPTIONS: list[dict[str, str]] = [
    {
        "key": "strand_orientation",
        "text": "Reads are in the original strand orientation after basecalling",
        "why": "Dorado outputs reads in template orientation by default",
    },
    {
        "key": "barcode_proximity",
        "text": "Barcodes are located within the search window of read ends",
        "why": "Standard library prep places barcodes adjacent to adapters",
    },
    {
        "key": "single_insert",
        "text": "Each read contains at most one target insert",
        "why": "Chimeric reads are rare and classified separately",
    },
    {
        "key": "adapter_integrity",
        "text": "Adapter sequences are intact and detectable by alignment",
        "why": "Degraded adapters would prevent proper read classification",
    },
]


class ExperimentConfig(BaseModel):
    """Top-level experiment configuration for SMA-seq analysis.

    Aggregates all sub-configurations: construct, demultiplexing,
    classification, quality metrics, truncation rules, and documented
    assumptions.
    """

    description: str = ""
    construct: ConstructConfig = Field(default_factory=ConstructConfig)
    demultiplexing: DemultiplexingConfig = Field(
        default_factory=DemultiplexingConfig
    )
    classification: ClassificationConfig = Field(
        default_factory=ClassificationConfig
    )
    quality: QualityMetrics = Field(default_factory=QualityMetrics)
    truncation: TruncationConfig = Field(default_factory=TruncationConfig)
    assumptions: list[Assumption] = Field(
        default_factory=lambda: [Assumption(**a) for a in _DEFAULT_ASSUMPTIONS],
    )
