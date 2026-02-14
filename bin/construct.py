"""Parse and validate extended Dorado TOML arrangement files for SMA-seq.

Reads a construct TOML file containing standard Dorado [arrangement] and
[scoring] sections plus an SMA-specific [sma] extension.  Returns a
``ConstructConfig`` dataclass with helper methods for barcode pairing and
flank access.

Usage::

    from construct import parse_construct_toml

    cfg = parse_construct_toml("construct.toml")
    mapping = cfg.barcode_pair_to_alias()   # dict[(str,str), str]
    ids = cfg.used_barcode_ids()            # set[str]
"""

from __future__ import annotations

import re
import tomllib
from dataclasses import dataclass, field
from pathlib import Path


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class ValidationError(Exception):
    """Raised when a construct TOML file fails validation."""


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_VALID_MODES = {"dual_independent", "start_only"}
_DNA_RE = re.compile(r"^[ACGTacgt]+$")
_BARCODE_ID_RE = re.compile(r"^NB(\d{2})$")


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class ArrangementConfig:
    """Standard Dorado [arrangement] section."""

    name: str
    kit: str
    mask1_front: str
    mask1_rear: str
    mask2_front: str | None = None
    mask2_rear: str | None = None
    barcode1_pattern: str = "NB%02i"
    barcode2_pattern: str | None = None
    first_index: int = 1
    last_index: int = 96


@dataclass
class ScoringConfig:
    """Standard Dorado [scoring] section with defaults."""

    max_barcode_penalty: int = 11
    min_barcode_penalty_dist: int = 3
    front_barcode_window: int = 100
    rear_barcode_window: int = 100
    min_flank_score: float = 0.5


@dataclass
class ConfidenceConfig:
    """SMA confidence thresholds ([sma.confidence])."""

    full_length_threshold: float = 0.75
    start_barcode_min: float = 0.6
    flank_max_error_rate: float = 0.5


@dataclass
class TruncationConfig:
    """SMA truncation settings ([sma.truncation])."""

    auto_generate_refs: bool = True
    min_target_length: int = 20


@dataclass
class TargetEntry:
    """A single [[sma.targets]] entry."""

    barcode1: str
    barcode2: str | None
    alias: str
    reference: str


@dataclass
class ConstructConfig:
    """Complete parsed and validated construct configuration."""

    arrangement: ArrangementConfig
    scoring: ScoringConfig
    mode: str
    barcode_fasta: str | None
    confidence: ConfidenceConfig
    truncation: TruncationConfig
    targets: list[TargetEntry] = field(default_factory=list)

    def barcode_pair_to_alias(self) -> dict[tuple[str, str], str]:
        """Return barcode-pair-to-alias mapping compatible with sample_sheet.

        Returns
        -------
        dict[tuple[str, str], str]
            Mapping of ``(upstream_barcode, downstream_barcode)`` tuples to
            the alias string.  Barcode IDs are normalised to lowercase ``nbNN``
            format.
        """
        mapping: dict[tuple[str, str], str] = {}
        for t in self.targets:
            bc1 = _normalise_barcode_id(t.barcode1)
            if self.mode == "start_only":
                bc2 = bc1
            else:
                bc2 = _normalise_barcode_id(t.barcode2)  # type: ignore[arg-type]
            mapping[(bc1, bc2)] = t.alias
        return mapping

    def used_barcode_ids(self) -> set[str]:
        """Return all barcode IDs referenced in targets (normalised to nbNN)."""
        ids: set[str] = set()
        for t in self.targets:
            ids.add(_normalise_barcode_id(t.barcode1))
            if t.barcode2 is not None:
                ids.add(_normalise_barcode_id(t.barcode2))
        return ids

    @property
    def flank_front(self) -> str:
        """Trailing flank after barcode1 (mask1_rear)."""
        return self.arrangement.mask1_rear

    @property
    def flank_rear(self) -> str:
        """Leading flank before RC(barcode2) (mask2_front).

        For start_only mode this returns an empty string.
        """
        return self.arrangement.mask2_front or ""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _normalise_barcode_id(barcode_id: str) -> str:
    """Normalise ``NB05`` to ``nb05``."""
    return barcode_id.lower()


def _validate_barcode_id(barcode_id: str) -> None:
    """Validate that a barcode ID matches NB01-NB96 format."""
    m = _BARCODE_ID_RE.match(barcode_id.upper())
    if not m:
        raise ValidationError(
            f"Invalid barcode ID '{barcode_id}': must match NB01-NB96"
        )
    num = int(m.group(1))
    if num < 1 or num > 96:
        raise ValidationError(
            f"Invalid barcode ID '{barcode_id}': number must be 01-96"
        )


def _validate_dna(seq: str, field_name: str) -> None:
    """Validate that a sequence contains only DNA bases."""
    if not _DNA_RE.match(seq):
        raise ValidationError(
            f"'{field_name}' contains non-DNA characters: '{seq}'"
        )


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


def parse_construct_toml(path: Path | str) -> ConstructConfig:
    """Parse and validate a construct TOML file.

    Parameters
    ----------
    path : Path or str
        Path to the construct TOML file.

    Returns
    -------
    ConstructConfig
        Validated configuration object.

    Raises
    ------
    ValidationError
        If required sections or fields are missing, or if values are invalid.
    """
    path = Path(path)
    with path.open("rb") as fh:
        data = tomllib.load(fh)

    # --- Validate required top-level sections ---
    if "arrangement" not in data:
        raise ValidationError(
            "Missing required [arrangement] section in construct TOML"
        )
    if "sma" not in data:
        raise ValidationError(
            "Missing required [sma] section in construct TOML"
        )

    # --- Parse [arrangement] ---
    arr_data = data["arrangement"]
    _require_fields(arr_data, ["name", "kit", "mask1_front", "mask1_rear",
                                "barcode1_pattern", "first_index", "last_index"],
                    "arrangement")

    # Validate DNA flanking sequences
    for flank_field in ["mask1_front", "mask1_rear"]:
        _validate_dna(arr_data[flank_field], flank_field)
    for flank_field in ["mask2_front", "mask2_rear"]:
        if flank_field in arr_data:
            _validate_dna(arr_data[flank_field], flank_field)

    arrangement = ArrangementConfig(
        name=arr_data["name"],
        kit=arr_data["kit"],
        mask1_front=arr_data["mask1_front"],
        mask1_rear=arr_data["mask1_rear"],
        mask2_front=arr_data.get("mask2_front"),
        mask2_rear=arr_data.get("mask2_rear"),
        barcode1_pattern=arr_data.get("barcode1_pattern", "NB%02i"),
        barcode2_pattern=arr_data.get("barcode2_pattern"),
        first_index=arr_data.get("first_index", 1),
        last_index=arr_data.get("last_index", 96),
    )

    # --- Parse [scoring] (all optional with defaults) ---
    scr_data = data.get("scoring", {})
    scoring = ScoringConfig(
        max_barcode_penalty=scr_data.get("max_barcode_penalty", 11),
        min_barcode_penalty_dist=scr_data.get("min_barcode_penalty_dist", 3),
        front_barcode_window=scr_data.get("front_barcode_window", 100),
        rear_barcode_window=scr_data.get("rear_barcode_window", 100),
        min_flank_score=scr_data.get("min_flank_score", 0.5),
    )

    # --- Parse [sma] ---
    sma_data = data["sma"]
    mode = sma_data.get("mode", "dual_independent")
    if mode not in _VALID_MODES:
        raise ValidationError(
            f"Invalid mode '{mode}': must be one of {sorted(_VALID_MODES)}"
        )

    barcode_fasta = sma_data.get("barcode_fasta")

    # --- Parse [sma.confidence] (all optional with defaults) ---
    conf_data = sma_data.get("confidence", {})
    confidence = ConfidenceConfig(
        full_length_threshold=conf_data.get("full_length_threshold", 0.75),
        start_barcode_min=conf_data.get("start_barcode_min", 0.6),
        flank_max_error_rate=conf_data.get("flank_max_error_rate", 0.5),
    )

    # --- Parse [sma.truncation] (all optional with defaults) ---
    trunc_data = sma_data.get("truncation", {})
    truncation = TruncationConfig(
        auto_generate_refs=trunc_data.get("auto_generate_refs", True),
        min_target_length=trunc_data.get("min_target_length", 20),
    )

    # --- Parse [[sma.targets]] ---
    targets_data = sma_data.get("targets", [])
    targets: list[TargetEntry] = []
    for i, tgt in enumerate(targets_data):
        _require_fields(tgt, ["barcode1", "alias", "reference"],
                        f"sma.targets[{i}]")

        bc1 = tgt["barcode1"]
        _validate_barcode_id(bc1)

        bc2 = tgt.get("barcode2")
        if mode == "dual_independent" and bc2 is None:
            raise ValidationError(
                f"sma.targets[{i}]: barcode2 is required in "
                f"dual_independent mode"
            )
        if bc2 is not None:
            _validate_barcode_id(bc2)

        targets.append(TargetEntry(
            barcode1=bc1,
            barcode2=bc2,
            alias=tgt["alias"],
            reference=tgt["reference"],
        ))

    return ConstructConfig(
        arrangement=arrangement,
        scoring=scoring,
        mode=mode,
        barcode_fasta=barcode_fasta,
        confidence=confidence,
        truncation=truncation,
        targets=targets,
    )


def _require_fields(
    data: dict,
    required: list[str],
    section: str,
) -> None:
    """Check that all required fields are present in a TOML section."""
    for field_name in required:
        if field_name not in data:
            raise ValidationError(
                f"Missing required field '{field_name}' in [{section}]"
            )
