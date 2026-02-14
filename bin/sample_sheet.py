"""MinKNOW sample sheet parser for duplexed barcode experiments.

Parses CSV sample sheets produced by MinKNOW, extracts duplexed barcode pairs,
maps them to target aliases, and detects barcode ambiguity.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_BARCODE_RE = re.compile(r"barcode(\d+)")
_DUPLEXED_RE = re.compile(r"^(barcode\d+)--(barcode\d+)$")

# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


def parse_sample_sheet(path: Path | str) -> dict[tuple[str, str], str]:
    """Parse a MinKNOW sample sheet CSV and return barcode-pair-to-alias mapping.

    Parameters
    ----------
    path : Path or str
        Path to the sample sheet CSV file.  Must contain at least ``barcode``
        and ``alias`` columns.

    Returns
    -------
    dict[tuple[str, str], str]
        Mapping of ``(upstream_barcode, downstream_barcode)`` tuples to the
        alias string.  Barcode names are normalised from ``barcodeNN`` format
        to ``nbNN`` (zero-padded to 2 digits).

    Raises
    ------
    ValueError
        If required columns (``barcode``, ``alias``) are missing or if a
        barcode entry is not in duplexed ``barcodeNN--barcodeNN`` format.
    """
    path = Path(path)
    mapping: dict[tuple[str, str], str] = {}

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []

        # Validate required columns
        if "barcode" not in fieldnames:
            raise ValueError(
                f"Sample sheet is missing required 'barcode' column. "
                f"Found columns: {fieldnames}"
            )
        if "alias" not in fieldnames:
            raise ValueError(
                f"Sample sheet is missing required 'alias' column. "
                f"Found columns: {fieldnames}"
            )

        for row_num, row in enumerate(reader, start=2):  # start=2 for 1-indexed + header
            barcode_field = row["barcode"].strip()
            alias = row["alias"].strip()

            # Parse duplexed format
            match = _DUPLEXED_RE.match(barcode_field)
            if not match:
                raise ValueError(
                    f"Row {row_num}: barcode '{barcode_field}' is not in duplexed "
                    f"'barcodeNN--barcodeNN' format. Only duplexed barcodes are supported."
                )

            upstream = _normalize_barcode_name(match.group(1))
            downstream = _normalize_barcode_name(match.group(2))
            mapping[(upstream, downstream)] = alias

    return mapping


def detect_barcode_ambiguity(mapping: dict[tuple[str, str], str]) -> bool:
    """Check if any single barcode appears as upstream in pairs with different aliases.

    Parameters
    ----------
    mapping : dict[tuple[str, str], str]
        Barcode-pair-to-alias mapping as returned by :func:`parse_sample_sheet`.

    Returns
    -------
    bool
        ``True`` if ambiguity is detected (a barcode is upstream for multiple
        different aliases), ``False`` otherwise.
    """
    upstream_aliases: dict[str, set[str]] = {}

    for (upstream, _downstream), alias in mapping.items():
        if upstream not in upstream_aliases:
            upstream_aliases[upstream] = set()
        upstream_aliases[upstream].add(alias)

    return any(len(aliases) > 1 for aliases in upstream_aliases.values())


def _normalize_barcode_name(name: str) -> str:
    """Normalise a MinKNOW barcode name to the ``nbNN`` short format.

    Parameters
    ----------
    name : str
        A barcode name like ``barcode05``, ``barcode5``, or ``barcode10``.

    Returns
    -------
    str
        Normalised name in ``nbNN`` format (zero-padded to 2 digits).

    Raises
    ------
    ValueError
        If the name does not match ``barcode(\\d+)`` or the number is outside
        the valid range 1-96.
    """
    match = _BARCODE_RE.fullmatch(name)
    if not match:
        raise ValueError(
            f"Invalid barcode name '{name}': expected format 'barcodeNN'"
        )

    num = int(match.group(1))
    if num < 1 or num > 96:
        raise ValueError(
            f"Barcode number {num} out of valid range 1-96"
        )

    return f"nb{num:02d}"
