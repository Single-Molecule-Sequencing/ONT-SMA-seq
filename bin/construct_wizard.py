"""Interactive CLI wizard for creating SMA-seq construct TOML files.

Provides non-interactive helper functions for reference discovery, construct
diagram formatting, and TOML generation, as well as an interactive wizard
that guides users through the 7-step construct creation process.

Usage::

    # Interactive wizard
    python construct_wizard.py -o construct.toml -rd references/

    # Edit an existing construct
    python construct_wizard.py -o construct.toml -i existing.toml

Non-interactive helpers (importable)::

    from construct_wizard import (
        KNOWN_KITS,
        discover_references,
        format_construct_diagram,
        build_construct_toml,
    )
"""

from __future__ import annotations

import argparse
from pathlib import Path

import tomli_w

from construct import parse_construct_toml

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

KNOWN_KITS: dict[str, dict] = {
    "SQK-NBD114-96": {"count": 96, "length": 24, "desc": "Native Barcoding Kit 96"},
    "SQK-NBD114-24": {"count": 24, "length": 24, "desc": "Native Barcoding Kit 24"},
    "SQK-RBK114-96": {"count": 96, "length": 24, "desc": "Rapid Barcoding Kit 96"},
    "SQK-RBK114-24": {"count": 24, "length": 24, "desc": "Rapid Barcoding Kit 24"},
}


# ---------------------------------------------------------------------------
# Non-interactive helpers
# ---------------------------------------------------------------------------


def discover_references(ref_dir: Path) -> list[dict]:
    """Find FASTA files in a directory and return metadata.

    Scans *ref_dir* for files with ``.fasta`` or ``.fa`` extensions, reads
    each file to determine the total sequence length (excluding header lines),
    and returns a list of dicts.

    Parameters
    ----------
    ref_dir : Path
        Directory to scan for FASTA files.

    Returns
    -------
    list[dict]
        Each dict has keys ``name`` (file stem), ``path`` (str), and
        ``length`` (total bases as int).
    """
    ref_dir = Path(ref_dir)
    results: list[dict] = []
    for ext in ("*.fasta", "*.fa"):
        for fasta_path in sorted(ref_dir.glob(ext)):
            seq_length = _fasta_sequence_length(fasta_path)
            results.append({
                "name": fasta_path.stem,
                "path": str(fasta_path),
                "length": seq_length,
            })
    return results


def _fasta_sequence_length(fasta_path: Path) -> int:
    """Read a FASTA file and return the total number of sequence bases."""
    total = 0
    for line in fasta_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        total += len(line)
    return total


def format_construct_diagram(
    mode: str,
    mask1_front: str,
    mask1_rear: str,
    mask2_front: str | None = None,
    mask2_rear: str | None = None,
) -> str:
    """Generate an ASCII diagram of the construct structure.

    Parameters
    ----------
    mode : str
        ``"dual_independent"`` or ``"start_only"``.
    mask1_front : str
        Flanking sequence before barcode 1.
    mask1_rear : str
        Flanking sequence after barcode 1.
    mask2_front : str, optional
        Flanking sequence before RC(barcode 2) (dual mode only).
    mask2_rear : str, optional
        Flanking sequence after RC(barcode 2) (dual mode only).

    Returns
    -------
    str
        Multi-line ASCII diagram of the construct.
    """
    lines: list[str] = []
    lines.append("Construct Diagram")
    lines.append("=" * 60)

    if mode == "dual_independent":
        lines.append("")
        lines.append(
            f"5'--[{mask1_front}]--[BC1]--[{mask1_rear}]"
            f"--TARGET--"
            f"[{mask2_front}]--[RC(BC2)]--[{mask2_rear}]--3'"
        )
        lines.append("")
        lines.append("Components:")
        lines.append(f"  mask1_front : {mask1_front}")
        lines.append(f"  BC1         : barcode 1")
        lines.append(f"  mask1_rear  : {mask1_rear}")
        lines.append(f"  TARGET      : insert sequence")
        lines.append(f"  mask2_front : {mask2_front}")
        lines.append(f"  RC(BC2)     : reverse-complement barcode 2")
        lines.append(f"  mask2_rear  : {mask2_rear}")
    else:
        # start_only
        lines.append("")
        lines.append(
            f"5'--[{mask1_front}]--[BC1]--[{mask1_rear}]"
            f"--TARGET--3'"
        )
        lines.append("")
        lines.append("Components:")
        lines.append(f"  mask1_front : {mask1_front}")
        lines.append(f"  BC1         : barcode 1")
        lines.append(f"  mask1_rear  : {mask1_rear}")
        lines.append(f"  TARGET      : insert sequence")

    lines.append("=" * 60)
    return "\n".join(lines)


def build_construct_toml(
    output_path: Path | str,
    name: str,
    kit: str,
    mode: str,
    mask1_front: str,
    mask1_rear: str,
    targets: list[dict],
    mask2_front: str | None = None,
    mask2_rear: str | None = None,
    barcode_fasta: str | None = None,
    confidence: dict | None = None,
    truncation: dict | None = None,
) -> None:
    """Build and write a valid construct TOML file.

    The output is compatible with :func:`construct.parse_construct_toml`.

    Parameters
    ----------
    output_path : Path or str
        Where to write the TOML file.
    name : str
        Construct name (e.g. ``"SMA_CYP2D6_dual"``).
    kit : str
        Barcoding kit identifier (e.g. ``"SQK-NBD114-96"``).
    mode : str
        ``"dual_independent"`` or ``"start_only"``.
    mask1_front : str
        DNA sequence flanking before barcode 1.
    mask1_rear : str
        DNA sequence flanking after barcode 1.
    targets : list[dict]
        List of target dicts, each with ``barcode1``, ``alias``,
        ``reference``, and optionally ``barcode2``.
    mask2_front : str, optional
        DNA sequence flanking before RC(barcode 2) (dual mode).
    mask2_rear : str, optional
        DNA sequence flanking after RC(barcode 2) (dual mode).
    barcode_fasta : str, optional
        Path to custom barcode FASTA file.
    confidence : dict, optional
        Override confidence thresholds. Keys: ``full_length_threshold``,
        ``start_barcode_min``, ``flank_max_error_rate``.
    truncation : dict, optional
        Override truncation settings. Keys: ``auto_generate_refs``,
        ``min_target_length``.
    """
    output_path = Path(output_path)

    # Determine kit barcode count for last_index
    kit_info = KNOWN_KITS.get(kit, {})
    last_index = kit_info.get("count", 96)

    # Build the arrangement section
    arrangement: dict = {
        "name": name,
        "kit": kit,
        "mask1_front": mask1_front,
        "mask1_rear": mask1_rear,
        "barcode1_pattern": "NB%02i",
        "first_index": 1,
        "last_index": last_index,
    }

    if mode == "dual_independent":
        arrangement["mask2_front"] = mask2_front or ""
        arrangement["mask2_rear"] = mask2_rear or ""
        arrangement["barcode2_pattern"] = "NB%02i"

    # Build the scoring section with defaults
    scoring: dict = {
        "max_barcode_penalty": 11,
        "min_barcode_penalty_dist": 3,
        "front_barcode_window": 100,
        "rear_barcode_window": 100,
        "min_flank_score": 0.5,
    }

    # Build the sma section
    sma: dict = {
        "mode": mode,
    }

    if barcode_fasta is not None:
        sma["barcode_fasta"] = barcode_fasta

    # Confidence sub-section
    conf = {
        "full_length_threshold": 0.75,
        "start_barcode_min": 0.6,
        "flank_max_error_rate": 0.5,
    }
    if confidence:
        conf.update(confidence)
    sma["confidence"] = conf

    # Truncation sub-section
    trunc = {
        "auto_generate_refs": True,
        "min_target_length": 20,
    }
    if truncation:
        trunc.update(truncation)
    sma["truncation"] = trunc

    # Targets sub-section
    sma["targets"] = []
    for tgt in targets:
        entry: dict = {
            "barcode1": tgt["barcode1"],
            "alias": tgt["alias"],
            "reference": tgt["reference"],
        }
        if "barcode2" in tgt and tgt["barcode2"] is not None:
            entry["barcode2"] = tgt["barcode2"]
        sma["targets"].append(entry)

    # Assemble full document
    doc: dict = {
        "arrangement": arrangement,
        "scoring": scoring,
        "sma": sma,
    }

    # Write using tomli_w
    with output_path.open("wb") as fh:
        tomli_w.dump(doc, fh)


# ---------------------------------------------------------------------------
# Interactive wizard
# ---------------------------------------------------------------------------


def run_wizard(
    output_path: Path | str,
    ref_dir: Path | str | None = None,
    input_toml: Path | str | None = None,
) -> None:
    """Interactive 7-step wizard for construct TOML creation.

    Steps:
        1. Mode selection (dual_independent or start_only)
        2. Kit selection
        3. Flanking sequences
        4. Target definitions
        5. Confidence thresholds
        6. Truncation settings
        7. Review and save

    Parameters
    ----------
    output_path : Path or str
        Where to write the output TOML file.
    ref_dir : Path or str, optional
        Directory containing reference FASTA files for auto-discovery.
    input_toml : Path or str, optional
        Existing TOML file to use as starting values (edit mode).
    """
    output_path = Path(output_path)

    # Load defaults from existing TOML if provided
    defaults: dict = {}
    if input_toml is not None:
        cfg = parse_construct_toml(input_toml)
        defaults = {
            "name": cfg.arrangement.name,
            "kit": cfg.arrangement.kit,
            "mode": cfg.mode,
            "mask1_front": cfg.arrangement.mask1_front,
            "mask1_rear": cfg.arrangement.mask1_rear,
            "mask2_front": cfg.arrangement.mask2_front,
            "mask2_rear": cfg.arrangement.mask2_rear,
        }

    # Discover references if ref_dir provided
    available_refs: list[dict] = []
    if ref_dir is not None:
        available_refs = discover_references(Path(ref_dir))

    # --- Step 1: Mode ---
    print("\n=== Step 1/7: Mode Selection ===")
    print("  1. dual_independent  (BC1 + RC(BC2) flanking target)")
    print("  2. start_only        (BC1 only, no second barcode)")
    default_mode = defaults.get("mode", "dual_independent")
    mode_input = input(f"Select mode [default: {default_mode}]: ").strip()
    if mode_input == "2":
        mode = "start_only"
    elif mode_input == "1" or not mode_input:
        mode = default_mode
    else:
        mode = mode_input

    # --- Step 2: Kit ---
    print("\n=== Step 2/7: Kit Selection ===")
    for i, (kit_id, info) in enumerate(KNOWN_KITS.items(), 1):
        print(f"  {i}. {kit_id} - {info['desc']} ({info['count']} barcodes)")
    default_kit = defaults.get("kit", "SQK-NBD114-96")
    kit_input = input(f"Select kit [default: {default_kit}]: ").strip()
    kit_list = list(KNOWN_KITS.keys())
    if kit_input.isdigit() and 1 <= int(kit_input) <= len(kit_list):
        kit = kit_list[int(kit_input) - 1]
    elif kit_input:
        kit = kit_input
    else:
        kit = default_kit

    # --- Step 3: Flanking sequences ---
    print("\n=== Step 3/7: Flanking Sequences ===")
    default_m1f = defaults.get("mask1_front", "")
    mask1_front = input(f"mask1_front [default: {default_m1f}]: ").strip() or default_m1f
    default_m1r = defaults.get("mask1_rear", "")
    mask1_rear = input(f"mask1_rear [default: {default_m1r}]: ").strip() or default_m1r

    mask2_front = None
    mask2_rear = None
    if mode == "dual_independent":
        default_m2f = defaults.get("mask2_front", "")
        mask2_front = input(f"mask2_front [default: {default_m2f}]: ").strip() or default_m2f
        default_m2r = defaults.get("mask2_rear", "")
        mask2_rear = input(f"mask2_rear [default: {default_m2r}]: ").strip() or default_m2r

    # --- Step 4: Targets ---
    print("\n=== Step 4/7: Target Definitions ===")
    if available_refs:
        print("Discovered references:")
        for i, ref in enumerate(available_refs, 1):
            print(f"  {i}. {ref['name']} ({ref['length']} bp) - {ref['path']}")

    name = input("Construct name: ").strip() or "SMA_construct"
    targets: list[dict] = []
    while True:
        print(f"\n  Adding target {len(targets) + 1}:")
        bc1 = input("    barcode1 (e.g. NB01): ").strip()
        if not bc1:
            break
        alias = input("    alias: ").strip()
        reference = input("    reference path: ").strip()
        tgt: dict = {"barcode1": bc1, "alias": alias, "reference": reference}
        if mode == "dual_independent":
            bc2 = input("    barcode2 (e.g. NB02): ").strip()
            tgt["barcode2"] = bc2
        targets.append(tgt)
        more = input("  Add another target? [y/N]: ").strip().lower()
        if more != "y":
            break

    # --- Step 5: Confidence ---
    print("\n=== Step 5/7: Confidence Thresholds ===")
    conf_input = input("Use defaults (0.75 / 0.6 / 0.5)? [Y/n]: ").strip().lower()
    confidence = None
    if conf_input == "n":
        confidence = {
            "full_length_threshold": float(input("  full_length_threshold: ").strip()),
            "start_barcode_min": float(input("  start_barcode_min: ").strip()),
            "flank_max_error_rate": float(input("  flank_max_error_rate: ").strip()),
        }

    # --- Step 6: Truncation ---
    print("\n=== Step 6/7: Truncation Settings ===")
    trunc_input = input("Use defaults (auto_generate_refs=true, min_target_length=20)? [Y/n]: ").strip().lower()
    truncation_settings = None
    if trunc_input == "n":
        truncation_settings = {
            "auto_generate_refs": input("  auto_generate_refs (true/false): ").strip().lower() == "true",
            "min_target_length": int(input("  min_target_length: ").strip()),
        }

    # --- Step 7: Review & Save ---
    print("\n=== Step 7/7: Review & Save ===")
    diagram = format_construct_diagram(
        mode=mode,
        mask1_front=mask1_front,
        mask1_rear=mask1_rear,
        mask2_front=mask2_front,
        mask2_rear=mask2_rear,
    )
    print(diagram)
    print(f"\nName: {name}")
    print(f"Kit:  {kit}")
    print(f"Mode: {mode}")
    print(f"Targets: {len(targets)}")
    for tgt in targets:
        print(f"  - {tgt['alias']}: {tgt['barcode1']}", end="")
        if "barcode2" in tgt:
            print(f" / {tgt['barcode2']}", end="")
        print(f" -> {tgt['reference']}")

    confirm = input(f"\nWrite to {output_path}? [Y/n]: ").strip().lower()
    if confirm == "n":
        print("Aborted.")
        return

    build_construct_toml(
        output_path=output_path,
        name=name,
        kit=kit,
        mode=mode,
        mask1_front=mask1_front,
        mask1_rear=mask1_rear,
        targets=targets,
        mask2_front=mask2_front,
        mask2_rear=mask2_rear,
        confidence=confidence,
        truncation=truncation_settings,
    )
    print(f"Construct TOML written to {output_path}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """CLI entry point for the construct wizard."""
    parser = argparse.ArgumentParser(
        description="Interactive wizard for creating SMA-seq construct TOML files.",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        help="Output path for the construct TOML file.",
    )
    parser.add_argument(
        "-rd", "--ref-dir",
        type=Path,
        default=None,
        help="Directory containing reference FASTA files for auto-discovery.",
    )
    parser.add_argument(
        "-i", "--input",
        type=Path,
        default=None,
        help="Existing construct TOML to use as starting values (edit mode).",
    )
    args = parser.parse_args()
    run_wizard(
        output_path=args.output,
        ref_dir=args.ref_dir,
        input_toml=args.input,
    )


if __name__ == "__main__":
    main()
