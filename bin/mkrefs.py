"""Generate truncated reference variants from a construct TOML and barcode FASTAs.

Reads a construct TOML file, loads barcode sequences, reads target reference
FASTAs, and generates truncated reference variants at three levels:

- ``full``: complete construct with both barcodes and flanks
- ``bc1_target``: barcode1 flanks + target only
- ``bc1_only``: barcode1 flanks only

Usage::

    python bin/mkrefs.py -c construct.toml -o references/

Output structure::

    references/truncated/
      {alias}_full.fasta
      {alias}_bc1_target.fasta
      {bc1_id}_bc_only.fasta
      manifest.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from barcodes import BARCODES, load_barcodes_from_fasta, reverse_complement
from construct import parse_construct_toml


# ---------------------------------------------------------------------------
# Core logic (pure, no I/O)
# ---------------------------------------------------------------------------


def generate_truncated_refs(
    alias: str,
    target_seq: str,
    bc1_seq: str,
    bc2_seq: str | None,
    mask1_front: str,
    mask1_rear: str,
    mask2_front: str | None,
    mask2_rear: str | None,
) -> dict[str, str]:
    """Generate truncated reference sequences at multiple levels.

    Parameters
    ----------
    alias : str
        Human-readable alias for the target (used in output naming).
    target_seq : str
        The target insert DNA sequence.
    bc1_seq : str
        Barcode 1 DNA sequence (forward orientation).
    bc2_seq : str or None
        Barcode 2 DNA sequence (forward orientation, will be reverse-
        complemented in the full construct).  ``None`` for start_only mode.
    mask1_front : str
        Flanking sequence before barcode 1.
    mask1_rear : str
        Flanking sequence after barcode 1.
    mask2_front : str or None
        Flanking sequence before RC(barcode 2).  ``None`` for start_only.
    mask2_rear : str or None
        Flanking sequence after RC(barcode 2).  ``None`` for start_only.

    Returns
    -------
    dict[str, str]
        Mapping of truncation level name to reference sequence.  Keys are a
        subset of ``{"full", "bc1_target", "bc1_only"}``.  In start_only mode
        (bc2_seq is None) the ``"full"`` key is omitted.
    """
    bc1_prefix = mask1_front + bc1_seq + mask1_rear

    refs: dict[str, str] = {}

    # Full construct (only when bc2 is provided)
    if bc2_seq is not None:
        bc2_suffix = (mask2_front or "") + reverse_complement(bc2_seq) + (mask2_rear or "")
        refs["full"] = bc1_prefix + target_seq + bc2_suffix

    # BC1 + target
    refs["bc1_target"] = bc1_prefix + target_seq

    # BC1 only
    refs["bc1_only"] = bc1_prefix

    return refs


# ---------------------------------------------------------------------------
# Manifest writer
# ---------------------------------------------------------------------------


def write_manifest(entries: list[dict], path: Path) -> None:
    """Write a TSV manifest of generated reference files.

    Parameters
    ----------
    entries : list[dict]
        Each dict must contain keys: ``alias``, ``level``, ``length``, ``path``.
    path : Path
        Output path for the TSV file.
    """
    columns = ["alias", "level", "length", "path"]
    with path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for entry in entries:
            row = [str(entry[col]) for col in columns]
            fh.write("\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# FASTA I/O helpers
# ---------------------------------------------------------------------------


def _read_fasta_sequence(path: Path) -> str:
    """Read a single-record FASTA and return the sequence as uppercase.

    Parameters
    ----------
    path : Path
        Path to a FASTA file containing exactly one record.

    Returns
    -------
    str
        The concatenated sequence lines in uppercase.

    Raises
    ------
    ValueError
        If the file is empty or contains no sequence data.
    """
    seq_parts: list[str] = []
    with path.open("r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)

    if not seq_parts:
        raise ValueError(f"No sequence found in {path}")

    return "".join(seq_parts).upper()


def _write_fasta(path: Path, header: str, sequence: str) -> None:
    """Write a single-record FASTA file.

    Parameters
    ----------
    path : Path
        Output file path.
    header : str
        FASTA header (without the leading ``>``).
    sequence : str
        DNA sequence to write (will be wrapped at 80 characters).
    """
    with path.open("w") as fh:
        fh.write(f">{header}\n")
        # Wrap sequence at 80 characters per line
        for i in range(0, len(sequence), 80):
            fh.write(sequence[i : i + 80] + "\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> None:
    """CLI entry point for truncated reference generation.

    Parameters
    ----------
    argv : list[str] or None
        Command-line arguments.  Defaults to ``sys.argv[1:]``.
    """
    parser = argparse.ArgumentParser(
        description="Generate truncated reference variants from a construct TOML.",
    )
    parser.add_argument(
        "-c",
        "--construct",
        required=True,
        type=Path,
        help="Path to construct TOML file.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for reference files.",
    )
    args = parser.parse_args(argv)

    # 1. Parse construct TOML
    cfg = parse_construct_toml(args.construct)

    # 2. Load barcode sequences
    if cfg.barcode_fasta:
        barcode_seqs = load_barcodes_from_fasta(cfg.barcode_fasta)
    else:
        barcode_seqs = BARCODES

    # 3. Prepare output directory
    trunc_dir = args.outdir / "truncated"
    trunc_dir.mkdir(parents=True, exist_ok=True)

    manifest_entries: list[dict] = []

    # 4. Process each target
    for target in cfg.targets:
        # Look up barcode sequences
        bc1_id = target.barcode1.lower()
        bc1_seq = barcode_seqs[bc1_id]

        bc2_seq: str | None = None
        if target.barcode2 is not None:
            bc2_id = target.barcode2.lower()
            bc2_seq = barcode_seqs[bc2_id]

        # Read target reference FASTA
        ref_path = Path(target.reference)
        if not ref_path.is_absolute():
            # Resolve relative to construct TOML directory
            ref_path = args.construct.parent / ref_path
        target_seq = _read_fasta_sequence(ref_path)

        # Generate truncated references
        refs = generate_truncated_refs(
            alias=target.alias,
            target_seq=target_seq,
            bc1_seq=bc1_seq,
            bc2_seq=bc2_seq,
            mask1_front=cfg.arrangement.mask1_front,
            mask1_rear=cfg.arrangement.mask1_rear,
            mask2_front=cfg.arrangement.mask2_front,
            mask2_rear=cfg.arrangement.mask2_rear,
        )

        # Write FASTA files for each truncation level
        for level, seq in refs.items():
            if level == "bc1_only":
                fname = f"{bc1_id}_bc_only.fasta"
            else:
                fname = f"{target.alias}_{level}.fasta"

            fasta_path = trunc_dir / fname
            header = f"{target.alias}_{level} len={len(seq)}"
            _write_fasta(fasta_path, header, seq)

            manifest_entries.append({
                "alias": target.alias,
                "level": level,
                "length": len(seq),
                "path": f"truncated/{fname}",
            })

    # 5. Write manifest
    manifest_path = trunc_dir / "manifest.tsv"
    write_manifest(manifest_entries, manifest_path)

    # 6. Print summary
    print(f"Generated {len(manifest_entries)} truncated references")
    print(f"Output directory: {trunc_dir}")
    print(f"Manifest: {manifest_path}")
    for entry in manifest_entries:
        print(f"  {entry['alias']:30s} {entry['level']:15s} {entry['length']:>6d}bp  {entry['path']}")


if __name__ == "__main__":
    main()
