#!/usr/bin/env python3
"""
inputInit.py - Input Standardization Script

Standardizes input file paths and extracts metadata from the incoming BAM filename.
Creates symlinks in a standardized Input/ directory structure.

Usage:
    python inputInit.py --bam <path_to_bam> --pod5 <path_to_pod5_dir> --ref <path_to_ref_fasta>

Naming Convention:
    BAM files must follow: {exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam

Output:
    Input/
    ├── {exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam -> source BAM
    ├── {exp_id}_pod5/ -> source Pod5 directory
    └── {exp_id}.fa -> source reference FASTA
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple


# Modification name to bitflag mapping
MOD_BITFLAGS = {
    "non": 0,
    "6mA": 1,
    "5mCG_5hmCG": 2,
    "5mC_5hmC": 4,
    "4mC_5mC": 8,
    "5mC": 16,
}


def parse_bam_filename(bam_path: str) -> Optional[Dict]:
    """
    Parse metadata from BAM filename.

    Expected format: {exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam

    Args:
        bam_path: Path to the BAM file

    Returns:
        Dictionary with parsed metadata or None if parsing fails
    """
    filename = os.path.basename(bam_path)

    # Remove .bam extension
    if not filename.endswith(".bam"):
        print(f"Error: File must have .bam extension: {filename}", file=sys.stderr)
        return None

    base = filename[:-4]  # Remove .bam

    # Pattern: {exp_id}_{model_tier}_v{model_ver}_{trim}_{mods}
    # exp_id can contain underscores, so we parse from the end
    # model_tier is single char: s, h, or f
    # model_ver is like 5.2.0
    # trim is 0 or 1
    # mods is the modification string

    # Try regex pattern
    # This pattern handles exp_id with underscores by using non-greedy match
    pattern = r"^(.+)_([shf])_v([\d.]+)_([01])_(.+)$"
    match = re.match(pattern, base)

    if not match:
        print(f"Error: BAM filename does not match expected pattern: {filename}", file=sys.stderr)
        print("Expected: {exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam")
        return None

    exp_id, model_tier, model_ver, trim, mods_str = match.groups()

    # Parse modification string to bitflag
    # Modifications can be combined with '+' (e.g., "6mA+5mC_5hmC")
    mod_bitflag = 0
    if mods_str != "non":
        mod_parts = mods_str.split("+")
        for mod in mod_parts:
            mod = mod.strip()
            if mod in MOD_BITFLAGS:
                mod_bitflag += MOD_BITFLAGS[mod]
            else:
                print(f"Warning: Unknown modification '{mod}', treating as 0", file=sys.stderr)

    return {
        "exp_id": exp_id,
        "model_tier": model_tier,
        "model_ver": model_ver,
        "trim": int(trim),
        "mods_str": mods_str,
        "mod_bitflag": mod_bitflag,
        "original_filename": filename,
    }


def create_symlink(source: str, target: str, force: bool = False) -> bool:
    """
    Create a symbolic link from target to source.

    Args:
        source: Source file/directory path
        target: Target symlink path
        force: If True, remove existing target before creating

    Returns:
        True if successful, False otherwise
    """
    source_path = Path(source).resolve()
    target_path = Path(target)

    if not source_path.exists():
        print(f"Error: Source does not exist: {source}", file=sys.stderr)
        return False

    if target_path.exists() or target_path.is_symlink():
        if force:
            target_path.unlink()
        else:
            print(f"Warning: Target already exists: {target}")
            return True

    try:
        target_path.symlink_to(source_path)
        return True
    except OSError as e:
        print(f"Error creating symlink: {e}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Standardize input file paths and extract BAM metadata",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python inputInit.py --bam data/EXP001_h_v5.2.0_1_6mA.bam --pod5 data/pod5/ --ref data/ref.fa
    python inputInit.py -b raw/sample.bam -p raw/pod5 -r raw/reference.fasta --force
        """
    )
    parser.add_argument(
        "--bam", "-b",
        type=str,
        required=True,
        help="Path to raw uBAM file (must follow naming convention)"
    )
    parser.add_argument(
        "--pod5", "-p",
        type=str,
        required=True,
        help="Path to raw Pod5 directory"
    )
    parser.add_argument(
        "--ref", "-r",
        type=str,
        required=True,
        help="Path to reference FASTA file"
    )
    parser.add_argument(
        "--force", "-f",
        action="store_true",
        help="Force overwrite existing symlinks"
    )
    parser.add_argument(
        "--output-json", "-o",
        type=str,
        help="Optional: Save parsed metadata to JSON file"
    )

    args = parser.parse_args()

    # Validate inputs exist
    if not os.path.exists(args.bam):
        print(f"Error: BAM file not found: {args.bam}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(args.pod5):
        print(f"Error: Pod5 directory not found: {args.pod5}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.ref):
        print(f"Error: Reference FASTA not found: {args.ref}", file=sys.stderr)
        sys.exit(1)

    # Parse BAM filename for metadata
    metadata = parse_bam_filename(args.bam)
    if metadata is None:
        sys.exit(1)

    exp_id = metadata["exp_id"]
    print(f"Parsed experiment ID: {exp_id}")
    print(f"  Model: {metadata['model_tier']}_v{metadata['model_ver']}")
    print(f"  Trim: {metadata['trim']}")
    print(f"  Modifications: {metadata['mods_str']} (bitflag: {metadata['mod_bitflag']})")

    # Create Input/ directory
    input_dir = Path("Input")
    input_dir.mkdir(exist_ok=True)

    # Create symlinks
    success = True

    # BAM symlink (keep original filename for metadata preservation)
    bam_target = input_dir / metadata["original_filename"]
    if create_symlink(args.bam, str(bam_target), args.force):
        print(f"  Created: {bam_target}")
    else:
        success = False

    # Pod5 directory symlink
    pod5_target = input_dir / f"{exp_id}_pod5"
    if create_symlink(args.pod5, str(pod5_target), args.force):
        print(f"  Created: {pod5_target}")
    else:
        success = False

    # Reference FASTA symlink
    ref_target = input_dir / f"{exp_id}.fa"
    if create_symlink(args.ref, str(ref_target), args.force):
        print(f"  Created: {ref_target}")
    else:
        success = False

    # Optionally save metadata to JSON
    if args.output_json:
        metadata["input_dir"] = str(input_dir.resolve())
        metadata["bam_path"] = str(bam_target.resolve())
        metadata["pod5_path"] = str(pod5_target.resolve())
        metadata["ref_path"] = str(ref_target.resolve())

        with open(args.output_json, "w") as f:
            json.dump(metadata, f, indent=2)
        print(f"  Metadata saved: {args.output_json}")

    if success:
        print("\nInput standardization complete!")
        print(f"Next step: python ingest.py {exp_id}")
    else:
        print("\nInput standardization completed with errors.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
