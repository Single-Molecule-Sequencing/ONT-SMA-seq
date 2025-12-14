#!/usr/bin/env python3
"""
Input Standardization Script for ONT-SMA-seq Pipeline

Standardizes input file paths and extracts metadata from BAM filename.
"""

import argparse
import os
import re
import sys
from pathlib import Path


def parse_bam_filename(bam_filename):
    """
    Parse metadata from BAM filename.
    
    Expected format: {exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam
    
    Args:
        bam_filename: Name of the BAM file (without path)
    
    Returns:
        Dictionary with parsed metadata
    """
    # Remove .bam extension
    name = bam_filename
    if name.endswith('.bam'):
        name = name[:-4]
    
    # Parse using regex pattern
    # Pattern: {exp_id}_{model_tier}_v{model_ver}_{trim}_{modifications}
    pattern = r'^(.+?)_([shf])_v([\d.]+)_([01])_(.+)$'
    match = re.match(pattern, name)
    
    if not match:
        raise ValueError(
            f'BAM filename does not match expected format: '
            f'{{exp_id}}_{{bc_model_type}}_v{{bc_model_version}}_{{trim}}_{{modifications}}.bam\n'
            f'Got: {bam_filename}'
        )
    
    exp_id, model_tier, model_ver, trim, modifications = match.groups()
    
    return {
        'exp_id': exp_id,
        'model_tier': model_tier,
        'model_ver': model_ver,
        'trim': trim,
        'modifications': modifications
    }


def create_symlink(source, target, force=False):
    """
    Create a symbolic link.
    
    Args:
        source: Source file/directory path
        target: Target symlink path
        force: If True, overwrite existing symlink
    """
    source_path = Path(source).resolve()
    target_path = Path(target)
    
    if not source_path.exists():
        raise FileNotFoundError(f'Source does not exist: {source}')
    
    # Remove existing symlink if force is True
    if target_path.exists() or target_path.is_symlink():
        if force:
            target_path.unlink()
        else:
            raise FileExistsError(f'Target already exists: {target}')
    
    # Create symlink
    target_path.symlink_to(source_path)


def main():
    """Main entry point for input standardization."""
    parser = argparse.ArgumentParser(
        description='Standardize input files for ONT-SMA-seq pipeline'
    )
    parser.add_argument(
        'bam',
        type=str,
        help='Path to raw uBAM file (must follow naming convention)'
    )
    parser.add_argument(
        'pod5_dir',
        type=str,
        help='Path to raw Pod5 directory'
    )
    parser.add_argument(
        'reference',
        type=str,
        help='Path to reference FASTA file'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='Input',
        help='Output directory for standardized inputs (default: Input)'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite existing symlinks'
    )
    
    args = parser.parse_args()
    
    # Get BAM filename (without path)
    bam_filename = Path(args.bam).name
    
    try:
        # Parse metadata from BAM filename
        metadata = parse_bam_filename(bam_filename)
        exp_id = metadata['exp_id']
        
        print(f'Parsed metadata from BAM filename:')
        print(f'  Experiment ID: {exp_id}')
        print(f'  Model Tier: {metadata["model_tier"]}')
        print(f'  Model Version: {metadata["model_ver"]}')
        print(f'  Trim: {metadata["trim"]}')
        print(f'  Modifications: {metadata["modifications"]}')
        
        # Create Input directory if it doesn't exist
        input_dir = Path(args.output_dir)
        input_dir.mkdir(parents=True, exist_ok=True)
        print(f'\nCreated directory: {input_dir}')
        
        # Create symlink for BAM file
        bam_target = input_dir / bam_filename
        create_symlink(args.bam, bam_target, force=args.force)
        print(f'Created symlink: {bam_target} -> {args.bam}')
        
        # Create symlink for Pod5 directory
        pod5_target = input_dir / f'{exp_id}_pod5'
        create_symlink(args.pod5_dir, pod5_target, force=args.force)
        print(f'Created symlink: {pod5_target} -> {args.pod5_dir}')
        
        # Create symlink for reference FASTA
        ref_target = input_dir / f'{exp_id}.fa'
        create_symlink(args.reference, ref_target, force=args.force)
        print(f'Created symlink: {ref_target} -> {args.reference}')
        
        print(f'\nSuccessfully standardized inputs in {input_dir}')
        return 0
        
    except ValueError as e:
        print(f'Error parsing filename: {e}', file=sys.stderr)
        return 1
    except FileNotFoundError as e:
        print(f'File not found: {e}', file=sys.stderr)
        return 1
    except FileExistsError as e:
        print(f'File exists: {e}', file=sys.stderr)
        print('Use --force to overwrite existing symlinks', file=sys.stderr)
        return 1
    except Exception as e:
        print(f'Error: {e}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
