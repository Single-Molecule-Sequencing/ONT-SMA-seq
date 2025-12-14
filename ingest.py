#!/usr/bin/env python3
"""
Main Processing Script for ONT-SMA-seq Pipeline

Parses inputs, calculates metrics, tags BAMs, and populates the database.
"""

import argparse
import hashlib
import math
import re
import sqlite3
import sys
from pathlib import Path

try:
    import pysam
except ImportError:
    print('Error: pysam is not installed. Install with: pip install pysam', file=sys.stderr)
    sys.exit(1)

try:
    import edlib
except ImportError:
    print('Error: edlib is not installed. Install with: pip install edlib', file=sys.stderr)
    sys.exit(1)

try:
    import pod5 as p5
except ImportError:
    print('Error: pod5 is not installed. Install with: pip install pod5', file=sys.stderr)
    sys.exit(1)


def parse_reference_fasta(fasta_path):
    """
    Parse reference FASTA file.
    
    Expected to contain exactly 2 sequences (long and short).
    Returns dictionary: {defline: [seq, len, range]}
    """
    references = {}
    
    with pysam.FastaFile(str(fasta_path)) as fasta:
        for ref_name in fasta.references:
            seq = fasta.fetch(ref_name)
            seq_len = len(seq)
            # Calculate refseqrange: (Len_seq - 150, Len_seq + 150)
            seq_range = (seq_len - 150, seq_len + 150)
            references[ref_name] = {
                'seq': seq,
                'len': seq_len,
                'range': seq_range
            }
    
    if len(references) != 2:
        raise ValueError(f'Expected exactly 2 reference sequences, found {len(references)}')
    
    return references


def parse_pod5_metadata(pod5_dir):
    """
    Parse Pod5 files to extract end reason metadata.
    
    Returns dictionary: {read_id: end_reason_string}
    """
    end_reasons = {}
    pod5_path = Path(pod5_dir)
    
    if not pod5_path.exists():
        raise FileNotFoundError(f'Pod5 directory not found: {pod5_dir}')
    
    # Find all .pod5 files in directory
    pod5_files = list(pod5_path.glob('*.pod5'))
    
    if not pod5_files:
        raise FileNotFoundError(f'No .pod5 files found in {pod5_dir}')
    
    for pod5_file in pod5_files:
        with p5.Reader(str(pod5_file)) as reader:
            for read in reader.reads():
                end_reason = read.end_reason.name if read.end_reason else 'unknown'
                end_reasons[str(read.read_id)] = end_reason
    
    return end_reasons


def calculate_q_bc(qualities):
    """
    Calculate probability-averaged basecall quality (q_bc).
    
    q_bc = -10 * log10(sum(10^(-Q_base/10)) / n)
    """
    if not qualities or len(qualities) == 0:
        return 0.0
    
    # Convert quality scores to error probabilities
    error_probs = [10 ** (-q / 10.0) for q in qualities]
    # Calculate average error probability
    avg_error_prob = sum(error_probs) / len(error_probs)
    # Convert back to quality score
    if avg_error_prob > 0:
        q_bc = -10 * math.log10(avg_error_prob)
    else:
        q_bc = 60.0  # Maximum quality
    
    return q_bc


def calculate_levenshtein(read_seq, ref_seq):
    """
    Calculate Levenshtein distance between read and reference.
    """
    result = edlib.align(read_seq, ref_seq, task='distance')
    return result['editDistance']


def calculate_q_ld(edit_distance, ref_len):
    """
    Calculate Levenshtein quality (q_ld).
    
    q_ld = -10 * log10(min(max(1/L^2, ed/L), 1))
    """
    if ref_len == 0:
        return 0.0
    
    # Calculate error rate
    error_rate = edit_distance / ref_len
    # Apply min/max bounds
    min_error = 1.0 / (ref_len ** 2)
    bounded_error = min(max(min_error, error_rate), 1.0)
    # Convert to quality score
    if bounded_error > 0:
        q_ld = -10 * math.log10(bounded_error)
    else:
        q_ld = 60.0  # Maximum quality
    
    return q_ld


def match_reference(read_len, references):
    """
    Match read to reference based on length range.
    
    Returns refseq_id if match found, None otherwise.
    """
    for ref_id, ref_data in references.items():
        min_len, max_len = ref_data['range']
        if min_len <= read_len <= max_len:
            return ref_id
    return None


def generate_uniq_id(exp_id, model_tier, model_ver, trim, mod_flag, read_id):
    """
    Generate unique ID for read.
    
    Format: {exp_id}{tier}{ver}t{trim}m{mod_flag}_{read_hash}
    """
    # Create a short hash of the read_id
    hash_obj = hashlib.md5(read_id.encode())
    read_hash = hash_obj.hexdigest()[:8]
    
    # Format version without dots
    ver_str = model_ver.replace('.', '')
    
    uniq_id = f'{exp_id}{model_tier}{ver_str}t{trim}m{mod_flag}_{read_hash}'
    return uniq_id


def parse_bam_filename(bam_filename):
    """Parse metadata from BAM filename."""
    name = bam_filename
    if name.endswith('.bam'):
        name = name[:-4]
    
    pattern = r'^(.+?)_([shf])_v([\d.]+)_([01])_(.+)$'
    match = re.match(pattern, name)
    
    if not match:
        raise ValueError(f'BAM filename does not match expected format: {bam_filename}')
    
    exp_id, model_tier, model_ver, trim, modifications = match.groups()
    
    return {
        'exp_id': exp_id,
        'model_tier': model_tier,
        'model_ver': model_ver,
        'trim': int(trim),
        'modifications': modifications
    }


def get_mod_bitflag(modifications):
    """
    Convert modification string to bitflag.
    
    Simple mapping for now - can be extended.
    """
    mod_map = {
        'non': 0,
        '6mA': 1,
        '5mCG_5hmCG': 2,
        '5mC_5hmC': 4,
        '4mC_5mC': 8,
        '5mC': 16,
    }
    
    # Handle combinations (e.g., "6mA+5mC_5hmC")
    if '+' in modifications:
        parts = modifications.split('+')
        bitflag = 0
        for part in parts:
            bitflag += mod_map.get(part, 0)
        return bitflag
    
    return mod_map.get(modifications, 0)


def main():
    """Main entry point for ingestion."""
    parser = argparse.ArgumentParser(
        description='Process reads and populate ONT-SMA-seq database'
    )
    parser.add_argument(
        'exp_id',
        type=str,
        help='Experiment ID'
    )
    parser.add_argument(
        '--input-dir',
        type=str,
        default='Input',
        help='Input directory with standardized files (default: Input)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='Output',
        help='Output directory for tagged BAM (default: Output)'
    )
    parser.add_argument(
        '--db-dir',
        type=str,
        default='.',
        help='Directory containing database file (default: current directory)'
    )
    
    args = parser.parse_args()
    
    exp_id = args.exp_id
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    db_path = Path(args.db_dir) / f'SMA_{exp_id}.db'
    
    # Check database exists
    if not db_path.exists():
        print(f'Error: Database not found: {db_path}', file=sys.stderr)
        print('Run mkdb.py first to create the database', file=sys.stderr)
        return 1
    
    try:
        # Find input BAM file
        bam_files = list(input_dir.glob(f'{exp_id}_*.bam'))
        if not bam_files:
            raise FileNotFoundError(f'No BAM file found for experiment {exp_id} in {input_dir}')
        bam_path = bam_files[0]
        
        # Parse metadata from BAM filename
        metadata = parse_bam_filename(bam_path.name)
        model_tier = metadata['model_tier']
        model_ver = metadata['model_ver']
        trim = metadata['trim']
        mod_bitflag = get_mod_bitflag(metadata['modifications'])
        
        print(f'Processing experiment: {exp_id}')
        print(f'Input BAM: {bam_path}')
        print(f'Model: {model_tier} v{model_ver}, Trim: {trim}, Modifications: {mod_bitflag}')
        
        # Parse reference FASTA
        ref_path = input_dir / f'{exp_id}.fa'
        if not ref_path.exists():
            raise FileNotFoundError(f'Reference FASTA not found: {ref_path}')
        
        print(f'Parsing reference: {ref_path}')
        references = parse_reference_fasta(ref_path)
        print(f'Found {len(references)} reference sequences')
        
        # Connect to database
        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()
        
        # Insert reference sequences into database
        for ref_id, ref_data in references.items():
            cursor.execute(
                'INSERT OR IGNORE INTO Refseq (refseq_id, refseq, reflen) VALUES (?, ?, ?)',
                (ref_id, ref_data['seq'], ref_data['len'])
            )
        conn.commit()
        print('Inserted reference sequences into database')
        
        # Parse Pod5 metadata
        pod5_dir = input_dir / f'{exp_id}_pod5'
        print(f'Parsing Pod5 metadata from: {pod5_dir}')
        end_reasons = parse_pod5_metadata(pod5_dir)
        print(f'Loaded end reasons for {len(end_reasons)} reads')
        
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        output_bam_path = output_dir / f'{exp_id}.bam'
        
        # Process reads
        print(f'Processing reads...')
        reads_processed = 0
        reads_matched = 0
        
        with pysam.AlignmentFile(str(bam_path), 'rb', check_sq=False) as inbam:
            with pysam.AlignmentFile(str(output_bam_path), 'wb', header=inbam.header) as outbam:
                for read in inbam:
                    read_id = read.query_name
                    read_seq = read.query_sequence
                    read_len = len(read_seq) if read_seq else 0
                    
                    if not read_seq:
                        continue  # Skip reads without sequence
                    
                    # Get end reason from Pod5 data
                    end_reason = end_reasons.get(read_id, 'unknown')
                    
                    # Add ER tag to read
                    read.set_tag('ER', end_reason, value_type='Z')
                    
                    # Write tagged read to output BAM
                    outbam.write(read)
                    
                    # Match to reference
                    refseq_id = match_reference(read_len, references)
                    
                    # Calculate quality metrics
                    qualities = read.query_qualities
                    if qualities is not None:
                        q_bc = calculate_q_bc(qualities)
                    else:
                        q_bc = 0.0
                    
                    # Calculate Levenshtein metrics if reference matched
                    ed = None
                    q_ld = None
                    if refseq_id:
                        ref_data = references[refseq_id]
                        ed = calculate_levenshtein(read_seq, ref_data['seq'])
                        q_ld = calculate_q_ld(ed, ref_data['len'])
                        reads_matched += 1
                    
                    # Generate unique ID
                    uniq_id = generate_uniq_id(exp_id, model_tier, model_ver, trim, mod_bitflag, read_id)
                    
                    # Extract additional metadata from read tags
                    forced = read.get_tag('pt') if read.has_tag('pt') else None
                    channel = read.get_tag('ch') if read.has_tag('ch') else None
                    well = read.get_tag('ws') if read.has_tag('ws') else None
                    pore_type = read.get_tag('rn') if read.has_tag('rn') else None
                    num_samples = read.get_tag('ns') if read.has_tag('ns') else None
                    start_sample = read.get_tag('ts') if read.has_tag('ts') else None
                    median_before = read.get_tag('sm') if read.has_tag('sm') else None
                    scale = read.get_tag('sd') if read.has_tag('sd') else None
                    offset = read.get_tag('sv') if read.has_tag('sv') else None
                    
                    # Insert into database
                    cursor.execute('''
                        INSERT INTO Reads (
                            uniq_id, exp_id, refseq_id, read_id, readseq, readlen,
                            model_tier, model_ver, trim, mod_bitflag,
                            ed, q_bc, q_ld, ER,
                            forced, channel, well, pore_type,
                            num_samples, start_sample, median_before, scale, offset
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ''', (
                        uniq_id, exp_id, refseq_id, read_id, read_seq, read_len,
                        model_tier, model_ver, trim, mod_bitflag,
                        ed, q_bc, q_ld, end_reason,
                        forced, channel, well, pore_type,
                        num_samples, start_sample, median_before, scale, offset
                    ))
                    
                    reads_processed += 1
                    if reads_processed % 1000 == 0:
                        conn.commit()
                        print(f'Processed {reads_processed} reads...')
        
        # Final commit
        conn.commit()
        conn.close()
        
        print(f'\nProcessing complete!')
        print(f'Total reads processed: {reads_processed}')
        print(f'Reads matched to reference: {reads_matched}')
        print(f'Output BAM: {output_bam_path}')
        print(f'Database updated: {db_path}')
        
        return 0
        
    except FileNotFoundError as e:
        print(f'File not found: {e}', file=sys.stderr)
        return 1
    except Exception as e:
        print(f'Error: {e}', file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
