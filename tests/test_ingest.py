"""Tests for ingest.py - Main Processing Script"""

import tempfile
from pathlib import Path
import sys
import os

# Add parent directory to path to import ingest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import ingest


def test_calculate_q_bc():
    """Test basecall quality calculation."""
    # Perfect quality scores (all 40)
    qualities = [40, 40, 40, 40]
    q_bc = ingest.calculate_q_bc(qualities)
    assert abs(q_bc - 40.0) < 0.1
    
    # Mixed quality scores
    qualities = [20, 30, 40]
    q_bc = ingest.calculate_q_bc(qualities)
    assert 20 < q_bc < 40
    
    # Empty list
    q_bc = ingest.calculate_q_bc([])
    assert q_bc == 0.0


def test_calculate_levenshtein():
    """Test Levenshtein distance calculation."""
    # Identical sequences
    ed = ingest.calculate_levenshtein('ACGT', 'ACGT')
    assert ed == 0
    
    # One substitution
    ed = ingest.calculate_levenshtein('ACGT', 'ACTT')
    assert ed == 1
    
    # One insertion
    ed = ingest.calculate_levenshtein('ACGT', 'ACGTT')
    assert ed == 1
    
    # One deletion
    ed = ingest.calculate_levenshtein('ACGT', 'ACT')
    assert ed == 1
    
    # Multiple differences
    ed = ingest.calculate_levenshtein('ACGT', 'TGCA')
    assert ed == 4


def test_calculate_q_ld():
    """Test Levenshtein quality calculation."""
    # Perfect match (ed=0) - bounded by 1/L^2
    q_ld = ingest.calculate_q_ld(0, 100)
    assert q_ld == 40.0  # -10 * log10(1/10000) = 40
    
    # 1 error in 100 bases
    q_ld = ingest.calculate_q_ld(1, 100)
    assert q_ld > 0
    assert q_ld < 40
    
    # High error rate
    q_ld = ingest.calculate_q_ld(50, 100)
    assert q_ld > 0
    assert q_ld < 10


def test_match_reference():
    """Test reference matching based on read length."""
    references = {
        'short': {
            'seq': 'A' * 500,
            'len': 500,
            'range': (350, 650)  # 500 +/- 150
        },
        'long': {
            'seq': 'A' * 1000,
            'len': 1000,
            'range': (850, 1150)  # 1000 +/- 150
        }
    }
    
    # Read matches short reference
    assert ingest.match_reference(500, references) == 'short'
    assert ingest.match_reference(400, references) == 'short'
    assert ingest.match_reference(600, references) == 'short'
    
    # Read matches long reference
    assert ingest.match_reference(1000, references) == 'long'
    assert ingest.match_reference(900, references) == 'long'
    assert ingest.match_reference(1100, references) == 'long'
    
    # Read doesn't match any reference
    assert ingest.match_reference(200, references) is None
    assert ingest.match_reference(2000, references) is None


def test_generate_uniq_id():
    """Test unique ID generation."""
    uniq_id = ingest.generate_uniq_id(
        'exp001', 's', '5.2.0', 1, 0, 'read_12345'
    )
    
    assert uniq_id.startswith('exp001')
    assert 's' in uniq_id
    assert '520' in uniq_id  # Version without dots
    assert 't1' in uniq_id  # Trim
    assert 'm0' in uniq_id  # Mod flag
    assert '_' in uniq_id
    
    # Different read IDs should produce different unique IDs
    uniq_id2 = ingest.generate_uniq_id(
        'exp001', 's', '5.2.0', 1, 0, 'read_67890'
    )
    assert uniq_id != uniq_id2


def test_parse_bam_filename():
    """Test BAM filename parsing."""
    metadata = ingest.parse_bam_filename('exp001_s_v5.2.0_1_6mA.bam')
    
    assert metadata['exp_id'] == 'exp001'
    assert metadata['model_tier'] == 's'
    assert metadata['model_ver'] == '5.2.0'
    assert metadata['trim'] == 1
    assert metadata['modifications'] == '6mA'


def test_get_mod_bitflag():
    """Test modification bitflag conversion."""
    assert ingest.get_mod_bitflag('non') == 0
    assert ingest.get_mod_bitflag('6mA') == 1
    assert ingest.get_mod_bitflag('5mCG_5hmCG') == 2
    assert ingest.get_mod_bitflag('5mC_5hmC') == 4
    assert ingest.get_mod_bitflag('4mC_5mC') == 8
    assert ingest.get_mod_bitflag('5mC') == 16
    
    # Test combination
    assert ingest.get_mod_bitflag('6mA+5mC_5hmC') == 5  # 1 + 4


def test_parse_reference_fasta():
    """Test reference FASTA parsing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a test FASTA file with 2 sequences
        fasta_path = Path(tmpdir) / 'test.fa'
        fasta_content = '>short\n' + 'A' * 500 + '\n>long\n' + 'C' * 1000 + '\n'
        fasta_path.write_text(fasta_content)
        
        references = ingest.parse_reference_fasta(fasta_path)
        
        assert len(references) == 2
        assert 'short' in references
        assert 'long' in references
        
        assert references['short']['len'] == 500
        assert references['short']['range'] == (350, 650)
        
        assert references['long']['len'] == 1000
        assert references['long']['range'] == (850, 1150)


def test_parse_reference_fasta_wrong_count():
    """Test reference FASTA with wrong number of sequences raises error."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a FASTA with only 1 sequence
        fasta_path = Path(tmpdir) / 'test.fa'
        fasta_content = '>only_one\nACGT\n'
        fasta_path.write_text(fasta_content)
        
        try:
            ingest.parse_reference_fasta(fasta_path)
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert 'Expected exactly 2 reference sequences' in str(e)


if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-v'])
