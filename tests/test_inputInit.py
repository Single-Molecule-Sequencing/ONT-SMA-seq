"""Tests for inputInit.py - Input Standardization Script"""

import tempfile
from pathlib import Path
import sys
import os

# Add parent directory to path to import inputInit
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import inputInit


def test_parse_bam_filename_valid():
    """Test parsing valid BAM filename."""
    filename = 'exp001_s_v5.2.0_1_6mA.bam'
    metadata = inputInit.parse_bam_filename(filename)
    
    assert metadata['exp_id'] == 'exp001'
    assert metadata['model_tier'] == 's'
    assert metadata['model_ver'] == '5.2.0'
    assert metadata['trim'] == '1'
    assert metadata['modifications'] == '6mA'


def test_parse_bam_filename_complex():
    """Test parsing BAM filename with complex experiment ID."""
    filename = 'my_complex_exp_h_v4.0.0_0_5mC_5hmC.bam'
    metadata = inputInit.parse_bam_filename(filename)
    
    assert metadata['exp_id'] == 'my_complex_exp'
    assert metadata['model_tier'] == 'h'
    assert metadata['model_ver'] == '4.0.0'
    assert metadata['trim'] == '0'
    assert metadata['modifications'] == '5mC_5hmC'


def test_parse_bam_filename_no_mods():
    """Test parsing BAM filename with no modifications."""
    filename = 'exp_f_v5.0.0_1_non.bam'
    metadata = inputInit.parse_bam_filename(filename)
    
    assert metadata['exp_id'] == 'exp'
    assert metadata['model_tier'] == 'f'
    assert metadata['model_ver'] == '5.0.0'
    assert metadata['trim'] == '1'
    assert metadata['modifications'] == 'non'


def test_parse_bam_filename_invalid():
    """Test parsing invalid BAM filename raises ValueError."""
    filename = 'invalid_filename.bam'
    
    try:
        inputInit.parse_bam_filename(filename)
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert 'does not match expected format' in str(e)


def test_create_symlink():
    """Test creating symbolic links."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create source file
        source = Path(tmpdir) / 'source.txt'
        source.write_text('test content')
        
        # Create symlink
        target = Path(tmpdir) / 'target.txt'
        inputInit.create_symlink(source, target)
        
        assert target.exists()
        assert target.is_symlink()
        assert target.read_text() == 'test content'


def test_create_symlink_directory():
    """Test creating symbolic link to directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create source directory
        source = Path(tmpdir) / 'source_dir'
        source.mkdir()
        (source / 'file.txt').write_text('test')
        
        # Create symlink
        target = Path(tmpdir) / 'target_dir'
        inputInit.create_symlink(source, target)
        
        assert target.exists()
        assert target.is_symlink()
        assert (target / 'file.txt').read_text() == 'test'


def test_create_symlink_source_not_exists():
    """Test creating symlink with non-existent source raises FileNotFoundError."""
    with tempfile.TemporaryDirectory() as tmpdir:
        source = Path(tmpdir) / 'nonexistent.txt'
        target = Path(tmpdir) / 'target.txt'
        
        try:
            inputInit.create_symlink(source, target)
            assert False, "Should have raised FileNotFoundError"
        except FileNotFoundError:
            pass


def test_create_symlink_target_exists():
    """Test creating symlink when target exists raises FileExistsError."""
    with tempfile.TemporaryDirectory() as tmpdir:
        source = Path(tmpdir) / 'source.txt'
        source.write_text('test')
        
        target = Path(tmpdir) / 'target.txt'
        target.write_text('existing')
        
        try:
            inputInit.create_symlink(source, target)
            assert False, "Should have raised FileExistsError"
        except FileExistsError:
            pass


def test_create_symlink_force_overwrite():
    """Test creating symlink with force=True overwrites existing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        source = Path(tmpdir) / 'source.txt'
        source.write_text('new content')
        
        target = Path(tmpdir) / 'target.txt'
        target.write_text('old content')
        
        inputInit.create_symlink(source, target, force=True)
        
        assert target.is_symlink()
        assert target.read_text() == 'new content'


def test_main_creates_symlinks():
    """Test main function creates all necessary symlinks."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create source files
        bam_source = Path(tmpdir) / 'exp001_s_v5.2.0_1_6mA.bam'
        bam_source.write_text('BAM content')
        
        pod5_source = Path(tmpdir) / 'pod5_data'
        pod5_source.mkdir()
        (pod5_source / 'data.pod5').write_text('POD5 content')
        
        ref_source = Path(tmpdir) / 'reference.fa'
        ref_source.write_text('>ref\nACGT')
        
        output_dir = Path(tmpdir) / 'Input'
        
        # Mock sys.argv
        sys.argv = [
            'inputInit.py',
            str(bam_source),
            str(pod5_source),
            str(ref_source),
            '--output-dir', str(output_dir)
        ]
        
        result = inputInit.main()
        
        assert result == 0
        assert output_dir.exists()
        
        # Check symlinks exist
        bam_link = output_dir / 'exp001_s_v5.2.0_1_6mA.bam'
        assert bam_link.exists()
        assert bam_link.is_symlink()
        
        pod5_link = output_dir / 'exp001_pod5'
        assert pod5_link.exists()
        assert pod5_link.is_symlink()
        
        ref_link = output_dir / 'exp001.fa'
        assert ref_link.exists()
        assert ref_link.is_symlink()


if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-v'])
