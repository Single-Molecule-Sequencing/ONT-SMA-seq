"""Tests for run grouping and BAM merge logic."""
from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from calibrate.discover import RunInfo
from calibrate.merge import (
    MergeGroup,
    ValidationError,
    find_bam_files,
    group_runs,
    merge_bams,
    validate_basecalling,
)


def _make_run_info(
    tmp_path: Path,
    *,
    flow_cell_id: str = "FBD69411",
    sample_id: str = "no_sample_id",
    started: str = "2025-12-28T22:19:00Z",
    suffix: str = "",
) -> RunInfo:
    """Create a RunInfo with a real directory on disk."""
    run_dir = tmp_path / f"run_{flow_cell_id}_{sample_id}_{started.replace(':', '')}_{suffix}"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "bam_pass").mkdir(exist_ok=True)
    fs_path = run_dir / "final_summary.txt"
    fs_path.write_text(
        f"flow_cell_id={flow_cell_id}\n"
        f"sample_id={sample_id}\n"
        f"started={started}\n"
        f"device_id=MD-100098\n"
        f"experiment_id=exp1\n"
        f"protocol_run_id=prid-{suffix or '0'}\n"
    )
    return RunInfo(
        run_dir=run_dir,
        flow_cell_id=flow_cell_id,
        device_id="MD-100098",
        sample_id=sample_id,
        experiment_id="exp1",
        started=started,
        protocol_run_id=f"prid-{suffix or '0'}",
        final_summary_path=fs_path,
        sequencing_summary_path=None,
        sample_sheet_path=None,
    )


def _write_tiny_bam(path: Path, model: str = "dna_r10.4.1_e8.2_400bps_sup@v4.3.0", trim: bool = False) -> None:
    """Write a minimal unaligned BAM with a @PG header."""
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unknown"},
        "PG": [
            {
                "ID": "basecaller",
                "PN": "dorado",
                "CL": f"dorado basecaller --model {model}"
                + (" --trim" if trim else ""),
            }
        ],
    })
    with pysam.AlignmentFile(str(path), "wb", header=header) as outf:
        seg = pysam.AlignedSegment(outf.header)
        seg.query_name = "read1"
        seg.query_sequence = "ACGT"
        seg.query_qualities = pysam.qualitystring_to_array("IIII")
        seg.flag = 4  # unmapped
        outf.write(seg)


# ---------------------------------------------------------------------------
# TestGroupRuns
# ---------------------------------------------------------------------------

class TestGroupRuns:
    def test_groups_same_flow_cell_and_sample(self, tmp_path):
        """Two runs same FC+sample within 24h -> 1 group with 2 runs."""
        r1 = _make_run_info(tmp_path, started="2025-12-28T22:00:00Z", suffix="a")
        r2 = _make_run_info(tmp_path, started="2025-12-29T10:00:00Z", suffix="b")
        groups = group_runs([r1, r2])
        assert len(groups) == 1
        assert len(groups[0].runs) == 2

    def test_separates_different_flow_cells(self, tmp_path):
        """Different flow_cell_id -> 2 groups."""
        r1 = _make_run_info(tmp_path, flow_cell_id="FC_AAA", started="2025-12-28T22:00:00Z", suffix="a")
        r2 = _make_run_info(tmp_path, flow_cell_id="FC_BBB", started="2025-12-28T23:00:00Z", suffix="b")
        groups = group_runs([r1, r2])
        assert len(groups) == 2
        fc_ids = {g.flow_cell_id for g in groups}
        assert fc_ids == {"FC_AAA", "FC_BBB"}

    def test_separates_runs_more_than_24h_apart(self, tmp_path):
        """Same FC but >24h gap -> 2 groups."""
        r1 = _make_run_info(tmp_path, started="2025-12-28T00:00:00Z", suffix="a")
        r2 = _make_run_info(tmp_path, started="2025-12-30T01:00:00Z", suffix="b")
        groups = group_runs([r1, r2])
        assert len(groups) == 2
        # Each group has 1 run
        assert all(len(g.runs) == 1 for g in groups)

    def test_separates_different_samples(self, tmp_path):
        """Different sample_id -> 2 groups."""
        r1 = _make_run_info(tmp_path, sample_id="sample_A", started="2025-12-28T22:00:00Z", suffix="a")
        r2 = _make_run_info(tmp_path, sample_id="sample_B", started="2025-12-28T23:00:00Z", suffix="b")
        groups = group_runs([r1, r2])
        assert len(groups) == 2
        sample_ids = {g.sample_id for g in groups}
        assert sample_ids == {"sample_A", "sample_B"}


# ---------------------------------------------------------------------------
# TestFindBamFiles
# ---------------------------------------------------------------------------

class TestFindBamFiles:
    def test_finds_flat_bams(self, tmp_path):
        """BAMs directly in bam_pass/."""
        run = _make_run_info(tmp_path, suffix="flat")
        bam_pass = run.run_dir / "bam_pass"
        _write_tiny_bam(bam_pass / "reads1.bam")
        _write_tiny_bam(bam_pass / "reads2.bam")
        bams = find_bam_files(run.run_dir)
        assert len(bams) == 2
        assert all(b.suffix == ".bam" for b in bams)

    def test_finds_per_barcode_bams(self, tmp_path):
        """BAMs in bam_pass/alias/ subdirs."""
        run = _make_run_info(tmp_path, suffix="barcoded")
        bam_pass = run.run_dir / "bam_pass"
        for bc in ("barcode01", "barcode02"):
            bc_dir = bam_pass / bc
            bc_dir.mkdir()
            _write_tiny_bam(bc_dir / f"{bc}_reads.bam")
        bams = find_bam_files(run.run_dir)
        assert len(bams) == 2

    def test_returns_empty_for_no_bams(self, tmp_path):
        """Empty bam_pass/ -> empty list."""
        run = _make_run_info(tmp_path, suffix="empty")
        bams = find_bam_files(run.run_dir)
        assert bams == []


# ---------------------------------------------------------------------------
# TestValidateBasecalling
# ---------------------------------------------------------------------------

class TestValidateBasecalling:
    def test_consistent_model_passes(self, tmp_path):
        model = "dna_r10.4.1_e8.2_400bps_sup@v4.3.0"
        bam1 = tmp_path / "a.bam"
        bam2 = tmp_path / "b.bam"
        _write_tiny_bam(bam1, model=model)
        _write_tiny_bam(bam2, model=model)
        result = validate_basecalling([bam1, bam2])
        assert result == model

    def test_different_models_raises(self, tmp_path):
        bam1 = tmp_path / "a.bam"
        bam2 = tmp_path / "b.bam"
        _write_tiny_bam(bam1, model="dna_r10.4.1_e8.2_400bps_sup@v4.3.0")
        _write_tiny_bam(bam2, model="dna_r10.4.1_e8.2_400bps_hac@v4.3.0")
        with pytest.raises(ValidationError, match="Multiple basecalling models"):
            validate_basecalling([bam1, bam2])

    def test_trim_without_no_trim_raises(self, tmp_path):
        bam = tmp_path / "trimmed.bam"
        _write_tiny_bam(bam, trim=True)
        with pytest.raises(ValidationError, match="--trim"):
            validate_basecalling([bam])


# ---------------------------------------------------------------------------
# TestMergeBams
# ---------------------------------------------------------------------------

class TestMergeBams:
    def test_merges_bams_from_group(self, tmp_path):
        r1 = _make_run_info(tmp_path, started="2025-12-28T22:00:00Z", suffix="m1")
        r2 = _make_run_info(tmp_path, started="2025-12-29T10:00:00Z", suffix="m2")
        _write_tiny_bam(r1.run_dir / "bam_pass" / "reads.bam")
        _write_tiny_bam(r2.run_dir / "bam_pass" / "reads.bam")

        group = MergeGroup(flow_cell_id="FBD69411", sample_id="no_sample_id", runs=[r1, r2])
        out = tmp_path / "merged" / "out.bam"
        prov = merge_bams(group, out)

        assert out.exists()
        assert prov["num_bams"] == 2
        assert prov["num_runs"] == 2
        # Verify we can open the merged BAM
        with pysam.AlignmentFile(str(out), "rb", check_sq=False) as af:
            reads = list(af)
            assert len(reads) == 2
