"""Tests for sma_merge data models."""
from __future__ import annotations


def test_run_info_creation():
    from sma_merge.models import RunInfo
    from pathlib import Path

    info = RunInfo(
        run_dir=Path("/tmp/run"),
        flow_cell_id="FBD66244",
        device_id="MD-101527",
        protocol_group_id="test_protocol",
        basecall_model="dna_r10.4.1_e8.2_400bps_hac@v5.2.0",
        sample_id="sample1",
        run_id="abc123",
        sample_rate=5000,
        pod5_dir=Path("/tmp/run/pod5_pass"),
        pod5_count=10,
        mod_base_models="",
    )
    assert info.flow_cell_id == "FBD66244"
    assert info.basecall_model == "dna_r10.4.1_e8.2_400bps_hac@v5.2.0"


def test_run_group_creation():
    from sma_merge.models import RunInfo, RunGroup
    from pathlib import Path

    run = RunInfo(
        run_dir=Path("/tmp/run"),
        flow_cell_id="FBD66244",
        device_id="MD-101527",
        protocol_group_id="test_protocol",
        basecall_model="dna_r10.4.1_e8.2_400bps_hac@v5.2.0",
        sample_id="",
        run_id="abc123",
        sample_rate=5000,
        pod5_dir=Path("/tmp/run/pod5_pass"),
        pod5_count=10,
        mod_base_models="",
    )
    group = RunGroup(
        flow_cell_id="FBD66244",
        runs=[run],
        basecall_model="dna_r10.4.1_e8.2_400bps_hac@v5.2.0",
        is_consistent=True,
        issues=[],
    )
    assert group.is_consistent
    assert len(group.runs) == 1


def test_merge_result_creation():
    from sma_merge.models import MergeResult
    from pathlib import Path

    result = MergeResult(
        merged_pod5=Path("/tmp/out/merged.pod5"),
        output_bam=Path("/tmp/out/merged.bam"),
        total_reads=5000,
        reads_tagged=4950,
    )
    assert result.total_reads == 5000
