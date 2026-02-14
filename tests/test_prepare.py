"""Tests for prepare.py — manifest, DB initialization, and CLI."""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest


class TestManifest:

    def test_create_and_save(self, tmp_path):
        from prepare import Manifest
        m = Manifest(exp_id="TEST_EXP", expdir="/fake/path")
        m.mark_stage("discover")
        path = tmp_path / "manifest.json"
        m.save(path)
        assert path.exists()
        data = json.loads(path.read_text())
        assert data["exp_id"] == "TEST_EXP"
        assert "discover" in data["stages_completed"]

    def test_load_existing(self, tmp_path):
        from prepare import Manifest
        m = Manifest(exp_id="TEST_EXP", expdir="/fake/path")
        m.mark_stage("discover")
        m.mark_stage("plan")
        path = tmp_path / "manifest.json"
        m.save(path)

        loaded = Manifest.load(path)
        assert loaded.exp_id == "TEST_EXP"
        assert "plan" in loaded.stages_completed

    def test_load_nonexistent_returns_none(self, tmp_path):
        from prepare import Manifest
        assert Manifest.load(tmp_path / "nope.json") is None


class TestInitDB:

    def test_creates_database_with_schema(self, tmp_path):
        from prepare import init_database
        db_path = tmp_path / "test.db"
        init_database(db_path, "TEST_EXP", "FC001", "S001", "alias1")
        assert db_path.exists()

        conn = sqlite3.connect(str(db_path))
        c = conn.cursor()
        # Check tables exist
        c.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = {row[0] for row in c.fetchall()}
        assert "Reads" in tables
        assert "Mods" in tables
        assert "Exp" in tables
        assert "Target" in tables
        assert "RunMetadata" in tables

        # Check Mods populated
        c.execute("SELECT COUNT(*) FROM Mods")
        assert c.fetchone()[0] == 10

        # Check Exp populated
        c.execute("SELECT exp_id FROM Exp")
        assert c.fetchone()[0] == "TEST_EXP"

        conn.close()

    def test_populate_run_metadata(self, tmp_path):
        from prepare import init_database, insert_run_metadata
        db_path = tmp_path / "test.db"
        init_database(db_path, "TEST_EXP", "FC001", "S001", "alias1")

        insert_run_metadata(db_path, {
            "run_id": "run1",
            "flow_cell_id": "FC001",
            "device_id": "MD-100",
            "sample_id": "S001",
            "experiment_id": "TEST_EXP",
            "kit": "SQK-NBD114-24",
            "basecall_model": "sup@v5.2.0",
            "source_bam_count": 100,
        })

        conn = sqlite3.connect(str(db_path))
        c = conn.cursor()
        c.execute("SELECT run_id, basecall_model FROM RunMetadata")
        row = c.fetchone()
        assert row[0] == "run1"
        assert row[1] == "sup@v5.2.0"
        conn.close()


class TestSymlinkPod5:

    def test_creates_symlinks(self, tmp_path):
        from prepare import symlink_pod5s
        from sma_merge.models import RunInfo

        # Create fake POD5 files
        run_dir = tmp_path / "run1" / "pod5_pass"
        run_dir.mkdir(parents=True)
        (run_dir / "chunk1.pod5").write_text("fake")
        (run_dir / "chunk2.pod5").write_text("fake")

        run = RunInfo(
            run_dir=run_dir.parent, flow_cell_id="FC1", device_id="D1",
            protocol_group_id="p1", basecall_model="sup", sample_id="",
            run_id="run1_abc", sample_rate=5000, pod5_dir=run_dir,
            pod5_count=2, mod_base_models="",
        )

        out = tmp_path / "output" / "pod5"
        symlink_pod5s([run], out)
        links = list(out.rglob("*.pod5"))
        assert len(links) == 2
        assert all(link.is_symlink() for link in links)


class TestCollectBams:

    def test_finds_all_bams(self, tmp_path):
        from prepare import collect_bam_files
        from sma_merge.models import RunInfo

        run_dir = tmp_path / "run1"
        bam_dir = run_dir / "bam_pass" / "barcode01"
        bam_dir.mkdir(parents=True)
        (bam_dir / "chunk1.bam").write_text("fake")
        (bam_dir / "chunk2.bam").write_text("fake")
        unclass = run_dir / "bam_pass" / "unclassified"
        unclass.mkdir(parents=True)
        (unclass / "chunk3.bam").write_text("fake")

        run = RunInfo(
            run_dir=run_dir, flow_cell_id="FC1", device_id="D1",
            protocol_group_id="p1", basecall_model="sup", sample_id="",
            run_id="run1", sample_rate=5000, pod5_dir=run_dir / "pod5_pass",
            pod5_count=0, mod_base_models="",
        )

        bams = collect_bam_files([run])
        assert len(bams) == 3


class TestFormatMergePlan:

    def test_formats_plan(self):
        from prepare import format_merge_plan
        from sma_merge.models import RunGroup, RunInfo
        from pathlib import Path

        run = RunInfo(
            run_dir=Path("/tmp"), flow_cell_id="FC1", device_id="D1",
            protocol_group_id="p1", basecall_model="sup@v5.2.0",
            sample_id="", run_id="run1", sample_rate=5000,
            pod5_dir=Path("/tmp/pod5"), pod5_count=10, mod_base_models="",
        )
        group = RunGroup(flow_cell_id="FC1", runs=[run],
                        basecall_model="sup@v5.2.0", is_consistent=True)

        text = format_merge_plan([group], "TEST_EXP")
        assert "FC1" in text
        assert "REUSE" in text or "sup" in text


class TestBuildParser:

    def test_all_required_args(self):
        from prepare import build_parser
        parser = build_parser()
        # Should fail without required args
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parses_valid_args(self):
        from prepare import build_parser
        parser = build_parser()
        args = parser.parse_args([
            "-d", "/fake/exp", "-e", "FC1_S1_A", "-r", "/fake/ref.fa",
        ])
        assert args.expdir == Path("/fake/exp")
        assert args.expid == "FC1_S1_A"
        assert args.ref == Path("/fake/ref.fa")
        assert args.outdir == Path("Output")  # default

    def test_optional_flags(self):
        from prepare import build_parser
        parser = build_parser()
        args = parser.parse_args([
            "-d", "/x", "-e", "E1", "-r", "/r.fa",
            "--force-rebasecall", "--dry-run", "-o", "/out",
        ])
        assert args.force_rebasecall is True
        assert args.dry_run is True
        assert args.outdir == Path("/out")
