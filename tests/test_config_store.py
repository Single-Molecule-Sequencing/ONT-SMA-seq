"""Tests for the ConfigStore class.

Tests cover:
- Initialisation: default TOML creation, loading existing TOML
- Sample sheet: read/write MinKNOW CSV
- Targets: read/write/delete reference FASTA files
- Arrangement: read default (None), write and read back
- Experiment config: save and reload round-trip
- Validation: consistent state, missing reference detection
- Versioning: backup creation on write operations
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from viz.config_store import ConfigStore
from viz.models import (
    ArrangementConfig,
    ExperimentConfig,
    SampleSheetEntry,
    TargetRef,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def experiment_dir(tmp_path: Path) -> Path:
    """Create a bare experiment directory."""
    exp = tmp_path / "my_experiment"
    exp.mkdir()
    return exp


@pytest.fixture()
def store(experiment_dir: Path) -> ConfigStore:
    """Return a ConfigStore pointed at the experiment directory."""
    return ConfigStore(experiment_dir)


@pytest.fixture()
def sample_entries() -> list[SampleSheetEntry]:
    """Two minimal sample sheet entries for testing."""
    return [
        SampleSheetEntry(
            flow_cell_id="FC001",
            kit="SQK-NBD114.96",
            barcode="barcode01--barcode02",
            alias="CYP2D6_fwd",
            type="test_sample",
        ),
        SampleSheetEntry(
            flow_cell_id="FC001",
            kit="SQK-NBD114.96",
            barcode="barcode02--barcode01",
            alias="CYP2D6_rev",
            type="test_sample",
        ),
    ]


# ---------------------------------------------------------------------------
# TestConfigStoreInit
# ---------------------------------------------------------------------------


class TestConfigStoreInit:
    """Verify ConfigStore initialisation behaviour."""

    def test_creates_default_toml_if_missing(self, experiment_dir: Path):
        """When no sma_experiment.toml exists, init should create one with defaults."""
        store = ConfigStore(experiment_dir)
        toml_path = experiment_dir / "sma_experiment.toml"
        assert toml_path.exists()
        # The loaded config should have default assumptions
        assert len(store.experiment_config.assumptions) == 4

    def test_loads_existing_toml(self, experiment_dir: Path):
        """When a sma_experiment.toml exists, init should load it."""
        # Create a store to generate the default TOML
        store1 = ConfigStore(experiment_dir)
        store1.experiment_config.description = "Test experiment"
        store1.save_experiment_config()

        # Create a new store and verify it loaded the description
        store2 = ConfigStore(experiment_dir)
        assert store2.experiment_config.description == "Test experiment"


# ---------------------------------------------------------------------------
# TestSampleSheet
# ---------------------------------------------------------------------------


class TestSampleSheet:
    """Verify sample sheet read/write round-trips."""

    def test_read_sample_sheet(self, experiment_dir: Path):
        """Write a CSV by hand, then read it back via ConfigStore."""
        csv_path = experiment_dir / "sample_sheet.csv"
        csv_path.write_text(
            "flow_cell_id,kit,barcode,alias,type\n"
            "FC001,SQK-NBD114.96,barcode05--barcode10,CYP2D6_fwd,test_sample\n"
        )
        store = ConfigStore(experiment_dir)
        entries = store.read_sample_sheet()
        assert len(entries) == 1
        assert entries[0].alias == "CYP2D6_fwd"
        assert entries[0].upstream_barcode == "nb05"
        assert entries[0].downstream_barcode == "nb10"

    def test_write_sample_sheet(
        self, store: ConfigStore, sample_entries: list[SampleSheetEntry]
    ):
        """Write entries via ConfigStore, then read the raw CSV back."""
        store.write_sample_sheet(sample_entries)
        csv_path = store.dir / "sample_sheet.csv"
        assert csv_path.exists()

        # Read back and verify
        entries = store.read_sample_sheet()
        assert len(entries) == 2
        aliases = {e.alias for e in entries}
        assert aliases == {"CYP2D6_fwd", "CYP2D6_rev"}


# ---------------------------------------------------------------------------
# TestTargets
# ---------------------------------------------------------------------------


class TestTargets:
    """Verify target reference FASTA read/write/delete."""

    def test_read_targets(self, store: ConfigStore):
        """Write a FASTA by hand, then read targets via ConfigStore."""
        refs_dir = store.dir / "references"
        refs_dir.mkdir(exist_ok=True)
        fasta = refs_dir / "CYP2D6_fwd.fasta"
        fasta.write_text(">CYP2D6_fwd\nACGTACGTACGT\n")
        targets = store.read_targets()
        assert len(targets) == 1
        assert targets[0].tgt_id == "CYP2D6_fwd"
        assert targets[0].sequence == "ACGTACGTACGT"
        assert targets[0].length == 12

    def test_write_target(self, store: ConfigStore):
        """Write a target via ConfigStore and verify the file is created."""
        store.write_target("GeneX", "AAACCCGGGTTT")
        refs_dir = store.dir / "references"
        fasta_path = refs_dir / "GeneX.fasta"
        assert fasta_path.exists()

        targets = store.read_targets()
        assert len(targets) == 1
        assert targets[0].tgt_id == "GeneX"
        assert targets[0].sequence == "AAACCCGGGTTT"

    def test_delete_target(self, store: ConfigStore):
        """Write then delete a target, verify file is removed."""
        store.write_target("GeneX", "AAACCCGGGTTT")
        assert len(store.read_targets()) == 1

        store.delete_target("GeneX")
        assert len(store.read_targets()) == 0
        assert not (store.dir / "references" / "GeneX.fasta").exists()


# ---------------------------------------------------------------------------
# TestArrangement
# ---------------------------------------------------------------------------


class TestArrangement:
    """Verify Dorado arrangement TOML read/write."""

    def test_read_arrangement_default(self, store: ConfigStore):
        """Without an arrangement.toml, read_arrangement returns None."""
        assert store.read_arrangement() is None

    def test_write_and_read_arrangement(self, store: ConfigStore):
        """Write an arrangement, then read it back and verify fields."""
        arr = ArrangementConfig(
            name="custom_kit",
            kit="SQK-NBD114.96",
            first_index=1,
            last_index=24,
        )
        store.write_arrangement(arr)
        assert (store.dir / "arrangement.toml").exists()

        loaded = store.read_arrangement()
        assert loaded is not None
        assert loaded.name == "custom_kit"
        assert loaded.kit == "SQK-NBD114.96"
        assert loaded.last_index == 24
        # Scoring defaults should be preserved
        assert loaded.scoring.max_barcode_penalty == 11


# ---------------------------------------------------------------------------
# TestExperimentConfig
# ---------------------------------------------------------------------------


class TestExperimentConfig:
    """Verify experiment config save/reload round-trip."""

    def test_save_and_reload(self, store: ConfigStore):
        """Modify config, save, create new store, verify changes persisted."""
        store.experiment_config.description = "My cool experiment"
        store.experiment_config.construct.insert_type = "genomic"
        store.save_experiment_config()

        store2 = ConfigStore(store.dir)
        assert store2.experiment_config.description == "My cool experiment"
        assert store2.experiment_config.construct.insert_type == "genomic"


# ---------------------------------------------------------------------------
# TestValidation
# ---------------------------------------------------------------------------


class TestValidation:
    """Verify cross-file consistency checks."""

    def test_validate_consistent(
        self, store: ConfigStore, sample_entries: list[SampleSheetEntry]
    ):
        """When every alias has a matching reference, validate returns no errors."""
        store.write_sample_sheet(sample_entries)
        for entry in sample_entries:
            store.write_target(entry.alias, "ACGTACGT" * 10)

        errors = store.validate()
        assert errors == []

    def test_validate_missing_reference(
        self, store: ConfigStore, sample_entries: list[SampleSheetEntry]
    ):
        """When a sample sheet alias has no reference, validate reports it."""
        store.write_sample_sheet(sample_entries)
        # Only write one reference, leaving the other missing
        store.write_target(sample_entries[0].alias, "ACGTACGT" * 10)

        errors = store.validate()
        assert len(errors) >= 1
        # The error should mention the missing alias
        missing_alias = sample_entries[1].alias
        assert any(missing_alias in err for err in errors)


# ---------------------------------------------------------------------------
# TestVersioning
# ---------------------------------------------------------------------------


class TestVersioning:
    """Verify file backup on writes."""

    def test_backup_created_on_write(self, store: ConfigStore):
        """Writing experiment config should create a backup of the old file."""
        # First save creates the file
        store.save_experiment_config()
        # Second save should back up the first
        store.experiment_config.description = "Updated"
        store.save_experiment_config()

        backups_dir = store.dir / ".backups"
        assert backups_dir.exists()
        backup_files = list(backups_dir.glob("sma_experiment.toml.*"))
        assert len(backup_files) >= 1
