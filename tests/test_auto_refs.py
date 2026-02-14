"""Tests for the /api/construct/auto-refs endpoint.

Validates that the auto-reference generation creates the correct truncated
reference FASTA variants for each sample sheet entry, using adapter and
barcode sequences from the experiment configuration.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from viz.app import app, set_experiment_dir
from viz.config_store import ConfigStore

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def experiment_with_adapter(tmp_path: Path) -> Path:
    """Build an experiment directory with adapter sequences configured.

    Creates:
    - sample_sheet.csv with one entry: barcode05--barcode10 -> target_A
    - references/target_A.fasta with sequence ACGTACGTACGTACGT
    - sma_experiment.toml with adapter_5prime="AATTCCGG", adapter_3prime="GGCCTTAA"
    """
    exp = tmp_path / "experiment"
    exp.mkdir()
    (exp / "references").mkdir()

    (exp / "sample_sheet.csv").write_text(
        "flow_cell_id,kit,barcode,alias,type\n"
        "FAL12345,SQK-NBD114-96,barcode05--barcode10,target_A,test_sample\n"
    )

    (exp / "references" / "target_A.fasta").write_text(
        ">target_A\nACGTACGTACGTACGT\n"
    )

    # Create config with adapter sequences via ConfigStore
    store = ConfigStore(exp)
    store.experiment_config.construct.adapter_5prime = "AATTCCGG"
    store.experiment_config.construct.adapter_3prime = "GGCCTTAA"
    store.save_experiment_config()

    return exp


@pytest.fixture()
def client(experiment_with_adapter: Path) -> TestClient:
    """Set experiment dir and return a TestClient for the app."""
    set_experiment_dir(experiment_with_adapter)
    return TestClient(app)


# ---------------------------------------------------------------------------
# TestAutoRefs
# ---------------------------------------------------------------------------


class TestAutoRefs:
    """Validate POST /api/construct/auto-refs."""

    def test_generates_four_variants(self, client: TestClient, experiment_with_adapter: Path):
        """POST should return 200 with exactly 4 generated variant names."""
        resp = client.post("/api/construct/auto-refs")
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] == "ok"
        assert len(data["generated"]) == 4
        assert "target_A_full" in data["generated"]
        assert "target_A_no_adapter" in data["generated"]
        assert "target_A_no_end_bc" in data["generated"]
        assert "target_A_bc_only" in data["generated"]

    def test_auto_ref_files_exist(self, client: TestClient, experiment_with_adapter: Path):
        """POST should create references/auto/ with 4 FASTA files."""
        client.post("/api/construct/auto-refs")
        auto_dir = experiment_with_adapter / "references" / "auto"
        assert auto_dir.is_dir()
        fasta_files = sorted(auto_dir.glob("*.fasta"))
        assert len(fasta_files) == 4
        names = {f.stem for f in fasta_files}
        assert names == {
            "target_A_full",
            "target_A_no_adapter",
            "target_A_no_end_bc",
            "target_A_bc_only",
        }

    def test_full_variant_contains_all_components(
        self, client: TestClient, experiment_with_adapter: Path
    ):
        """The _full variant should contain adapter + bc1 + target + rc_bc2 + rc_adapter."""
        client.post("/api/construct/auto-refs")
        full_fasta = experiment_with_adapter / "references" / "auto" / "target_A_full.fasta"
        content = full_fasta.read_text()
        assert ">target_A_full" in content

        # Expected components:
        # 5' adapter:  AATTCCGG
        # bc1 (nb05):  AAGGTTACACAAACCCTGGACAAG
        # target:      ACGTACGTACGTACGT
        # rc_bc2 (rc of nb10): AAGCGTTGAAACCTTTGTCCTCTC
        # rc 3' adapter: TTAAGGCC
        assert "AATTCCGG" in content
        assert "AAGGTTACACAAACCCTGGACAAG" in content
        assert "ACGTACGTACGTACGT" in content
        assert "AAGCGTTGAAACCTTTGTCCTCTC" in content
        assert "TTAAGGCC" in content

        # Verify full sequence is the concatenation
        sequence_line = [ln for ln in content.splitlines() if not ln.startswith(">")][0]
        expected = (
            "AATTCCGG"
            "AAGGTTACACAAACCCTGGACAAG"
            "ACGTACGTACGTACGT"
            "AAGCGTTGAAACCTTTGTCCTCTC"
            "TTAAGGCC"
        )
        assert sequence_line == expected

    def test_no_adapter_variant(self, client: TestClient, experiment_with_adapter: Path):
        """The _no_adapter variant should contain bc1 + target + rc_bc2 (no adapters)."""
        client.post("/api/construct/auto-refs")
        fasta = experiment_with_adapter / "references" / "auto" / "target_A_no_adapter.fasta"
        content = fasta.read_text()
        assert ">target_A_no_adapter" in content

        sequence_line = [ln for ln in content.splitlines() if not ln.startswith(">")][0]
        expected = (
            "AAGGTTACACAAACCCTGGACAAG"  # bc1 (nb05)
            "ACGTACGTACGTACGT"  # target
            "AAGCGTTGAAACCTTTGTCCTCTC"  # rc_bc2 (rc of nb10)
        )
        assert sequence_line == expected
        # Ensure adapters are absent
        assert "AATTCCGG" not in sequence_line
        assert "TTAAGGCC" not in sequence_line

    def test_no_end_bc_variant(self, client: TestClient, experiment_with_adapter: Path):
        """The _no_end_bc variant should contain adapter + bc1 + target (no end barcode)."""
        client.post("/api/construct/auto-refs")
        fasta = experiment_with_adapter / "references" / "auto" / "target_A_no_end_bc.fasta"
        content = fasta.read_text()
        assert ">target_A_no_end_bc" in content

        sequence_line = [ln for ln in content.splitlines() if not ln.startswith(">")][0]
        expected = (
            "AATTCCGG"  # 5' adapter
            "AAGGTTACACAAACCCTGGACAAG"  # bc1 (nb05)
            "ACGTACGTACGTACGT"  # target
        )
        assert sequence_line == expected

    def test_bc_only_variant(self, client: TestClient, experiment_with_adapter: Path):
        """The _bc_only variant should contain only bc1 + rc_bc2."""
        client.post("/api/construct/auto-refs")
        fasta = experiment_with_adapter / "references" / "auto" / "target_A_bc_only.fasta"
        content = fasta.read_text()
        assert ">target_A_bc_only" in content

        sequence_line = [ln for ln in content.splitlines() if not ln.startswith(">")][0]
        expected = (
            "AAGGTTACACAAACCCTGGACAAG"  # bc1 (nb05)
            "AAGCGTTGAAACCTTTGTCCTCTC"  # rc_bc2 (rc of nb10)
        )
        assert sequence_line == expected
        # No target, no adapters
        assert "ACGTACGTACGTACGT" not in sequence_line
        assert "AATTCCGG" not in sequence_line

    def test_idempotent_regeneration(self, client: TestClient, experiment_with_adapter: Path):
        """Calling auto-refs twice should overwrite cleanly without duplicates."""
        resp1 = client.post("/api/construct/auto-refs")
        assert resp1.status_code == 200
        resp2 = client.post("/api/construct/auto-refs")
        assert resp2.status_code == 200
        assert resp2.json()["generated"] == resp1.json()["generated"]

        auto_dir = experiment_with_adapter / "references" / "auto"
        fasta_files = list(auto_dir.glob("*.fasta"))
        assert len(fasta_files) == 4

    def test_no_match_skips_entry(self, client: TestClient, experiment_with_adapter: Path):
        """Entries without matching target references are silently skipped."""
        # Add a sample sheet entry with no matching reference
        csv_path = experiment_with_adapter / "sample_sheet.csv"
        csv_path.write_text(
            "flow_cell_id,kit,barcode,alias,type\n"
            "FAL12345,SQK-NBD114-96,barcode05--barcode10,target_A,test_sample\n"
            "FAL12345,SQK-NBD114-96,barcode01--barcode02,missing_target,test_sample\n"
        )
        # Re-init so the store picks up the new sample sheet
        set_experiment_dir(experiment_with_adapter)

        resp = client.post("/api/construct/auto-refs")
        assert resp.status_code == 200
        data = resp.json()
        # Only target_A variants should be generated (missing_target has no ref)
        assert len(data["generated"]) == 4
        assert all("target_A" in name for name in data["generated"])
