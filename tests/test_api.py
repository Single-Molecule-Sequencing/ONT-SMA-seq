"""Tests for the REST API endpoints (viz.api).

Test classes:
- TestSampleSheetAPI: GET/PUT sample sheet entries
- TestTargetsAPI: GET/POST/DELETE target references
- TestValidationAPI: GET validate endpoint
- TestExperimentAPI: GET/PUT experiment config
"""

from __future__ import annotations

from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from viz.app import app, set_experiment_dir


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def experiment_dir(tmp_path: Path) -> Path:
    """Create an experiment directory with a sample sheet and one reference."""
    exp = tmp_path / "test_experiment"
    exp.mkdir()

    # Write a minimal sample sheet CSV
    csv_path = exp / "sample_sheet.csv"
    csv_path.write_text(
        "flow_cell_id,kit,barcode,alias,experiment_id,sample_id,type\n"
        "FC001,SQK-NBD114.96,barcode01--barcode02,CYP2D6_fwd,EXP001,S001,test_sample\n"
    )

    # Write a reference FASTA
    refs_dir = exp / "references"
    refs_dir.mkdir()
    fasta = refs_dir / "CYP2D6_fwd.fasta"
    fasta.write_text(">CYP2D6_fwd\nACGTACGTACGTACGTACGT\n")

    return exp


@pytest.fixture()
def client(experiment_dir: Path) -> TestClient:
    """Set experiment dir and return a TestClient for the app."""
    set_experiment_dir(experiment_dir)
    return TestClient(app)


# ---------------------------------------------------------------------------
# TestSampleSheetAPI
# ---------------------------------------------------------------------------


class TestSampleSheetAPI:
    """Validate GET and PUT for /api/sample-sheet."""

    def test_get(self, client: TestClient):
        """GET /api/sample-sheet should return 1 entry with correct alias."""
        resp = client.get("/api/sample-sheet")
        assert resp.status_code == 200
        data = resp.json()
        assert isinstance(data, list)
        assert len(data) == 1
        assert data[0]["alias"] == "CYP2D6_fwd"

    def test_put(self, client: TestClient):
        """PUT /api/sample-sheet should replace entries; GET returns new count."""
        new_entries = [
            {
                "flow_cell_id": "FC001",
                "kit": "SQK-NBD114.96",
                "barcode": "barcode01--barcode02",
                "alias": "CYP2D6_fwd",
                "experiment_id": "EXP001",
                "sample_id": "S001",
                "type": "test_sample",
            },
            {
                "flow_cell_id": "FC001",
                "kit": "SQK-NBD114.96",
                "barcode": "barcode03--barcode04",
                "alias": "CYP2D6_rev",
                "experiment_id": "EXP001",
                "sample_id": "S002",
                "type": "test_sample",
            },
        ]
        resp = client.put("/api/sample-sheet", json=new_entries)
        assert resp.status_code == 200
        assert resp.json()["status"] == "ok"

        # Verify the new count
        resp2 = client.get("/api/sample-sheet")
        assert len(resp2.json()) == 2


# ---------------------------------------------------------------------------
# TestTargetsAPI
# ---------------------------------------------------------------------------


class TestTargetsAPI:
    """Validate GET/POST/DELETE for /api/targets."""

    def test_get(self, client: TestClient):
        """GET /api/targets should return 1 target."""
        resp = client.get("/api/targets")
        assert resp.status_code == 200
        data = resp.json()
        assert isinstance(data, list)
        assert len(data) == 1
        assert data[0]["tgt_id"] == "CYP2D6_fwd"

    def test_post(self, client: TestClient):
        """POST /api/targets should add a new target."""
        resp = client.post(
            "/api/targets",
            json={"tgt_id": "GeneX", "sequence": "AAACCCGGGTTT"},
        )
        assert resp.status_code == 200
        assert resp.json()["status"] == "ok"

        # Should now have 2 targets
        resp2 = client.get("/api/targets")
        assert len(resp2.json()) == 2

    def test_delete(self, client: TestClient):
        """DELETE /api/targets/{alias} should remove the target."""
        resp = client.delete("/api/targets/CYP2D6_fwd")
        assert resp.status_code == 200
        assert resp.json()["status"] == "ok"

        # Should now have 0 targets
        resp2 = client.get("/api/targets")
        assert len(resp2.json()) == 0


# ---------------------------------------------------------------------------
# TestValidationAPI
# ---------------------------------------------------------------------------


class TestValidationAPI:
    """Validate GET /api/validate."""

    def test_validate(self, client: TestClient):
        """GET /api/validate should return errors list and valid bool."""
        resp = client.get("/api/validate")
        assert resp.status_code == 200
        data = resp.json()
        assert "errors" in data
        assert "valid" in data
        assert isinstance(data["errors"], list)
        assert isinstance(data["valid"], bool)


# ---------------------------------------------------------------------------
# TestExperimentAPI
# ---------------------------------------------------------------------------


class TestExperimentAPI:
    """Validate GET and PUT for /api/experiment."""

    def test_get(self, client: TestClient):
        """GET /api/experiment should return experiment config dict."""
        resp = client.get("/api/experiment")
        assert resp.status_code == 200
        data = resp.json()
        assert "description" in data
        assert "construct" in data
        assert "assumptions" in data

    def test_put(self, client: TestClient):
        """PUT /api/experiment should update description and persist it."""
        resp = client.put(
            "/api/experiment",
            json={"description": "Updated via API"},
        )
        assert resp.status_code == 200
        assert resp.json()["status"] == "ok"

        # Verify the update took effect
        resp2 = client.get("/api/experiment")
        assert resp2.json()["description"] == "Updated via API"
