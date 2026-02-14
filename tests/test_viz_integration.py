"""Integration test: full workflow from config creation to export."""

from __future__ import annotations

from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from viz.app import app, set_experiment_dir


@pytest.fixture()
def empty_experiment(tmp_path: Path) -> Path:
    exp = tmp_path / "new_experiment"
    exp.mkdir()
    return exp


@pytest.fixture()
def client(empty_experiment: Path) -> TestClient:
    set_experiment_dir(empty_experiment)
    return TestClient(app)


class TestFullWorkflow:
    """Test creating an experiment config from scratch and exporting."""

    def test_create_and_export(self, client: TestClient, empty_experiment: Path):
        # 1. Set experiment description
        resp = client.put("/api/experiment", json={
            "description": "Integration test experiment",
        })
        assert resp.status_code == 200

        # 2. Add sample sheet entries
        resp = client.put("/api/sample-sheet", json=[
            {
                "flow_cell_id": "FAL99999",
                "kit": "SQK-NBD114-96",
                "barcode": "barcode05--barcode10",
                "alias": "target_A",
                "type": "test_sample",
            },
        ])
        assert resp.status_code == 200

        # 3. Add target reference
        resp = client.post("/api/targets", json={
            "tgt_id": "target_A",
            "sequence": "ACGTACGTACGTACGTACGTACGT",
        })
        assert resp.status_code == 200

        # 4. Validate - should be clean
        resp = client.get("/api/validate")
        assert resp.status_code == 200
        assert resp.json()["valid"] is True

        # 5. Export
        resp = client.post("/api/export")
        assert resp.status_code == 200
        export_path = Path(resp.json()["path"])
        assert export_path.is_dir()
        assert (export_path / "index.html").exists()
        assert (export_path / "sample_sheet.html").exists()
        assert (export_path / "construct.html").exists()

        # 6. Verify exported HTML contains our data
        index = (export_path / "index.html").read_text()
        assert "Integration test experiment" in index
        assert "target_A" in index

    def test_empty_experiment_validates(self, client: TestClient):
        """An empty experiment should validate without errors."""
        resp = client.get("/api/validate")
        assert resp.status_code == 200
        assert resp.json()["valid"] is True

    def test_validation_catches_missing_reference(self, client: TestClient):
        """Adding a sample sheet entry without a matching reference should fail validation."""
        client.put("/api/sample-sheet", json=[
            {
                "flow_cell_id": "FAL99999",
                "kit": "SQK-NBD114-96",
                "barcode": "barcode01--barcode02",
                "alias": "missing_target",
                "type": "test_sample",
            },
        ])
        resp = client.get("/api/validate")
        assert resp.status_code == 200
        assert resp.json()["valid"] is False
        assert any("missing_target" in e for e in resp.json()["errors"])
