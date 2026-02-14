"""REST API router for SMA-seq config visualizer.

Provides CRUD endpoints for all configuration elements: experiment config,
sample sheet, barcodes, construct, targets, assumptions, validation, and export.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import edlib
from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse, HTMLResponse

from viz.config_store import ConfigStore
from viz.models import (
    ArrangementConfig,
    Assumption,
    SampleSheetEntry,
)

# Try importing barcodes module (lives in bin/ alongside viz/)
import barcodes as _barcodes_mod

router = APIRouter(prefix="/api")

# ---------------------------------------------------------------------------
# Module-level store state
# ---------------------------------------------------------------------------

_store: ConfigStore | None = None


def set_store(store: ConfigStore) -> None:
    """Set the module-level ConfigStore instance."""
    global _store
    _store = store


def get_store() -> ConfigStore:
    """Return the module-level ConfigStore, raising 500 if not set."""
    if _store is None:
        raise HTTPException(
            status_code=500,
            detail="ConfigStore not initialised. Set an experiment directory first.",
        )
    return _store


# ---------------------------------------------------------------------------
# 1. Experiment
# ---------------------------------------------------------------------------


@router.get("/experiment")
async def get_experiment() -> dict[str, Any]:
    """Return the current experiment configuration."""
    store = get_store()
    return store.experiment_config.model_dump()


@router.put("/experiment")
async def put_experiment(payload: dict[str, Any]) -> dict[str, str]:
    """Update fields on the experiment configuration and save."""
    store = get_store()
    cfg = store.experiment_config

    # Update top-level simple fields
    for key, value in payload.items():
        if hasattr(cfg, key):
            setattr(cfg, key, value)

    store.save_experiment_config()
    return {"status": "ok"}


# ---------------------------------------------------------------------------
# 2. Sample Sheet
# ---------------------------------------------------------------------------


@router.get("/sample-sheet")
async def get_sample_sheet() -> list[dict[str, Any]]:
    """Return sample sheet entries as a list of dicts."""
    store = get_store()
    entries = store.read_sample_sheet()
    return [e.model_dump() for e in entries]


@router.put("/sample-sheet")
async def put_sample_sheet(payload: list[dict[str, Any]]) -> dict[str, str]:
    """Replace sample sheet entries."""
    store = get_store()
    entries = [SampleSheetEntry(**row) for row in payload]
    store.write_sample_sheet(entries)
    return {"status": "ok"}


@router.get("/sample-sheet/export")
async def export_sample_sheet() -> FileResponse:
    """Download sample_sheet.csv as a file."""
    store = get_store()
    csv_path = store.dir / "sample_sheet.csv"
    if not csv_path.exists():
        raise HTTPException(status_code=404, detail="sample_sheet.csv not found")
    return FileResponse(
        path=str(csv_path),
        filename="sample_sheet.csv",
        media_type="text/csv",
    )


# ---------------------------------------------------------------------------
# 3. Barcodes
# ---------------------------------------------------------------------------


@router.get("/barcodes")
async def get_barcodes() -> dict[str, Any]:
    """Return arrangement, custom sequences, and ONT built-in barcodes."""
    store = get_store()
    arrangement = store.read_arrangement()
    custom_seqs = store.read_barcode_sequences()
    return {
        "arrangement": arrangement.model_dump() if arrangement else None,
        "custom_sequences": custom_seqs,
        "ont_barcodes": _barcodes_mod.BARCODES,
    }


@router.put("/barcodes/arrangement")
async def put_barcode_arrangement(payload: dict[str, Any]) -> dict[str, str]:
    """Write a new ArrangementConfig."""
    store = get_store()
    arr = ArrangementConfig(**payload)
    store.write_arrangement(arr)
    return {"status": "ok"}


@router.put("/barcodes/sequences")
async def put_barcode_sequences(payload: dict[str, str]) -> dict[str, str]:
    """Write custom barcode sequences to barcodes.fasta."""
    store = get_store()
    store.write_barcode_sequences(payload)
    return {"status": "ok"}


@router.get("/barcodes/distances")
async def get_barcode_distances() -> dict[str, Any]:
    """Compute pairwise edit distances for barcodes used in the sample sheet."""
    store = get_store()
    entries = store.read_sample_sheet()

    # Collect unique barcode IDs referenced by sample sheet
    used_ids: set[str] = set()
    for entry in entries:
        if entry.upstream_barcode:
            used_ids.add(entry.upstream_barcode)
        if entry.downstream_barcode:
            used_ids.add(entry.downstream_barcode)

    # Resolve sequences (prefer custom, fall back to ONT built-in)
    custom_seqs = store.read_barcode_sequences()
    all_barcodes = {**_barcodes_mod.BARCODES, **custom_seqs}

    bc_seqs: dict[str, str] = {}
    for bc_id in sorted(used_ids):
        if bc_id in all_barcodes:
            bc_seqs[bc_id] = all_barcodes[bc_id]

    # Compute pairwise edit distances
    ids = sorted(bc_seqs.keys())
    distances: list[dict[str, Any]] = []
    for i, id_a in enumerate(ids):
        for j in range(i + 1, len(ids)):
            id_b = ids[j]
            result = edlib.align(bc_seqs[id_a], bc_seqs[id_b], mode="NW", task="distance")
            distances.append({
                "barcode_a": id_a,
                "barcode_b": id_b,
                "edit_distance": result["editDistance"],
            })

    return {"barcodes": ids, "distances": distances}


# ---------------------------------------------------------------------------
# 4. Construct
# ---------------------------------------------------------------------------


@router.get("/construct")
async def get_construct() -> dict[str, Any]:
    """Return construct, demultiplexing, truncation, and classification configs."""
    store = get_store()
    cfg = store.experiment_config
    return {
        "construct": cfg.construct.model_dump(),
        "demultiplexing": cfg.demultiplexing.model_dump(),
        "truncation": cfg.truncation.model_dump(),
        "classification": cfg.classification.model_dump(),
    }


@router.put("/construct")
async def put_construct(payload: dict[str, Any]) -> dict[str, str]:
    """Update nested construct-related configs."""
    store = get_store()
    cfg = store.experiment_config

    if "construct" in payload:
        for k, v in payload["construct"].items():
            if hasattr(cfg.construct, k):
                setattr(cfg.construct, k, v)

    if "demultiplexing" in payload:
        for k, v in payload["demultiplexing"].items():
            if hasattr(cfg.demultiplexing, k):
                setattr(cfg.demultiplexing, k, v)

    if "truncation" in payload:
        for k, v in payload["truncation"].items():
            if hasattr(cfg.truncation, k):
                setattr(cfg.truncation, k, v)

    if "classification" in payload:
        for k, v in payload["classification"].items():
            if hasattr(cfg.classification, k):
                setattr(cfg.classification, k, v)

    store.save_experiment_config()
    return {"status": "ok"}


@router.post("/construct/auto-refs")
async def auto_refs() -> dict[str, Any]:
    """Generate truncated reference FASTAs for each sample sheet entry.

    For each entry with a matching target, generates 4 variants:
    - {alias}_full: adapter + bc1 + target + rc_bc2 + rc_adapter
    - {alias}_no_adapter: bc1 + target + rc_bc2
    - {alias}_no_end_bc: adapter + bc1 + target
    - {alias}_bc_only: bc1 + rc_bc2
    """
    store = get_store()
    cfg = store.experiment_config
    entries = store.read_sample_sheet()
    targets = {t.tgt_id: t.sequence for t in store.read_targets()}

    adapter_5 = cfg.construct.adapter_5prime
    adapter_3 = cfg.construct.adapter_3prime
    rc_adapter_3 = _barcodes_mod.reverse_complement(adapter_3) if adapter_3 else ""

    # Merge custom and ONT barcodes
    custom_seqs = store.read_barcode_sequences()
    all_barcodes = {**_barcodes_mod.BARCODES, **custom_seqs}

    auto_dir = store.dir / "references" / "auto"
    auto_dir.mkdir(parents=True, exist_ok=True)

    generated: list[str] = []

    for entry in entries:
        alias = entry.alias
        if alias not in targets:
            continue

        target_seq = targets[alias]
        bc1_id = entry.upstream_barcode
        bc2_id = entry.downstream_barcode

        bc1_seq = all_barcodes.get(bc1_id, "")
        bc2_seq = all_barcodes.get(bc2_id, "")
        rc_bc2 = _barcodes_mod.reverse_complement(bc2_seq) if bc2_seq else ""

        variants = {
            f"{alias}_full": adapter_5 + bc1_seq + target_seq + rc_bc2 + rc_adapter_3,
            f"{alias}_no_adapter": bc1_seq + target_seq + rc_bc2,
            f"{alias}_no_end_bc": adapter_5 + bc1_seq + target_seq,
            f"{alias}_bc_only": bc1_seq + rc_bc2,
        }

        for variant_name, sequence in variants.items():
            fasta_path = auto_dir / f"{variant_name}.fasta"
            fasta_path.write_text(f">{variant_name}\n{sequence}\n")
            generated.append(variant_name)

    return {"status": "ok", "generated": generated}


# ---------------------------------------------------------------------------
# 5. Targets
# ---------------------------------------------------------------------------


@router.get("/targets")
async def get_targets() -> list[dict[str, Any]]:
    """Return all target references as a list of dicts."""
    store = get_store()
    targets = store.read_targets()
    return [t.model_dump() for t in targets]


@router.post("/targets")
async def post_target(payload: dict[str, str]) -> dict[str, str]:
    """Add a new target reference (tgt_id + sequence)."""
    store = get_store()
    tgt_id = payload.get("tgt_id", "")
    sequence = payload.get("sequence", "")
    if not tgt_id or not sequence:
        raise HTTPException(
            status_code=400,
            detail="Both 'tgt_id' and 'sequence' are required.",
        )
    store.write_target(tgt_id, sequence)
    return {"status": "ok"}


@router.put("/targets/{alias}")
async def put_target(alias: str, payload: dict[str, str]) -> dict[str, str]:
    """Update an existing target reference."""
    store = get_store()
    sequence = payload.get("sequence", "")
    if not sequence:
        raise HTTPException(status_code=400, detail="'sequence' is required.")
    store.write_target(alias, sequence)
    return {"status": "ok"}


@router.delete("/targets/{alias}")
async def delete_target(alias: str) -> dict[str, str]:
    """Delete a target reference."""
    store = get_store()
    store.delete_target(alias)
    return {"status": "ok"}


# ---------------------------------------------------------------------------
# 6. Assumptions
# ---------------------------------------------------------------------------


@router.get("/assumptions")
async def get_assumptions() -> list[dict[str, str]]:
    """Return the list of documented assumptions."""
    store = get_store()
    return [a.model_dump() for a in store.experiment_config.assumptions]


@router.put("/assumptions")
async def put_assumptions(payload: list[dict[str, str]]) -> dict[str, str]:
    """Replace all assumptions."""
    store = get_store()
    store.experiment_config.assumptions = [Assumption(**a) for a in payload]
    store.save_experiment_config()
    return {"status": "ok"}


# ---------------------------------------------------------------------------
# 7. Validation
# ---------------------------------------------------------------------------


@router.get("/validate")
async def validate() -> dict[str, Any]:
    """Cross-validate sample sheet against references."""
    store = get_store()
    errors = store.validate()
    return {"errors": errors, "valid": len(errors) == 0}


@router.get("/validate/badge")
async def validate_badge() -> HTMLResponse:
    """Return an HTML badge indicating validation status."""
    store = get_store()
    errors = store.validate()
    if errors:
        badge = '<span class="badge-error">&#x2717; Errors found</span>'
    else:
        badge = '<span class="badge-ok">&#x2713; Valid</span>'
    return HTMLResponse(content=badge)


# ---------------------------------------------------------------------------
# 8. Export
# ---------------------------------------------------------------------------


@router.post("/export")
async def export_config() -> dict[str, str]:
    """Placeholder for full config export."""
    return {"status": "not_implemented"}
