"""Generate SMA-seq construct TOML and sample sheet CSV files."""
from __future__ import annotations

import csv
from pathlib import Path

try:
    import tomli_w
    _HAS_TOMLI_W = True
except ImportError:
    _HAS_TOMLI_W = False


def generate_construct_toml(
    targets: list[dict],
    kit: str,
    mode: str,
    output_path: Path,
    mask1_front: str,
    mask1_rear: str,
    mask2_front: str | None = None,
    mask2_rear: str | None = None,
) -> None:
    """Write a TOML config file describing the construct arrangement.

    Parameters
    ----------
    targets : list[dict]
        Each dict has keys: barcode1, barcode2 (optional), alias.
    kit : str
        Barcoding kit name (e.g. "SQK-NBD114-96").
    mode : str
        Either "dual_independent" or "start_only".
    output_path : Path
        Where to write the TOML file.
    mask1_front, mask1_rear : str
        Mask sequences for the first barcode position.
    mask2_front, mask2_rear : str | None
        Mask sequences for the second barcode position (dual mode only).
    """
    arrangement: dict = {
        "kit": kit,
        "mask1_front": mask1_front,
        "mask1_rear": mask1_rear,
    }
    if mode == "dual_independent" and mask2_front is not None:
        arrangement["mask2_front"] = mask2_front
    if mode == "dual_independent" and mask2_rear is not None:
        arrangement["mask2_rear"] = mask2_rear

    data = {
        "arrangement": arrangement,
        "sma": {
            "mode": mode,
            "targets": targets,
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)

    if _HAS_TOMLI_W:
        with output_path.open("wb") as f:
            tomli_w.dump(data, f)
    else:
        _write_toml_manual(data, output_path)


def _write_toml_manual(data: dict, path: Path) -> None:
    """Fallback TOML writer when tomli_w is not available."""
    lines: list[str] = []

    # [arrangement] section
    lines.append("[arrangement]")
    for key, value in data["arrangement"].items():
        lines.append(f'{key} = "{value}"')
    lines.append("")

    # [sma] section
    lines.append("[sma]")
    lines.append(f'mode = "{data["sma"]["mode"]}"')

    # targets as array of inline tables
    target_strs: list[str] = []
    for t in data["sma"]["targets"]:
        parts = [f'{k} = "{v}"' for k, v in t.items()]
        target_strs.append("{" + ", ".join(parts) + "}")
    lines.append("targets = [")
    for i, ts in enumerate(target_strs):
        comma = "," if i < len(target_strs) - 1 else ""
        lines.append(f"    {ts}{comma}")
    lines.append("]")
    lines.append("")

    path.write_text("\n".join(lines))


def generate_sample_sheet(
    entries: list[dict],
    flow_cell_id: str,
    kit: str,
    experiment_id: str,
    output_path: Path,
) -> None:
    """Write a CSV sample sheet for MinKNOW / downstream tools.

    Parameters
    ----------
    entries : list[dict]
        Each dict has keys: barcode, alias, type.
    flow_cell_id : str
        Flow cell identifier.
    kit : str
        Barcoding kit name.
    experiment_id : str
        Experiment identifier.
    output_path : Path
        Where to write the CSV file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ["flow_cell_id", "kit", "barcode", "alias", "type", "experiment_id"]
    with output_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for entry in entries:
            writer.writerow({
                "flow_cell_id": flow_cell_id,
                "kit": kit,
                "barcode": entry["barcode"],
                "alias": entry["alias"],
                "type": entry["type"],
                "experiment_id": experiment_id,
            })
