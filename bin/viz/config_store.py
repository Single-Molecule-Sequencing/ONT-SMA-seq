"""ConfigStore: reads and writes ONT-native config files for SMA-seq experiments.

Manages the on-disk representation of an experiment directory, including:
- ``sma_experiment.toml`` (top-level experiment configuration)
- ``sample_sheet.csv`` (MinKNOW sample sheet)
- ``arrangement.toml`` (Dorado custom barcode arrangement)
- ``barcodes.fasta`` (barcode sequences)
- ``references/*.fasta`` (per-target reference sequences)

All file writes go through a backup mechanism that keeps timestamped copies
in a ``.backups/`` subdirectory.
"""

from __future__ import annotations

import csv
import shutil
import tomllib
from datetime import datetime, timezone
from pathlib import Path

import tomli_w

from viz.models import (
    ArrangementConfig,
    ExperimentConfig,
    SampleSheetEntry,
    TargetRef,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_EXPERIMENT_TOML = "sma_experiment.toml"
_SAMPLE_SHEET_CSV = "sample_sheet.csv"
_ARRANGEMENT_TOML = "arrangement.toml"
_BARCODES_FASTA = "barcodes.fasta"
_REFERENCES_DIR = "references"

_SAMPLE_SHEET_COLUMNS = [
    "flow_cell_id",
    "kit",
    "barcode",
    "alias",
    "experiment_id",
    "sample_id",
    "type",
]


# ---------------------------------------------------------------------------
# ConfigStore
# ---------------------------------------------------------------------------


class ConfigStore:
    """Manages reading and writing ONT-native config files for one experiment.

    Parameters
    ----------
    experiment_dir : Path
        Root directory of the experiment. Created if it does not exist.
    """

    MAX_BACKUPS = 10

    def __init__(self, experiment_dir: Path) -> None:
        self.dir = Path(experiment_dir).resolve()
        self.dir.mkdir(parents=True, exist_ok=True)
        self._backups_dir = self.dir / ".backups"
        self.experiment_config = self._load_or_create_experiment_config()

    # ------------------------------------------------------------------
    # sma_experiment.toml
    # ------------------------------------------------------------------

    def _load_or_create_experiment_config(self) -> ExperimentConfig:
        """Load ``sma_experiment.toml`` or create one with defaults."""
        toml_path = self.dir / _EXPERIMENT_TOML
        if toml_path.exists():
            return self._read_experiment_toml(toml_path)
        cfg = ExperimentConfig()
        self._write_experiment_toml(cfg)
        return cfg

    def _read_experiment_toml(self, path: Path) -> ExperimentConfig:
        """Parse an existing experiment TOML file into an ExperimentConfig."""
        raw = path.read_bytes()
        data = tomllib.loads(raw.decode())

        # The TOML is structured with [experiment] at the top level
        exp_data = data.get("experiment", {})

        return ExperimentConfig(**exp_data)

    def save_experiment_config(self) -> None:
        """Persist the current experiment config to disk (with backup)."""
        toml_path = self.dir / _EXPERIMENT_TOML
        if toml_path.exists():
            self._backup(_EXPERIMENT_TOML)
        self._write_experiment_toml(self.experiment_config)

    def _write_experiment_toml(self, cfg: ExperimentConfig) -> None:
        """Serialise an ExperimentConfig to ``sma_experiment.toml``."""
        # Build the dict structure: everything under [experiment]
        data = {"experiment": cfg.model_dump(mode="python")}
        toml_path = self.dir / _EXPERIMENT_TOML
        toml_path.write_bytes(tomli_w.dumps(data).encode())

    # ------------------------------------------------------------------
    # sample_sheet.csv  (MinKNOW CSV)
    # ------------------------------------------------------------------

    def read_sample_sheet(self) -> list[SampleSheetEntry]:
        """Read ``sample_sheet.csv`` and return a list of SampleSheetEntry."""
        csv_path = self.dir / _SAMPLE_SHEET_CSV
        if not csv_path.exists():
            return []

        entries: list[SampleSheetEntry] = []
        with csv_path.open(newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                # Only pass known fields to avoid pydantic complaints
                filtered = {
                    k: v.strip() for k, v in row.items() if k in _SAMPLE_SHEET_COLUMNS
                }
                entries.append(SampleSheetEntry(**filtered))
        return entries

    def write_sample_sheet(self, entries: list[SampleSheetEntry]) -> None:
        """Write sample sheet entries to ``sample_sheet.csv``."""
        csv_path = self.dir / _SAMPLE_SHEET_CSV
        if csv_path.exists():
            self._backup(_SAMPLE_SHEET_CSV)

        with csv_path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=_SAMPLE_SHEET_COLUMNS)
            writer.writeheader()
            for entry in entries:
                row = entry.model_dump(mode="python")
                writer.writerow({k: row.get(k, "") for k in _SAMPLE_SHEET_COLUMNS})

    # ------------------------------------------------------------------
    # arrangement.toml  (Dorado custom barcodes)
    # ------------------------------------------------------------------

    def read_arrangement(self) -> ArrangementConfig | None:
        """Read ``arrangement.toml`` or return ``None`` if it does not exist."""
        toml_path = self.dir / _ARRANGEMENT_TOML
        if not toml_path.exists():
            return None

        raw = toml_path.read_bytes()
        data = tomllib.loads(raw.decode())

        # scoring may be a nested table
        return ArrangementConfig(**data)

    def write_arrangement(self, arr: ArrangementConfig) -> None:
        """Write an ArrangementConfig to ``arrangement.toml``."""
        toml_path = self.dir / _ARRANGEMENT_TOML
        if toml_path.exists():
            self._backup(_ARRANGEMENT_TOML)

        data = arr.model_dump(mode="python")
        toml_path.write_bytes(tomli_w.dumps(data).encode())

    # ------------------------------------------------------------------
    # barcodes.fasta
    # ------------------------------------------------------------------

    def read_barcode_sequences(self) -> dict[str, str]:
        """Read ``barcodes.fasta`` and return id-to-sequence mapping."""
        fasta_path = self.dir / _BARCODES_FASTA
        if not fasta_path.exists():
            return {}
        return self._read_fasta(fasta_path)

    def write_barcode_sequences(self, sequences: dict[str, str]) -> None:
        """Write barcode sequences to ``barcodes.fasta``."""
        fasta_path = self.dir / _BARCODES_FASTA
        if fasta_path.exists():
            self._backup(_BARCODES_FASTA)
        self._write_fasta(fasta_path, sequences)

    # ------------------------------------------------------------------
    # references/*.fasta
    # ------------------------------------------------------------------

    def read_targets(self) -> list[TargetRef]:
        """Read all ``references/*.fasta`` files and return TargetRef list."""
        refs_dir = self.dir / _REFERENCES_DIR
        if not refs_dir.exists():
            return []

        targets: list[TargetRef] = []
        for fasta_path in sorted(refs_dir.glob("*.fasta")):
            seqs = self._read_fasta(fasta_path)
            for tgt_id, sequence in seqs.items():
                targets.append(
                    TargetRef(tgt_id=tgt_id, sequence=sequence, length=len(sequence))
                )
        return targets

    def write_target(self, tgt_id: str, sequence: str) -> None:
        """Write a single target reference to ``references/<tgt_id>.fasta``."""
        refs_dir = self.dir / _REFERENCES_DIR
        refs_dir.mkdir(exist_ok=True)
        fasta_path = refs_dir / f"{tgt_id}.fasta"
        if fasta_path.exists():
            self._backup(f"{_REFERENCES_DIR}/{tgt_id}.fasta")
        self._write_fasta(fasta_path, {tgt_id: sequence})

    def delete_target(self, tgt_id: str) -> None:
        """Delete a target reference FASTA file."""
        fasta_path = self.dir / _REFERENCES_DIR / f"{tgt_id}.fasta"
        if fasta_path.exists():
            self._backup(f"{_REFERENCES_DIR}/{tgt_id}.fasta")
            fasta_path.unlink()

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate(self) -> list[str]:
        """Cross-validate sample sheet aliases against reference FASTA files.

        Returns a list of error strings. An empty list means everything is
        consistent.

        Checks performed:
        - Every sample sheet alias has a matching ``references/<alias>.fasta``
        - Every reference FASTA is referenced by at least one sample sheet alias
        """
        errors: list[str] = []

        # Gather aliases from sample sheet
        entries = self.read_sample_sheet()
        aliases = {e.alias for e in entries}

        # Gather target IDs from references
        targets = self.read_targets()
        target_ids = {t.tgt_id for t in targets}

        # Check aliases without references
        for alias in sorted(aliases - target_ids):
            errors.append(
                f"Sample sheet alias '{alias}' has no matching reference FASTA"
            )

        # Check orphan references
        for tgt_id in sorted(target_ids - aliases):
            errors.append(
                f"Reference '{tgt_id}' is not referenced by any sample sheet entry"
            )

        return errors

    # ------------------------------------------------------------------
    # File versioning
    # ------------------------------------------------------------------

    def _backup(self, relative_path: str) -> None:
        """Copy a file to ``.backups/`` with a timestamp suffix.

        Prunes old backups to keep at most ``MAX_BACKUPS`` per file.
        """
        src = self.dir / relative_path
        if not src.exists():
            return

        self._backups_dir.mkdir(exist_ok=True)

        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%f")
        # Flatten the relative path for the backup filename
        flat_name = relative_path.replace("/", "__")
        dst = self._backups_dir / f"{flat_name}.{timestamp}"
        shutil.copy2(src, dst)

        # Prune old backups for this file
        pattern = f"{flat_name}.*"
        existing = sorted(self._backups_dir.glob(pattern))
        while len(existing) > self.MAX_BACKUPS:
            oldest = existing.pop(0)
            oldest.unlink()

    # ------------------------------------------------------------------
    # FASTA helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _read_fasta(path: Path) -> dict[str, str]:
        """Parse a FASTA file into an id-to-sequence mapping.

        Handles multi-line sequences (joins all non-header lines).
        """
        sequences: dict[str, str] = {}
        current_id: str | None = None
        current_seq: list[str] = []

        for line in path.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]  # take first word after >
                current_seq = []
            else:
                current_seq.append(line)

        if current_id is not None:
            sequences[current_id] = "".join(current_seq)

        return sequences

    @staticmethod
    def _write_fasta(path: Path, sequences: dict[str, str]) -> None:
        """Write sequences to a FASTA file, wrapping lines at 80 characters."""
        lines: list[str] = []
        for seq_id, seq in sequences.items():
            lines.append(f">{seq_id}")
            # Wrap sequence at 80 characters
            for i in range(0, len(seq), 80):
                lines.append(seq[i : i + 80])
        path.write_text("\n".join(lines) + "\n")
