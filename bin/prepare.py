"""SMA-seq experiment preparation: discover, plan, merge, init, align+qc."""
from __future__ import annotations

import argparse
import json
import sqlite3
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

from sma_merge.models import RunInfo, RunGroup
from sma_merge.validate import validate_runs


@dataclass
class Manifest:
    """Checkpoint / provenance manifest for the prepare pipeline."""

    exp_id: str
    expdir: str
    timestamp: str = ""
    runs: list[dict] = field(default_factory=list)
    merge_groups: list[list[str]] = field(default_factory=list)
    output_bam: str = ""
    output_bam_reads: int = 0
    database_path: str = ""
    ref_path: str = ""
    ref_targets: list[str] = field(default_factory=list)
    stages_completed: list[str] = field(default_factory=list)
    qc_plots: list[str] = field(default_factory=list)

    def __post_init__(self):
        if not self.timestamp:
            self.timestamp = time.strftime("%Y-%m-%dT%H:%M:%S")

    def mark_stage(self, stage: str) -> None:
        if stage not in self.stages_completed:
            self.stages_completed.append(stage)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.__dict__, f, indent=2, default=str)

    @classmethod
    def load(cls, path: Path) -> "Manifest | None":
        if not path.exists():
            return None
        with open(path) as f:
            data = json.load(f)
        m = cls(exp_id=data["exp_id"], expdir=data["expdir"])
        for k, v in data.items():
            if hasattr(m, k):
                setattr(m, k, v)
        return m


def init_database(
    db_path: Path, exp_id: str, flow_cell_id: str, sample_id: str, alias: str,
) -> None:
    """Create SMA-seq SQLite database with full schema."""
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    c = conn.cursor()

    c.execute('''CREATE TABLE Reads (
        uniq_id TEXT PRIMARY KEY, exp_id TEXT, tgt_id TEXT, read_id TEXT,
        readseq TEXT, readlen INTEGER, model_tier TEXT, model_ver TEXT,
        trim INTEGER, mod_bitflag INTEGER, ed INTEGER, q_bc REAL, q_ld REAL,
        ER TEXT, bc_start_id TEXT, bc_start_ed INTEGER, bc_start_conf REAL,
        bc_end_id TEXT, bc_end_ed INTEGER, bc_end_conf REAL,
        trunc_level TEXT, signal_duration_s REAL, mean_qscore REAL,
        FOREIGN KEY(tgt_id) REFERENCES Target(tgt_id),
        FOREIGN KEY(mod_bitflag) REFERENCES Mods(mod_bitflag)
    )''')

    c.execute('''CREATE TABLE Mods (
        mod_bitflag INTEGER PRIMARY KEY, mods TEXT
    )''')

    c.execute('''CREATE TABLE Exp (
        exp_id TEXT PRIMARY KEY, flow_cell_id TEXT,
        sample_id TEXT, alias TEXT, exp_desc TEXT
    )''')

    c.execute('''CREATE TABLE Target (
        tgt_id TEXT PRIMARY KEY, tgt_refseq TEXT, tgt_reflen INTEGER
    )''')

    c.execute('''CREATE TABLE IF NOT EXISTS RunMetadata (
        run_id TEXT PRIMARY KEY, flow_cell_id TEXT, device_id TEXT,
        sample_id TEXT, experiment_id TEXT, kit TEXT,
        protocol_run_id TEXT, start_time TEXT, basecall_model TEXT,
        source_bam_count INTEGER, source_bam_paths TEXT, merge_timestamp TEXT
    )''')

    c.execute('''CREATE TABLE IF NOT EXISTS ReadRun (
        read_id TEXT PRIMARY KEY, run_id TEXT,
        FOREIGN KEY (run_id) REFERENCES RunMetadata(run_id)
    )''')

    mods_data = [
        (0, "non"), (1, "6mA"), (2, "5mCG_5hmCG"), (4, "5mC_5hmC"),
        (8, "4mC_5mC"), (16, "5mC"), (3, "6mA,5mCG_5hmCG"),
        (5, "6mA,5mC_5hmC"), (9, "6mA,4mC_5mC"), (17, "6mA,5mC"),
    ]
    c.executemany("INSERT OR IGNORE INTO Mods (mod_bitflag, mods) VALUES (?, ?)", mods_data)

    c.execute(
        "INSERT INTO Exp (exp_id, flow_cell_id, sample_id, alias, exp_desc) VALUES (?, ?, ?, ?, ?)",
        (exp_id, flow_cell_id, sample_id, alias, "Initialized via prepare"),
    )

    conn.commit()
    conn.close()


def insert_run_metadata(db_path: Path, run_meta: dict) -> None:
    """Insert a single RunMetadata row."""
    conn = sqlite3.connect(str(db_path))
    c = conn.cursor()
    c.execute(
        """INSERT OR REPLACE INTO RunMetadata
        (run_id, flow_cell_id, device_id, sample_id, experiment_id, kit,
         protocol_run_id, start_time, basecall_model, source_bam_count,
         source_bam_paths, merge_timestamp)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (
            run_meta.get("run_id", ""),
            run_meta.get("flow_cell_id", ""),
            run_meta.get("device_id", ""),
            run_meta.get("sample_id", ""),
            run_meta.get("experiment_id", ""),
            run_meta.get("kit", ""),
            run_meta.get("protocol_run_id", ""),
            run_meta.get("start_time", ""),
            run_meta.get("basecall_model", ""),
            run_meta.get("source_bam_count", 0),
            run_meta.get("source_bam_paths", ""),
            time.strftime("%Y-%m-%dT%H:%M:%S"),
        ),
    )
    conn.commit()
    conn.close()


def symlink_pod5s(runs: list[RunInfo], output_dir: Path) -> None:
    """Create organized symlinks to POD5 files from each run."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for run in runs:
        run_subdir = output_dir / run.run_id
        run_subdir.mkdir(exist_ok=True)
        for pod5_file in sorted(run.pod5_dir.rglob("*.pod5")):
            link = run_subdir / pod5_file.name
            if link.exists() or link.is_symlink():
                link.unlink()
            link.symlink_to(pod5_file.resolve())


def collect_bam_files(runs: list[RunInfo]) -> list[Path]:
    """Collect all BAM files from bam_pass/ across runs."""
    bams: list[Path] = []
    for run in runs:
        bam_pass = run.run_dir / "bam_pass"
        if bam_pass.is_dir():
            bams.extend(sorted(bam_pass.rglob("*.bam")))
    return bams


def format_merge_plan(groups: list[RunGroup], exp_id: str) -> str:
    """Format the merge plan for user display."""
    lines = [
        f"[prepare] Merge Plan for {exp_id}",
        "-" * (28 + len(exp_id)),
    ]
    for g in groups:
        for i, run in enumerate(g.runs, 1):
            pod5_size_gb = 0
            for f in run.pod5_dir.rglob("*.pod5"):
                try:
                    pod5_size_gb += f.stat().st_size
                except OSError:
                    pass
            pod5_size_gb /= 1e9

            model_short = run.basecall_model.split("@")[-1] if "@" in run.basecall_model else run.basecall_model
            lines.append(f"Run {i}: {run.run_id}")
            lines.append(f"  Flow Cell: {g.flow_cell_id} | POD5: {run.pod5_count} files ({pod5_size_gb:.0f} GB) | Model: {model_short}")
            action = "REUSE existing BAMs" if g.is_consistent else "REBASECALL from POD5"
            lines.append(f"  Action: {action}")
            lines.append("")

        if len(g.runs) > 1:
            lines.append(f"Merge group: [{' + '.join(f'Run {i+1}' for i in range(len(g.runs)))}] -> same flow cell")
        for issue in g.issues:
            lines.append(f"  NOTE: {issue}")
        lines.append("")

    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="prepare",
        description="SMA-seq experiment preparation: discover, plan, merge, init, align+qc.",
    )
    parser.add_argument("-d", "--expdir", type=Path, required=True,
                        help="Path to experiment directory")
    parser.add_argument("-e", "--expid", required=True,
                        help="Experiment ID (FlowCell_Sample_Alias)")
    parser.add_argument("-r", "--ref", type=Path, required=True,
                        help="Multi-sequence FASTA with target sequences")
    parser.add_argument("-o", "--outdir", type=Path, default=Path("Output"),
                        help="Output directory (default: Output)")
    parser.add_argument("--force-rebasecall", action="store_true",
                        help="Skip smart detection, always re-basecall from POD5")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print merge plan without executing")
    return parser


def run_prepare(args: argparse.Namespace) -> None:
    """Execute the full prepare pipeline."""
    from sma_merge.discover import discover_runs
    from align import parse_fasta, process_bam
    from qc import run_qc

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / "prepare_manifest.json"

    # Check for existing manifest (resume)
    existing = Manifest.load(manifest_path)
    if existing and existing.exp_id == args.expid:
        print(f"[prepare] Found existing manifest. Stages completed: {existing.stages_completed}")
        manifest = existing
    else:
        manifest = Manifest(exp_id=args.expid, expdir=str(args.expdir))

    # --- Stage 1: DISCOVER ---
    if "discover" not in manifest.stages_completed:
        print(f"[prepare] Stage 1: Discovering runs in {args.expdir}...")
        runs = discover_runs(args.expdir)
        if not runs:
            sys.exit("[prepare] No MinKNOW runs found.")
        print(f"[prepare] Found {len(runs)} run(s).")
        manifest.runs = [
            {"run_id": r.run_id, "flow_cell_id": r.flow_cell_id,
             "device_id": r.device_id, "pod5_count": r.pod5_count,
             "basecall_model": r.basecall_model}
            for r in runs
        ]
        manifest.mark_stage("discover")
        manifest.save(manifest_path)
    else:
        print("[prepare] Stage 1: DISCOVER (cached)")
        runs = discover_runs(args.expdir)

    # --- Stage 2: PLAN ---
    if "plan" not in manifest.stages_completed:
        groups = validate_runs(runs)
        plan_text = format_merge_plan(groups, args.expid)
        print(plan_text)

        if args.dry_run:
            print("[prepare] Dry run — stopping after plan.")
            return

        response = input("Proceed? [Y/n] ").strip().lower()
        if response and response != "y":
            print("[prepare] Aborted by user.")
            return

        manifest.mark_stage("plan")
        manifest.save(manifest_path)
    else:
        print("[prepare] Stage 2: PLAN (cached)")
        groups = validate_runs(runs)

    # --- Stage 3: MERGE ---
    if "merge" not in manifest.stages_completed:
        print("[prepare] Stage 3: Merging...")
        for g in groups:
            # 3a: POD5 symlinks
            pod5_out = outdir / "pod5"
            symlink_pod5s(g.runs, pod5_out)

            # 3b: BAM merge
            if g.is_consistent and not args.force_rebasecall:
                import pysam
                bams = collect_bam_files(g.runs)
                if not bams:
                    print(f"[prepare] Warning: No BAMs found for {g.flow_cell_id}. Need re-basecall.")
                    continue
                merged_bam = outdir / "merged.bam"
                print(f"[prepare] Merging {len(bams)} BAM files...")
                pysam.merge("-f", "-n", str(merged_bam), *[str(b) for b in bams])
                print("[prepare] Sorting...")
                sorted_bam = outdir / "merged_sorted.bam"
                pysam.sort("-o", str(sorted_bam), str(merged_bam))
                merged_bam.unlink()
                sorted_bam.rename(merged_bam)
                print("[prepare] Indexing...")
                pysam.index(str(merged_bam))
                manifest.output_bam = str(merged_bam)
            else:
                from sma_merge.basecall import merge_pod5s, basecall
                print("[prepare] Re-basecalling from POD5...")
                merged_pod5 = outdir / f"{g.flow_cell_id}_merged.pod5"
                merge_pod5s([r.pod5_dir for r in g.runs], merged_pod5)
                merged_bam = outdir / "merged.bam"
                basecall(merged_pod5, merged_bam, g.basecall_model)
                manifest.output_bam = str(merged_bam)

        manifest.mark_stage("merge")
        manifest.save(manifest_path)
    else:
        print("[prepare] Stage 3: MERGE (cached)")

    # --- Stage 4: INIT ---
    if "init" not in manifest.stages_completed:
        print("[prepare] Stage 4: Initializing database...")
        exp_parts = args.expid.split("_")
        fc = exp_parts[0] if len(exp_parts) >= 1 else ""
        sid = exp_parts[1] if len(exp_parts) >= 2 else ""
        alias = exp_parts[-1] if len(exp_parts) >= 3 else ""

        db_path = outdir / f"SMA_{args.expid}.db"
        init_database(db_path, args.expid, fc, sid, alias)

        for g in groups:
            for run in g.runs:
                insert_run_metadata(db_path, {
                    "run_id": run.run_id,
                    "flow_cell_id": run.flow_cell_id,
                    "device_id": run.device_id,
                    "sample_id": run.sample_id,
                    "experiment_id": args.expid,
                    "basecall_model": run.basecall_model,
                    "source_bam_count": run.pod5_count,
                })

        manifest.database_path = str(db_path)
        manifest.mark_stage("init")
        manifest.save(manifest_path)
    else:
        print("[prepare] Stage 4: INIT (cached)")

    # --- Stage 5: ALIGN + QC ---
    if "align" not in manifest.stages_completed:
        print("[prepare] Stage 5a: Aligning all reads against all targets...")
        refs = parse_fasta(args.ref)
        manifest.ref_path = str(args.ref)
        manifest.ref_targets = list(refs.keys())

        merged_bam = Path(manifest.output_bam)
        align_tsv = outdir / "alignments.tsv"
        class_tsv = outdir / "classification.tsv"
        n = process_bam(merged_bam, refs, align_tsv, class_tsv)
        print(f"[prepare] Processed {n} reads.")
        manifest.output_bam_reads = n
        manifest.mark_stage("align")
        manifest.save(manifest_path)
    else:
        print("[prepare] Stage 5a: ALIGN (cached)")

    if "qc" not in manifest.stages_completed:
        print("[prepare] Stage 5b: Generating QC plots...")
        qc_dir = outdir / "qc"
        align_tsv = outdir / "alignments.tsv"
        class_tsv = outdir / "classification.tsv"
        ref_lengths = {k: len(v) for k, v in parse_fasta(args.ref).items()} if args.ref.exists() else None
        plots = run_qc(align_tsv, class_tsv, qc_dir, ref_lengths=ref_lengths)
        manifest.qc_plots = [str(p) for p in plots]
        manifest.mark_stage("qc")
        manifest.save(manifest_path)
    else:
        print("[prepare] Stage 5b: QC (cached)")

    print(f"\n[prepare] Complete. Outputs in {outdir}/")
    print(f"  - Merged BAM: {manifest.output_bam}")
    print(f"  - Database:   {manifest.database_path}")
    print(f"  - Alignments: {outdir}/alignments.tsv")
    print(f"  - Classified: {outdir}/classification.tsv")
    print(f"  - QC Plots:   {outdir}/qc/")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    run_prepare(args)


if __name__ == "__main__":
    main()
