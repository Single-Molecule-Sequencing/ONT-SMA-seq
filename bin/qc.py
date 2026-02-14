"""QC visualization module for SMA-seq alignment classification."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

# Plotting defaults
STYLE = "seaborn-v0_8-whitegrid"
FIGSIZE = (8, 5)
DPI = 300
SAVGOL_WINDOW = 51
SAVGOL_POLY = 3

FLOAT_COLS_CLASS = ["best_ned", "second_ned", "margin"]
INT_COLS_CLASS = ["read_len"]
FLOAT_COLS_ALIGN = [
    "ned", "identity", "ref_coverage", "read_to_ref_ratio",
    "seg5_identity", "seg5_contiguity", "segM_identity", "segM_contiguity",
    "seg3_identity", "seg3_contiguity", "five_prime_identity_20", "three_prime_identity_20",
    "margin",
]
INT_COLS_ALIGN = [
    "read_len", "ref_len", "ed", "max_ins", "max_del", "n_sig_indels",
    "five_prime_offset", "three_prime_offset", "rank",
]


def load_classification(tsv: Path) -> pd.DataFrame:
    """Load classification.tsv into a DataFrame with correct dtypes."""
    df = pd.read_csv(tsv, sep="\t")
    for col in FLOAT_COLS_CLASS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in INT_COLS_CLASS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    return df


def load_alignments(tsv: Path) -> pd.DataFrame:
    """Load alignments.tsv into a DataFrame with correct dtypes."""
    df = pd.read_csv(tsv, sep="\t")
    for col in FLOAT_COLS_ALIGN:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in INT_COLS_ALIGN:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    return df


def _kde_smooth(values: np.ndarray, bw: float = 0.02, n_points: int = 500) -> tuple[np.ndarray, np.ndarray]:
    """Compute Gaussian KDE with Savitzky-Golay smoothing."""
    from scipy.stats import gaussian_kde
    vmin, vmax = values.min(), values.max()
    pad = (vmax - vmin) * 0.1
    x = np.linspace(vmin - pad, vmax + pad, n_points)
    try:
        kde = gaussian_kde(values, bw_method=bw)
        y = kde(x)
    except (np.linalg.LinAlgError, ValueError):
        y = np.zeros_like(x)
    window = min(SAVGOL_WINDOW, len(y) - 1)
    if window % 2 == 0:
        window -= 1
    if window >= SAVGOL_POLY + 2:
        y = savgol_filter(y, window, SAVGOL_POLY)
        y = np.maximum(y, 0)
    return x, y


def _save(fig: plt.Figure, outdir: Path, name: str) -> None:
    """Save figure as PNG and PDF."""
    fig.savefig(outdir / f"{name}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(outdir / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# --- Plot functions ---

def plot_margin_kde(df: pd.DataFrame, outdir: Path) -> None:
    """Plot 2: Score margin distribution KDE."""
    plt.style.use(STYLE)
    fig, ax = plt.subplots(figsize=FIGSIZE)

    margins = df["margin"].dropna().values
    x, y = _kde_smooth(margins)
    ax.fill_between(x, y, alpha=0.3)
    ax.plot(x, y, linewidth=1.5)

    for thresh in [0.05, 0.1, 0.2]:
        ax.axvline(thresh, color="red", linestyle="--", alpha=0.5, label=f"threshold={thresh}")

    ax.set_xlabel("Classification Margin (NED_2nd - NED_best)")
    ax.set_ylabel("Density")
    ax.set_title("Score Margin Distribution")
    ax.legend(fontsize=8)
    _save(fig, outdir, "margin_kde")


def plot_ned_by_target(df_align: pd.DataFrame, df_class: pd.DataFrame, outdir: Path) -> None:
    """Plot 3: NED distribution by assigned target (best + second-best)."""
    plt.style.use(STYLE)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    assign = df_class[["read_id", "assigned_ref"]].drop_duplicates()
    merged = df_align.merge(assign, on="read_id", how="left")

    targets = sorted(merged["assigned_ref"].dropna().unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(targets), 1)))

    for i, tgt in enumerate(targets):
        best = merged[(merged["assigned_ref"] == tgt) & (merged["rank"] == 1)]
        if len(best) > 2:
            x, y = _kde_smooth(best["ned"].values)
            ax1.plot(x, y, color=colors[i], label=tgt, linewidth=1.5)
            ax1.fill_between(x, y, alpha=0.2, color=colors[i])

        second = merged[(merged["assigned_ref"] == tgt) & (merged["rank"] == 2)]
        if len(second) > 2:
            x, y = _kde_smooth(second["ned"].values)
            ax2.plot(x, y, color=colors[i], label=tgt, linewidth=1.5)
            ax2.fill_between(x, y, alpha=0.2, color=colors[i])

    ax1.set_ylabel("Density")
    ax1.set_title("NED to Assigned Target (Rank 1)")
    ax1.legend(fontsize=8)
    ax2.set_xlabel("Normalized Edit Distance")
    ax2.set_ylabel("Density")
    ax2.set_title("NED to Second-Best Target (Rank 2)")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    _save(fig, outdir, "ned_by_target")


def plot_length_vs_ned(df: pd.DataFrame, outdir: Path, ref_lengths: dict[str, int] | None = None) -> None:
    """Plot 5: Read length vs NED scatter."""
    plt.style.use(STYLE)
    fig, ax = plt.subplots(figsize=FIGSIZE)

    targets = sorted(df["assigned_ref"].dropna().unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(targets), 1)))

    for i, tgt in enumerate(targets):
        sub = df[df["assigned_ref"] == tgt]
        ax.scatter(sub["read_len"], sub["best_ned"], s=5, alpha=0.4, color=colors[i], label=tgt)

    for thresh in [0.1, 0.3, 0.5]:
        ax.axhline(thresh, color="gray", linestyle="--", alpha=0.5, linewidth=0.8)

    if ref_lengths:
        for name, length in ref_lengths.items():
            ax.axvline(length, color="red", linestyle=":", alpha=0.5, linewidth=0.8)

    ax.set_xlabel("Read Length (bp)")
    ax.set_ylabel("NED to Assigned Target")
    ax.set_title("Read Length vs. Normalized Edit Distance")
    ax.legend(fontsize=8, markerscale=3)
    _save(fig, outdir, "length_vs_ned_scatter")


def plot_segmented_violins(df_align: pd.DataFrame, df_class: pd.DataFrame, outdir: Path) -> None:
    """Plot 4: Segmented identity violin plots per target."""
    plt.style.use(STYLE)

    assign = df_class[["read_id", "assigned_ref"]].drop_duplicates()
    best = df_align[df_align["rank"] == 1].merge(assign, on="read_id", how="left")

    targets = sorted(best["assigned_ref"].dropna().unique())
    n_targets = max(len(targets), 1)
    fig, axes = plt.subplots(1, n_targets, figsize=(4 * n_targets, 5), sharey=True, squeeze=False)

    segments = ["seg5_identity", "segM_identity", "seg3_identity"]
    labels = ["5'", "Mid", "3'"]

    for idx, tgt in enumerate(targets):
        ax = axes[0, idx]
        sub = best[best["assigned_ref"] == tgt]
        data = [sub[seg].dropna().values for seg in segments]
        if all(len(d) > 0 for d in data):
            vp = ax.violinplot(data, showmedians=True)
            ax.set_xticks([1, 2, 3])
            ax.set_xticklabels(labels)
        ax.set_title(tgt, fontsize=10)
        if idx == 0:
            ax.set_ylabel("Identity")

    fig.suptitle("Segmented Alignment Identity (5' / Mid / 3')")
    fig.tight_layout()
    _save(fig, outdir, "segmented_identity_violins")


def plot_five_prime_quality(df_align: pd.DataFrame, df_class: pd.DataFrame, outdir: Path) -> None:
    """Plot 6: 5' alignment quality (offset + identity KDEs)."""
    plt.style.use(STYLE)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

    assign = df_class[["read_id", "assigned_ref"]].drop_duplicates()
    best = df_align[df_align["rank"] == 1].merge(assign, on="read_id", how="left")
    targets = sorted(best["assigned_ref"].dropna().unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(targets), 1)))

    for i, tgt in enumerate(targets):
        sub = best[best["assigned_ref"] == tgt]

        offsets = sub["five_prime_offset"].dropna().values.astype(float)
        if len(offsets) > 2:
            x, y = _kde_smooth(offsets, bw=0.5)
            ax1.plot(x, y, color=colors[i], label=tgt)

        ident = sub["five_prime_identity_20"].dropna().values
        if len(ident) > 2:
            x, y = _kde_smooth(ident)
            ax2.plot(x, y, color=colors[i], label=tgt)

    ax1.set_xlabel("5' Offset (ref bases uncovered)")
    ax1.set_ylabel("Density")
    ax1.set_title("5' Alignment Offset")
    ax1.legend(fontsize=8)

    ax2.set_xlabel("5' Identity (first 20bp)")
    ax2.set_ylabel("Density")
    ax2.set_title("5' Alignment Identity")
    ax2.legend(fontsize=8)

    fig.tight_layout()
    _save(fig, outdir, "five_prime_quality")


def plot_indel_profile(df_align: pd.DataFrame, df_class: pd.DataFrame, outdir: Path) -> None:
    """Plot 7: Significant indel profile (size KDE + position histogram)."""
    plt.style.use(STYLE)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

    assign = df_class[["read_id", "assigned_ref"]].drop_duplicates()
    best = df_align[df_align["rank"] == 1].merge(assign, on="read_id", how="left")
    targets = sorted(best["assigned_ref"].dropna().unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(targets), 1)))

    for i, tgt in enumerate(targets):
        sub = best[best["assigned_ref"] == tgt]
        ins = sub["max_ins"].dropna().values.astype(float)
        dels = sub["max_del"].dropna().values.astype(float)
        if len(ins) > 2:
            x, y = _kde_smooth(ins, bw=0.5)
            ax1.plot(x, y, color=colors[i], linestyle="-", label=f"{tgt} ins")
        if len(dels) > 2:
            x, y = _kde_smooth(dels, bw=0.5)
            ax1.plot(x, y, color=colors[i], linestyle="--", label=f"{tgt} del")

    ax1.set_xlabel("Indel Size (bp)")
    ax1.set_ylabel("Density")
    ax1.set_title("Largest Indel Size Distribution")
    ax1.legend(fontsize=7)

    # Position histogram: bin indel counts by 5'/mid/3' segment
    seg_counts = {"5'": 0, "Mid": 0, "3'": 0}
    for _, row in best.iterrows():
        ref_len = row.get("ref_len", 0)
        if ref_len is None or ref_len <= 0:
            continue
        offset = row.get("five_prime_offset", 0) or 0
        if offset > ref_len / 3:
            seg_counts["5'"] += 1
        three_offset = row.get("three_prime_offset", 0) or 0
        if three_offset > ref_len / 3:
            seg_counts["3'"] += 1

    n_sig = best["n_sig_indels"].dropna()
    mid_count = max(0, len(n_sig[n_sig > 0]) - seg_counts["5'"] - seg_counts["3'"])
    ax2.bar(["5'", "Mid", "3'"], [seg_counts["5'"], mid_count, seg_counts["3'"]])
    ax2.set_xlabel("Reference Segment")
    ax2.set_ylabel("Reads with Significant Indels")
    ax2.set_title("Indel Hotspot by Segment")

    fig.tight_layout()
    _save(fig, outdir, "indel_profile")


def plot_threshold_sweep(df: pd.DataFrame, outdir: Path) -> None:
    """Plot 8: Misclassification threshold sweep."""
    plt.style.use(STYLE)
    fig, ax1 = plt.subplots(figsize=FIGSIZE)
    ax2 = ax1.twinx()

    thresholds = np.arange(0.1, 0.52, 0.02)
    frac_classified = []
    misclass_rate = []
    n_total = len(df)

    for t in thresholds:
        classified = df[df["best_ned"] <= t]
        frac = len(classified) / n_total if n_total > 0 else 0
        frac_classified.append(frac)

        if len(classified) > 0:
            ambiguous = classified[classified["margin"] <= 0.05]
            misclass_rate.append(len(ambiguous) / len(classified))
        else:
            misclass_rate.append(0)

    ax1.plot(thresholds, frac_classified, "b-", linewidth=2, label="Fraction Classified")
    ax2.plot(thresholds, misclass_rate, "r--", linewidth=2, label="Est. Misclassification Rate")

    ax1.set_xlabel("NED Threshold")
    ax1.set_ylabel("Fraction of Reads Classified", color="blue")
    ax2.set_ylabel("Est. Misclassification Rate", color="red")
    ax1.set_title("Classification Threshold Sweep")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=8)

    _save(fig, outdir, "threshold_sweep")


def plot_classification_heatmap(df_align: pd.DataFrame, df_class: pd.DataFrame, outdir: Path, max_reads: int = 5000) -> None:
    """Plot 1: Classification landscape heatmap (reads x refs, color=NED)."""
    plt.style.use(STYLE)

    assign = df_class[["read_id", "assigned_ref"]].drop_duplicates()

    # Subsample if too many reads
    read_ids = df_align["read_id"].unique()
    if len(read_ids) > max_reads:
        rng = np.random.default_rng(42)
        read_ids = rng.choice(read_ids, max_reads, replace=False)

    sub = df_align[df_align["read_id"].isin(read_ids)].merge(assign, on="read_id", how="left")
    refs = sorted(sub["ref_id"].unique())

    pivot = sub.pivot_table(index="read_id", columns="ref_id", values="ned", aggfunc="first")
    pivot = pivot.reindex(columns=refs)

    # Sort by assigned target, then by NED to that target
    pivot = pivot.merge(assign.set_index("read_id"), left_index=True, right_index=True)
    pivot = pivot.sort_values(["assigned_ref", refs[0] if refs else ""])
    pivot = pivot.drop(columns=["assigned_ref"])

    fig, ax = plt.subplots(figsize=(max(6, len(refs) * 2), 8))
    im = ax.imshow(pivot.values, aspect="auto", cmap="viridis_r", vmin=0, vmax=0.8)
    ax.set_xticks(range(len(refs)))
    ax.set_xticklabels(refs, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(f"Reads (n={len(pivot)})")
    ax.set_yticks([])
    ax.set_title("Classification Landscape (color = NED)")
    fig.colorbar(im, ax=ax, label="Normalized Edit Distance")
    fig.tight_layout()
    _save(fig, outdir, "classification_heatmap")


def plot_per_read_profiles(df_align: pd.DataFrame, outdir: Path, n_sample: int = 50) -> None:
    """Plot 9: Per-read bar charts showing NED to all targets (sampled)."""
    plt.style.use(STYLE)
    read_ids = df_align["read_id"].unique()
    rng = np.random.default_rng(42)
    sampled = rng.choice(read_ids, min(n_sample, len(read_ids)), replace=False)

    n = len(sampled)
    ncols = min(5, n)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 2.5 * nrows), squeeze=False)

    for idx, rid in enumerate(sampled):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        sub = df_align[df_align["read_id"] == rid].sort_values("ned")
        colors = ["green" if r == 1 else "gray" for r in sub["rank"]]
        ax.barh(sub["ref_id"], sub["ned"], color=colors)
        ax.set_xlim(0, 1)
        ax.set_title(rid[:8], fontsize=7)
        ax.tick_params(labelsize=6)

    # Hide empty subplots
    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    fig.suptitle(f"Per-Read NED Profiles (n={n} sampled)", fontsize=11)
    fig.tight_layout()
    _save(fig, outdir, "per_read_profiles_sample")


def plot_end_reason_facets(df_align: pd.DataFrame, df_class: pd.DataFrame, outdir: Path) -> None:
    """Plot 10: Margin KDE faceted by end reason."""
    plt.style.use(STYLE)

    ers = df_class["end_reason"].dropna().unique()
    if len(ers) < 2:
        fig, ax = plt.subplots(figsize=FIGSIZE)
        ax.text(0.5, 0.5, "Insufficient end reason diversity for faceting",
                ha="center", va="center", transform=ax.transAxes)
        _save(fig, outdir, "margin_by_endreason")
        return

    fig, axes = plt.subplots(1, len(ers), figsize=(5 * len(ers), 5), sharey=True, squeeze=False)

    for idx, er in enumerate(sorted(ers)):
        ax = axes[0, idx]
        sub = df_class[df_class["end_reason"] == er]
        margins = sub["margin"].dropna().values
        if len(margins) > 2:
            x, y = _kde_smooth(margins)
            ax.fill_between(x, y, alpha=0.3)
            ax.plot(x, y, linewidth=1.5)
        ax.set_title(er.replace("_", "\n"), fontsize=8)
        ax.set_xlabel("Margin")
        if idx == 0:
            ax.set_ylabel("Density")

    fig.suptitle("Classification Margin by End Reason")
    fig.tight_layout()
    _save(fig, outdir, "margin_by_endreason")


# --- Orchestrator ---

def run_qc(
    alignments_tsv: Path,
    classification_tsv: Path,
    outdir: Path,
    ref_lengths: dict[str, int] | None = None,
) -> list[Path]:
    """Generate all QC plots from alignment and classification TSVs.

    Returns list of generated plot file paths.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    df_align = load_alignments(alignments_tsv)
    df_class = load_classification(classification_tsv)

    plots = []

    plot_classification_heatmap(df_align, df_class, outdir)
    plots.append(outdir / "classification_heatmap.png")

    plot_margin_kde(df_class, outdir)
    plots.append(outdir / "margin_kde.png")

    plot_ned_by_target(df_align, df_class, outdir)
    plots.append(outdir / "ned_by_target.png")

    plot_segmented_violins(df_align, df_class, outdir)
    plots.append(outdir / "segmented_identity_violins.png")

    plot_length_vs_ned(df_class, outdir, ref_lengths=ref_lengths)
    plots.append(outdir / "length_vs_ned_scatter.png")

    plot_five_prime_quality(df_align, df_class, outdir)
    plots.append(outdir / "five_prime_quality.png")

    plot_indel_profile(df_align, df_class, outdir)
    plots.append(outdir / "indel_profile.png")

    plot_threshold_sweep(df_class, outdir)
    plots.append(outdir / "threshold_sweep.png")

    plot_per_read_profiles(df_align, outdir)
    plots.append(outdir / "per_read_profiles_sample.png")

    plot_end_reason_facets(df_align, df_class, outdir)
    plots.append(outdir / "margin_by_endreason.png")

    return plots
