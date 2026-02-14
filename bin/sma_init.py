#!/usr/bin/env python3
"""sma_init.py — Create SMA-seq database from merged BAM and generate QC report."""

import json as _json
import math
import sqlite3
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


# ---------------------------------------------------------------------------
# Database creation
# ---------------------------------------------------------------------------


def _mean_qscore(quals: list[int]) -> float:
    """Compute probability-averaged Q-score from quality scores."""
    if not quals:
        return 0.0
    probs = [10 ** (-q / 10) for q in quals]
    mean_prob = sum(probs) / len(probs)
    if mean_prob <= 0:
        return 40.0
    return -10 * math.log10(mean_prob)


def create_database(db_path: Path, exp_id: str) -> None:
    """Create the SMA-seq database with standard schema."""
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    c = conn.cursor()

    c.execute("""
        CREATE TABLE IF NOT EXISTS Reads (
            uniq_id TEXT PRIMARY KEY,
            exp_id TEXT,
            tgt_id TEXT,
            read_id TEXT,
            readseq TEXT,
            readlen INTEGER,
            model_tier TEXT,
            model_ver TEXT,
            trim INTEGER,
            mod_bitflag INTEGER,
            ed INTEGER,
            q_bc REAL,
            q_ld REAL,
            ER TEXT,
            bc_start_id TEXT,
            bc_start_ed INTEGER,
            bc_start_conf REAL,
            bc_end_id TEXT,
            bc_end_ed INTEGER,
            bc_end_conf REAL,
            bc_start_ed2 INTEGER,
            bc_start_conf2 REAL,
            bc_end_ed2 INTEGER,
            bc_end_conf2 REAL,
            trunc_level TEXT,
            signal_duration_s REAL,
            mean_qscore REAL
        )
    """)

    c.execute("""
        CREATE TABLE IF NOT EXISTS Exp (
            exp_id TEXT PRIMARY KEY,
            flow_cell_id TEXT,
            sample_id TEXT,
            alias TEXT,
            exp_desc TEXT
        )
    """)

    c.execute("""
        CREATE TABLE IF NOT EXISTS Target (
            tgt_id TEXT PRIMARY KEY,
            tgt_refseq TEXT,
            tgt_reflen INTEGER
        )
    """)

    c.execute("""
        CREATE TABLE IF NOT EXISTS RunMetadata (
            run_id TEXT PRIMARY KEY,
            flow_cell_id TEXT,
            device_id TEXT,
            sample_id TEXT,
            experiment_id TEXT,
            kit TEXT,
            protocol_run_id TEXT,
            start_time TEXT,
            basecall_model TEXT,
            source_bam_count INTEGER,
            source_bam_paths TEXT,
            merge_timestamp TEXT
        )
    """)

    c.execute("""
        CREATE TABLE IF NOT EXISTS ReadRun (
            read_id TEXT PRIMARY KEY,
            run_id TEXT,
            FOREIGN KEY (run_id) REFERENCES RunMetadata(run_id)
        )
    """)

    c.execute("""
        CREATE TABLE IF NOT EXISTS Mods (
            mod_bitflag INTEGER PRIMARY KEY,
            mods TEXT
        )
    """)

    conn.commit()
    conn.close()


def ingest_bam_to_db(
    bam_path: Path,
    db_path: Path,
    exp_id: str,
    summary_map: dict[str, dict] | None = None,
) -> int:
    """Ingest reads from BAM into database. Returns read count."""
    import pysam

    # Create DB if it doesn't exist
    if not db_path.exists():
        create_database(db_path, exp_id)

    conn = sqlite3.connect(str(db_path))
    c = conn.cursor()

    # Extract basecaller info from BAM header
    model_tier = None
    model_ver = None
    with pysam.AlignmentFile(str(bam_path), check_sq=False) as bam:
        header = bam.header.to_dict()
        for pg in header.get("PG", []):
            if pg.get("PN") in ("dorado", "basecaller"):
                model_ver = pg.get("VN", "")
                cl = pg.get("CL", "")
                if "sup" in cl:
                    model_tier = "sup"
                elif "hac" in cl:
                    model_tier = "hac"
                elif "fast" in cl:
                    model_tier = "fast"

    # Ingest reads
    count = 0
    if summary_map is None:
        summary_map = {}

    with pysam.AlignmentFile(str(bam_path), check_sq=False) as bam:
        for read in bam:
            read_id = read.query_name
            seq = read.query_sequence or ""
            quals = list(read.query_qualities) if read.query_qualities is not None else []
            q_bc = _mean_qscore(quals) if quals else None

            # Summary data
            sm = summary_map.get(read_id, {})
            end_reason = sm.get("end_reason")
            duration = sm.get("duration")
            mean_qscore = sm.get("mean_qscore")

            uniq_id = f"{exp_id}_{read_id}"

            c.execute(
                """INSERT OR IGNORE INTO Reads
                   (uniq_id, exp_id, read_id, readseq, readlen, q_bc,
                    model_tier, model_ver, ER, signal_duration_s, mean_qscore)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (uniq_id, exp_id, read_id, seq, len(seq), q_bc,
                 model_tier, model_ver, end_reason, duration, mean_qscore),
            )
            count += 1

    # Insert experiment record
    c.execute(
        "INSERT OR IGNORE INTO Exp (exp_id) VALUES (?)",
        (exp_id,),
    )

    conn.commit()
    conn.close()
    return count


# ---------------------------------------------------------------------------
# QC metrics
# ---------------------------------------------------------------------------


def compute_qc_metrics(db_path: Path) -> dict:
    """Compute pre-classification QC metrics from the database."""
    conn = sqlite3.connect(str(db_path))

    # Total read count
    total = conn.execute("SELECT COUNT(*) FROM Reads").fetchone()[0]

    # Read length distribution
    read_lengths = [
        r[0] for r in conn.execute("SELECT readlen FROM Reads WHERE readlen IS NOT NULL")
    ]

    # Q-score distribution
    qscores = [
        r[0] for r in conn.execute("SELECT q_bc FROM Reads WHERE q_bc IS NOT NULL")
    ]

    # End reason counts
    end_reason_rows = conn.execute(
        "SELECT ER, COUNT(*) FROM Reads GROUP BY ER ORDER BY COUNT(*) DESC"
    ).fetchall()
    end_reason_counts = {r[0] or "unknown": r[1] for r in end_reason_rows}

    # Signal duration
    durations = [
        r[0] for r in conn.execute(
            "SELECT signal_duration_s FROM Reads WHERE signal_duration_s IS NOT NULL"
        )
    ]

    # Duration vs length (for scatter)
    dur_vs_len = conn.execute(
        "SELECT readlen, signal_duration_s FROM Reads "
        "WHERE readlen IS NOT NULL AND signal_duration_s IS NOT NULL"
    ).fetchall()

    # Mean Q-score
    mean_qbc = None
    if qscores:
        probs = [10 ** (-q / 10) for q in qscores]
        mean_prob = sum(probs) / len(probs)
        mean_qbc = -10 * math.log10(mean_prob) if mean_prob > 0 else 0

    # Basic stats
    mean_len = sum(read_lengths) / len(read_lengths) if read_lengths else 0
    median_len = sorted(read_lengths)[len(read_lengths) // 2] if read_lengths else 0
    total_bases = sum(read_lengths)

    conn.close()

    return {
        "total_reads": total,
        "total_bases": total_bases,
        "mean_length": round(mean_len, 1),
        "median_length": median_len,
        "mean_qbc": round(mean_qbc, 2) if mean_qbc else None,
        "read_lengths": read_lengths,
        "qscores": qscores,
        "end_reason_counts": end_reason_counts,
        "durations": durations,
        "dur_vs_len": dur_vs_len,
    }


# ---------------------------------------------------------------------------
# HTML QC report
# ---------------------------------------------------------------------------


def generate_qc_html(metrics: dict, manifest: dict) -> str:
    """Generate interactive HTML QC report."""
    exp_id = manifest.get("experiment_id", "Unknown")
    kit = manifest.get("kit", "Unknown")
    fc_product = manifest.get("flow_cell_product", "Unknown")
    software = manifest.get("software", {})

    total = metrics["total_reads"]
    total_bases = metrics["total_bases"]
    mean_len = metrics["mean_length"]
    median_len = metrics["median_length"]
    mean_q = metrics["mean_qbc"]

    # Flow cell info
    flow_cells = manifest.get("flow_cells", [])
    n_runs = sum(len(fc.get("runs", [])) for fc in flow_cells)
    fc_ids = [fc.get("flow_cell_id", "?") for fc in flow_cells]

    # End reason table rows
    er_counts = metrics.get("end_reason_counts", {})
    er_rows = ""
    for er, cnt in sorted(er_counts.items(), key=lambda x: -x[1]):
        pct = cnt / total * 100 if total else 0
        er_rows += f"<tr><td>{er}</td><td>{cnt:,}</td><td>{pct:.1f}%</td></tr>\n"

    # Read length histogram data (bin into 50bp bins)
    lengths = metrics.get("read_lengths", [])
    if lengths:
        min_l = min(lengths)
        max_l = max(lengths)
        bin_size = max(50, (max_l - min_l) // 100)
        bins: dict[int, int] = {}
        for l in lengths:
            b = (l // bin_size) * bin_size
            bins[b] = bins.get(b, 0) + 1
        len_labels = _json.dumps(sorted(bins.keys()))
        len_values = _json.dumps([bins[k] for k in sorted(bins.keys())])
    else:
        len_labels = "[]"
        len_values = "[]"

    # Q-score histogram data (bin by 0.5)
    qscores = metrics.get("qscores", [])
    if qscores:
        q_bins: dict[float, int] = {}
        for q in qscores:
            b = round(q * 2) / 2  # 0.5 bins
            q_bins[b] = q_bins.get(b, 0) + 1
        q_labels = _json.dumps(sorted(q_bins.keys()))
        q_values = _json.dumps([q_bins[k] for k in sorted(q_bins.keys())])
    else:
        q_labels = "[]"
        q_values = "[]"

    # Duration histogram
    durations = metrics.get("durations", [])
    if durations:
        dur_max = min(max(durations), 100)  # cap at 100s
        dur_bin_size = max(1, dur_max // 50)
        dur_bins: dict[int, int] = {}
        for d in durations:
            b = int(min(d, dur_max) // dur_bin_size) * dur_bin_size
            dur_bins[b] = dur_bins.get(b, 0) + 1
        dur_labels = _json.dumps(sorted(dur_bins.keys()))
        dur_values = _json.dumps([dur_bins[k] for k in sorted(dur_bins.keys())])
    else:
        dur_labels = "[]"
        dur_values = "[]"

    # Duration vs length scatter (subsample for performance)
    dvl = metrics.get("dur_vs_len", [])
    import random
    if len(dvl) > 5000:
        rng = random.Random(42)
        dvl = rng.sample(dvl, 5000)
    scatter_data = _json.dumps([{"x": r[0], "y": r[1]} for r in dvl])

    # End reason pie data
    er_labels = _json.dumps(list(er_counts.keys()))
    er_values = _json.dumps(list(er_counts.values()))

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>SMA-seq QC Report — {exp_id}</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
       background: #f5f5f5; color: #333; }}
.header {{ background: linear-gradient(135deg, #1a5276, #2e86c1); color: white;
           padding: 24px 32px; }}
.header h1 {{ font-size: 24px; margin-bottom: 4px; }}
.header .subtitle {{ opacity: 0.85; font-size: 14px; }}
.tabs {{ display: flex; background: #fff; border-bottom: 2px solid #e0e0e0;
         padding: 0 16px; overflow-x: auto; }}
.tab-btn {{ padding: 12px 20px; border: none; background: none; cursor: pointer;
            font-size: 14px; color: #666; border-bottom: 3px solid transparent;
            white-space: nowrap; }}
.tab-btn:hover {{ color: #333; }}
.tab-btn.active {{ color: #1a5276; border-bottom-color: #2e86c1; font-weight: 600; }}
.tab-content {{ display: none; padding: 24px 32px; max-width: 1200px; margin: 0 auto; }}
.tab-content.active {{ display: block; }}
.card {{ background: #fff; border-radius: 8px; padding: 20px; margin-bottom: 16px;
         box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
.card h3 {{ margin-bottom: 12px; color: #1a5276; }}
.stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
               gap: 16px; margin-bottom: 24px; }}
.stat-card {{ background: #fff; border-radius: 8px; padding: 16px; text-align: center;
              box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
.stat-card .value {{ font-size: 28px; font-weight: 700; color: #1a5276; }}
.stat-card .label {{ font-size: 12px; color: #888; margin-top: 4px; }}
table {{ width: 100%; border-collapse: collapse; }}
th, td {{ padding: 8px 12px; text-align: left; border-bottom: 1px solid #eee; }}
th {{ background: #f8f9fa; font-weight: 600; color: #555; }}
.chart-container {{ position: relative; height: 350px; }}
</style>
</head>
<body>

<div class="header">
  <h1>SMA-seq QC Report</h1>
  <div class="subtitle">{exp_id} | {kit} | {fc_product} | {', '.join(fc_ids)}</div>
</div>

<div class="tabs">
  <button class="tab-btn active" data-tab="overview">Overview</button>
  <button class="tab-btn" data-tab="lengths">Read Lengths</button>
  <button class="tab-btn" data-tab="quality">Quality</button>
  <button class="tab-btn" data-tab="endreasons">End Reasons</button>
  <button class="tab-btn" data-tab="duration">Signal Duration</button>
</div>

<!-- Overview Tab -->
<div id="overview" class="tab-content active">
  <div class="stats-grid">
    <div class="stat-card"><div class="value">{total:,}</div><div class="label">Total Reads</div></div>
    <div class="stat-card"><div class="value">{total_bases / 1e6:.1f}M</div><div class="label">Total Bases</div></div>
    <div class="stat-card"><div class="value">{mean_len:.0f}</div><div class="label">Mean Length (bp)</div></div>
    <div class="stat-card"><div class="value">{median_len:,}</div><div class="label">Median Length (bp)</div></div>
    <div class="stat-card"><div class="value">{mean_q if mean_q else 'N/A'}</div><div class="label">Mean Q-score</div></div>
    <div class="stat-card"><div class="value">{n_runs}</div><div class="label">Sequencing Runs</div></div>
  </div>
  <div class="card">
    <h3>Experiment Details</h3>
    <table>
      <tr><th>Experiment ID</th><td>{exp_id}</td></tr>
      <tr><th>Kit</th><td>{kit}</td></tr>
      <tr><th>Flow Cell</th><td>{fc_product}</td></tr>
      <tr><th>Flow Cell IDs</th><td>{', '.join(fc_ids)}</td></tr>
      <tr><th>Runs</th><td>{n_runs}</td></tr>
      <tr><th>Basecaller</th><td>{software.get('basecaller_version', 'N/A')}</td></tr>
      <tr><th>MinKNOW</th><td>{software.get('distribution_version', 'N/A')}</td></tr>
    </table>
  </div>
</div>

<!-- Read Lengths Tab -->
<div id="lengths" class="tab-content">
  <div class="card">
    <h3>Read Length Distribution</h3>
    <div class="chart-container">
      <canvas id="lengthChart"></canvas>
    </div>
  </div>
</div>

<!-- Quality Tab -->
<div id="quality" class="tab-content">
  <div class="card">
    <h3>Q-Score Distribution</h3>
    <div class="chart-container">
      <canvas id="qscoreChart"></canvas>
    </div>
  </div>
</div>

<!-- End Reasons Tab -->
<div id="endreasons" class="tab-content">
  <div class="card">
    <h3>End Reason Breakdown</h3>
    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 24px;">
      <div>
        <table>
          <thead><tr><th>End Reason</th><th>Count</th><th>%</th></tr></thead>
          <tbody>{er_rows}</tbody>
        </table>
      </div>
      <div class="chart-container">
        <canvas id="erPieChart"></canvas>
      </div>
    </div>
  </div>
</div>

<!-- Duration Tab -->
<div id="duration" class="tab-content">
  <div class="card">
    <h3>Signal Duration Distribution</h3>
    <div class="chart-container">
      <canvas id="durationChart"></canvas>
    </div>
  </div>
  <div class="card">
    <h3>Duration vs Read Length</h3>
    <div class="chart-container">
      <canvas id="scatterChart"></canvas>
    </div>
  </div>
</div>

<script>
// Tab switching
document.querySelectorAll('.tab-btn').forEach(btn => {{
  btn.addEventListener('click', () => {{
    document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
    document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
    btn.classList.add('active');
    document.getElementById(btn.dataset.tab).classList.add('active');
  }});
}});

// Read Length Chart
new Chart(document.getElementById('lengthChart'), {{
  type: 'bar',
  data: {{
    labels: {len_labels},
    datasets: [{{ label: 'Read Count', data: {len_values},
      backgroundColor: 'rgba(26,82,118,0.6)', borderColor: 'rgba(26,82,118,1)', borderWidth: 1 }}]
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    scales: {{ x: {{ title: {{ display: true, text: 'Read Length (bp)' }} }},
               y: {{ title: {{ display: true, text: 'Count' }} }} }},
    plugins: {{ legend: {{ display: false }} }}
  }}
}});

// Q-Score Chart
new Chart(document.getElementById('qscoreChart'), {{
  type: 'bar',
  data: {{
    labels: {q_labels},
    datasets: [{{ label: 'Count', data: {q_values},
      backgroundColor: 'rgba(46,134,193,0.6)', borderColor: 'rgba(46,134,193,1)', borderWidth: 1 }}]
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    scales: {{ x: {{ title: {{ display: true, text: 'Q-Score' }} }},
               y: {{ title: {{ display: true, text: 'Count' }} }} }},
    plugins: {{ legend: {{ display: false }} }}
  }}
}});

// End Reason Pie
new Chart(document.getElementById('erPieChart'), {{
  type: 'doughnut',
  data: {{
    labels: {er_labels},
    datasets: [{{ data: {er_values},
      backgroundColor: ['#1a5276','#2e86c1','#3498db','#5dade2','#85c1e9','#aed6f1','#d4e6f1'] }}]
  }},
  options: {{ responsive: true, maintainAspectRatio: false }}
}});

// Duration Chart
new Chart(document.getElementById('durationChart'), {{
  type: 'bar',
  data: {{
    labels: {dur_labels},
    datasets: [{{ label: 'Count', data: {dur_values},
      backgroundColor: 'rgba(39,174,96,0.6)', borderColor: 'rgba(39,174,96,1)', borderWidth: 1 }}]
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    scales: {{ x: {{ title: {{ display: true, text: 'Duration (s)' }} }},
               y: {{ title: {{ display: true, text: 'Count' }} }} }},
    plugins: {{ legend: {{ display: false }} }}
  }}
}});

// Scatter Chart
new Chart(document.getElementById('scatterChart'), {{
  type: 'scatter',
  data: {{
    datasets: [{{ label: 'Reads', data: {scatter_data},
      backgroundColor: 'rgba(26,82,118,0.3)', borderColor: 'rgba(26,82,118,0.5)',
      pointRadius: 2 }}]
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    scales: {{ x: {{ title: {{ display: true, text: 'Read Length (bp)' }} }},
               y: {{ title: {{ display: true, text: 'Duration (s)' }} }} }},
    plugins: {{ legend: {{ display: false }} }}
  }}
}});
</script>
</body>
</html>"""
    return html


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Create SMA-seq database from merged BAM and generate QC report"
    )
    parser.add_argument("--bam", type=Path, required=True,
                        help="Merged BAM file")
    parser.add_argument("--sequencing-summary", type=Path, default=None,
                        help="Merged sequencing summary TSV")
    parser.add_argument("--manifest", type=Path, required=True,
                        help="Manifest JSON from sma_scan.py")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Output directory")
    args = parser.parse_args()

    # Load manifest
    manifest = _json.loads(args.manifest.read_text())
    exp_id = manifest.get("experiment_id", "experiment")

    # Parse sequencing summary
    summary_map: dict[str, dict] = {}
    if args.sequencing_summary and args.sequencing_summary.exists():
        print(f"[sma_init] Parsing sequencing summary: {args.sequencing_summary}")
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from sma_basecall import parse_sequencing_summary
        summary_map = parse_sequencing_summary(str(args.sequencing_summary))
        print(f"[sma_init] Loaded {len(summary_map)} summary records")

    # Create database
    db_path = args.output / f"{exp_id}.db"
    print(f"[sma_init] Creating database: {db_path}")
    count = ingest_bam_to_db(args.bam, db_path, exp_id, summary_map)
    print(f"[sma_init] Ingested {count} reads")

    # Compute QC metrics
    print("[sma_init] Computing QC metrics")
    metrics = compute_qc_metrics(db_path)

    # Generate HTML report
    print("[sma_init] Generating QC report")
    html = generate_qc_html(metrics, manifest)
    report_path = args.output / f"{exp_id}_qc_report.html"
    report_path.write_text(html)
    print(f"[sma_init] QC report: {report_path} ({len(html) / 1024:.0f} KB)")

    # Write provenance
    provenance = {
        "bam": str(args.bam),
        "db": str(db_path),
        "report": str(report_path),
        "total_reads": count,
        "total_bases": metrics["total_bases"],
        "mean_qbc": metrics["mean_qbc"],
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    prov_path = args.output / f"{exp_id}_init.json"
    prov_path.write_text(_json.dumps(provenance, indent=2))
    print(f"[sma_init] Done. {count} reads, {metrics['total_bases'] / 1e6:.1f}M bases")


if __name__ == "__main__":
    main()
