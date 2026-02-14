#!/usr/bin/env python3
"""Generate interactive HTML report for barcode classification results.

Reads classification data from an SMA database (produced by ingest.py) and
generates a self-contained HTML report with interactive tables, barcode pair
matrices, confidence histograms, and annotated read sequence viewers.
"""

import argparse
import sqlite3
import sys
from pathlib import Path

from sample_sheet import parse_sample_sheet
from report_analysis import analyze_classification
from report_template import generate_html


# ---------------------------------------------------------------------------
# argparse
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Generate barcode classification HTML report from SMA database"
)
parser.add_argument("-d", "--database", required=True,
    help="SMA SQLite database (from mkdb.py + ingest.py)")
parser.add_argument("-ss", "--sample-sheet", required=True,
    help="MinKNOW sample sheet CSV (duplexed barcode format)")
parser.add_argument("-o", "--output", required=True,
    help="Output HTML file path")
parser.add_argument("--max-detail-reads", type=int, default=16,
    help="Max reads shown in Read Details tab [%(default)s]")
parser.add_argument("--flank-front", default=None,
    help="Forward front flanking sequence (e.g., AAGGTTAA)")
parser.add_argument("--flank-rear", default=None,
    help="Forward rear flanking sequence (e.g., CAGCACCT)")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# Load data from database
# ---------------------------------------------------------------------------

DB_PATH = Path(args.database)
SS_PATH = Path(args.sample_sheet)
OUTPUT_PATH = Path(args.output)

if not DB_PATH.exists():
    sys.exit(f"[report] Error: Database not found: {DB_PATH}")
if not SS_PATH.exists():
    sys.exit(f"[report] Error: Sample sheet not found: {SS_PATH}")

print(f"[report] Loading data from {DB_PATH}")

conn = sqlite3.connect(str(DB_PATH))
conn.row_factory = sqlite3.Row

# Load experiment metadata
exp_metadata = {"exp_id": "", "flow_cell_id": "", "sample_id": ""}
try:
    exp_row = conn.execute("SELECT * FROM Exp LIMIT 1").fetchone()
    if exp_row:
        exp_metadata = {
            "exp_id": exp_row["exp_id"],
            "flow_cell_id": exp_row["flow_cell_id"],
            "sample_id": exp_row["sample_id"],
        }
except Exception:
    # Exp table may not exist or be empty; continue with defaults
    pass

# Load reads with barcode data
db_reads = []
for row in conn.execute(
    "SELECT read_id, readseq, readlen, tgt_id, ed, q_bc, q_ld, ER, "
    "bc_start_id, bc_start_ed, bc_start_conf, "
    "bc_end_id, bc_end_ed, bc_end_conf FROM Reads"
):
    db_reads.append(dict(row))

print(f"[report] Loaded {len(db_reads)} reads")

# Load references from Target table
references = {}
try:
    for row in conn.execute("SELECT tgt_id, tgt_refseq, tgt_reflen FROM Target"):
        references[row["tgt_id"]] = (row["tgt_refseq"], row["tgt_reflen"])
    print(f"[report] Loaded {len(references)} target references")
except Exception:
    pass

conn.close()

if not db_reads:
    sys.exit("[report] Error: No reads found in database.")


# ---------------------------------------------------------------------------
# Parse sample sheet
# ---------------------------------------------------------------------------

print(f"[report] Parsing sample sheet: {SS_PATH}")
barcode_pair_to_alias = parse_sample_sheet(SS_PATH)
print(f"[report] {len(barcode_pair_to_alias)} barcode pairs loaded")


# ---------------------------------------------------------------------------
# Analyze
# ---------------------------------------------------------------------------

print("[report] Analyzing classification results...")
analysis = analyze_classification(
    db_reads=db_reads,
    barcode_pair_to_alias=barcode_pair_to_alias,
    max_detail_reads=args.max_detail_reads,
    flank_front=args.flank_front,
    flank_rear=args.flank_rear,
    references=references if references else None,
)

summary = analysis["summary"]
print(f"  Total:       {summary['total']}")
print(f"  Full-length: {summary['full_length']}")
print(f"  Truncated:   {summary['truncated']}")
print(f"  Matched:     {summary['matched']}")
print(f"  Unmatched:   {summary['unmatched']}")


# ---------------------------------------------------------------------------
# Generate HTML
# ---------------------------------------------------------------------------

print("[report] Generating HTML report...")
html = generate_html(analysis, exp_metadata)

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
OUTPUT_PATH.write_text(html)

size_kb = len(html) / 1024
print(f"[report] Written to {OUTPUT_PATH} ({size_kb:.0f} KB)")
