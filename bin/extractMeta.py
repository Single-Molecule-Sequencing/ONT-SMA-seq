#!/usr/bin/env python3
# extractMeta.py
# Extract ER metadata from Pod5 files using pod5 view CLI for fast parsing

import argparse
import sys
import subprocess
from pathlib import Path


############
# argparse #
############

parser = argparse.ArgumentParser(description="Extract Metadata from Pod5 files")
parser.add_argument("-i", "--input", required=True, 
	help="Path to input directory containing .pod5 files")
parser.add_argument("-o", "--output", default="./pod5_ER_summary.tsv",
	help="Path to output summary TSV file")
args = parser.parse_args()


########
# path #
########

pod5_dir   = Path(args.input)
output_tsv = Path(args.output)

# Input Validation
if not pod5_dir.exists():
	sys.exit(f"[extractMeta] Error: Input directory {pod5_dir} does not exist.")

if not pod5_dir.is_dir():
	sys.exit(f"[extractMeta] Error: {pod5_dir} is not a directory.")

# Output Validation
if not output_tsv.parent.exists():
	print(f"[extractMeta] Warning: Output directory {output_tsv.parent} does not exist. Creating it.")
	output_tsv.parent.mkdir(parents=True, exist_ok=True)


#############
# execution #
#############

print(f"[extractMeta] Input Directory: {pod5_dir}")
print(f"[extractMeta] Output TSV:      {output_tsv}")

# Columns to extract
columns = "read_id,end_reason"

# Construct the command
# Note: Using --recursive to avoid shell globbing (*.pod5)
cmd = [
	"pod5", "view",
	str(pod5_dir),          # Input directory
	"--recursive",          # Find all .pod5 files inside
	"--include", columns,   # Strict column filtering
	"--output", str(output_tsv),
	"--force-overwrite"     # Overwrite if exists
]

# Join command for display
cmd_str = " ".join(cmd)
print(f"[extractMeta] Running command:\n{cmd_str}")

try:
	subprocess.run(cmd, shell=False, check=True)
except subprocess.CalledProcessError as e:
	sys.exit(f"[extractMeta] Error: pod5 view failed with exit code {e.returncode}")
except Exception as e:
	sys.exit(f"[extractMeta] Unexpected error: {e}")

print(f"[extractMeta] Summary saved to {output_tsv}")
