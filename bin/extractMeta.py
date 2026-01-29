#!/usr/bin/env python3
# extractMeta.py
# Extract Read ID and End Reason from Pod5 files using pod5 view CLI

import argparse
import sys
import subprocess
from pathlib import Path


############
# argparse #
############

parser = argparse.ArgumentParser(description="Extract end reason from Pod5 files")
parser.add_argument("-i", "--input", default="Input/pod5",
	help="Path to input directory containing .pod5 files [%(default)s]")
parser.add_argument("-o", "--output", default="Input/summary.tsv",
	help="Path to output summary TSV file [%(default)s]")
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
cmd = [
	"pod5", "view",
	str(pod5_dir),          # Input directory
	"--recursive",          # Walk directory tree
	"--include", columns,   # Strict column filtering
	"--output", str(output_tsv),
	"--force-overwrite"     # Overwrite if exists
]

# Join command for display
cmd_str = " ".join(cmd)
print(f"[extractMeta] Running command:\n{cmd_str}")

try:
	subprocess.run(cmd, check=True)
	print(f"[extractMeta] Success: Summary written to {output_tsv}")
except subprocess.CalledProcessError as e:
	sys.exit(f"[extractMeta] Error: pod5 view failed with exit code {e.returncode}")
except FileNotFoundError:
	sys.exit("[extractMeta] Error: 'pod5' executable not found.")
