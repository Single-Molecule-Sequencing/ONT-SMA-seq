#!/usr/bin/env python3
# inputInit.py
# parses BAM filename for metadata and standardize symlinks

import argparse
import os
import sys
from pathlib import Path


#########
# funcs #
#########

def safe_symlink(src, dst):
	# Resolve absolute paths to avoid broken links if CWD changes
	src = os.path.abspath(src)
	dst = os.path.abspath(dst)
	
	# Check if dst exists (link or file)
	if os.path.islink(dst) or os.path.exists(dst):
		os.remove(dst) # Clean overwrite
		
	os.symlink(src, dst)
	print(f"Symlinked: {dst} from {src}")


############
# argparse #
############

parser = argparse.ArgumentParser(description="Standardize Inputs and Symlink")
parser.add_argument("-b", "--bam", required=True, 
	help="Path to raw uBAM file")
parser.add_argument("-p", "--pod5_dir", required=True,
	help="Path to raw Pod5 directory")
parser.add_argument("-r", "--ref", required=True, 
	help="Path to Reference FASTA")
args = parser.parse_args()


##############
# inputcheck #
##############

if not os.path.isfile(args.bam):
	sys.exit(f"[inputInit] Error: BAM file {args.bam} does not exist.")
if not os.path.isdir(args.pod5_dir):
	sys.exit(f"[inputInit] Error: Pod5 directory {args.pod5_dir} does not exist.")
if not os.path.isfile(args.ref):
	sys.exit(f"[inputInit] Error: Reference FASTA {args.ref} does not exist.")


############
# metadata #
############

bam_filename = os.path.basename(args.bam)
filename_root, ext = os.path.splitext(bam_filename)
parts = filename_root.split("_")

# Format: {exp_id}_{bc_model_type}_v{bc_model_version}_{trim}_{modifications}.bam
if len(parts) < 5:
	sys.exit(f"[inputInit] Error: BAM filename {bam_filename} does not conform to expected format.")

# Parsing based on your strict convention
exp_id = parts[0]
tier = parts[1]
ver = parts[2].lstrip("v")
trim = parts[3][-1] # Extract '0' from 'trim0' or '1' from 'trim1'
mods = parts[4]

# --- Validation ---
if tier not in ["fast", "hac", "sup"]:
	sys.exit(f"[inputInit] Error: Invalid model tier '{tier}' in BAM filename.")

if not all(x.isdigit() for x in ver.split(".")):
	sys.exit(f"[inputInit] Error: Invalid model version '{ver}' in BAM filename.")

if trim not in ["0", "1"]:
	sys.exit(f"[inputInit] Error: Invalid trim status '{trim}' (must be 0 or 1).")

valid_mod_flags = {"0", "1", "2", "4", "8", "16", "3", "5", "9", "17"}
if mods not in valid_mod_flags:
	sys.exit(f"[inputInit] Error: Invalid modifications '{mods}' in BAM filename.")

print(f"[inputInit] Metadata Parsed: Exp={exp_id}, Tier={tier}, Ver={ver}, Trim={trim}, Mods={mods}")


###########
# linking #
###########

input_dir = Path("Input")
input_dir.mkdir(exist_ok=True)

# Symlink BAM (Keep original name so ingest.py can re-parse it if needed)
target_bam = input_dir / bam_filename
safe_symlink(args.bam, target_bam)

# Symlink Pod5 Dir (Standardized name: {exp_id}_pod5)
target_pod5 = input_dir / f"{exp_id}_pod5"
safe_symlink(args.pod5_dir, target_pod5)

# Symlink RefSeq (Standardized name: {exp_id}.fa)
target_ref = input_dir / f"{exp_id}.fa"
safe_symlink(args.ref, target_ref)

print("[inputInit] Input initialization complete.")
