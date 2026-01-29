#!/usr/bin/env python3
# inputInit.py
# Standardize Inputs: Symlink BAM/Pod5 and Sanitize Reference FASTA

import argparse
import os
import sys
from pathlib import Path


#########
# funcs #
#########

def safe_symlink(src, dst):
	"""
	Creates a symbolic link from src to dst.
	Removes existing file/link at dst if it exists.
	"""
	src = os.path.abspath(src)
	dst = os.path.abspath(dst)
	
	# Check if dst exists (link or file)
	if os.path.islink(dst) or os.path.exists(dst):
		os.remove(dst) # Clean overwrite
		
	os.symlink(src, dst)
	print(f"[inputInit] Symlinked: {dst} -> {src}")


def extract_first_sequence(src, dst):
	"""
	Reads the source FASTA and writes ONLY the first sequence to dst.
	Ignores any subsequent sequences. Raises error if no sequence found.
	"""
	found_header = False

	with open(src, 'r') as f_in, open(dst, 'w') as f_out:
		for line in f_in:
			line_strip = line.strip()
			if not line_strip: continue

			if line.startswith('>'):
				if found_header:
					print(f"[inputInit] Warning: Multiple sequences detected in source. Using only the first one.")
					break
				found_header = True
				f_out.write(line)
			elif found_header:
				f_out.write(line)

	if not found_header:
		sys.exit(f"[inputInit] Error: No sequence found in Reference FASTA {src}.")

	print(f"[inputInit] Extracted target sequence to {dst}")


############
# argparse #
############

parser = argparse.ArgumentParser(description="Standardize Inputs, Symlink BAM/Pod5, Sanitize Ref")
parser.add_argument("-b", "--bam", required=True,
	help="Path to raw uBAM file")
parser.add_argument("-p", "--pod5_dir", required=True,
	help="Path to raw Pod5 directory")
parser.add_argument("-r", "--ref", required=True,
	help="Path to Target Reference FASTA")
args = parser.parse_args()


##############
# inputcheck #
##############

bam_path = Path(args.bam)
pod5_dir = Path(args.pod5_dir)
ref_path = Path(args.ref)

if not bam_path.is_file():
	sys.exit(f"[inputInit] Error: BAM file {bam_path} does not exist.")
if not pod5_dir.is_dir():
	sys.exit(f"[inputInit] Error: Pod5 directory {pod5_dir} does not exist.")
if not ref_path.is_file():
	sys.exit(f"[inputInit] Error: Reference FASTA {ref_path} does not exist.")


############
# metadata #
############

parts = bam_path.stem.split('_')

if len(parts) < 5:
	sys.exit(f"[inputInit] Error: BAM filename '{bam_path.name}' does not conform to format: "
			 "{exp_id}_{tier}_v{ver}_{trim}_{mods}.bam")

mods = parts[-1]
trim_part = parts[-2]
ver_str = parts[-3]
tier = parts[-4]
exp_id = "_".join(parts[:-4])

# --- Meta Validation ---

# 1. Tier
if tier not in ["fast", "hac", "sup"]:
	sys.exit(f"[inputInit] Error: Invalid model tier '{tier}' in BAM filename.")

# 2. Version
ver = ver_str.lstrip("v")
if not all(x.isdigit() for x in ver.split(".")):
	sys.exit(f"[inputInit] Error: Invalid model version '{ver_str}' in BAM filename.")

# 3. Trim
trim = trim_part[-1]
if trim not in ["0", "1"]:
	sys.exit(f"[inputInit] Error: Invalid trim status '{trim_part}' (must end in 0 or 1).")

# 4. Mods
valid_mod_flags = {"0", "1", "2", "4", "8", "16", "3", "5", "9", "17"}
if mods not in valid_mod_flags:
	sys.exit(f"[inputInit] Error: Invalid modifications '{mods}' in BAM filename.")

print(f"[inputInit] Metadata Parsed:")
print(f"  - Exp ID: {exp_id}")
print(f"  - Tier:   {tier}")
print(f"  - Ver:    {ver}")
print(f"  - Trim:   {trim}")
print(f"  - Mods:   {mods}")


###########
# linking #
###########

input_dir = Path("Input")
input_dir.mkdir(exist_ok=True)

# 1. BAM -> Input/reads.bam
safe_symlink(bam_path, input_dir / "reads.bam")

# 2. Pod5 Dir -> Input/pod5/
safe_symlink(pod5_dir, input_dir / "pod5")

# 3. Ref FASTA -> Input/target.fa
extract_first_sequence(ref_path, input_dir / "target.fa")

print("[inputInit] Initialization complete.")
