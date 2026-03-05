# meta.py
# Extracts read_id and end_reason metadata from Pod5 files.

import shutil
import subprocess
from pathlib import Path


def run(pod5_input, tsv_output):
	"""
	Wraps 'pod5 view' to extract read_id and end_reason into a TSV.
	"""
	pod5_dir = Path(pod5_input)
	output_tsv = Path(tsv_output)

	print(f"[meta] Input Directory: {pod5_dir}")
	print(f"[meta] Output TSV:      {output_tsv}")

	# --- Validation ---
	if not pod5_dir.exists():
		raise FileNotFoundError(f"[meta] Input directory does not exist: {pod5_dir}")

	if not pod5_dir.is_dir():
		raise NotADirectoryError(f"[meta] Not a directory: {pod5_dir}")

	if shutil.which("pod5") is None:
		raise RuntimeError("[meta] 'pod5' command not found. Please install pod5-file-format.")

	if not output_tsv.parent.exists():
		print(f"[meta] Creating output directory: {output_tsv.parent}")
		output_tsv.parent.mkdir(parents=True, exist_ok=True)

	if not any(pod5_dir.rglob("*.pod5")):
		raise RuntimeError(f"[meta] No .pod5 files found in {pod5_dir}")

	columns = "read_id,end_reason"

	cmd = [
		"pod5", "view",
		str(pod5_dir),
		"--recursive",
		"--include", columns,   # limit to required columns
		"--output", str(output_tsv),
		"--force-overwrite"
	]

	print(f"[meta] Executing: {' '.join(cmd)}")

	try:
		subprocess.run(cmd, check=True)
		print(f"[meta] Success: Metadata extracted to {output_tsv}")
	except subprocess.CalledProcessError as e:
		raise RuntimeError(
			f"[meta] pod5 view failed with exit code {e.returncode}.") from e
	except Exception as e:
		raise RuntimeError(f"[meta] Critical Error: {e}") from e
