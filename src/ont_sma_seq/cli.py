# cli.py
# Primary command-line interface for ONT-SMA-Seq.

import argparse
import sys
import yaml
from ont_sma_seq import mkdb
from ont_sma_seq import init
from ont_sma_seq import meta
from ont_sma_seq import ingest
from ont_sma_seq import merge

_REQUIRED_KEYS = ["bam", "pod5_dir", "ref", "outdir", "summary_tsv"]


def _resolve_metadata(cfg):
	"""
	Resolves experiment metadata (exp_id, model_tier, model_ver, trim, mods_bitflag)
	from config and/or BAM filename.

	Raises:
		ValueError: If any required metadata field is missing or invalid.
	"""
	from pathlib import Path

	bam_path = Path(cfg["bam"])

	# Attempt filename parse unconditionally — needed for any blank fields.
	parsed = None
	try:
		parsed = ingest.parse_bam_filename(bam_path)
	except ValueError:
		pass

	p_exp_id, p_tier, p_ver, p_trim, p_mods = parsed if parsed else (None, None, None, None, None)

	def _pick(cfg_key, parsed_val):
		"""Return config value if non-blank, else parsed value, else raise."""
		val = cfg.get(cfg_key)
		if val is not None and str(val).strip():
			return str(val).strip()
		if parsed_val is not None:
			return str(parsed_val)
		raise ValueError(
			f"'{cfg_key}' is not set in config.yml and could not be derived from "
			f"the BAM filename. Either use the standard naming scheme or set "
			f"'{cfg_key}' explicitly in config.yml."
		)

	exp_id = _pick("exp_id",      p_exp_id)
	tier   = _pick("model_tier",  p_tier)
	ver    = _pick("model_ver",   p_ver)

	# trim: config accepts "yes"/"no"/"1"/"0"; BAM filename gives 0/1.
	trim_raw = cfg.get("trim")
	if trim_raw is not None and str(trim_raw).strip():
		t = str(trim_raw).strip().lower()
		if t in ("yes", "true", "1"):
			trim = 1
		elif t in ("no", "false", "0"):
			trim = 0
		else:
			raise ValueError(
				f"'trim' value '{trim_raw}' is invalid. Use 'yes'/'no' or 1/0."
			)
	elif p_trim is not None:
		trim = int(p_trim)
	else:
		raise ValueError(
			"'trim' is not set in config.yml and could not be derived from the BAM filename."
		)

	# mods_bitflag: normalize to int.
	mods_raw = cfg.get("mods_bitflag")
	if mods_raw is not None and str(mods_raw).strip():
		try:
			mods = int(str(mods_raw).strip())
		except ValueError:
			raise ValueError(
				f"'mods_bitflag' value '{mods_raw}' is not a valid integer."
			)
	elif p_mods is not None:
		mods = int(p_mods)
	else:
		raise ValueError(
			"'mods_bitflag' is not set in config.yml and could not be derived from the BAM filename."
		)

	return exp_id, tier, ver, trim, mods


def _step(name, fn):
	"""Run a pipeline step, re-raising any exception labelled with the step name."""
	print(f"\n[run] --- {name} ---")
	try:
		return fn()
	except Exception as e:
		raise RuntimeError(f"Step '{name}' failed: {e}") from e


def run_pipeline(config_path):
	"""
	Loads config_path (YAML) and chains mkdb → init → meta → ingest.

	Raises:
		FileNotFoundError: If config_path does not exist.
		ValueError: If required config keys are missing.
		RuntimeError: Wraps any step-level failure with its step name.
	"""
	with open(config_path, 'r') as f:
		cfg = yaml.safe_load(f)

	missing = [k for k in _REQUIRED_KEYS if k not in cfg]
	if missing:
		raise ValueError(f"Config missing required keys: {missing}")

	# Resolve experiment metadata early — exits immediately if any field is
	# unresolvable before any DB work begins.
	exp_id, tier, ver, trim, mods = _resolve_metadata(cfg)

	print(f"[run] Pipeline: {exp_id}")
	print(f"[run] Config:   {config_path}")
	print(f"[run] Metadata: tier={tier}  ver={ver}  trim={trim}  mods={mods}")

	db_path = _step("mkdb",   lambda: mkdb.run(exp_id=exp_id, out_dir=cfg["outdir"]))
	_step("init",   lambda: init.run(db_path=db_path, ref_fasta=cfg["ref"]))
	_step("meta",   lambda: meta.run(pod5_input=cfg["pod5_dir"], tsv_output=cfg["summary_tsv"]))
	_step("ingest", lambda: ingest.run(
		bam_path=cfg["bam"], db_path=db_path, meta_path=cfg["summary_tsv"],
		exp_id=exp_id, tier=tier, ver=ver, trim=trim, mods=mods,
		len_min_mult=cfg.get("len_min_mult", ingest.LENGTH_MIN_MULT),
		len_max_mult=cfg.get("len_max_mult", ingest.LENGTH_MAX_MULT),
	))

	print(f"\n[run] Pipeline complete.")
	print(f"  DB: {db_path}")


def main():
	"""Primary entry point for the ont-sma CLI."""
	parser = argparse.ArgumentParser(
		description="ONT Single-Molecule-Accuracy Sequencing (ONT-SMA-Seq)"
	)
	subparsers = parser.add_subparsers(dest='command', help='Available commands')

	# --- mkdb ---
	parser_mkdb = subparsers.add_parser('mkdb', help='Initialize the SQLite database schema')
	parser_mkdb.add_argument("-e", "--expid", required=True,
		help="Experiment ID (Format: FlowCell_Sample_Alias)")
	parser_mkdb.add_argument("-o", "--outdir", default="Output",
		help="Output directory [%(default)s]")

	# --- init ---
	parser_init = subparsers.add_parser('init', help='Lock reference FASTA into the database')
	parser_init.add_argument("-d", "--db", required=True,
		help="Path to existing SQLite database")
	parser_init.add_argument("-r", "--ref", required=True,
		help="Path to reference FASTA")

	# --- meta ---
	parser_meta = subparsers.add_parser('meta', help='Extract end_reason metadata from Pod5 files')
	parser_meta.add_argument("-i", "--input", required=True,
		help="Input directory containing .pod5 files")
	parser_meta.add_argument("-o", "--output", default="Output/summary.tsv",
		help="Output TSV path [%(default)s]")

	# --- ingest ---
	parser_ingest = subparsers.add_parser('ingest', help='Stream BAM, calculate metrics, insert to DB')
	parser_ingest.add_argument("-b", "--bam", required=True, help="Input BAM file")
	parser_ingest.add_argument("-d", "--db", required=True, help="Target database")
	parser_ingest.add_argument("-m", "--meta", required=True, help="Metadata TSV")
	parser_ingest.add_argument("--len-min", type=float, default=ingest.LENGTH_MIN_MULT,
		metavar="MULT",
		help="Min read length as multiplier of tgt_reflen [%(default)s]")
	parser_ingest.add_argument("--len-max", type=float, default=ingest.LENGTH_MAX_MULT,
		metavar="MULT",
		help="Max read length as multiplier of tgt_reflen [%(default)s]")

	# --- merge ---
	parser_merge = subparsers.add_parser('merge', help='Merge per-run databases into a master DB')
	parser_merge.add_argument("-o", "--output", required=True,
		help="Path to master output database")
	parser_merge.add_argument("inputs", nargs="+",
		help="One or more source .db files to merge")

	# --- run ---
	parser_run = subparsers.add_parser('run', help='Run full pipeline from a config.yml')
	parser_run.add_argument("-c", "--config", default="config.yml",
		help="Path to config file [%(default)s]")

	args = parser.parse_args()

	if args.command is None:
		parser.print_help()
		return

	try:
		if args.command == 'mkdb':
			mkdb.run(exp_id=args.expid, out_dir=args.outdir)
		elif args.command == 'init':
			init.run(db_path=args.db, ref_fasta=args.ref)
		elif args.command == 'meta':
			meta.run(pod5_input=args.input, tsv_output=args.output)
		elif args.command == 'ingest':
			ingest.run(bam_path=args.bam, db_path=args.db, meta_path=args.meta,
			           len_min_mult=args.len_min, len_max_mult=args.len_max)
		elif args.command == 'merge':
			merge.run(output_db=args.output, input_dbs=args.inputs)
		elif args.command == 'run':
			run_pipeline(args.config)
	except Exception as e:
		print(f"\n[{args.command}] Error: {e}", file=sys.stderr)
		sys.exit(1)


if __name__ == "__main__":
	main()
