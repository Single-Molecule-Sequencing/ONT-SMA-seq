"""Alignment engine: all-vs-all forward-strand alignment with segmented metrics."""
from __future__ import annotations

import re
from dataclasses import dataclass, field

CIGAR_RE = re.compile(r"(\d+)([=XIDM])")
SIG_INDEL_THRESHOLD = 5


def parse_cigar(cigar: str) -> list[tuple[str, int]]:
    """Parse an edlib extended CIGAR string into (op, length) tuples."""
    return [(op, int(length)) for length, op in CIGAR_RE.findall(cigar)]


def compute_cigar_metrics(cigar: str, ref_len: int, read_len: int) -> dict:
    """Compute whole-alignment metrics from a CIGAR string.

    Returns dict with: matches, mismatches, insertions, deletions,
    identity, ref_coverage, max_ins, max_del, n_sig_indels, alignment_length.
    """
    ops = parse_cigar(cigar)
    matches = 0
    mismatches = 0
    ins_total = 0
    del_total = 0
    max_ins = 0
    max_del = 0
    n_sig_indels = 0
    ref_bases_covered = 0

    for op, length in ops:
        if op == "=" or op == "M":
            matches += length
            ref_bases_covered += length
        elif op == "X":
            mismatches += length
            ref_bases_covered += length
        elif op == "I":
            ins_total += length
            max_ins = max(max_ins, length)
            if length >= SIG_INDEL_THRESHOLD:
                n_sig_indels += 1
        elif op == "D":
            del_total += length
            ref_bases_covered += length
            max_del = max(max_del, length)
            if length >= SIG_INDEL_THRESHOLD:
                n_sig_indels += 1

    alignment_length = matches + mismatches + ins_total + del_total
    identity = matches / alignment_length if alignment_length > 0 else 0.0
    ref_coverage = (matches + mismatches) / ref_len if ref_len > 0 else 0.0

    return {
        "matches": matches,
        "mismatches": mismatches,
        "insertions": ins_total,
        "deletions": del_total,
        "alignment_length": alignment_length,
        "identity": identity,
        "ref_coverage": ref_coverage,
        "max_ins": max_ins,
        "max_del": max_del,
        "n_sig_indels": n_sig_indels,
    }
