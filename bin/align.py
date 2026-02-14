"""Alignment engine: all-vs-all forward-strand alignment with segmented metrics."""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import edlib

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


def compute_segmented_metrics(cigar: str, ref_len: int) -> dict:
    """Compute identity/contiguity/gap_count within 5', middle, 3' reference thirds.

    Walks the CIGAR, tracking the current reference position to assign each
    operation to the appropriate segment.
    """
    seg1_end = ref_len // 3
    seg2_end = 2 * ref_len // 3
    seg_len = [seg1_end, seg2_end - seg1_end, ref_len - seg2_end]

    # Per-segment accumulators
    seg_matches = [0, 0, 0]
    seg_mismatches = [0, 0, 0]
    seg_gap_count = [0, 0, 0]
    seg_contiguity = [0, 0, 0]
    seg_current_run = [0, 0, 0]

    def _seg_index(ref_pos: int) -> int:
        if ref_pos < seg1_end:
            return 0
        elif ref_pos < seg2_end:
            return 1
        else:
            return 2

    ref_pos = 0
    ops = parse_cigar(cigar)

    for op, length in ops:
        if op in ("=", "M"):
            for _ in range(length):
                if ref_pos < ref_len:
                    s = _seg_index(ref_pos)
                    seg_matches[s] += 1
                    seg_current_run[s] += 1
                    seg_contiguity[s] = max(seg_contiguity[s], seg_current_run[s])
                ref_pos += 1
        elif op == "X":
            for _ in range(length):
                if ref_pos < ref_len:
                    s = _seg_index(ref_pos)
                    seg_mismatches[s] += 1
                    seg_current_run[s] = 0
                ref_pos += 1
        elif op == "I":
            # Insertion: ref_pos doesn't advance. Assign to current segment.
            if ref_pos < ref_len:
                s = _seg_index(ref_pos)
            else:
                s = 2
            seg_gap_count[s] += 1
            seg_current_run[s] = 0
        elif op == "D":
            for i in range(length):
                if ref_pos < ref_len:
                    s = _seg_index(ref_pos)
                    if i == 0:  # count deletion event once per segment
                        seg_gap_count[s] += 1
                    seg_current_run[s] = 0
                ref_pos += 1

    names = ["five_prime", "middle", "three_prime"]
    result = {}
    for i, name in enumerate(names):
        result[f"{name}_identity"] = seg_matches[i] / seg_len[i] if seg_len[i] > 0 else 0.0
        result[f"{name}_contiguity"] = seg_contiguity[i]
        result[f"{name}_gap_count"] = seg_gap_count[i]
    return result


def compute_terminal_metrics(cigar: str, ref_len: int, window: int = 20) -> dict:
    """Compute 5' and 3' alignment quality metrics.

    five_prime_offset: number of ref bases at 5' end covered by deletions
    five_prime_identity_20: identity within first `window` ref bases
    five_prime_first_match: ref position of first match/mismatch
    three_prime_offset: ref bases at 3' end covered by deletions
    three_prime_identity_20: identity within last `window` ref bases
    """
    ops = parse_cigar(cigar)
    win = min(window, ref_len)

    # Build per-ref-position array: 1=match, 0=mismatch, -1=deletion
    ref_status = [-1] * ref_len  # default: uncovered (deletion)
    ref_pos = 0
    for op, length in ops:
        if op in ("=", "M"):
            for _ in range(length):
                if ref_pos < ref_len:
                    ref_status[ref_pos] = 1
                ref_pos += 1
        elif op == "X":
            for _ in range(length):
                if ref_pos < ref_len:
                    ref_status[ref_pos] = 0
                ref_pos += 1
        elif op == "I":
            pass  # no ref consumed
        elif op == "D":
            for _ in range(length):
                if ref_pos < ref_len:
                    ref_status[ref_pos] = -1
                ref_pos += 1

    # 5' offset: leading deletions
    five_prime_offset = 0
    for s in ref_status:
        if s == -1:
            five_prime_offset += 1
        else:
            break

    # 3' offset: trailing deletions
    three_prime_offset = 0
    for s in reversed(ref_status):
        if s == -1:
            three_prime_offset += 1
        else:
            break

    # 5' first match position
    five_prime_first_match = ref_len
    for i, s in enumerate(ref_status):
        if s >= 0:
            five_prime_first_match = i
            break

    # 5' identity within window
    five_window = ref_status[:win]
    five_prime_identity = sum(1 for s in five_window if s == 1) / len(five_window) if five_window else 0.0

    # 3' identity within window
    three_window = ref_status[-win:] if win > 0 else []
    three_prime_identity = sum(1 for s in three_window if s == 1) / len(three_window) if three_window else 0.0

    return {
        "five_prime_offset": five_prime_offset,
        "three_prime_offset": three_prime_offset,
        "five_prime_identity_20": five_prime_identity,
        "three_prime_identity_20": three_prime_identity,
        "five_prime_first_match": five_prime_first_match,
    }


def compute_indel_positions(cigar: str, min_size: int = 3) -> list[tuple[str, int, int]]:
    """Return list of (type, size, ref_position) for indels >= min_size."""
    ops = parse_cigar(cigar)
    result = []
    ref_pos = 0

    for op, length in ops:
        if op in ("=", "M", "X"):
            ref_pos += length
        elif op == "I":
            if length >= min_size:
                result.append(("I", length, ref_pos))
        elif op == "D":
            if length >= min_size:
                result.append(("D", length, ref_pos))
            ref_pos += length

    return result


METRIC_COLUMNS = [
    "read_id", "ref_id", "read_len", "ref_len",
    "ed", "ned", "identity", "ref_coverage", "read_to_ref_ratio",
    "seg5_identity", "seg5_contiguity",
    "segM_identity", "segM_contiguity",
    "seg3_identity", "seg3_contiguity",
    "max_ins", "max_del", "n_sig_indels",
    "five_prime_offset", "five_prime_identity_20",
    "three_prime_offset", "three_prime_identity_20",
    "rank", "margin",
]


def align_read_to_ref(read_seq: str, ref_seq: str, ref_id: str) -> dict:
    """Align a single read to a single reference (forward-strand, NW mode).

    Returns a dict with all metric columns except rank/margin (added during classification).
    """
    result = edlib.align(read_seq, ref_seq, mode="NW", task="path")
    cigar = result.get("cigar", "")
    ed = result.get("editDistance", len(ref_seq))

    ref_len = len(ref_seq)
    read_len = len(read_seq)

    cigar_m = compute_cigar_metrics(cigar, ref_len, read_len)
    seg_m = compute_segmented_metrics(cigar, ref_len)
    term_m = compute_terminal_metrics(cigar, ref_len)

    return {
        "read_id": "",
        "ref_id": ref_id,
        "read_len": read_len,
        "ref_len": ref_len,
        "ed": ed,
        "ned": ed / ref_len if ref_len > 0 else float("inf"),
        "identity": cigar_m["identity"],
        "ref_coverage": cigar_m["ref_coverage"],
        "read_to_ref_ratio": read_len / ref_len if ref_len > 0 else 0.0,
        "seg5_identity": seg_m["five_prime_identity"],
        "seg5_contiguity": seg_m["five_prime_contiguity"],
        "segM_identity": seg_m["middle_identity"],
        "segM_contiguity": seg_m["middle_contiguity"],
        "seg3_identity": seg_m["three_prime_identity"],
        "seg3_contiguity": seg_m["three_prime_contiguity"],
        "max_ins": cigar_m["max_ins"],
        "max_del": cigar_m["max_del"],
        "n_sig_indels": cigar_m["n_sig_indels"],
        "five_prime_offset": term_m["five_prime_offset"],
        "five_prime_identity_20": term_m["five_prime_identity_20"],
        "three_prime_offset": term_m["three_prime_offset"],
        "three_prime_identity_20": term_m["three_prime_identity_20"],
        "rank": 0,
        "margin": 0.0,
    }


def classify_read(
    read_seq: str, read_id: str, refs: dict[str, str],
) -> dict:
    """Align read against all references and classify by best NED.

    Returns dict with: read_id, read_len, assigned_ref, best_ned, second_ned,
    margin, confidence_flag, pairwise (list of per-ref metric dicts).
    """
    pairwise = []
    for ref_id, ref_seq in refs.items():
        m = align_read_to_ref(read_seq, ref_seq, ref_id)
        m["read_id"] = read_id
        pairwise.append(m)

    pairwise.sort(key=lambda x: x["ned"])
    for i, m in enumerate(pairwise):
        m["rank"] = i + 1

    best_ned = pairwise[0]["ned"]
    second_ned = pairwise[1]["ned"] if len(pairwise) > 1 else float("inf")
    margin = second_ned - best_ned

    for m in pairwise:
        m["margin"] = margin

    if best_ned > 0.5:
        flag = "POOR"
    elif margin <= 0.1:
        flag = "LOW"
    else:
        flag = "HIGH"

    return {
        "read_id": read_id,
        "read_len": len(read_seq),
        "assigned_ref": pairwise[0]["ref_id"],
        "best_ned": best_ned,
        "second_ned": second_ned,
        "margin": margin,
        "confidence_flag": flag,
        "pairwise": pairwise,
    }
