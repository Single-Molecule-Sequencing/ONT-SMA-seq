"""Pure analysis module for barcode classification reports.

Processes SQLite DB results from ingest.py into structured data suitable for
HTML report generation. No file I/O -- all data is received as arguments and
returned as structured dicts.
"""

from __future__ import annotations

import collections
from typing import Any

import edlib

from barcodes import BARCODES, BARCODE_LENGTH, reverse_complement

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_FULL_LENGTH_THRESHOLD = 0.75
_FLANK_SEARCH_WINDOW = 50
_READ_ID_TRUNCATE = 12
_CONFIDENCE_THRESHOLDS = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _compute_n50(lengths: list[int]) -> int:
    """Compute N50 from a list of read lengths.

    Sort descending, walk cumulative sum; N50 is the length at which the
    cumulative sum reaches 50% of the total bases.
    """
    if not lengths:
        return 0
    sorted_desc = sorted(lengths, reverse=True)
    total_bases = sum(sorted_desc)
    half = total_bases / 2.0
    cumulative = 0
    for length in sorted_desc:
        cumulative += length
        if cumulative >= half:
            return length
    return sorted_desc[-1]


def _compute_median(values: list[float | int]) -> float:
    """Compute median by sorting and taking the middle value."""
    if not values:
        return 0.0
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2 == 1:
        return float(s[mid])
    return (s[mid - 1] + s[mid]) / 2.0


def classify_barcode_detailed(
    segment: str,
    expected_barcodes: dict[str, str],
) -> list[dict[str, Any]]:
    """Classify a read segment against all expected barcodes with full detail.

    Like ``barcodes.classify_barcode`` but returns results for ALL barcodes,
    sorted by edit distance (ascending), with alignment positions.

    Parameters
    ----------
    segment : str
        The read segment to classify.
    expected_barcodes : dict[str, str]
        Mapping of barcode_id -> sequence to search against.

    Returns
    -------
    list[dict]
        Sorted list (best first) of dicts with keys:
        ``bc_id``, ``edit_distance``, ``confidence``, ``start``, ``end``.
    """
    results: list[dict[str, Any]] = []

    for bc_id, bc_seq in expected_barcodes.items():
        aln = edlib.align(bc_seq, segment, mode="HW", task="locations")
        ed: int = aln["editDistance"]
        confidence = max(0.0, 1.0 - ed / BARCODE_LENGTH)

        # Extract best alignment location
        locations = aln.get("locations", [])
        if locations:
            start, end = locations[0]
        else:
            start, end = 0, 0

        results.append({
            "bc_id": bc_id,
            "edit_distance": ed,
            "confidence": confidence,
            "start": start,
            "end": end,
        })

    results.sort(key=lambda r: (r["edit_distance"], r["bc_id"]))
    return results


def find_flank_position(
    read_seq: str,
    flank_seq: str,
    region_start: int,
    region_end: int,
) -> dict[str, int] | None:
    """Find a flanking sequence via edlib HW alignment within a region.

    Parameters
    ----------
    read_seq : str
        Full read sequence.
    flank_seq : str
        Flanking sequence to search for.
    region_start : int
        Start position of the search region (inclusive).
    region_end : int
        End position of the search region (exclusive).

    Returns
    -------
    dict or None
        Dict with ``start``, ``end``, ``edit_distance`` if found, else None.
    """
    region = read_seq[region_start:region_end]
    if not region or not flank_seq:
        return None

    aln = edlib.align(flank_seq, region, mode="HW", task="locations")
    ed: int = aln["editDistance"]

    # Reject poor matches (more than half the flank is errors)
    if ed > len(flank_seq) // 2:
        return None

    locations = aln.get("locations", [])
    if not locations:
        return None

    start, end = locations[0]
    return {
        "start": region_start + start,
        "end": region_start + end + 1,
        "edit_distance": ed,
    }


# ---------------------------------------------------------------------------
# Main analysis function
# ---------------------------------------------------------------------------


def analyze_classification(
    db_reads: list[dict],
    barcode_pair_to_alias: dict[tuple[str, str], str],
    max_detail_reads: int = 16,
    flank_front: str | None = None,
    flank_rear: str | None = None,
    segment_len: int = 100,
    references: dict[str, tuple[str, int]] | None = None,
) -> dict[str, Any]:
    """Analyze barcode classification results from the SQLite database.

    Parameters
    ----------
    db_reads : list[dict]
        List of read dicts from the Reads table. Each dict has keys:
        read_id, readseq, readlen, tgt_id, ed, q_bc, q_ld, ER,
        bc_start_id, bc_start_ed, bc_start_conf,
        bc_end_id, bc_end_ed, bc_end_conf.
    barcode_pair_to_alias : dict[tuple[str, str], str]
        Mapping from (upstream_barcode, downstream_barcode) to alias.
    max_detail_reads : int
        Maximum number of detailed read records to include.
    flank_front : str or None
        Front flanking sequence for region annotation.
    flank_rear : str or None
        Rear flanking sequence for region annotation.
    segment_len : int
        Length of segment to search for barcodes at read start/end.
    references : dict[str, tuple[str, int]] or None
        Mapping of tgt_id -> (sequence, length). If None, omitted from output.

    Returns
    -------
    dict
        Analysis results with keys: summary, per_target_stats,
        confidence_distributions, pair_matrix, reads_table,
        detailed_reads, barcode_info, pairing_table, references.
    """
    result: dict[str, Any] = {}

    # Build expected barcodes from sample sheet
    used_bc_ids: set[str] = set()
    for bc_up, bc_down in barcode_pair_to_alias:
        used_bc_ids.add(bc_up)
        used_bc_ids.add(bc_down)

    expected_barcodes = {
        bc_id: BARCODES[bc_id] for bc_id in used_bc_ids if bc_id in BARCODES
    }
    expected_barcodes_rc = {
        bc_id: reverse_complement(seq) for bc_id, seq in expected_barcodes.items()
    }

    # ------------------------------------------------------------------
    # summary
    # ------------------------------------------------------------------
    total = len(db_reads)
    full_length = 0
    matched = 0
    unmatched = 0

    for rd in db_reads:
        bc_end_conf = rd.get("bc_end_conf")
        if bc_end_conf is not None and bc_end_conf >= _FULL_LENGTH_THRESHOLD:
            full_length += 1

        tgt_id = rd.get("tgt_id", "")
        if tgt_id.startswith("unmatched_"):
            unmatched += 1
        else:
            matched += 1

    truncated = total - full_length

    result["summary"] = {
        "total": total,
        "full_length": full_length,
        "truncated": truncated,
        "matched": matched,
        "unmatched": unmatched,
    }

    # ------------------------------------------------------------------
    # truncation_counts
    # ------------------------------------------------------------------
    truncation_counts: dict[str, int] = {}
    for r in db_reads:
        level = r.get("trunc_level")
        if level:
            truncation_counts[level] = truncation_counts.get(level, 0) + 1

    result["truncation_counts"] = truncation_counts

    # ------------------------------------------------------------------
    # per_target_stats
    # ------------------------------------------------------------------
    target_groups: dict[str, list[dict]] = collections.defaultdict(list)

    for rd in db_reads:
        tgt_id = rd.get("tgt_id", "")
        if tgt_id.startswith("unmatched_"):
            target_groups["unmatched"].append(rd)
        else:
            target_groups[tgt_id].append(rd)

    per_target_stats: dict[str, dict[str, Any]] = {}

    for tgt_id, reads in target_groups.items():
        count = len(reads)
        lengths = [r["readlen"] for r in reads]
        start_confs = [
            r["bc_start_conf"] for r in reads if r.get("bc_start_conf") is not None
        ]
        end_confs = [
            r["bc_end_conf"] for r in reads if r.get("bc_end_conf") is not None
        ]
        eds = [r["ed"] for r in reads if r.get("ed") is not None]
        q_lds = [r["q_ld"] for r in reads if r.get("q_ld") is not None]

        fl_count = sum(
            1 for r in reads
            if r.get("bc_end_conf") is not None
            and r["bc_end_conf"] >= _FULL_LENGTH_THRESHOLD
        )

        per_target_stats[tgt_id] = {
            "count": count,
            "avg_length": sum(lengths) / count if count else 0.0,
            "avg_start_conf": (
                sum(start_confs) / len(start_confs) if start_confs else 0.0
            ),
            "avg_end_conf": (
                sum(end_confs) / len(end_confs) if end_confs else 0.0
            ),
            "full_length_count": fl_count,
            "full_length_pct": (fl_count / count * 100.0) if count else 0.0,
            "avg_ed": sum(eds) / len(eds) if eds else 0.0,
            "avg_q_ld": sum(q_lds) / len(q_lds) if q_lds else 0.0,
            "min_length": min(lengths) if lengths else 0,
            "max_length": max(lengths) if lengths else 0,
            "median_length": _compute_median(lengths),
            "lengths": lengths,
            "start_confs": start_confs,
            "end_confs": end_confs,
        }

    result["per_target_stats"] = per_target_stats

    # ------------------------------------------------------------------
    # confidence_distributions
    # ------------------------------------------------------------------
    start_confs_all: list[float] = []
    end_confs_all: list[float] = []

    for rd in db_reads:
        if rd.get("bc_start_conf") is not None:
            start_confs_all.append(rd["bc_start_conf"])
        if rd.get("bc_end_conf") is not None:
            end_confs_all.append(rd["bc_end_conf"])

    result["confidence_distributions"] = {
        "start_confs": start_confs_all,
        "end_confs": end_confs_all,
    }

    # ------------------------------------------------------------------
    # pair_matrix
    # ------------------------------------------------------------------
    pair_counter: collections.Counter[tuple[str, str]] = collections.Counter()
    all_barcodes_seen: set[str] = set()

    for rd in db_reads:
        bc_s = rd.get("bc_start_id")
        bc_e = rd.get("bc_end_id")
        if bc_s is not None and bc_e is not None:
            pair_counter[(bc_s, bc_e)] += 1
            all_barcodes_seen.add(bc_s)
            all_barcodes_seen.add(bc_e)

    result["pair_matrix"] = {
        "counts": pair_counter,
        "used_barcodes": sorted(all_barcodes_seen),
    }

    # ------------------------------------------------------------------
    # reads_table
    # ------------------------------------------------------------------
    reads_table: list[dict[str, Any]] = []

    for rd in db_reads:
        read_id_full: str = rd["read_id"]
        bc_start = rd.get("bc_start_id", "")
        bc_end = rd.get("bc_end_id", "")
        bc_end_conf = rd.get("bc_end_conf")
        is_fl = (
            bc_end_conf is not None and bc_end_conf >= _FULL_LENGTH_THRESHOLD
        )

        q_ld_raw = rd.get("q_ld")
        q_ld_rounded = round(q_ld_raw, 1) if q_ld_raw is not None else None

        bc_start_conf_raw = rd.get("bc_start_conf")
        bc_start_conf_rounded = (
            round(bc_start_conf_raw, 3) if bc_start_conf_raw is not None else None
        )
        bc_end_conf_rounded = (
            round(bc_end_conf, 3) if bc_end_conf is not None else None
        )

        truncated_id = (
            read_id_full[:_READ_ID_TRUNCATE] + "..."
            if len(read_id_full) > _READ_ID_TRUNCATE
            else read_id_full
        )

        reads_table.append({
            "read_id": truncated_id,
            "read_id_full": read_id_full,
            "length": rd["readlen"],
            "bc_start": bc_start,
            "bc_start_ed": rd.get("bc_start_ed"),
            "bc_start_conf": bc_start_conf_rounded,
            "bc_end": bc_end,
            "bc_end_ed": rd.get("bc_end_ed"),
            "bc_end_conf": bc_end_conf_rounded,
            "pair": f"{bc_start}--{bc_end}",
            "tgt_id": rd.get("tgt_id", ""),
            "is_full_length": is_fl,
            "ed": rd.get("ed"),
            "q_ld": q_ld_rounded,
        })

    result["reads_table"] = reads_table

    # ------------------------------------------------------------------
    # detailed_reads
    # ------------------------------------------------------------------
    result["detailed_reads"] = _build_detailed_reads(
        db_reads=db_reads,
        target_groups=target_groups,
        expected_barcodes=expected_barcodes,
        expected_barcodes_rc=expected_barcodes_rc,
        max_detail_reads=max_detail_reads,
        flank_front=flank_front,
        flank_rear=flank_rear,
        segment_len=segment_len,
    )

    # ------------------------------------------------------------------
    # barcode_info
    # ------------------------------------------------------------------
    barcode_info: dict[str, dict[str, str]] = {}
    for bc_id in sorted(used_bc_ids):
        if bc_id in BARCODES:
            fwd = BARCODES[bc_id]
            barcode_info[bc_id] = {
                "fwd": fwd,
                "rc": reverse_complement(fwd),
            }

    result["barcode_info"] = barcode_info

    # ------------------------------------------------------------------
    # pairing_table
    # ------------------------------------------------------------------
    pairing_table: list[dict[str, str]] = []
    for (upstream, downstream), alias in sorted(barcode_pair_to_alias.items()):
        pairing_table.append({
            "upstream": upstream,
            "downstream": downstream,
            "alias": alias,
        })

    result["pairing_table"] = pairing_table

    # ------------------------------------------------------------------
    # references (optional)
    # ------------------------------------------------------------------
    if references is not None:
        refs_out: dict[str, dict[str, Any]] = {}
        for tgt_id, (seq, length) in references.items():
            refs_out[tgt_id] = {"seq": seq, "length": length}
        result["references"] = refs_out

    # ------------------------------------------------------------------
    # length_stats (new)
    # ------------------------------------------------------------------
    all_lengths: list[int] = [rd["readlen"] for rd in db_reads]
    lengths_full: list[int] = []
    lengths_trunc: list[int] = []

    for rd in db_reads:
        bc_end_conf = rd.get("bc_end_conf")
        is_fl = bc_end_conf is not None and bc_end_conf >= _FULL_LENGTH_THRESHOLD
        if is_fl:
            lengths_full.append(rd["readlen"])
        else:
            lengths_trunc.append(rd["readlen"])

    result["length_stats"] = {
        "min": min(all_lengths) if all_lengths else 0,
        "max": max(all_lengths) if all_lengths else 0,
        "mean": sum(all_lengths) / len(all_lengths) if all_lengths else 0.0,
        "median": _compute_median(all_lengths),
        "n50": _compute_n50(all_lengths),
        "lengths": all_lengths,
        "lengths_full": lengths_full,
        "lengths_trunc": lengths_trunc,
    }

    # ------------------------------------------------------------------
    # quality_stats (new)
    # ------------------------------------------------------------------
    q_bc_values: list[float] = []
    q_ld_values: list[float] = []
    ed_values: list[int] = []

    for rd in db_reads:
        q_bc = rd.get("q_bc")
        if q_bc is not None:
            q_bc_values.append(q_bc)

        tgt_id = rd.get("tgt_id", "")
        is_matched = not tgt_id.startswith("unmatched_")
        if is_matched:
            q_ld = rd.get("q_ld")
            if q_ld is not None:
                q_ld_values.append(q_ld)
            ed = rd.get("ed")
            if ed is not None:
                ed_values.append(ed)

    result["quality_stats"] = {
        "q_bc_values": q_bc_values,
        "q_ld_values": q_ld_values,
        "ed_values": ed_values,
        "mean_q_bc": sum(q_bc_values) / len(q_bc_values) if q_bc_values else 0.0,
        "mean_q_ld": sum(q_ld_values) / len(q_ld_values) if q_ld_values else 0.0,
        "mean_ed": sum(ed_values) / len(ed_values) if ed_values else 0.0,
    }

    # ------------------------------------------------------------------
    # confidence_scatter (new)
    # ------------------------------------------------------------------
    confidence_scatter: list[dict[str, Any]] = []

    for rd in db_reads:
        start_conf = rd.get("bc_start_conf")
        end_conf = rd.get("bc_end_conf")
        if start_conf is not None and end_conf is not None:
            tgt_id = rd.get("tgt_id", "")
            is_fl = end_conf >= _FULL_LENGTH_THRESHOLD
            confidence_scatter.append({
                "x": start_conf,
                "y": end_conf,
                "tgt": tgt_id,
                "fl": is_fl,
                "len": rd["readlen"],
            })

    result["confidence_scatter"] = confidence_scatter

    # ------------------------------------------------------------------
    # end_reason_counts (new)
    # ------------------------------------------------------------------
    er_counter: collections.Counter[str] = collections.Counter()
    for rd in db_reads:
        er = rd.get("ER")
        if er is not None:
            er_counter[er] += 1

    result["end_reason_counts"] = dict(er_counter)

    # ------------------------------------------------------------------
    # classification_by_confidence (new)
    # ------------------------------------------------------------------
    classification_by_confidence: list[dict[str, Any]] = []

    for threshold in _CONFIDENCE_THRESHOLDS:
        count_at_threshold = 0
        for rd in db_reads:
            start_conf = rd.get("bc_start_conf")
            end_conf = rd.get("bc_end_conf")
            if (
                start_conf is not None
                and end_conf is not None
                and start_conf >= threshold
                and end_conf >= threshold
            ):
                count_at_threshold += 1

        pct = (count_at_threshold / total * 100.0) if total else 0.0
        classification_by_confidence.append({
            "threshold": threshold,
            "matched": count_at_threshold,
            "pct": round(pct, 1),
        })

    result["classification_by_confidence"] = classification_by_confidence

    return result


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _build_detailed_reads(
    db_reads: list[dict],
    target_groups: dict[str, list[dict]],
    expected_barcodes: dict[str, str],
    expected_barcodes_rc: dict[str, str],
    max_detail_reads: int,
    flank_front: str | None,
    flank_rear: str | None,
    segment_len: int,
) -> list[dict[str, Any]]:
    """Select representative reads and build detailed annotations.

    Selects up to 1 full-length and 1 truncated read per target group,
    limited to *max_detail_reads* total.

    Parameters
    ----------
    db_reads : list[dict]
        All reads from the database.
    target_groups : dict[str, list[dict]]
        Reads grouped by target ID (with "unmatched" as aggregate key).
    expected_barcodes : dict[str, str]
        Forward barcode sequences keyed by ID.
    expected_barcodes_rc : dict[str, str]
        Reverse-complement barcode sequences keyed by ID.
    max_detail_reads : int
        Maximum number of detailed reads to return.
    flank_front : str or None
        Front flanking sequence.
    flank_rear : str or None
        Rear flanking sequence.
    segment_len : int
        Length of segment window at read start/end.

    Returns
    -------
    list[dict]
        Detailed read annotations.
    """
    selected: list[dict] = []

    for tgt_id in sorted(target_groups.keys()):
        if len(selected) >= max_detail_reads:
            break

        reads = target_groups[tgt_id]
        fl_read: dict | None = None
        tr_read: dict | None = None

        for rd in reads:
            bc_end_conf = rd.get("bc_end_conf")
            is_fl = (
                bc_end_conf is not None
                and bc_end_conf >= _FULL_LENGTH_THRESHOLD
            )
            if is_fl and fl_read is None:
                fl_read = rd
            elif not is_fl and tr_read is None:
                tr_read = rd

            if fl_read is not None and tr_read is not None:
                break

        for rd in (fl_read, tr_read):
            if rd is not None and len(selected) < max_detail_reads:
                selected.append(rd)

    # Build detailed annotation for each selected read
    detailed: list[dict[str, Any]] = []

    for rd in selected:
        read_seq: str = rd["readseq"]
        read_len: int = rd["readlen"]
        bc_end_conf = rd.get("bc_end_conf")
        is_fl = (
            bc_end_conf is not None and bc_end_conf >= _FULL_LENGTH_THRESHOLD
        )

        bc_start = rd.get("bc_start_id", "")
        bc_end = rd.get("bc_end_id", "")

        # --- Competitor tables ---
        start_segment = read_seq[:segment_len]
        start_competitors = classify_barcode_detailed(
            start_segment, expected_barcodes
        )

        end_segment = read_seq[-segment_len:] if read_len >= segment_len else read_seq
        end_competitors = classify_barcode_detailed(
            end_segment, expected_barcodes_rc
        )

        # --- Region annotations ---
        regions: list[dict[str, Any]] = []

        # Start barcode position from competitor table
        start_bc_info = _find_bc_in_competitors(bc_start, start_competitors)
        if start_bc_info is not None:
            regions.append({
                "start": start_bc_info["start"],
                "end": start_bc_info["end"] + 1,
                "label": f"bc_start ({bc_start})",
                "type": "barcode_start",
                "ed": start_bc_info["edit_distance"],
            })

        # End barcode position from competitor table (offset by read position)
        end_offset = max(0, read_len - segment_len)
        end_bc_info = _find_bc_in_competitors(bc_end, end_competitors)
        if end_bc_info is not None:
            regions.append({
                "start": end_offset + end_bc_info["start"],
                "end": end_offset + end_bc_info["end"] + 1,
                "label": f"bc_end ({bc_end})",
                "type": "barcode_end",
                "ed": end_bc_info["edit_distance"],
            })

        # Flank annotations
        if flank_front is not None:
            # Front flank in first 50bp
            ff_result = find_flank_position(
                read_seq, flank_front, 0, _FLANK_SEARCH_WINDOW
            )
            if ff_result is not None:
                regions.append({
                    "start": ff_result["start"],
                    "end": ff_result["end"],
                    "label": "flank_front",
                    "type": "flank",
                    "ed": ff_result["edit_distance"],
                })

        if flank_rear is not None and start_bc_info is not None:
            # Rear flank after start barcode
            bc_end_pos = start_bc_info["end"] + 1
            fr_result = find_flank_position(
                read_seq, flank_rear,
                bc_end_pos,
                min(bc_end_pos + _FLANK_SEARCH_WINDOW, read_len),
            )
            if fr_result is not None:
                regions.append({
                    "start": fr_result["start"],
                    "end": fr_result["end"],
                    "label": "flank_rear",
                    "type": "flank",
                    "ed": fr_result["edit_distance"],
                })

            # RC(flank_rear) near end of read
            rc_flank_rear = reverse_complement(flank_rear)
            rc_fr_result = find_flank_position(
                read_seq, rc_flank_rear,
                max(0, read_len - segment_len),
                read_len,
            )
            if rc_fr_result is not None:
                regions.append({
                    "start": rc_fr_result["start"],
                    "end": rc_fr_result["end"],
                    "label": "rc_flank_rear",
                    "type": "flank",
                    "ed": rc_fr_result["edit_distance"],
                })

        if flank_front is not None:
            # RC(flank_front) near end of read
            rc_flank_front = reverse_complement(flank_front)
            rc_ff_result = find_flank_position(
                read_seq, rc_flank_front,
                max(0, read_len - segment_len),
                read_len,
            )
            if rc_ff_result is not None:
                regions.append({
                    "start": rc_ff_result["start"],
                    "end": rc_ff_result["end"],
                    "label": "rc_flank_front",
                    "type": "flank",
                    "ed": rc_ff_result["edit_distance"],
                })

        # Sort regions by start position
        regions.sort(key=lambda r: r["start"])

        detailed.append({
            "read_id": rd["read_id"],
            "length": read_len,
            "tgt_id": rd.get("tgt_id", ""),
            "is_full_length": is_fl,
            "pair": f"{bc_start}--{bc_end}",
            "seq": read_seq,
            "regions": regions,
            "start_competitors": start_competitors,
            "end_competitors": end_competitors,
            "db_ed": rd.get("ed"),
            "db_q_ld": rd.get("q_ld"),
        })

    return detailed


def _find_bc_in_competitors(
    bc_id: str,
    competitors: list[dict[str, Any]],
) -> dict[str, Any] | None:
    """Find a specific barcode in the competitors list.

    Parameters
    ----------
    bc_id : str
        Barcode ID to find.
    competitors : list[dict]
        Output from :func:`classify_barcode_detailed`.

    Returns
    -------
    dict or None
        The matching competitor dict, or None if not found.
    """
    for comp in competitors:
        if comp["bc_id"] == bc_id:
            return comp
    return None
