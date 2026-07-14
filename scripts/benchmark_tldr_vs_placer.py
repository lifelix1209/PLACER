#!/usr/bin/env python3
"""Benchmark PLACER final calls against a TLDR table.

The comparison is a coordinate concordance benchmark, not a truth-set claim:
TLDR can be used as a reference callset, but the resulting precision/recall
labels mean agreement with TLDR under the chosen matching radius.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)


@dataclass(frozen=True)
class PlacerCall:
    index: int
    chrom: str
    pos: int
    bp_left: int
    bp_right: int
    family: str
    subfamily: str
    qc: str
    posterior_qc: str
    lfdr_qc: str
    conformal_qc: str
    support_reads: int
    insert_len: int
    best_te_identity: float
    best_te_query_coverage: float
    raw: dict[str, str]


@dataclass(frozen=True)
class TldrCall:
    index: int
    uuid: str
    chrom: str
    start: int
    end: int
    strand: str
    family: str
    subfamily: str
    length_ins: int
    te_match: float
    used_reads: int
    span_reads: int
    remappable: str
    filt: str
    raw: dict[str, str]

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2.0


@dataclass(frozen=True)
class Match:
    distance: float
    placer_index: int
    tldr_index: int


def parse_int(value: str, default: int = 0) -> int:
    value = (value or "").strip()
    if not value or value.upper() == "NA":
        return default
    try:
        return int(float(value))
    except ValueError:
        return default


def parse_float(value: str, default: float = 0.0) -> float:
    value = (value or "").strip()
    if not value or value.upper() == "NA":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def norm_missing(value: str) -> str | None:
    value = (value or "").strip()
    if not value:
        return None
    if value.upper() in {"NA", "N/A", "NONE", "UNKNOWN", "."}:
        return None
    return value.upper()


def read_placer_scientific(path: Path) -> tuple[list[PlacerCall], dict[str, str]]:
    summary: dict[str, str] = {}
    calls: list[PlacerCall] = []
    header: list[str] | None = None

    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0] == "#chrom":
                header = row[:]
                header[0] = "chrom"
                break
            if row[0].startswith("#"):
                continue
            if len(row) >= 2:
                summary[row[0]] = row[1]

        if header is None:
            raise ValueError(f"Could not find #chrom table header in {path}")

        idx = {name: i for i, name in enumerate(header)}
        required = [
            "chrom",
            "pos",
            "bp_left",
            "bp_right",
            "family",
            "subfamily",
            "qc",
            "posterior_qc",
            "lfdr_qc",
            "conformal_qc",
            "support_reads",
            "insert_len",
            "best_te_identity",
            "best_te_query_coverage",
        ]
        missing = [name for name in required if name not in idx]
        if missing:
            raise ValueError(f"PLACER table missing columns: {', '.join(missing)}")

        for row in reader:
            if not row or len(row) < len(header):
                continue
            raw = {name: row[i] for name, i in idx.items()}
            calls.append(
                PlacerCall(
                    index=len(calls),
                    chrom=raw["chrom"],
                    pos=parse_int(raw["pos"]),
                    bp_left=parse_int(raw["bp_left"]),
                    bp_right=parse_int(raw["bp_right"]),
                    family=raw["family"],
                    subfamily=raw["subfamily"],
                    qc=raw["qc"],
                    posterior_qc=raw["posterior_qc"],
                    lfdr_qc=raw["lfdr_qc"],
                    conformal_qc=raw["conformal_qc"],
                    support_reads=parse_int(raw["support_reads"]),
                    insert_len=parse_int(raw["insert_len"]),
                    best_te_identity=parse_float(raw["best_te_identity"]),
                    best_te_query_coverage=parse_float(raw["best_te_query_coverage"]),
                    raw=raw,
                )
            )

    return calls, summary


def read_tldr_table(path: Path) -> list[TldrCall]:
    calls: list[TldrCall] = []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        idx = {name: i for i, name in enumerate(header)}
        required = [
            "UUID",
            "Chrom",
            "Start",
            "End",
            "Strand",
            "Family",
            "Subfamily",
            "LengthIns",
            "TEMatch",
            "UsedReads",
            "SpanReads",
            "Remappable",
            "Filter",
        ]
        missing = [name for name in required if name not in idx]
        if missing:
            raise ValueError(f"TLDR table missing columns: {', '.join(missing)}")

        for row in reader:
            if not row or len(row) <= max(idx.values()):
                continue
            raw = {name: row[i] for name, i in idx.items() if i < len(row)}
            calls.append(
                TldrCall(
                    index=len(calls),
                    uuid=raw["UUID"],
                    chrom=raw["Chrom"],
                    start=parse_int(raw["Start"]),
                    end=parse_int(raw["End"]),
                    strand=raw["Strand"],
                    family=raw["Family"],
                    subfamily=raw["Subfamily"],
                    length_ins=parse_int(raw["LengthIns"]),
                    te_match=parse_float(raw["TEMatch"]),
                    used_reads=parse_int(raw["UsedReads"]),
                    span_reads=parse_int(raw["SpanReads"]),
                    remappable=raw["Remappable"],
                    filt=raw["Filter"],
                    raw=raw,
                )
            )
    return calls


def index_tldr_by_chrom(calls: list[TldrCall]) -> dict[str, tuple[list[float], list[int]]]:
    grouped: dict[str, list[tuple[float, int]]] = defaultdict(list)
    for call in calls:
        grouped[call.chrom].append((call.midpoint, call.index))
    indexed: dict[str, tuple[list[float], list[int]]] = {}
    for chrom, values in grouped.items():
        values.sort()
        indexed[chrom] = ([point for point, _ in values], [index for _, index in values])
    return indexed


def greedy_point_match(
    placer_calls: list[PlacerCall],
    tldr_calls: list[TldrCall],
    radius: int,
) -> list[Match]:
    tldr_index = index_tldr_by_chrom(tldr_calls)
    candidates: list[tuple[float, int, int]] = []

    for placer in placer_calls:
        indexed = tldr_index.get(placer.chrom)
        if indexed is None:
            continue
        points, indices = indexed
        lo = bisect.bisect_left(points, placer.pos - radius)
        hi = bisect.bisect_right(points, placer.pos + radius)
        for offset in range(lo, hi):
            tldr_call = tldr_calls[indices[offset]]
            distance = abs(float(placer.pos) - tldr_call.midpoint)
            candidates.append((distance, placer.index, tldr_call.index))

    candidates.sort()
    used_placer: set[int] = set()
    used_tldr: set[int] = set()
    matches: list[Match] = []
    for distance, placer_index, tldr_index_value in candidates:
        if placer_index in used_placer or tldr_index_value in used_tldr:
            continue
        used_placer.add(placer_index)
        used_tldr.add(tldr_index_value)
        matches.append(Match(distance, placer_index, tldr_index_value))
    return matches


def percentile(values: list[float], pct: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    rank = (len(ordered) - 1) * pct
    low = math.floor(rank)
    high = math.ceil(rank)
    if low == high:
        return ordered[int(rank)]
    weight = rank - low
    return ordered[low] * (1 - weight) + ordered[high] * weight


def family_stats(
    matches: list[Match],
    placer_calls: list[PlacerCall],
    tldr_calls: list[TldrCall],
) -> dict[str, int | float]:
    family_known = 0
    family_match = 0
    subfamily_known = 0
    subfamily_match = 0
    for match in matches:
        placer = placer_calls[match.placer_index]
        tldr = tldr_calls[match.tldr_index]
        pfam = norm_missing(placer.family)
        tfam = norm_missing(tldr.family)
        psub = norm_missing(placer.subfamily)
        tsub = norm_missing(tldr.subfamily)
        if pfam is not None and tfam is not None:
            family_known += 1
            family_match += int(pfam == tfam)
        if psub is not None and tsub is not None:
            subfamily_known += 1
            subfamily_match += int(psub == tsub)
    return {
        "family_known_pairs": family_known,
        "family_exact_matches": family_match,
        "family_exact_rate": family_match / family_known if family_known else math.nan,
        "subfamily_known_pairs": subfamily_known,
        "subfamily_exact_matches": subfamily_match,
        "subfamily_exact_rate": subfamily_match / subfamily_known if subfamily_known else math.nan,
    }


def interval_coverage_stats(
    label: str,
    placer_calls: list[PlacerCall],
    tldr_calls: list[TldrCall],
    point_matches: list[Match] | None = None,
) -> tuple[dict[str, object], set[int]]:
    """Count TLDR calls whose midpoint is covered by any PLACER interval.

    This is a diagnostic rather than the main concordance metric because broad
    PLACER event intervals can cover multiple TLDR calls.
    """

    import heapq

    intervals_by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    tldr_by_chrom: dict[str, list[TldrCall]] = defaultdict(list)
    for placer in placer_calls:
        lo = min(placer.bp_left, placer.bp_right, placer.pos)
        hi = max(placer.bp_left, placer.bp_right, placer.pos)
        intervals_by_chrom[placer.chrom].append((lo, hi, placer.family))
    for call in tldr_calls:
        tldr_by_chrom[call.chrom].append(call)

    covered_ids: set[int] = set()
    multiplicity = Counter()
    family_known = 0
    family_exact = 0

    for chrom, calls in tldr_by_chrom.items():
        intervals = sorted(intervals_by_chrom.get(chrom, []))
        if not intervals:
            continue
        calls = sorted(calls, key=lambda call: call.midpoint)
        interval_cursor = 0
        active_heap: list[tuple[int, int, str | None]] = []
        active_count = 0
        active_known_family_count = 0
        active_family = Counter()
        serial = 0

        for call in calls:
            point = call.midpoint
            while interval_cursor < len(intervals) and intervals[interval_cursor][0] <= point:
                _, hi, family = intervals[interval_cursor]
                norm_family = norm_missing(family)
                heapq.heappush(active_heap, (hi, serial, norm_family))
                serial += 1
                active_count += 1
                if norm_family is not None:
                    active_known_family_count += 1
                    active_family[norm_family] += 1
                interval_cursor += 1

            while active_heap and active_heap[0][0] < point:
                _, _, norm_family = heapq.heappop(active_heap)
                active_count -= 1
                if norm_family is not None:
                    active_known_family_count -= 1
                    active_family[norm_family] -= 1
                    if active_family[norm_family] <= 0:
                        del active_family[norm_family]

            if active_count <= 0:
                continue
            covered_ids.add(call.index)
            multiplicity[active_count] += 1
            tldr_family = norm_missing(call.family)
            if tldr_family is not None and active_known_family_count > 0:
                family_known += 1
                family_exact += int(active_family.get(tldr_family, 0) > 0)

    point_matched_ids = {match.tldr_index for match in point_matches} if point_matches is not None else set()
    row: dict[str, object] = {
        "comparison": label,
        "tldr_calls": len(tldr_calls),
        "covered_by_any_placer_interval": len(covered_ids),
        "tldr_interval_coverage": len(covered_ids) / len(tldr_calls) if tldr_calls else math.nan,
        "not_interval_covered": len(tldr_calls) - len(covered_ids),
        "point_matched": len(point_matched_ids) if point_matches is not None else "NA",
        "interval_covered_not_point_matched": len(covered_ids - point_matched_ids)
        if point_matches is not None
        else "NA",
        "family_known_interval_covered": family_known,
        "family_exact_any_cover": family_exact,
        "family_exact_any_cover_rate": family_exact / family_known if family_known else math.nan,
        "covered_once": multiplicity.get(1, 0),
        "covered_multiple": len(covered_ids) - multiplicity.get(1, 0),
    }
    return row, covered_ids


def metrics_for_matches(
    label: str,
    radius: int,
    matches: list[Match],
    placer_total: int,
    tldr_total: int,
    placer_calls: list[PlacerCall],
    tldr_calls: list[TldrCall],
) -> dict[str, str | int | float]:
    distances = [match.distance for match in matches]
    matched = len(matches)
    placer_agreement = matched / placer_total if placer_total else math.nan
    tldr_recovery = matched / tldr_total if tldr_total else math.nan
    f1 = (
        2 * placer_agreement * tldr_recovery / (placer_agreement + tldr_recovery)
        if placer_agreement + tldr_recovery
        else math.nan
    )
    stats = family_stats(matches, placer_calls, tldr_calls)
    return {
        "comparison": label,
        "radius_bp": radius,
        "placer_calls": placer_total,
        "tldr_calls": tldr_total,
        "matched_pairs": matched,
        "placer_agreement": placer_agreement,
        "tldr_recovery": tldr_recovery,
        "f1_vs_tldr": f1,
        "median_abs_distance_bp": statistics.median(distances) if distances else math.nan,
        "p90_abs_distance_bp": percentile(distances, 0.90),
        "p99_abs_distance_bp": percentile(distances, 0.99),
        **stats,
    }


def format_value(value: object) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        if abs(value) >= 1000 or (0 < abs(value) < 0.001):
            return f"{value:.6g}"
        return f"{value:.4f}"
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_value(row.get(key, "")) for key in fields})


def match_pair_rows(
    matches: list[Match],
    placer_calls: list[PlacerCall],
    tldr_calls: list[TldrCall],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for match in sorted(matches, key=lambda item: (placer_calls[item.placer_index].chrom, placer_calls[item.placer_index].pos)):
        placer = placer_calls[match.placer_index]
        tldr = tldr_calls[match.tldr_index]
        rows.append(
            {
                "distance_bp": match.distance,
                "placer_chrom": placer.chrom,
                "placer_pos": placer.pos,
                "placer_bp_left": placer.bp_left,
                "placer_bp_right": placer.bp_right,
                "placer_family": placer.family,
                "placer_subfamily": placer.subfamily,
                "placer_qc": placer.qc,
                "placer_support_reads": placer.support_reads,
                "placer_insert_len": placer.insert_len,
                "tldr_uuid": tldr.uuid,
                "tldr_chrom": tldr.chrom,
                "tldr_start": tldr.start,
                "tldr_end": tldr.end,
                "tldr_family": tldr.family,
                "tldr_subfamily": tldr.subfamily,
                "tldr_filter": tldr.filt,
                "tldr_used_reads": tldr.used_reads,
                "tldr_span_reads": tldr.span_reads,
                "tldr_length_ins": tldr.length_ins,
            }
        )
    return rows


def unmatched_placer_rows(
    matches: list[Match],
    placer_calls: list[PlacerCall],
    limit: int = 200,
) -> list[dict[str, object]]:
    matched = {match.placer_index for match in matches}
    remaining = [call for call in placer_calls if call.index not in matched]
    remaining.sort(key=lambda call: (-call.support_reads, call.chrom, call.pos))
    rows: list[dict[str, object]] = []
    for call in remaining[:limit]:
        rows.append(
            {
                "chrom": call.chrom,
                "pos": call.pos,
                "bp_left": call.bp_left,
                "bp_right": call.bp_right,
                "family": call.family,
                "subfamily": call.subfamily,
                "qc": call.qc,
                "posterior_qc": call.posterior_qc,
                "lfdr_qc": call.lfdr_qc,
                "conformal_qc": call.conformal_qc,
                "support_reads": call.support_reads,
                "insert_len": call.insert_len,
                "best_te_identity": call.best_te_identity,
                "best_te_query_coverage": call.best_te_query_coverage,
            }
        )
    return rows


def unmatched_tldr_rows(
    matches: list[Match],
    tldr_calls: list[TldrCall],
    limit: int = 200,
) -> list[dict[str, object]]:
    matched = {match.tldr_index for match in matches}
    remaining = [call for call in tldr_calls if call.index not in matched]
    remaining.sort(key=lambda call: (-call.used_reads, call.chrom, call.start, call.end))
    rows: list[dict[str, object]] = []
    for call in remaining[:limit]:
        rows.append(
            {
                "uuid": call.uuid,
                "chrom": call.chrom,
                "start": call.start,
                "end": call.end,
                "strand": call.strand,
                "family": call.family,
                "subfamily": call.subfamily,
                "length_ins": call.length_ins,
                "te_match": call.te_match,
                "used_reads": call.used_reads,
                "span_reads": call.span_reads,
                "remappable": call.remappable,
                "filter": call.filt,
            }
        )
    return rows


def counter_table(counter: Counter[str], limit: int = 20) -> list[tuple[str, int]]:
    return counter.most_common(limit)


def markdown_table(rows: list[dict[str, object]], fields: list[str]) -> str:
    lines = ["|" + "|".join(fields) + "|", "|" + "|".join(["---"] * len(fields)) + "|"]
    for row in rows:
        lines.append("|" + "|".join(format_value(row.get(field, "")) for field in fields) + "|")
    return "\n".join(lines)


def write_report(
    path: Path,
    placer_path: Path,
    tldr_path: Path,
    summary: dict[str, str],
    placer_calls: list[PlacerCall],
    tldr_calls: list[TldrCall],
    tldr_pass: list[TldrCall],
    metric_rows: list[dict[str, object]],
    interval_rows: list[dict[str, object]],
    report_radius: int,
) -> None:
    selected = [
        row
        for row in metric_rows
        if row["radius_bp"] == report_radius and row["comparison"] in {"TLDR all", "TLDR PASS"}
    ]
    tldr_filters = Counter(call.filt for call in tldr_calls)
    placer_qc = Counter(call.qc for call in placer_calls)
    placer_lfdr = Counter(call.lfdr_qc for call in placer_calls)
    placer_conformal = Counter(call.conformal_qc for call in placer_calls)

    fields = [
        "comparison",
        "radius_bp",
        "placer_calls",
        "tldr_calls",
        "matched_pairs",
        "placer_agreement",
        "tldr_recovery",
        "f1_vs_tldr",
        "median_abs_distance_bp",
        "p90_abs_distance_bp",
        "family_exact_rate",
        "family_known_pairs",
    ]

    lines = [
        "# PLACER vs TLDR Benchmark",
        "",
        "## Inputs",
        f"- PLACER final call table: `{placer_path}`",
        f"- TLDR table: `{tldr_path}`",
        "",
        "## Method",
        "- PLACER calls are read from the final-call table embedded in `scientific.txt`.",
        "- TLDR is evaluated in two views: all rows and rows with `Filter == PASS`.",
        "- Loci are represented by PLACER `pos` and TLDR midpoint `(Start + End) / 2`.",
        "- Matching is same-chromosome, one-to-one greedy assignment by smallest absolute coordinate distance.",
        "- Radii are reported as a sweep, so the conclusion does not depend on a single hand-picked tolerance.",
        "- Agreement metrics treat TLDR as the comparison callset, not as biological ground truth.",
        "",
        "## Callset Sizes",
        f"- PLACER final calls: {len(placer_calls)}",
        f"- TLDR all rows: {len(tldr_calls)}",
        f"- TLDR PASS rows: {len(tldr_pass)}",
    ]
    if summary:
        for key in [
            "total_reads",
            "gate1_passed",
            "processed_bins",
            "components",
            "event_consensus_calls",
            "genotype_calls",
            "final_pass_calls",
            "performance_pipeline_wall_seconds",
        ]:
            if key in summary:
                lines.append(f"- PLACER `{key}`: {summary[key]}")

    lines.extend(
        [
            "",
            f"## Main Checkpoint ({report_radius} bp)",
            markdown_table(selected, fields),
            "",
            "## Radius Sweep",
            markdown_table(metric_rows, fields),
            "",
            "## Interval Coverage Diagnostic",
            "This diagnostic counts a TLDR call as covered when its midpoint falls inside any PLACER `[bp_left, bp_right]` interval. It is not one-to-one and can overestimate agreement for broad PLACER clusters, but it separates coordinate-offset effects from absence of any PLACER event.",
            markdown_table(
                interval_rows,
                [
                    "comparison",
                    "tldr_calls",
                    "covered_by_any_placer_interval",
                    "tldr_interval_coverage",
                    "not_interval_covered",
                    "point_matched",
                    "interval_covered_not_point_matched",
                    "family_exact_any_cover_rate",
                    "family_known_interval_covered",
                    "covered_once",
                    "covered_multiple",
                ],
            ),
            "",
            "## TLDR Filter Distribution",
            "|filter|count|",
            "|---|---|",
        ]
    )
    lines.extend(f"|{name}|{count}|" for name, count in counter_table(tldr_filters, 25))

    lines.extend(
        [
            "",
            "## PLACER QC Distribution",
            "|qc|count|",
            "|---|---|",
        ]
    )
    lines.extend(f"|{name}|{count}|" for name, count in counter_table(placer_qc, 25))

    lines.extend(
        [
            "",
            "## PLACER LFDR QC",
            "|lfdr_qc|count|",
            "|---|---|",
        ]
    )
    lines.extend(f"|{name}|{count}|" for name, count in counter_table(placer_lfdr, 25))

    lines.extend(
        [
            "",
            "## PLACER Conformal QC",
            "|conformal_qc|count|",
            "|---|---|",
        ]
    )
    lines.extend(f"|{name}|{count}|" for name, count in counter_table(placer_conformal, 25))

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def remap_tldr_indices(calls: list[TldrCall]) -> list[TldrCall]:
    remapped: list[TldrCall] = []
    for index, call in enumerate(calls):
        remapped.append(
            TldrCall(
                index=index,
                uuid=call.uuid,
                chrom=call.chrom,
                start=call.start,
                end=call.end,
                strand=call.strand,
                family=call.family,
                subfamily=call.subfamily,
                length_ins=call.length_ins,
                te_match=call.te_match,
                used_reads=call.used_reads,
                span_reads=call.span_reads,
                remappable=call.remappable,
                filt=call.filt,
                raw=call.raw,
            )
        )
    return remapped


def parse_radii(values: str) -> list[int]:
    radii = sorted({int(value.strip()) for value in values.split(",") if value.strip()})
    if not radii:
        raise argparse.ArgumentTypeError("at least one radius is required")
    return radii


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--placer", required=True, type=Path, help="PLACER scientific.txt path")
    parser.add_argument("--tldr", required=True, type=Path, help="TLDR .table.txt path")
    parser.add_argument("--out-dir", required=True, type=Path, help="Directory for benchmark outputs")
    parser.add_argument(
        "--radii",
        type=parse_radii,
        default=parse_radii("0,50,100,250,500,1000,5000"),
        help="Comma-separated coordinate radii in bp",
    )
    parser.add_argument(
        "--report-radius",
        type=int,
        default=500,
        help="Radius used for detailed matched/unmatched output files",
    )
    parser.add_argument(
        "--unmatched-limit",
        type=int,
        default=200,
        help="Maximum unmatched rows to write for each side",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    placer_calls, summary = read_placer_scientific(args.placer)
    tldr_calls = read_tldr_table(args.tldr)
    tldr_pass = remap_tldr_indices([call for call in tldr_calls if call.filt == "PASS"])

    metric_rows: list[dict[str, object]] = []
    match_cache: dict[tuple[str, int], list[Match]] = {}
    for radius in args.radii:
        all_matches = greedy_point_match(placer_calls, tldr_calls, radius)
        pass_matches = greedy_point_match(placer_calls, tldr_pass, radius)
        match_cache[("all", radius)] = all_matches
        match_cache[("pass", radius)] = pass_matches
        metric_rows.append(
            metrics_for_matches(
                "TLDR all",
                radius,
                all_matches,
                len(placer_calls),
                len(tldr_calls),
                placer_calls,
                tldr_calls,
            )
        )
        metric_rows.append(
            metrics_for_matches(
                "TLDR PASS",
                radius,
                pass_matches,
                len(placer_calls),
                len(tldr_pass),
                placer_calls,
                tldr_pass,
            )
        )

    write_tsv(args.out_dir / "benchmark_metrics.tsv", metric_rows)

    if args.report_radius not in args.radii:
        all_matches = greedy_point_match(placer_calls, tldr_calls, args.report_radius)
        pass_matches = greedy_point_match(placer_calls, tldr_pass, args.report_radius)
    else:
        all_matches = match_cache[("all", args.report_radius)]
        pass_matches = match_cache[("pass", args.report_radius)]

    write_tsv(
        args.out_dir / f"matched_pairs_tldr_all_{args.report_radius}bp.tsv",
        match_pair_rows(all_matches, placer_calls, tldr_calls),
    )
    write_tsv(
        args.out_dir / f"matched_pairs_tldr_pass_{args.report_radius}bp.tsv",
        match_pair_rows(pass_matches, placer_calls, tldr_pass),
    )
    write_tsv(
        args.out_dir / f"unmatched_placer_vs_tldr_pass_{args.report_radius}bp_top.tsv",
        unmatched_placer_rows(pass_matches, placer_calls, args.unmatched_limit),
    )
    write_tsv(
        args.out_dir / f"unmatched_tldr_pass_{args.report_radius}bp_top.tsv",
        unmatched_tldr_rows(pass_matches, tldr_pass, args.unmatched_limit),
    )
    interval_rows = [
        interval_coverage_stats(f"TLDR all; point radius {args.report_radius} bp", placer_calls, tldr_calls, all_matches)[0],
        interval_coverage_stats(f"TLDR PASS; point radius {args.report_radius} bp", placer_calls, tldr_pass, pass_matches)[0],
    ]
    write_tsv(args.out_dir / "interval_coverage.tsv", interval_rows)
    write_report(
        args.out_dir / "benchmark_report.md",
        args.placer,
        args.tldr,
        summary,
        placer_calls,
        tldr_calls,
        tldr_pass,
        metric_rows,
        interval_rows,
        args.report_radius,
    )

    print(f"PLACER final calls: {len(placer_calls)}")
    print(f"TLDR all rows: {len(tldr_calls)}")
    print(f"TLDR PASS rows: {len(tldr_pass)}")
    print(f"Wrote: {args.out_dir / 'benchmark_report.md'}")
    print(f"Wrote: {args.out_dir / 'benchmark_metrics.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
