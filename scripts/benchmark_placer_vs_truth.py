#!/usr/bin/env python3
"""Step 4/4: benchmark PLACER TE calls against the GraffiTE assembly truth set.

Unlike the TLDR comparison, the GraffiTE truth is caller-independent (derived from
an assembly), so precision/recall here approximate real error rates -- bounded only
by assembly completeness (short-read-N50 fragmentation understates recall).

PLACER TE calls are taken from the evidence ledger (emit_te_call == 1 by default)
and positioned by their Sniffles candidate breakpoint from the PASS_INS FASTA.
Each call is matched to the nearest truth breakpoint within --radius bp.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.calibrate_target_q_vs_tldr import parse_ins_fa, short_id  # noqa: E402


def read_truth(path: Path) -> tuple[dict[str, list[int]], dict[tuple[str, int], str]]:
    """chrom -> sorted breakpoints, and (chrom,pos) -> family."""
    by_chrom: dict[str, list[int]] = defaultdict(list)
    family: dict[tuple[str, int], str] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chrom, start = row["chrom"], int(row["start"])
            by_chrom[chrom].append(start)
            family[(chrom, start)] = row.get("family", "")
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom, family


def read_placer_te_calls(ledger: Path, ins_fa: Path, select: str) -> list[tuple[str, int, str]]:
    """Return (chrom, pos, family) for emitted TE calls, positioned via the ins FASTA."""
    ins = parse_ins_fa(ins_fa)
    calls: list[tuple[str, int, str]] = []
    with ledger.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get(select, "0") not in ("1", "True", "true", "yes"):
                continue
            cid = short_id(row["candidate_id"])
            chrom, pos, _ = ins.get(cid, ("", 0, ""))
            if not chrom:
                continue
            calls.append((chrom, pos, row.get("family", "")))
    return calls


def nearest_within(sorted_pos: list[int], pos: int, radius: int) -> int | None:
    if not sorted_pos:
        return None
    i = bisect.bisect_left(sorted_pos, pos)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(sorted_pos):
            d = abs(sorted_pos[j] - pos)
            if d <= radius and (best is None or d < abs(best - pos)):
                best = sorted_pos[j]
    return best


def norm_family(value: str) -> str:
    return (value or "").strip().lower().split("/")[0].split("-")[0]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ledger", required=True, type=Path, help="evidence_ledger_py.tsv")
    ap.add_argument("--ins-fa", required=True, type=Path, help="D2.sniffles.PASS_INS.fa")
    ap.add_argument("--truth", required=True, type=Path, help="truth TSV from graffite_vcf_to_truth.py")
    ap.add_argument("--radius", type=int, default=100, help="breakpoint match radius bp (default 100)")
    ap.add_argument("--select", default="emit_te_call", help="ledger column that is 1 for an emitted TE call")
    ap.add_argument("--out-matched", type=Path, help="optional TSV of matched pairs")
    args = ap.parse_args()

    for f in (args.ledger, args.ins_fa, args.truth):
        if not f.exists():
            print(f"ERROR: not found: {f}", file=sys.stderr)
            return 1

    truth_pos, truth_fam = read_truth(args.truth)
    total_truth = sum(len(v) for v in truth_pos.values())
    calls = read_placer_te_calls(args.ledger, args.ins_fa, args.select)

    matched = 0
    matched_truth: set[tuple[str, int]] = set()
    fam_known = 0
    fam_agree = 0
    rows_out: list[tuple] = []
    for chrom, pos, fam in calls:
        hit = nearest_within(truth_pos.get(chrom, []), pos, args.radius)
        if hit is None:
            continue
        matched += 1
        matched_truth.add((chrom, hit))
        tfam = truth_fam.get((chrom, hit), "")
        if fam and tfam:
            fam_known += 1
            if norm_family(fam) == norm_family(tfam):
                fam_agree += 1
        rows_out.append((chrom, pos, hit, fam, tfam))

    emitted = len(calls)
    precision = matched / emitted if emitted else 0.0
    recall = len(matched_truth) / total_truth if total_truth else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    fam_conc = fam_agree / fam_known if fam_known else 0.0

    print(f"PLACER TE calls emitted (select={args.select}): {emitted}")
    print(f"GraffiTE truth insertions:                     {total_truth}")
    print(f"match radius:                                  {args.radius} bp\n")
    print(f"matched calls:        {matched}")
    print(f"precision:            {precision:.3f}   (matched / emitted)")
    print(f"recall:               {recall:.3f}   (distinct truth hit / total truth)")
    print(f"f1:                   {f1:.3f}")
    print(f"family concordance:   {fam_conc:.3f}   (over {fam_known} matched pairs with both families)")

    if args.out_matched:
        with args.out_matched.open("w") as out:
            out.write("chrom\tplacer_pos\ttruth_pos\tplacer_family\ttruth_family\n")
            for r in rows_out:
                out.write("\t".join(map(str, r)) + "\n")
        print(f"\nwrote matched pairs -> {args.out_matched}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
