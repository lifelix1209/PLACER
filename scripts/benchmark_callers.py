#!/usr/bin/env python3
"""Score multiple TE / SV callers against the simulation truth and emit one results table.

Each caller is supplied as `name=calls.tsv` where calls.tsv has a header with at
least `chrom` and `pos` (and optionally `family`). All callers are scored the same
way against the planted truth from simulate_te_insertions.py: a call is a true
positive when its breakpoint is within --radius bp of an unclaimed truth insertion.

Outputs (to --out-tsv) per caller: emitted, tp, fp, fn, precision, recall, f1,
family committed/correct/wrong counts, family concordance, family commit rate,
and median breakpoint error. PLACER's explicit family_status distinguishes an
abstention from a legitimate library family named Unknown.
This is the machine-readable input to the benchmark figure
(plot_caller_benchmark.py).
"""

from __future__ import annotations

import argparse
import bisect
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def read_truth(path: Path):
    by_chrom: dict[str, list[int]] = defaultdict(list)
    fam: dict[tuple[str, int], str] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chrom, pos = row["chrom"], int(row["pos"])
            by_chrom[chrom].append(pos)
            fam[(chrom, pos)] = row.get("family", "")
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom, fam


def read_calls(path: Path):
    calls = []
    with path.open() as handle:
        fields = None
        for line in handle:
            candidate = line.rstrip("\n")
            if candidate.startswith("#"):
                candidate = candidate[1:]
            candidate_fields = candidate.split("\t")
            if "chrom" in candidate_fields and "pos" in candidate_fields:
                fields = candidate_fields
                break
        if fields is None:
            return calls
        for row in csv.DictReader(handle, delimiter="\t", fieldnames=fields):
            try:
                status = (row.get("family_status") or "").strip().upper()
                committed = status == "COMMITTED" if status else None
                calls.append((
                    row["chrom"],
                    int(row["pos"]),
                    (row.get("family") or "").strip(),
                    committed,
                ))
            except (KeyError, ValueError):
                continue
    return calls


def norm_family(value: str) -> str:
    return (value or "").strip().lower().split("/")[0].split(":")[0].split("-")[0]


def is_committed_family(value: str, explicit_status: bool | None = None) -> bool:
    if explicit_status is not None:
        return explicit_status
    # Legacy files lack family_status. Preserve case so the PLACER abstention
    # sentinel UNKNOWN is distinct from a real library family named Unknown.
    return (value or "").strip() not in {"", "UNKNOWN", "NA", "NONE", "."}


def nearest_unused(sorted_pos, used, pos, radius):
    if not sorted_pos:
        return None
    i = bisect.bisect_left(sorted_pos, pos)
    best, best_d = None, radius + 1
    for j in range(max(0, i - 2), min(len(sorted_pos), i + 2)):
        p = sorted_pos[j]
        d = abs(p - pos)
        if p not in used and d <= radius and d < best_d:
            best, best_d = p, d
    return best


def score(calls, truth_pos, truth_fam, radius):
    total_truth = sum(len(v) for v in truth_pos.values())
    used = defaultdict(set)
    tp = fam_known = fam_agree = 0
    bp_errors = []
    # order calls nearest-first is unnecessary; greedy by input order is fine at this density
    for call in calls:
        chrom, pos, fam = call[:3]
        explicit_status = call[3] if len(call) > 3 else None
        hit = nearest_unused(truth_pos.get(chrom, []), used[chrom], pos, radius)
        if hit is None:
            continue
        used[chrom].add(hit)
        tp += 1
        bp_errors.append(abs(hit - pos))
        tfam = truth_fam.get((chrom, hit), "")
        if is_committed_family(fam, explicit_status) and tfam:
            fam_known += 1
            if norm_family(fam) == norm_family(tfam):
                fam_agree += 1
    emitted = len(calls)
    precision = tp / emitted if emitted else 0.0
    recall = tp / total_truth if total_truth else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return {
        "emitted": emitted,
        "tp": tp,
        "fp": emitted - tp,
        "fn": total_truth - tp,
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1": round(f1, 4),
        "family_committed": fam_known,
        "family_correct": fam_agree,
        "family_wrong": fam_known - fam_agree,
        "family_concordance": round(fam_agree / fam_known, 4) if fam_known else "",
        "family_commit_rate": round(fam_known / tp, 4) if tp else "",
        "median_bp_error": round(statistics.median(bp_errors), 1) if bp_errors else "",
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth", required=True, type=Path)
    ap.add_argument("--radius", type=int, default=100)
    ap.add_argument("--out-tsv", required=True, type=Path)
    ap.add_argument("caller", nargs="+", help="name=calls.tsv (order preserved in output)")
    args = ap.parse_args()

    truth_pos, truth_fam = read_truth(args.truth)
    total_truth = sum(len(v) for v in truth_pos.values())

    fields = ["caller", "emitted", "tp", "fp", "fn", "precision", "recall", "f1",
              "family_committed", "family_correct", "family_wrong",
              "family_concordance", "family_commit_rate", "median_bp_error"]
    results = []
    for spec in args.caller:
        if "=" not in spec:
            print(f"ERROR: bad caller spec (need name=path): {spec}", file=sys.stderr)
            return 1
        name, path = spec.split("=", 1)
        p = Path(path)
        if not p.exists():
            print(f"WARN: skip {name}: not found {p}", file=sys.stderr)
            continue
        row = {"caller": name, **score(read_calls(p), truth_pos, truth_fam, args.radius)}
        results.append(row)

    with args.out_tsv.open("w") as out:
        w = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    width = max((len(r["caller"]) for r in results), default=10)
    print(f"truth insertions: {total_truth}   match radius: {args.radius} bp\n")
    print(f"{'caller':<{width}}  {'emit':>5} {'TP':>4} {'FP':>4} {'FN':>4} "
          f"{'prec':>6} {'recall':>6} {'F1':>6} {'famcon':>6} {'commit':>6} {'bperr':>6}")
    for r in results:
        fam = f"{r['family_concordance']:.3f}" if r["family_concordance"] != "" else "  -  "
        commit = f"{r['family_commit_rate']:.3f}" if r["family_commit_rate"] != "" else "  -  "
        bpe = f"{r['median_bp_error']:.0f}" if r["median_bp_error"] != "" else "  -  "
        print(f"{r['caller']:<{width}}  {r['emitted']:>5} {r['tp']:>4} {r['fp']:>4} {r['fn']:>4} "
              f"{r['precision']:>6.3f} {r['recall']:>6.3f} {r['f1']:>6.3f} "
              f"{fam:>6} {commit:>6} {bpe:>6}")
    print(f"\nwrote -> {args.out_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
