#!/usr/bin/env python3
"""Step 3/3: score PLACER calls against the simulation truth from simulate_te_insertions.py.

Because the truth is planted (exact breakpoint, family, strand, length) this gives
real precision AND recall -- no assembly/caller bias. PLACER final calls are read
from scientific.txt; a call matches a truth insertion when their breakpoints are
within --radius bp. Each truth insertion can be matched at most once.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


def read_truth(path: Path):
    by_chrom: dict[str, list[int]] = defaultdict(list)
    meta: dict[tuple[str, int], tuple[str, int]] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chrom, pos = row["chrom"], int(row["pos"])
            by_chrom[chrom].append(pos)
            meta[(chrom, pos)] = (row["family"], int(row["te_len"]))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom, meta


def read_scientific(path: Path, qc_prefix: str):
    """Yield (chrom, pos, family, insert_len, qc) for final calls whose qc starts with qc_prefix."""
    header: list[str] | None = None
    calls = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#chrom\t"):
                header = line[1:].split("\t")
                continue
            if header is None or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < len(header):
                continue
            rec = dict(zip(header, fields))
            qc = rec.get("qc", "")
            if qc_prefix and not qc.startswith(qc_prefix):
                continue
            try:
                calls.append((rec["chrom"], int(rec["pos"]), rec.get("family", ""),
                              int(rec.get("insert_len", "0") or 0), qc))
            except ValueError:
                continue
    return calls


def nearest_unused(sorted_pos: list[int], used: set[int], pos: int, radius: int):
    if not sorted_pos:
        return None
    i = bisect.bisect_left(sorted_pos, pos)
    best = None
    best_d = radius + 1
    for j in range(max(0, i - 2), min(len(sorted_pos), i + 2)):
        p = sorted_pos[j]
        d = abs(p - pos)
        if p not in used and d <= radius and d < best_d:
            best, best_d = p, d
    return best


def norm_family(value: str) -> str:
    return (value or "").strip().lower().split("/")[0].split(":")[0].split("-")[0]


def is_abstained(family: str) -> bool:
    """A call that did not commit to a specific family (UNKNOWN / FAMILY_ONLY / blank)."""
    return norm_family(family) in ("", "unknown", "na")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scientific", required=True, type=Path, help="PLACER scientific.txt")
    ap.add_argument("--truth", required=True, type=Path, help="truth TSV from simulate_te_insertions.py")
    ap.add_argument("--radius", type=int, default=100, help="breakpoint match radius bp (default 100)")
    ap.add_argument("--qc-prefix", default="PASS_TE", help="final qc prefix to count as a TE call (default PASS_TE; use '' for all)")
    ap.add_argument("--out-matched", type=Path, help="optional TSV of matched/unmatched detail")
    ap.add_argument("--out-family", type=Path, help="optional TSV of the (truth_family, placer_family) confusion matrix")
    args = ap.parse_args()

    for f in (args.scientific, args.truth):
        if not f.exists():
            print(f"ERROR: not found: {f}", file=sys.stderr)
            return 1

    truth_pos, truth_meta = read_truth(args.truth)
    total_truth = sum(len(v) for v in truth_pos.values())
    calls = read_scientific(args.scientific, args.qc_prefix)

    used: dict[str, set[int]] = defaultdict(set)
    matched = 0
    # family outcome over matched (true-positive) calls
    fam_correct = fam_wrong = fam_abstain = 0
    confusion: Counter[tuple[str, str]] = Counter()  # (truth_family, placer_family) for wrong calls
    detail = []
    for chrom, pos, fam, ilen, qc in calls:
        hit = nearest_unused(truth_pos.get(chrom, []), used[chrom], pos, args.radius)
        if hit is None:
            detail.append(("FP", chrom, pos, fam, "", qc))
            continue
        used[chrom].add(hit)
        matched += 1
        tfam, _ = truth_meta[(chrom, hit)]
        if is_abstained(fam):
            fam_abstain += 1
        elif norm_family(fam) == norm_family(tfam):
            fam_correct += 1
        else:
            fam_wrong += 1
            confusion[(norm_family(tfam), norm_family(fam))] += 1
        detail.append(("TP", chrom, pos, fam, tfam, qc))

    emitted = len(calls)
    precision = matched / emitted if emitted else 0.0
    recall = matched / total_truth if total_truth else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    committed = fam_correct + fam_wrong
    # accuracy among calls that committed to a family (the number A/B optimisation targets)
    fam_accuracy = fam_correct / committed if committed else 0.0
    commit_rate = committed / matched if matched else 0.0
    wrong_rate = fam_wrong / matched if matched else 0.0

    print(f"qc filter:                 {args.qc_prefix or '(all final calls)'}")
    print(f"PLACER calls emitted:      {emitted}")
    print(f"planted truth insertions:  {total_truth}")
    print(f"match radius:              {args.radius} bp\n")
    print(f"true positives (matched):  {matched}")
    print(f"false positives:           {emitted - matched}")
    print(f"false negatives (missed):  {total_truth - matched}")
    print(f"precision:                 {precision:.3f}")
    print(f"recall:                    {recall:.3f}")
    print(f"f1:                        {f1:.3f}\n")
    print("family assignment (over matched calls):")
    print(f"  correct:                 {fam_correct}")
    print(f"  wrong (confident):       {fam_wrong}")
    print(f"  abstained (UNKNOWN/fam): {fam_abstain}")
    print(f"  committed family acc.:   {fam_accuracy:.3f}  (correct / committed, over {committed})")
    print(f"  commit rate:             {commit_rate:.3f}  (committed / matched)")
    print(f"  wrong-family rate:       {wrong_rate:.3f}  (wrong / matched)")
    if confusion:
        print("  top family confusions (truth -> placer : count):")
        for (tf, pf), n in confusion.most_common(10):
            print(f"    {tf or '(blank)'} -> {pf or '(blank)'} : {n}")

    if args.out_family:
        with args.out_family.open("w") as out:
            out.write("truth_family\tplacer_family\tcount\n")
            for (tf, pf), n in confusion.most_common():
                out.write(f"{tf}\t{pf}\t{n}\n")
        print(f"\nwrote family confusion -> {args.out_family}")

    if args.out_matched:
        with args.out_matched.open("w") as out:
            out.write("status\tchrom\tplacer_pos\tplacer_family\ttruth_family\tqc\n")
            for r in detail:
                out.write("\t".join(map(str, r)) + "\n")
            # record the missed truth insertions too
            for chrom, positions in truth_pos.items():
                for p in positions:
                    if p not in used[chrom]:
                        tfam, _ = truth_meta[(chrom, p)]
                        out.write(f"FN\t{chrom}\t{p}\t\t{tfam}\t\n")
        print(f"\nwrote detail -> {args.out_matched}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
