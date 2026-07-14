#!/usr/bin/env python3
"""Step 3/4: normalize a GraffiTE VCF into a TE-insertion truth table.

Keeps non-reference TE insertions (SVTYPE=INS) that RepeatMasker annotated as a
transposable element, and emits a flat TSV consumed by benchmark_placer_vs_truth.py:

    chrom  start  end  svlen  family  te_class  n_hits

GraffiTE INFO fields used: SVTYPE, SVLEN, n_hits, repeat_ids, matching_classes,
mam_filter_1 (the >=80%-coverage single-family filter; "None" means no family
explains the SV well enough to be a confident pME).
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def parse_info(field: str) -> dict[str, str]:
    info: dict[str, str] = {}
    for item in field.split(";"):
        if not item:
            continue
        key, _, value = item.partition("=")
        info[key] = value
    return info


def primary_family(repeat_ids: str) -> str:
    """First repeat id; GraffiTE lists comma-separated hits, best first."""
    if not repeat_ids:
        return ""
    return repeat_ids.split(",")[0].strip()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True, type=Path, help="GraffiTE pangenome VCF (.vcf or .vcf.gz)")
    ap.add_argument("--out", required=True, type=Path, help="output truth TSV")
    ap.add_argument("--min-svlen", type=int, default=50, help="minimum |SVLEN| to keep (default 50)")
    ap.add_argument("--require-mam1", action="store_true", default=True,
                    help="require mam_filter_1 != None (a single TE family covers >=80%; default on)")
    ap.add_argument("--no-require-mam1", dest="require_mam1", action="store_false")
    args = ap.parse_args()

    if not args.vcf.exists():
        print(f"ERROR: VCF not found: {args.vcf}", file=sys.stderr)
        return 1

    kept = 0
    seen = 0
    with _open(args.vcf) as handle, args.out.open("w") as out:
        out.write("chrom\tstart\tend\tsvlen\tfamily\tte_class\tn_hits\n")
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            seen += 1
            chrom, pos = cols[0], cols[1]
            info = parse_info(cols[7])

            if info.get("SVTYPE") != "INS":
                continue
            try:
                svlen = abs(int(info.get("SVLEN", "0")))
                n_hits = int(info.get("n_hits", "0"))
                start = int(pos)
            except ValueError:
                continue
            if n_hits < 1 or svlen < args.min_svlen:
                continue
            if args.require_mam1 and info.get("mam_filter_1", "None") in ("", "None"):
                continue

            family = primary_family(info.get("repeat_ids", ""))
            te_class = primary_family(info.get("matching_classes", ""))
            if not family:
                continue
            out.write(f"{chrom}\t{start}\t{start + 1}\t{svlen}\t{family}\t{te_class}\t{n_hits}\n")
            kept += 1

    print(f"parsed {seen} VCF records -> kept {kept} TE-insertion truth rows -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
