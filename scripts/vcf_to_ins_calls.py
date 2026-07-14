#!/usr/bin/env python3
"""Normalize any SV-caller VCF into an insertion-call table for the caller benchmark.

Emits `chrom  pos  svlen  family` (family blank for non-TE-aware callers) for every
SVTYPE=INS record above --min-svlen. Works for Sniffles2 / cuteSV / SVIM style VCFs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pysam


def info_get(rec, key, default=None):
    try:
        val = rec.info.get(key, default)
    except (KeyError, ValueError):
        return default
    if isinstance(val, tuple):
        return val[0] if val else default
    return val


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-svlen", type=int, default=50)
    ap.add_argument("--min-qual", type=float, default=None,
                    help="drop records with QUAL below this (SVIM recommends 10)")
    ap.add_argument("--pass-only", action="store_true", help="keep only FILTER=PASS/. records")
    args = ap.parse_args()

    kept = 0
    with pysam.VariantFile(str(args.vcf)) as vcf, args.out.open("w") as out:
        out.write("chrom\tpos\tsvlen\tfamily\n")
        for rec in vcf.fetch():
            if info_get(rec, "SVTYPE") != "INS":
                continue
            if args.pass_only:
                keys = list(rec.filter.keys())
                if keys and keys != ["PASS"]:
                    continue
            if args.min_qual is not None and rec.qual is not None and rec.qual < args.min_qual:
                continue
            try:
                svlen = abs(int(info_get(rec, "SVLEN", 0) or 0))
            except (TypeError, ValueError):
                svlen = 0
            if svlen and svlen < args.min_svlen:
                continue
            out.write(f"{rec.chrom}\t{int(rec.pos)}\t{svlen}\t\n")
            kept += 1
    print(f"{args.vcf.name}: wrote {kept} INS calls -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
