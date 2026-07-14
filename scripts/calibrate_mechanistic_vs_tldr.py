#!/usr/bin/env python3
"""Validate the mechanistic (TPRT-hallmark) gate against TLDR PASS.

Computes the mechanistic score per candidate from:
  * endonuclease target motif (reference flanks around the breakpoint),
  * TSD (ledger tsd columns),
  * poly(A) tail (from the inserted sequence),
  * TE body (minimap2 besthit identity/coverage + core length),
then sweeps the hallmark-count / log-LR threshold and reports precision / recall
against TLDR PASS (a reference callset, not biological ground truth).
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from placer_py.evidence.segmental_explanation import build_segmental_explanation_from_alignment
from placer_py.model.mechanistic import endonuclease_motif_score, mechanistic_te_score
from scripts.calibrate_target_q_vs_tldr import (  # noqa: E402
    _f,
    _i,
    match_precision_recall,
    parse_besthits,
    parse_ins_fa,
    read_tldr_pass,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ledger", required=True, type=Path)
    ap.add_argument("--ins-fa", required=True, type=Path)
    ap.add_argument("--besthits", required=True, type=Path)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--tldr", required=True, type=Path)
    ap.add_argument("--radius", type=int, default=100)
    ap.add_argument("--flank", type=int, default=30)
    args = ap.parse_args()

    ins = parse_ins_fa(args.ins_fa)
    besthits = parse_besthits(args.besthits)
    tldr = read_tldr_pass(args.tldr)
    ref = pysam.FastaFile(str(args.reference))
    ref_contigs = set(ref.references)
    rows = list(csv.DictReader(args.ledger.open(), delimiter="\t"))

    # Per-candidate mechanistic score.
    scored = []  # (chrom, pos, hallmark_count, loglr, tldr_matchable)
    for row in rows:
        cid = row["candidate_id"]
        chrom, pos, seq = ins.get(cid, ("", 0, ""))
        if not chrom:
            continue
        bh = besthits.get(cid)
        insert_len = _i(row["sequence_len"]) or (len(seq) if seq else 0)

        # Poly(A) + TE core from the alignment-driven segmental explanation.
        if bh is not None and seq:
            seg = build_segmental_explanation_from_alignment(
                seq, identity=bh["identity"], query_coverage=bh["query_cov"],
                aln_query_start=bh["qstart"], aln_query_end=bh["qend"], mapq=bh["mapq"],
            )
            te_core = seg.te_core_fraction
            polya = seg.polya_fraction
            identity = bh["identity"]
        else:
            te_core, polya, identity = 0.0, 0.0, 0.0
        core_len = te_core * insert_len

        # Endonuclease motif from reference flanks.
        en_score, en_eval = 0.0, False
        if chrom in ref_contigs:
            length = ref.get_reference_length(chrom)
            left = ref.fetch(chrom, max(0, pos - args.flank), max(0, pos)).upper()
            right = ref.fetch(chrom, min(length, pos), min(length, pos + args.flank)).upper()
            en_score, en_eval = endonuclease_motif_score(left, right)

        tsd_len = _i(row.get("tsd_len"))
        tsd_sig = row.get("tsd_qc", "") == "PASS_TSD_DUP"
        tsd_bg_p = _f(row.get("tsd_bg_p"), 1.0)

        mech = mechanistic_te_score(
            endonuclease=en_score, endonuclease_evaluated=en_eval,
            tsd_len=tsd_len, tsd_significant=tsd_sig, tsd_bg_p=tsd_bg_p,
            polya_fraction=polya, insert_len=insert_len,
            te_core_fraction=te_core, identity=identity, core_len=core_len,
        )
        scored.append((chrom, pos, mech.hallmark_count, mech.loglr))

    total_pass = sum(len(v) for v in tldr.values())
    print(f"scored candidates={len(scored)}  TLDR PASS rows={total_pass}  radius={args.radius}bp\n")

    print("== hallmark-count gate ==")
    for hmin in (1, 2, 3, 4):
        pts = [(c, p) for (c, p, h, lr) in scored if h >= hmin]
        e, m, prec, rec, f1 = match_precision_recall(pts, tldr, args.radius)
        print(f"  hallmarks>={hmin}  emitted={e:6d} matched={m:5d} precision={prec:.3f} recall={rec:.3f} f1={f1:.3f}")

    print("\n== hallmarks>=2 AND log-LR gate ==")
    for thr in (2.0, 3.0, 4.0, 5.0, 6.0):
        pts = [(c, p) for (c, p, h, lr) in scored if h >= 2 and lr >= thr]
        e, m, prec, rec, f1 = match_precision_recall(pts, tldr, args.radius)
        print(f"  h>=2 & loglr>={thr}  emitted={e:6d} matched={m:5d} precision={prec:.3f} recall={rec:.3f} f1={f1:.3f}")


if __name__ == "__main__":
    main()
