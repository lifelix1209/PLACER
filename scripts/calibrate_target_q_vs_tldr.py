#!/usr/bin/env python3
"""Calibrate the precision-first target_q against a TLDR PASS reference.

For the D2 sample we have:
  * the Python evidence ledger (BAM + reference/boundary evidence per candidate),
  * the real inserted sequences (PASS_INS.fa, keyed by candidate id + chrom:pos),
  * the real minimap2 TE alignments (besthits.tsv: identity, coverage, span),
  * a TLDR table (D1.table.txt) whose PASS rows are the reference truth.

This sweeps target_q and reports coordinate concordance (precision = fraction of
emitted TE calls matching a TLDR PASS row; recall = fraction of TLDR PASS rows
matched) for three sequence-evidence modes:
  * alignment  -- segmental explanation driven by the minimap2 alignment,
  * kmer       -- segmental explanation from exact 17-mers,
  * fallback   -- TeHit coverage fallback (no raw sequence).

TLDR is a reference callset, not biological ground truth; metrics mean agreement
with TLDR under the matching radius.
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

from placer_py.evidence.reference_boundary import ReferenceEvidence
from placer_py.evidence.segmental_explanation import (
    SegmentalExplanation,
    build_segmental_explanation,
    build_segmental_explanation_from_alignment,
)
from placer_py.evidence.sequence import build_te_kmer_index, read_fasta_records
from placer_py.model.placer_discriminator import build_discriminator_evidence, evaluate_joint_hypotheses
from placer_py.models import BamEvidence, InsertionCandidate, SequenceFeatures, TeHit

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)


def _f(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def _i(x, d=0):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return d


def short_id(full: str) -> str:
    return full.split("|", 1)[0]


def parse_ins_fa(path: Path) -> dict[str, tuple[str, int, str]]:
    """short_id -> (chrom, pos, sequence)."""
    out: dict[str, tuple[str, int, str]] = {}
    name = None
    chrom = ""
    pos = 0
    parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out[name] = (chrom, pos, "".join(parts))
                header = line[1:]
                name = short_id(header)
                chrom, pos = "", 0
                for field in header.split("|"):
                    if ":" in field and field.split(":", 1)[0].startswith("chr"):
                        c, p = field.split(":", 1)
                        chrom, pos = c, _i(p)
                parts = []
            elif name is not None:
                parts.append(line.strip().upper())
    if name is not None:
        out[name] = (chrom, pos, "".join(parts))
    return out


def parse_besthits(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out[short_id(row["query"])] = {
                "identity": _f(row.get("identity")),
                "query_cov": _f(row.get("query_cov")),
                "qstart": _i(row.get("qstart")),
                "qend": _i(row.get("qend")),
                "mapq": _i(row.get("mapq")),
                "te_hit": row.get("te_hit", ""),
            }
    return out


def read_tldr_pass(path: Path) -> dict[str, list[float]]:
    """chrom -> sorted list of PASS midpoints."""
    by_chrom: dict[str, list[float]] = defaultdict(list)
    with path.open(errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("Filter") != "PASS":
                continue
            by_chrom[row["Chrom"]].append((_i(row["Start"]) + _i(row["End"])) / 2.0)
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


def match_precision_recall(points: list[tuple[str, int]], tldr: dict[str, list[float]], radius: int):
    """Greedy one-to-one match of emitted (chrom,pos) to TLDR PASS midpoints."""
    total_pass = sum(len(v) for v in tldr.values())
    used: dict[str, set[int]] = defaultdict(set)
    # Build (distance, point_i, chrom, tldr_i) sorted by distance for one-to-one.
    cand = []
    for pi, (chrom, pos) in enumerate(points):
        arr = tldr.get(chrom)
        if not arr:
            continue
        lo = bisect.bisect_left(arr, pos - radius)
        hi = bisect.bisect_right(arr, pos + radius)
        for ti in range(lo, hi):
            cand.append((abs(pos - arr[ti]), pi, chrom, ti))
    cand.sort()
    used_pt: set[int] = set()
    matched = 0
    for _, pi, chrom, ti in cand:
        if pi in used_pt or ti in used[chrom]:
            continue
        used_pt.add(pi)
        used[chrom].add(ti)
        matched += 1
    emitted = len(points)
    precision = matched / emitted if emitted else 0.0
    recall = matched / total_pass if total_pass else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return emitted, matched, precision, recall, f1


def build_evidence(row, seq, besthit, mode):
    sid = row["candidate_id"]
    length = _i(row["sequence_len"])
    cand = InsertionCandidate(
        "chr1", 101, 101, sid, length or 1, "PASS", _i(row["raw_cigar_insert_reads"]),
        alt_sequence=(seq or None),
    )
    bam = BamEvidence(
        raw_cigar_insert_reads=_i(row["raw_cigar_insert_reads"]),
        split_read_support=_i(row["split_read_support"]),
        left_clip_reads=_i(row["left_clip_reads"]),
        right_clip_reads=_i(row["right_clip_reads"]),
        ref_span_reads=_i(row["ref_span_reads"]),
        low_mapq_ref_span_reads=_i(row["low_mapq_ref_span_reads"]),
        local_depth=_i(row["local_depth"]),
    )
    feat = SequenceFeatures(
        length=length, gc_fraction=_f(row["sequence_gc"]),
        entropy=_f(row["sequence_entropy"]), low_complexity_fraction=_f(row["sequence_low_complexity"]),
    )
    hit = TeHit(
        family=row.get("family", ""), subfamily=row.get("subfamily", ""),
        identity=_f(row["te_identity"]), query_coverage=_f(row["te_query_coverage"]),
        orientation=row.get("te_orientation", ""),
    )
    ref_ok = row.get("reference_qc", "") not in ("", "REFERENCE_NOT_EVALUATED")
    ref_ev = ReferenceEvidence(
        reference_available=ref_ok,
        reference_qc=row.get("reference_qc", "REFERENCE_NOT_EVALUATED"),
        left_anchor_len=_i(row.get("left_ref_anchor_len")),
        right_anchor_len=_i(row.get("right_ref_anchor_len")),
        insert_core_sequence=("X" * _i(row.get("insert_core_len")) if _i(row.get("insert_core_len")) else None),
        remap_qc=row.get("remap_qc", "REFERENCE_REMAP_NOT_EVALUATED"),
        boundary_qc=row.get("boundary_qc", "BOUNDARY_NOT_EVALUATED"),
        boundary_score=1.0 if row.get("boundary_qc", "").startswith("PASS_BOUNDARY_") else -0.25,
    )
    segmental = None
    if mode == "alignment":
        if besthit is not None and seq:
            fam = besthit["te_hit"].split(":", 1)
            segmental = build_segmental_explanation_from_alignment(
                seq, identity=besthit["identity"], query_coverage=besthit["query_cov"],
                aln_query_start=besthit["qstart"], aln_query_end=besthit["qend"], mapq=besthit["mapq"],
                family=fam[0], subfamily=fam[1] if len(fam) > 1 else fam[0],
            )
        else:
            # No minimap2 TE alignment -> no TE sequence evidence (do not inherit
            # the degenerate exact-k-mer coverage via the TeHit fallback).
            segmental = SegmentalExplanation()
    elif mode == "kmer" and seq:
        segmental = build_segmental_explanation(seq, REFS, INDEX)
    # mode == "fallback" -> segmental None (derived from TeHit inside builder)
    return build_discriminator_evidence(cand, bam, feat, hit, reference_evidence=(ref_ev if ref_ok else None), segmental=segmental)


REFS: dict[str, str] = {}
INDEX: dict[str, set[str]] = {}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ledger", required=True, type=Path)
    ap.add_argument("--ins-fa", required=True, type=Path)
    ap.add_argument("--besthits", required=True, type=Path)
    ap.add_argument("--te-fasta", required=True, type=Path)
    ap.add_argument("--tldr", required=True, type=Path)
    ap.add_argument("--radius", type=int, default=100)
    ap.add_argument("--modes", default="alignment,kmer")
    ap.add_argument("--target-qs", default="0.02,0.05,0.10,0.20,0.30")
    args = ap.parse_args()

    global REFS, INDEX
    ins = parse_ins_fa(args.ins_fa)
    besthits = parse_besthits(args.besthits)
    tldr = read_tldr_pass(args.tldr)
    rows = list(csv.DictReader(args.ledger.open(), delimiter="\t"))
    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    if "kmer" in modes:
        REFS = read_fasta_records(args.te_fasta)
        INDEX = build_te_kmer_index(REFS)

    total_pass = sum(len(v) for v in tldr.values())
    print(f"ledger rows={len(rows)}  TLDR PASS rows={total_pass}  radius={args.radius}bp")

    # OLD baseline from recorded ledger emit_te_call.
    old_points = []
    for row in rows:
        if _i(row.get("emit_te_call")):
            chrom, pos, _ = ins.get(row["candidate_id"], ("", 0, ""))
            if chrom:
                old_points.append((chrom, pos))
    e, m, p, r, f1 = match_precision_recall(old_points, tldr, args.radius)
    print(f"\nOLD (recorded emit)            emitted={e:6d} matched={m:5d} precision={p:.3f} recall={r:.3f} f1={f1:.3f}")

    qs = [float(x) for x in args.target_qs.split(",")]
    for mode in modes:
        print(f"\n=== mode={mode} ===")
        # Precompute evidence once per row (independent of target_q); decision depends on q.
        cache = []
        for row in rows:
            chrom, pos, seq = ins.get(row["candidate_id"], ("", 0, ""))
            if not chrom:
                continue
            bh = besthits.get(row["candidate_id"])
            ev = build_evidence(row, seq, bh, mode)
            cache.append((chrom, pos, ev))
        for q in qs:
            pts = [(c, p) for (c, p, ev) in cache if evaluate_joint_hypotheses(ev, target_q=q).emit_te_call]
            e, m, prec, rec, f1 = match_precision_recall(pts, tldr, args.radius)
            print(f"  target_q={q:<5} emitted={e:6d} matched={m:5d} precision={prec:.3f} recall={rec:.3f} f1={f1:.3f}")


if __name__ == "__main__":
    main()
