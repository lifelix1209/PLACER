#!/usr/bin/env python3
"""Analyze missed truth loci on the sTELLeR supplementary dataset.

This script explains likely miss reasons by combining:
1) truth set and final calls
2) BAM-level insertion evidence near each truth locus
3) fragment extraction outputs (ins_fragments.fasta / ins_fragment_hits.tsv)

Focus is on cross-bin evidence loss: PLACER bins reads by read start position,
while evidence windows are clipped to current bin boundaries in component build.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pysam


@dataclass
class TruthSite:
    chrom: str
    pos: int
    family: str


def load_truth(path: Path) -> List[TruthSite]:
    out: List[TruthSite] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            out.append(TruthSite(chrom=fields[0], pos=int(fields[1]), family=fields[2]))
    return out


def load_calls(path: Path) -> List[Tuple[str, int, str]]:
    out: List[Tuple[str, int, str]] = []
    if not path.exists():
        return out
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            out.append((fields[0], int(fields[1]), fields[2]))
    return out


def truth_hit_status(
    truth: List[TruthSite], calls: List[Tuple[str, int, str]], match_window: int
) -> Dict[Tuple[str, int], bool]:
    status: Dict[Tuple[str, int], bool] = {}
    for site in truth:
        hit = False
        for c_chrom, c_pos, _ in calls:
            if c_chrom == site.chrom and abs(c_pos - site.pos) <= match_window:
                hit = True
                break
        status[(site.chrom, site.pos)] = hit
    return status


def parse_fragment_fasta(path: Path) -> Dict[Tuple[str, int], List[str]]:
    by_anchor: Dict[Tuple[str, int], List[str]] = defaultdict(list)
    if not path.exists():
        return by_anchor
    header_re = re.compile(r"^>([^:]+):(\d+)\|")
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            line = line.strip()
            m = header_re.match(line)
            if not m:
                continue
            chrom = m.group(1)
            pos = int(m.group(2))
            by_anchor[(chrom, pos)].append(line[1:])
    return by_anchor


def parse_fragment_hits(path: Path) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    if not path.exists():
        return out
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            frag_id = row.get("fragment_id")
            if not frag_id:
                continue
            out[frag_id] = row
    return out


def pick_best_hit(
    site: TruthSite,
    fragments_by_anchor: Dict[Tuple[str, int], List[str]],
    hit_rows: Dict[str, Dict[str, str]],
    match_window: int,
) -> Tuple[int, str, float]:
    best_te = "NA"
    best_support = 0.0
    n_frag = 0
    for (chrom, pos), frag_ids in fragments_by_anchor.items():
        if chrom != site.chrom or abs(pos - site.pos) > match_window:
            continue
        for frag_id in frag_ids:
            n_frag += 1
            row = hit_rows.get(frag_id)
            if not row:
                continue
            support = float(row.get("kmer_support", "0") or "0")
            if support > best_support:
                best_support = support
                best_te = row.get("te", "NA") or "NA"
    return n_frag, best_te, best_support


def collect_bam_support(
    bam: pysam.AlignmentFile,
    site: TruthSite,
    event_window: int,
    fetch_window: int,
    bin_size: int,
) -> Dict[str, int]:
    reads_with_ins50_start: Dict[str, int] = {}
    total_ins50_reads: Set[str] = set()
    same_bin_reads: Set[str] = set()
    prev_bin_reads: Set[str] = set()
    next_bin_reads: Set[str] = set()

    locus_bin = (site.pos - 1) // bin_size
    start = max(0, site.pos - fetch_window)
    end = site.pos + fetch_window
    for aln in bam.fetch(site.chrom, start, end):
        if aln.is_unmapped or aln.is_secondary:
            continue
        if aln.cigartuples is None:
            continue
        ref = aln.reference_start + 1
        has_ins50 = False
        for op, length in aln.cigartuples:
            if op == 1 and length >= 50 and abs(ref - site.pos) <= event_window:
                has_ins50 = True
                break
            if op in (0, 2, 3, 7, 8):
                ref += length
        if not has_ins50:
            continue

        qname = aln.query_name
        total_ins50_reads.add(qname)
        read_start = aln.reference_start + 1
        if qname not in reads_with_ins50_start:
            reads_with_ins50_start[qname] = read_start

    for qname, read_start in reads_with_ins50_start.items():
        read_bin = (read_start - 1) // bin_size
        if read_bin == locus_bin:
            same_bin_reads.add(qname)
        elif read_bin < locus_bin:
            prev_bin_reads.add(qname)
        else:
            next_bin_reads.add(qname)

    return {
        "ins50_total": len(total_ins50_reads),
        "ins50_same_bin_start": len(same_bin_reads),
        "ins50_prev_bin_start": len(prev_bin_reads),
        "ins50_next_bin_start": len(next_bin_reads),
    }


def infer_reason(site_hit: bool, stats: Dict[str, int], n_frag: int) -> str:
    if site_hit:
        return "CALLED"
    if stats["ins50_total"] == 0:
        return "NO_STRONG_INSERTION_SIGNAL_IN_BAM"
    if stats["ins50_same_bin_start"] == 0:
        return "ALL_SUPPORT_READ_STARTS_OUTSIDE_LOCUS_BIN (cross-bin signal clipping)"
    if n_frag <= 1:
        return "TOO_FEW_FRAGMENTS_FOR_ASSEMBLY (requires >=2)"
    return "FILTERED_DOWNSTREAM (check te/evidence gates)"


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    default_data_root = repo_root.parent / "sTELLeR_supplementary"
    default_out_root = repo_root / "benchmark_steller_out"

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        default=str(default_data_root / "testdata.bam"),
    )
    parser.add_argument(
        "--truth",
        default=str(default_data_root / "testbamTRUTH.txt"),
    )
    parser.add_argument(
        "--calls",
        default=str(default_out_root / "calls.best_run1.tsv"),
    )
    parser.add_argument(
        "--fragments-fasta",
        default=str(repo_root / "ins_fragments.fasta"),
    )
    parser.add_argument(
        "--fragment-hits",
        default=str(repo_root / "ins_fragment_hits.tsv"),
    )
    parser.add_argument(
        "--out",
        default=str(default_out_root / "miss_reason.tsv"),
    )
    parser.add_argument("--match-window", type=int, default=200)
    parser.add_argument("--event-window", type=int, default=120)
    parser.add_argument("--fetch-window", type=int, default=25000)
    parser.add_argument("--bin-size", type=int, default=10000)
    args = parser.parse_args()

    truth = load_truth(Path(args.truth))
    calls = load_calls(Path(args.calls))
    hit_status = truth_hit_status(truth, calls, args.match_window)
    fragments_by_anchor = parse_fragment_fasta(Path(args.fragments_fasta))
    fragment_hits = parse_fragment_hits(Path(args.fragment_hits))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with pysam.AlignmentFile(args.bam, "rb") as bam, out_path.open("w", encoding="utf-8") as out:
        header = [
            "chrom",
            "pos",
            "family",
            "status",
            "ins50_total",
            "ins50_same_bin_start",
            "ins50_prev_bin_start",
            "ins50_next_bin_start",
            "n_fragments_near_site",
            "best_fragment_te",
            "best_fragment_kmer_support",
            "inferred_reason",
        ]
        out.write("\t".join(header) + "\n")

        for site in truth:
            key = (site.chrom, site.pos)
            status = "hit" if hit_status.get(key, False) else "miss"
            stats = collect_bam_support(
                bam=bam,
                site=site,
                event_window=args.event_window,
                fetch_window=args.fetch_window,
                bin_size=args.bin_size,
            )
            n_frag, best_te, best_support = pick_best_hit(
                site=site,
                fragments_by_anchor=fragments_by_anchor,
                hit_rows=fragment_hits,
                match_window=args.match_window,
            )
            reason = infer_reason(hit_status.get(key, False), stats, n_frag)
            out.write(
                "\t".join(
                    [
                        site.chrom,
                        str(site.pos),
                        site.family,
                        status,
                        str(stats["ins50_total"]),
                        str(stats["ins50_same_bin_start"]),
                        str(stats["ins50_prev_bin_start"]),
                        str(stats["ins50_next_bin_start"]),
                        str(n_frag),
                        best_te,
                        f"{best_support:.4f}",
                        reason,
                    ]
                )
                + "\n"
            )

    print(f"[analyze-steller] wrote: {out_path}")


if __name__ == "__main__":
    main()
