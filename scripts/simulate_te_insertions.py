#!/usr/bin/env python3
"""Plant known TE insertions into a reference sequence to build a simulation truth set.

This is the *controlled* truth for benchmarking PLACER: unlike an assembly- or
caller-derived truth, every planted insertion has an exact breakpoint, family,
strand and length, so precision AND recall are exact (no assembly fragmentation
to understate recall).

Pipeline position (step 1/3):
    simulate_te_insertions.py  -> modified.fa + truth.tsv     (this script)
    badread + minimap2         -> reads aligned to the ORIGINAL reference
    benchmark_sim_truth.py     -> precision / recall vs truth.tsv

The modified FASTA is emitted in ORIGINAL reference coordinates: we splice TE
sequence into a copy of the chosen contig(s) and keep the contig name, so reads
simulated from it and aligned back to the untouched reference land at exactly the
recorded breakpoint. Truth positions are therefore reference coordinates.

No third-party deps (stdlib only) so it runs unchanged on a login node or a
compute node without a conda env.
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    """Load a FASTA into {short_name: sequence}. Short name = header up to first whitespace."""
    seqs: dict[str, list[str]] = {}
    name = None
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def te_family(header: str) -> str:
    """`family:copy` -> `family`; fall back to the whole header."""
    return header.split(":", 1)[0] if ":" in header else header


def valid_te_copies(te: dict[str, str], min_len: int, max_len: int) -> list[str]:
    return [name for name, seq in te.items() if min_len <= len(seq) <= max_len and "N" not in seq[:50]]


def n_run_mask(seq: str, window: int) -> list[bool]:
    """True where a +/-window neighbourhood is free of N (safe to place an insertion)."""
    # prefix count of N so any window's N-count is O(1)
    prefix = [0] * (len(seq) + 1)
    for i, base in enumerate(seq):
        prefix[i + 1] = prefix[i] + (1 if base == "N" else 0)
    ok = [False] * len(seq)
    for i in range(window, len(seq) - window):
        if prefix[i + window] - prefix[i - window] == 0:
            ok[i] = True
    return ok


def choose_sites(seq: str, count: int, flank: int, min_gap: int, rng: random.Random) -> list[int]:
    """Pick `count` insertion breakpoints with clean flanks and >= min_gap spacing."""
    ok = n_run_mask(seq, flank)
    candidates = [i for i in range(flank, len(seq) - flank) if ok[i]]
    if not candidates:
        raise SystemExit("ERROR: no N-free sites with the requested flank; lower --flank or pick another contig")
    rng.shuffle(candidates)
    chosen: list[int] = []
    for pos in candidates:
        if all(abs(pos - c) >= min_gap for c in chosen):
            chosen.append(pos)
            if len(chosen) == count:
                break
    chosen.sort()
    return chosen


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, type=Path, help="reference FASTA (untouched target of alignment)")
    ap.add_argument("--te", required=True, type=Path, help="TE library FASTA (family:copy headers)")
    ap.add_argument("--out-fa", required=True, type=Path, help="output modified FASTA to simulate reads from")
    ap.add_argument("--out-truth", required=True, type=Path, help="output truth TSV")
    ap.add_argument("--contigs", default="", help="comma-separated contigs to modify (default: all in --ref)")
    ap.add_argument("--num", type=int, default=120, help="insertions per contig (default 120)")
    ap.add_argument("--min-te-len", type=int, default=250, help="min planted TE length bp (default 250)")
    ap.add_argument("--max-te-len", type=int, default=8000, help="max planted TE length bp (default 8000)")
    ap.add_argument("--flank", type=int, default=2000, help="min N-free flank around each site (default 2000)")
    ap.add_argument("--min-gap", type=int, default=50000, help="min spacing between insertions bp (default 50000)")
    ap.add_argument("--seed", type=int, default=1, help="RNG seed (default 1)")
    args = ap.parse_args()

    for f in (args.ref, args.te):
        if not f.exists():
            print(f"ERROR: not found: {f}", file=sys.stderr)
            return 1

    rng = random.Random(args.seed)
    ref = read_fasta(args.ref)
    te = read_fasta(args.te)
    copies = valid_te_copies(te, args.min_te_len, args.max_te_len)
    if not copies:
        print("ERROR: no TE copies in the requested length range", file=sys.stderr)
        return 1

    contigs = [c.strip() for c in args.contigs.split(",") if c.strip()] or list(ref)
    missing = [c for c in contigs if c not in ref]
    if missing:
        print(f"ERROR: contigs not in reference: {missing}", file=sys.stderr)
        return 1

    truth_rows: list[tuple] = []
    with args.out_fa.open("w") as fa:
        for contig in contigs:
            seq = ref[contig]
            sites = choose_sites(seq, args.num, args.flank, args.min_gap, rng)
            pieces: list[str] = []
            prev = 0
            for pos in sites:
                copy = rng.choice(copies)
                strand = rng.choice(["+", "-"])
                ins = te[copy] if strand == "+" else revcomp(te[copy])
                pieces.append(seq[prev:pos])
                pieces.append(ins)
                prev = pos
                truth_rows.append((contig, pos, pos + 1, te_family(copy), copy, strand, len(ins)))
            pieces.append(seq[prev:])
            fa.write(f">{contig}\n")
            modified = "".join(pieces)
            for i in range(0, len(modified), 60):
                fa.write(modified[i:i + 60] + "\n")

    with args.out_truth.open("w") as out:
        out.write("chrom\tpos0\tpos\tfamily\tte_copy\tstrand\tte_len\n")
        for row in truth_rows:
            out.write("\t".join(map(str, row)) + "\n")

    print(f"planted {len(truth_rows)} insertions across {len(contigs)} contig(s)")
    print(f"  modified FASTA -> {args.out_fa}")
    print(f"  truth TSV      -> {args.out_truth}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
