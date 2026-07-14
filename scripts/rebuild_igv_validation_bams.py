#!/usr/bin/env python3
"""Randomly sample truth loci and build local BAMs for IGV inspection.

This is the inspection-only companion to random_truth_interval_eval.py. It uses
the same TLDR truth parsing, deterministic sampling, window calculation, and
read-complete subset BAM extraction, but it does not run PLACER.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from random_truth_interval_eval import (  # noqa: E402
    PreparedSubset,
    TruthRow,
    WindowSpec,
    build_window,
    fail,
    load_contig_lengths,
    load_sampled_truth_rows,
    load_truth_rows,
    prepare_sample_subset_bams,
    require_file,
    resolve_samtools_threads,
    sample_truth_rows,
    write_tsv,
)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Randomly sample TLDR truth loci and create read-complete local "
            "subset BAMs for IGV review."
        )
    )
    parser.add_argument("--ground-truth", required=True, type=Path, help="TLDR table path.")
    parser.add_argument("--bam", required=True, type=Path, help="Input BAM path.")
    parser.add_argument(
        "--ref",
        required=True,
        type=Path,
        help="Reference FASTA path. Its .fai is used to clip windows.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help=(
            "Output directory. Default: "
            "placer_out/igv_validation_bams_<timestamp>_seed<seed>."
        ),
    )
    parser.add_argument(
        "--sampled-truth-tsv",
        type=Path,
        default=None,
        help="Replay an existing sampled_truth.tsv instead of drawing a new sample.",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=10,
        help="Number of truth rows to sample per run.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed. Default: generated automatically and recorded in summary.json.",
    )
    parser.add_argument(
        "--truth-filter",
        default="PASS",
        help="Only sample rows whose Filter column equals this value.",
    )
    parser.add_argument(
        "--min-length-ins",
        type=int,
        default=None,
        help="Only sample truth rows whose LengthIns is >= this threshold.",
    )
    parser.add_argument(
        "--window-flank-bp",
        type=int,
        default=2000,
        help="Flank size added to each side of the truth interval.",
    )
    parser.add_argument(
        "--samtools",
        default="samtools",
        help="samtools executable path.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="samtools threads used for view/index. Default: all detected CPUs.",
    )
    return parser.parse_args(argv)


def resolve_outdir(outdir: Optional[Path], seed: int) -> Path:
    if outdir is not None:
        resolved = outdir.expanduser().resolve()
    else:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        resolved = Path("placer_out") / f"igv_validation_bams_{stamp}_seed{seed}"
        resolved = resolved.resolve()
    if resolved.exists():
        fail(f"output directory already exists: {resolved}")
    return resolved


def build_sample_manifest(
    sampled_rows: Sequence[TruthRow],
    contig_lengths: Dict[str, int],
    flank_bp: int,
) -> tuple[List[Dict[str, str]], Dict[str, WindowSpec]]:
    rows: List[Dict[str, str]] = []
    windows_by_sample_id: Dict[str, WindowSpec] = {}
    for row in sampled_rows:
        window = build_window(row, flank_bp, contig_lengths)
        windows_by_sample_id[row.sample_id] = window
        rows.append(
            {
                "sample_id": row.sample_id,
                "uuid": row.uuid,
                "chrom": row.chrom,
                "truth_start": str(row.start),
                "truth_end": str(row.end),
                "truth_strand": row.strand,
                "truth_family": row.family,
                "truth_subfamily": row.subfamily,
                "truth_filter": row.filter_value,
                "truth_length_ins": str(row.length_ins),
                "window_start": str(window.window_start),
                "window_end": str(window.window_end),
                "window_region": window.region,
            }
        )
    return rows, windows_by_sample_id


def build_bam_manifest(
    sampled_rows: Sequence[TruthRow],
    windows_by_sample_id: Dict[str, WindowSpec],
    prepared_subsets: Dict[str, PreparedSubset],
) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for row in sampled_rows:
        subset = prepared_subsets[row.sample_id]
        window = windows_by_sample_id[row.sample_id]
        rows.append(
            {
                "sample_id": row.sample_id,
                "uuid": row.uuid,
                "chrom": row.chrom,
                "truth_start": str(row.start),
                "truth_end": str(row.end),
                "window_region": window.region,
                "qname_count": str(subset.qname_count),
                "subset_bam": str(subset.subset_bam),
                "subset_bai": str(Path(str(subset.subset_bam) + ".bai")),
                "window_qnames_txt": str(subset.read_names_path),
                "run_log": str(subset.log_path),
            }
        )
    return rows


def build_summary(
    *,
    seed: Optional[int],
    sampled_rows: Sequence[TruthRow],
    truth_pool_size: int,
    args: argparse.Namespace,
    outdir: Path,
    total_elapsed_s: float,
) -> Dict[str, object]:
    summary: Dict[str, object] = {
        "seed": seed,
        "sample_size": len(sampled_rows),
        "truth_filter": args.truth_filter,
        "min_length_ins": args.min_length_ins,
        "truth_pool_size": truth_pool_size,
        "window_flank_bp": args.window_flank_bp,
        "samtools_threads": args.threads,
        "total_elapsed_s": total_elapsed_s,
        "ground_truth": str(args.ground_truth),
        "bam": str(args.bam),
        "ref": str(args.ref),
        "outdir": str(outdir),
        "sampled_truth_tsv": str(outdir / "sampled_truth.tsv"),
        "igv_bams_tsv": str(outdir / "igv_bams.tsv"),
    }
    if args.sampled_truth_tsv is not None:
        summary["input_sampled_truth_tsv"] = str(args.sampled_truth_tsv)
    return summary


def main(argv: Sequence[str]) -> int:
    started = time.time()
    args = parse_args(argv)
    args.threads = resolve_samtools_threads(args.threads)
    if args.sampled_truth_tsv is not None:
        args.sampled_truth_tsv = require_file(args.sampled_truth_tsv, "sampled truth TSV")
    args.ground_truth = require_file(args.ground_truth, "ground truth table")
    args.bam = require_file(args.bam, "BAM")
    args.ref = require_file(args.ref, "reference FASTA")
    if args.min_length_ins is not None and args.min_length_ins < 0:
        fail("--min-length-ins must be >= 0")
    if args.window_flank_bp < 0:
        fail("--window-flank-bp must be >= 0")
    if args.threads <= 0:
        fail("--threads must be > 0")

    if args.sampled_truth_tsv is None:
        seed = args.seed if args.seed is not None else random.SystemRandom().randrange(1, 2**31)
    else:
        seed = args.seed
    if seed is None and args.outdir is None:
        fail("--sampled-truth-tsv without --seed requires an explicit --outdir")

    outdir = resolve_outdir(args.outdir, seed if seed is not None else 0)
    outdir.mkdir(parents=True, exist_ok=False)

    contig_lengths = load_contig_lengths(args.ref)
    truth_pool = load_truth_rows(
        args.ground_truth,
        args.truth_filter,
        min_length_ins=args.min_length_ins,
    )
    if args.sampled_truth_tsv is not None:
        sampled_rows = load_sampled_truth_rows(args.sampled_truth_tsv)
    else:
        sampled_rows = sample_truth_rows(truth_pool, args.sample_size, seed)

    sampled_manifest, windows_by_sample_id = build_sample_manifest(
        sampled_rows,
        contig_lengths,
        args.window_flank_bp,
    )
    prepared_subsets = prepare_sample_subset_bams(
        sampled_rows=sampled_rows,
        windows_by_sample_id=windows_by_sample_id,
        args=args,
        outdir=outdir,
    )
    bam_manifest = build_bam_manifest(
        sampled_rows,
        windows_by_sample_id,
        prepared_subsets,
    )

    write_tsv(outdir / "sampled_truth.tsv", sampled_manifest)
    write_tsv(outdir / "igv_bams.tsv", bam_manifest)
    summary = build_summary(
        seed=seed,
        sampled_rows=sampled_rows,
        truth_pool_size=len(truth_pool),
        args=args,
        outdir=outdir,
        total_elapsed_s=time.time() - started,
    )
    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"sampled_truth_tsv\t{outdir / 'sampled_truth.tsv'}")
    print(f"igv_bams_tsv\t{outdir / 'igv_bams.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
