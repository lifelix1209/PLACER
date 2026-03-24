#!/usr/bin/env python3
"""Randomly sample truth loci and evaluate whether PLACER rediscovers them.

Evaluation semantics:
- Truth pool defaults to rows with `Filter == PASS` from the TLDR table.
- Each sampled locus is evaluated independently on a read-complete subset BAM.
  The subset keeps every alignment record for reads that overlap the local
  truth window, so supplementary alignments outside the window are preserved.
- A PLACER call is counted as detected when it is on the same chromosome and
  the minimum distance between `[bp_left, bp_right]` and `[Start, End]` is
  less than or equal to `--match-distance-bp`.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


@dataclass(frozen=True)
class TruthRow:
    sample_id: str
    uuid: str
    chrom: str
    start: int
    end: int
    strand: str
    family: str
    subfamily: str
    filter_value: str
    raw: Dict[str, str]


@dataclass(frozen=True)
class WindowSpec:
    chrom: str
    truth_start: int
    truth_end: int
    window_start: int
    window_end: int

    @property
    def region(self) -> str:
        return f"{self.chrom}:{self.window_start}-{self.window_end}"


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Randomly sample TLDR truth loci and check whether PLACER rediscovers them."
    )
    parser.add_argument("--ground-truth", required=True, type=Path, help="TLDR table path.")
    parser.add_argument("--bam", required=True, type=Path, help="Input BAM path.")
    parser.add_argument("--ref", required=True, type=Path, help="Reference FASTA path.")
    parser.add_argument("--te", required=True, type=Path, help="TE FASTA path.")
    parser.add_argument(
        "--placer",
        type=Path,
        default=Path("build/placer"),
        help="PLACER executable path.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory. Default: placer_out/random_truth_eval_<timestamp>_seed<seed>",
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
        help="Random seed. Default: generated automatically and recorded in the report.",
    )
    parser.add_argument(
        "--truth-filter",
        default="PASS",
        help="Only sample rows whose Filter column equals this value.",
    )
    parser.add_argument(
        "--window-flank-bp",
        type=int,
        default=2000,
        help="Flank size added to each side of the truth interval when extracting the local BAM.",
    )
    parser.add_argument(
        "--match-distance-bp",
        type=int,
        default=1000,
        help="Max interval distance for counting a PLACER call as detected.",
    )
    parser.add_argument(
        "--samtools",
        default="samtools",
        help="samtools executable path.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="samtools threads used for view/index.",
    )
    return parser.parse_args(argv)


def fail(message: str) -> "NoReturn":
    raise SystemExit(message)


def require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        fail(f"{label} not found: {resolved}")
    return resolved


def require_executable(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        fail(f"{label} not found: {resolved}")
    if not os.access(resolved, os.X_OK):
        fail(f"{label} is not executable: {resolved}")
    return resolved


def resolve_outdir(outdir: Optional[Path], seed: int) -> Path:
    if outdir is not None:
        resolved = outdir.expanduser().resolve()
    else:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        resolved = Path("placer_out") / f"random_truth_eval_{stamp}_seed{seed}"
        resolved = resolved.resolve()
    if resolved.exists():
        fail(f"output directory already exists: {resolved}")
    return resolved


def load_contig_lengths(reference_fasta: Path) -> Dict[str, int]:
    fai_path = Path(str(reference_fasta) + ".fai")
    if not fai_path.is_file():
        fail(f"reference FASTA index not found: {fai_path}")

    contig_lengths: Dict[str, int] = {}
    with fai_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            chrom = fields[0]
            try:
                contig_lengths[chrom] = int(fields[1])
            except ValueError as exc:
                fail(f"invalid contig length in {fai_path}: {line.rstrip()}")
    if not contig_lengths:
        fail(f"no contig lengths loaded from {fai_path}")
    return contig_lengths


def load_truth_rows(path: Path, truth_filter: str) -> List[TruthRow]:
    rows: List[TruthRow] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"missing header in truth table: {path}")
        required = {"UUID", "Chrom", "Start", "End", "Strand", "Family", "Subfamily", "Filter"}
        missing = required.difference(reader.fieldnames)
        if missing:
            fail(f"truth table missing columns {sorted(missing)}: {path}")
        for row_index, row in enumerate(reader, start=1):
            if row["Filter"] != truth_filter:
                continue
            try:
                start = int(row["Start"])
                end = int(row["End"])
            except ValueError:
                fail(f"invalid Start/End at truth row {row_index + 1}")
            if end < start:
                start, end = end, start
            rows.append(
                TruthRow(
                    sample_id="",
                    uuid=row["UUID"],
                    chrom=row["Chrom"],
                    start=start,
                    end=end,
                    strand=row["Strand"],
                    family=row["Family"],
                    subfamily=row["Subfamily"],
                    filter_value=row["Filter"],
                    raw=dict(row),
                )
            )
    if not rows:
        fail(f"no truth rows matched Filter == {truth_filter!r} in {path}")
    return rows


def sample_truth_rows(rows: Sequence[TruthRow], sample_size: int, seed: int) -> List[TruthRow]:
    if sample_size <= 0:
        fail("--sample-size must be > 0")
    if sample_size > len(rows):
        fail(
            f"--sample-size={sample_size} exceeds truth pool size {len(rows)} "
            f"for the selected filter"
        )
    rng = random.Random(seed)
    sampled = rng.sample(list(rows), sample_size)
    labelled: List[TruthRow] = []
    for index, row in enumerate(sampled, start=1):
        labelled.append(
            TruthRow(
                sample_id=f"S{index:02d}",
                uuid=row.uuid,
                chrom=row.chrom,
                start=row.start,
                end=row.end,
                strand=row.strand,
                family=row.family,
                subfamily=row.subfamily,
                filter_value=row.filter_value,
                raw=row.raw,
            )
        )
    return labelled


def build_window(row: TruthRow, flank_bp: int, contig_lengths: Dict[str, int]) -> WindowSpec:
    if flank_bp < 0:
        fail("--window-flank-bp must be >= 0")
    contig_len = contig_lengths.get(row.chrom)
    if contig_len is None:
        fail(f"truth contig not present in reference index: {row.chrom}")
    window_start = max(1, row.start - flank_bp)
    window_end = min(contig_len, row.end + flank_bp)
    if window_end < window_start:
        fail(
            f"invalid extraction window for {row.sample_id} {row.chrom}:{row.start}-{row.end}"
        )
    return WindowSpec(
        chrom=row.chrom,
        truth_start=row.start,
        truth_end=row.end,
        window_start=window_start,
        window_end=window_end,
    )


def run_command(
    cmd: Sequence[str],
    *,
    cwd: Optional[Path] = None,
    log_path: Optional[Path] = None,
) -> None:
    completed = subprocess.run(
        list(cmd),
        cwd=str(cwd) if cwd else None,
        check=False,
        text=True,
        capture_output=True,
    )
    if log_path is not None:
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write("$ " + " ".join(cmd) + "\n")
            if completed.stdout:
                handle.write(completed.stdout)
                if not completed.stdout.endswith("\n"):
                    handle.write("\n")
            if completed.stderr:
                handle.write(completed.stderr)
                if not completed.stderr.endswith("\n"):
                    handle.write("\n")
    if completed.returncode != 0:
        fail(
            f"command failed with exit code {completed.returncode}: "
            f"{' '.join(cmd)}"
        )


def extract_subset_bam(
    *,
    samtools: str,
    threads: int,
    source_bam: Path,
    region: str,
    read_names_path: Path,
    subset_bam: Path,
    log_path: Path,
) -> int:
    region_fetch = subprocess.run(
        [
            samtools,
            "view",
            "-@",
            str(threads),
            str(source_bam),
            region,
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("$ " + " ".join([samtools, "view", "-@", str(threads), str(source_bam), region]) + "\n")
        if region_fetch.stdout:
            handle.write(f"[window_records_stdout_lines]\t{len(region_fetch.stdout.splitlines())}\n")
        if region_fetch.stderr:
            handle.write(region_fetch.stderr)
            if not region_fetch.stderr.endswith("\n"):
                handle.write("\n")
    if region_fetch.returncode != 0:
        fail(
            f"command failed with exit code {region_fetch.returncode}: "
            f"{samtools} view -@ {threads} {source_bam} {region}"
        )

    qnames: List[str] = []
    seen: set[str] = set()
    for line in region_fetch.stdout.splitlines():
        if not line:
            continue
        fields = line.split("\t", 1)
        if not fields:
            continue
        qname = fields[0]
        if qname not in seen:
            seen.add(qname)
            qnames.append(qname)
    if not qnames:
        fail(f"no overlapping read names found for region {region}")
    read_names_path.write_text("".join(name + "\n" for name in qnames), encoding="utf-8")

    run_command(
        [
            samtools,
            "view",
            "-@",
            str(threads),
            "-N",
            str(read_names_path),
            "-b",
            "-o",
            str(subset_bam),
            str(source_bam),
        ],
        log_path=log_path,
    )
    run_command(
        [samtools, "index", "-@", str(threads), str(subset_bam)],
        log_path=log_path,
    )
    return len(qnames)


def parse_scientific_txt(path: Path) -> List[Dict[str, str]]:
    if not path.is_file():
        return []

    header: Optional[List[str]] = None
    rows: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#chrom\t"):
                header = line[1:].split("\t")
                continue
            if line.startswith("#"):
                continue
            if header is None:
                continue
            fields = line.split("\t")
            if len(fields) != len(header):
                fail(f"malformed scientific row in {path}: {line}")
            rows.append(dict(zip(header, fields)))
    if header is None:
        fail(f"scientific header not found in {path}")
    return rows


def interval_distance(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    if a_end < a_start:
        a_start, a_end = a_end, a_start
    if b_end < b_start:
        b_start, b_end = b_end, b_start
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0


def choose_best_call(
    truth: TruthRow,
    calls: Iterable[Dict[str, str]],
) -> Optional[Dict[str, str]]:
    best_row: Optional[Dict[str, str]] = None
    best_distance: Optional[int] = None
    for call in calls:
        if call.get("chrom") != truth.chrom:
            continue
        try:
            bp_left = int(call["bp_left"])
            bp_right = int(call["bp_right"])
        except (KeyError, ValueError):
            fail(f"scientific.txt missing valid bp_left/bp_right: {call}")
        distance = interval_distance(truth.start, truth.end, bp_left, bp_right)
        if best_distance is None or distance < best_distance:
            best_distance = distance
            best_row = dict(call)
            best_row["interval_distance_bp"] = str(distance)
    return best_row


def run_placer_sample(
    *,
    row: TruthRow,
    window: WindowSpec,
    args: argparse.Namespace,
    outdir: Path,
) -> Dict[str, str]:
    sample_dir = outdir / row.sample_id
    sample_dir.mkdir(parents=True, exist_ok=False)
    log_path = sample_dir / "run.log"
    read_names_path = sample_dir / "window_qnames.txt"
    subset_bam = sample_dir / "subset.bam"

    qname_count = extract_subset_bam(
        samtools=args.samtools,
        threads=args.threads,
        source_bam=args.bam,
        region=window.region,
        read_names_path=read_names_path,
        subset_bam=subset_bam,
        log_path=log_path,
    )

    completed = subprocess.run(
        [str(args.placer), str(subset_bam), str(args.ref), str(args.te)],
        cwd=str(sample_dir),
        check=False,
        text=True,
        capture_output=True,
    )
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(
            "$ " + " ".join([str(args.placer), str(subset_bam), str(args.ref), str(args.te)]) + "\n"
        )
        if completed.stdout:
            handle.write(completed.stdout)
            if not completed.stdout.endswith("\n"):
                handle.write("\n")
        if completed.stderr:
            handle.write(completed.stderr)
            if not completed.stderr.endswith("\n"):
                handle.write("\n")

    scientific_path = sample_dir / "scientific.txt"
    calls = parse_scientific_txt(scientific_path) if completed.returncode == 0 else []
    best_call = choose_best_call(row, calls)
    detected = (
        best_call is not None
        and int(best_call["interval_distance_bp"]) <= args.match_distance_bp
    )

    result = {
        "sample_id": row.sample_id,
        "uuid": row.uuid,
        "chrom": row.chrom,
        "truth_start": str(row.start),
        "truth_end": str(row.end),
        "truth_strand": row.strand,
        "truth_family": row.family,
        "truth_subfamily": row.subfamily,
        "truth_filter": row.filter_value,
        "window_start": str(window.window_start),
        "window_end": str(window.window_end),
        "window_region": window.region,
        "window_qname_count": str(qname_count),
        "window_qnames_txt": str(read_names_path),
        "subset_bam": str(subset_bam),
        "scientific_txt": str(scientific_path),
        "run_dir": str(sample_dir),
        "placer_exit_code": str(completed.returncode),
        "placer_call_count": str(len(calls)),
        "detected": "1" if detected else "0",
        "nearest_call_distance_bp": best_call["interval_distance_bp"] if best_call else "NA",
        "nearest_call_pos": best_call["pos"] if best_call else "NA",
        "nearest_call_bp_left": best_call["bp_left"] if best_call else "NA",
        "nearest_call_bp_right": best_call["bp_right"] if best_call else "NA",
        "nearest_call_te": best_call["te"] if best_call else "NA",
        "nearest_call_family": best_call["family"] if best_call else "NA",
        "nearest_call_subfamily": best_call["subfamily"] if best_call else "NA",
        "nearest_call_gt": best_call["gt"] if best_call else "NA",
        "nearest_call_af": best_call["af"] if best_call else "NA",
        "nearest_call_gq": best_call["gq"] if best_call else "NA",
        "nearest_call_qc": best_call["qc"] if best_call else "NA",
    }
    return result


def write_tsv(path: Path, rows: Sequence[Dict[str, str]]) -> None:
    if not rows:
        fail(f"no rows available for TSV output: {path}")
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    args.ground_truth = require_file(args.ground_truth, "ground truth table")
    args.bam = require_file(args.bam, "BAM")
    args.ref = require_file(args.ref, "reference FASTA")
    args.te = require_file(args.te, "TE FASTA")
    args.placer = require_executable(args.placer, "PLACER executable")
    if args.match_distance_bp < 0:
        fail("--match-distance-bp must be >= 0")
    if args.threads <= 0:
        fail("--threads must be > 0")

    seed = args.seed if args.seed is not None else random.SystemRandom().randrange(1, 2**31)
    outdir = resolve_outdir(args.outdir, seed)
    outdir.mkdir(parents=True, exist_ok=False)

    contig_lengths = load_contig_lengths(args.ref)
    truth_pool = load_truth_rows(args.ground_truth, args.truth_filter)
    sampled_rows = sample_truth_rows(truth_pool, args.sample_size, seed)

    sampled_manifest: List[Dict[str, str]] = []
    results: List[Dict[str, str]] = []
    for row in sampled_rows:
        window = build_window(row, args.window_flank_bp, contig_lengths)
        sampled_manifest.append(
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
                "window_start": str(window.window_start),
                "window_end": str(window.window_end),
                "window_region": window.region,
            }
        )
        results.append(run_placer_sample(row=row, window=window, args=args, outdir=outdir))

    write_tsv(outdir / "sampled_truth.tsv", sampled_manifest)
    write_tsv(outdir / "evaluation.tsv", results)

    detected_count = sum(1 for row in results if row["detected"] == "1")
    summary = {
        "seed": seed,
        "sample_size": args.sample_size,
        "truth_filter": args.truth_filter,
        "truth_pool_size": len(truth_pool),
        "window_flank_bp": args.window_flank_bp,
        "match_distance_bp": args.match_distance_bp,
        "detected_count": detected_count,
        "undetected_count": len(results) - detected_count,
        "detection_rate": detected_count / len(results),
        "ground_truth": str(args.ground_truth),
        "bam": str(args.bam),
        "ref": str(args.ref),
        "te": str(args.te),
        "placer": str(args.placer),
        "outdir": str(outdir),
    }
    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"evaluation_tsv\t{outdir / 'evaluation.tsv'}")
    print(f"sampled_truth_tsv\t{outdir / 'sampled_truth.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
