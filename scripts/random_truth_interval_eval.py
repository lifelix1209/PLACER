#!/usr/bin/env python3
"""Randomly sample truth loci and evaluate whether PLACER rediscovers them.

Evaluation semantics:
- Truth pool defaults to rows with `Filter == PASS` from the TLDR table.
- Each sampled locus is evaluated independently on a read-complete subset BAM.
  The subset keeps every alignment record for reads whose primary alignment
  overlaps the local truth window, so supplementary alignments outside the
  window are preserved without pulling in unrelated remote primary loci.
- A PLACER call is counted as detected when it is on the same chromosome and
  the minimum distance between `[bp_left, bp_right]` and `[Start, End]` is
  less than or equal to `--match-distance-bp`.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import os
import random
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence


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
    length_ins: int = 0


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
        help=(
            "Output directory. Default: placer_out/random_truth_eval_<timestamp>_seed<seed>. "
            "Must already exist when --reuse-existing-samples is set."
        ),
    )
    parser.add_argument(
        "--sampled-truth-tsv",
        type=Path,
        default=None,
        help=(
            "Replay an existing sampled cohort from sampled_truth.tsv instead of drawing a new "
            "random sample."
        ),
    )
    parser.add_argument(
        "--reuse-existing-samples",
        action="store_true",
        help=(
            "Do not rerun samtools/PLACER. Re-evaluate the existing SXX sample directories under "
            "--outdir using --sampled-truth-tsv and overwrite the top-level reports."
        ),
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
        "--min-length-ins",
        type=int,
        default=None,
        help="Only sample truth rows whose LengthIns is >= this threshold.",
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
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of sampled loci evaluated concurrently.",
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


def load_truth_rows(
    path: Path,
    truth_filter: str,
    min_length_ins: Optional[int] = None,
) -> List[TruthRow]:
    rows: List[TruthRow] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"missing header in truth table: {path}")
        required = {"UUID", "Chrom", "Start", "End", "Strand", "Family", "Subfamily", "Filter"}
        if min_length_ins is not None:
            required.add("LengthIns")
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

            length_ins = 0
            if min_length_ins is not None:
                try:
                    length_ins = int(row["LengthIns"])
                except ValueError:
                    fail(f"invalid LengthIns at truth row {row_index + 1}")
                if length_ins < min_length_ins:
                    continue

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
                    length_ins=length_ins,
                )
            )
    if not rows:
        fail(f"no truth rows matched Filter == {truth_filter!r} in {path}")
    return rows


def load_sampled_truth_rows(path: Path) -> List[TruthRow]:
    rows: List[TruthRow] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"missing header in sampled truth TSV: {path}")
        required = {
            "sample_id",
            "uuid",
            "chrom",
            "truth_start",
            "truth_end",
            "truth_strand",
            "truth_family",
            "truth_subfamily",
            "truth_filter",
            "truth_length_ins",
        }
        missing = required.difference(reader.fieldnames)
        if missing:
            fail(f"sampled truth TSV missing columns {sorted(missing)}: {path}")
        for row_index, row in enumerate(reader, start=1):
            sample_id = row["sample_id"].strip()
            if not sample_id:
                fail(f"empty sample_id at sampled truth row {row_index + 1}")
            try:
                start = int(row["truth_start"])
                end = int(row["truth_end"])
                length_ins = int(row["truth_length_ins"])
            except ValueError:
                fail(f"invalid sampled truth numeric field at row {row_index + 1}")
            if end < start:
                start, end = end, start
            rows.append(
                TruthRow(
                    sample_id=sample_id,
                    uuid=row["uuid"],
                    chrom=row["chrom"],
                    start=start,
                    end=end,
                    strand=row["truth_strand"],
                    family=row["truth_family"],
                    subfamily=row["truth_subfamily"],
                    filter_value=row["truth_filter"],
                    raw=dict(row),
                    length_ins=length_ins,
                )
            )
    if not rows:
        fail(f"no sampled rows loaded from {path}")
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
                length_ins=row.length_ins,
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
        fields = line.split("\t")
        if len(fields) < 2:
            fail(f"malformed samtools view row for region {region}: {line}")
        qname = fields[0]
        try:
            flag = int(fields[1])
        except ValueError:
            fail(f"invalid FLAG in samtools view row for region {region}: {line}")
        if flag & (0x4 | 0x100 | 0x800):
            continue
        if qname not in seen:
            seen.add(qname)
            qnames.append(qname)
    if not qnames:
        fail(f"no primary overlapping read names found for region {region}")
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


def parse_joint_component_metrics(log_path: Path) -> Dict[str, str]:
    metrics = {
        "joint_best_kind": "NA",
        "joint_runner_up_kind": "NA",
        "joint_final_qc": "NA",
        "joint_emit_te": "NA",
        "joint_emit_unknown": "NA",
        "consensus_qc": "NA",
        "consensus_len": "NA",
        "seg_qc": "NA",
        "seg_insert_len": "NA",
        "te_qc": "NA",
        "te_pass": "NA",
    }
    if not log_path.is_file():
        return metrics

    for raw_line in log_path.read_text(encoding="utf-8").splitlines():
        if "[Pipeline][joint]" not in raw_line:
            continue
        parsed: Dict[str, str] = {}
        for token in raw_line.split():
            if "=" not in token:
                continue
            key, value = token.split("=", 1)
            parsed[key] = value
        metrics = {
            "joint_best_kind": parsed.get("best_kind", "NA"),
            "joint_runner_up_kind": parsed.get("runner_kind", "NA"),
            "joint_final_qc": parsed.get("final_qc", "NA"),
            "joint_emit_te": parsed.get("emit_te", "NA"),
            "joint_emit_unknown": parsed.get("emit_unknown", "NA"),
            "consensus_qc": parsed.get("consensus_qc", "NA"),
            "consensus_len": parsed.get("consensus_len", "NA"),
            "seg_qc": parsed.get("seg_qc", "NA"),
            "seg_insert_len": parsed.get("seg_insert_len", "NA"),
            "te_qc": parsed.get("te_qc", "NA"),
            "te_pass": parsed.get("te_pass", "NA"),
        }
    return metrics


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


def count_nonempty_lines(path: Path) -> int:
    if not path.is_file():
        fail(f"required file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def classify_failure_stage(*, detected: bool, joint_metrics: Dict[str, str]) -> str:
    if detected:
        return "DETECTED"

    if joint_metrics["joint_final_qc"] == "NA":
        return "NO_EVENT_SIGNAL"

    if joint_metrics["consensus_qc"] not in {"NA", "PASS_EVENT_CONSENSUS"}:
        return "EVENT_BUT_NO_CONSENSUS"

    if joint_metrics["seg_qc"] not in {"NA", "PASS_EVENT_SEGMENTATION"}:
        return "CONSENSUS_BUT_NO_SEGMENTATION"

    seg_insert_len = joint_metrics["seg_insert_len"]
    if (
        seg_insert_len != "NA"
        and int(seg_insert_len) > 0
        and joint_metrics["te_qc"] not in {
            "PASS_INSERT_TE_ALIGNMENT",
            "PASS_INSERT_TE_ALIGNMENT_UNKNOWN",
        }
    ):
        return "SEGMENTED_BUT_NO_TE_ALIGNMENT"

    if (
        joint_metrics["joint_best_kind"] == "NON_TE_INSERTION"
        and (
            joint_metrics["te_qc"] in {
                "PASS_INSERT_TE_ALIGNMENT",
                "PASS_INSERT_TE_ALIGNMENT_UNKNOWN",
            }
            or joint_metrics["joint_runner_up_kind"] in {"TE_UNKNOWN", "TE_RESOLVED"}
        )
    ):
        return "TE_EVIDENCE_BUT_NON_TE_WINS"

    return "UNCLASSIFIED"


def build_result_row(
    *,
    row: TruthRow,
    window: WindowSpec,
    sample_dir: Path,
    read_names_path: Path,
    subset_bam: Path,
    scientific_path: Path,
    qname_count: int,
    placer_exit_code: str,
    calls: Sequence[Dict[str, str]],
    best_call: Optional[Dict[str, str]],
    detected: bool,
    extract_elapsed_s: float,
    placer_elapsed_s: float,
    sample_elapsed_s: float,
    joint_metrics: Dict[str, str],
) -> Dict[str, str]:
    normalized_joint_metrics = {
        "joint_best_kind": joint_metrics.get("joint_best_kind", "NA"),
        "joint_runner_up_kind": joint_metrics.get("joint_runner_up_kind", "NA"),
        "joint_final_qc": joint_metrics.get("joint_final_qc", "NA"),
        "joint_emit_te": joint_metrics.get("joint_emit_te", "NA"),
        "joint_emit_unknown": joint_metrics.get("joint_emit_unknown", "NA"),
        "consensus_qc": joint_metrics.get("consensus_qc", "NA"),
        "consensus_len": joint_metrics.get("consensus_len", "NA"),
        "seg_qc": joint_metrics.get("seg_qc", "NA"),
        "seg_insert_len": joint_metrics.get("seg_insert_len", "NA"),
        "te_qc": joint_metrics.get("te_qc", "NA"),
        "te_pass": joint_metrics.get("te_pass", "NA"),
    }
    failure_stage = classify_failure_stage(
        detected=detected,
        joint_metrics=normalized_joint_metrics,
    )
    return {
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
        "window_qname_count": str(qname_count),
        "window_qnames_txt": str(read_names_path),
        "subset_bam": str(subset_bam),
        "scientific_txt": str(scientific_path),
        "run_dir": str(sample_dir),
        "placer_exit_code": placer_exit_code,
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
        "extract_elapsed_s": f"{extract_elapsed_s:.6f}",
        "placer_elapsed_s": f"{placer_elapsed_s:.6f}",
        "sample_elapsed_s": f"{sample_elapsed_s:.6f}",
        "joint_best_kind": normalized_joint_metrics["joint_best_kind"],
        "joint_runner_up_kind": normalized_joint_metrics["joint_runner_up_kind"],
        "joint_final_qc": normalized_joint_metrics["joint_final_qc"],
        "joint_emit_te": normalized_joint_metrics["joint_emit_te"],
        "joint_emit_unknown": normalized_joint_metrics["joint_emit_unknown"],
        "consensus_qc": normalized_joint_metrics["consensus_qc"],
        "consensus_len": normalized_joint_metrics["consensus_len"],
        "seg_qc": normalized_joint_metrics["seg_qc"],
        "seg_insert_len": normalized_joint_metrics["seg_insert_len"],
        "te_qc": normalized_joint_metrics["te_qc"],
        "te_pass": normalized_joint_metrics["te_pass"],
        "failure_stage": failure_stage,
    }


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

    sample_started = time.time()
    extract_started = time.time()
    qname_count = extract_subset_bam(
        samtools=args.samtools,
        threads=args.threads,
        source_bam=args.bam,
        region=window.region,
        read_names_path=read_names_path,
        subset_bam=subset_bam,
        log_path=log_path,
    )
    extract_elapsed_s = time.time() - extract_started

    placer_started = time.time()
    completed = subprocess.run(
        [str(args.placer), str(subset_bam), str(args.ref), str(args.te)],
        cwd=str(sample_dir),
        check=False,
        text=True,
        capture_output=True,
    )
    placer_elapsed_s = time.time() - placer_started
    sample_elapsed_s = time.time() - sample_started
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
    joint_metrics = parse_joint_component_metrics(log_path)
    return build_result_row(
        row=row,
        window=window,
        sample_dir=sample_dir,
        read_names_path=read_names_path,
        subset_bam=subset_bam,
        scientific_path=scientific_path,
        qname_count=qname_count,
        placer_exit_code=str(completed.returncode),
        calls=calls,
        best_call=best_call,
        detected=detected,
        extract_elapsed_s=extract_elapsed_s,
        placer_elapsed_s=placer_elapsed_s,
        sample_elapsed_s=sample_elapsed_s,
        joint_metrics=joint_metrics,
    )


def evaluate_existing_sample(
    *,
    row: TruthRow,
    window: WindowSpec,
    args: argparse.Namespace,
    outdir: Path,
) -> Dict[str, str]:
    sample_dir = outdir / row.sample_id
    if not sample_dir.is_dir():
        fail(f"existing sample directory not found: {sample_dir}")
    read_names_path = sample_dir / "window_qnames.txt"
    subset_bam = sample_dir / "subset.bam"
    scientific_path = sample_dir / "scientific.txt"
    log_path = sample_dir / "run.log"
    if not subset_bam.is_file():
        fail(f"existing subset BAM not found: {subset_bam}")
    if not scientific_path.is_file():
        fail(f"existing scientific.txt not found: {scientific_path}")

    qname_count = count_nonempty_lines(read_names_path)
    calls = parse_scientific_txt(scientific_path)
    best_call = choose_best_call(row, calls)
    detected = (
        best_call is not None
        and int(best_call["interval_distance_bp"]) <= args.match_distance_bp
    )
    joint_metrics = parse_joint_component_metrics(log_path)
    return build_result_row(
        row=row,
        window=window,
        sample_dir=sample_dir,
        read_names_path=read_names_path,
        subset_bam=subset_bam,
        scientific_path=scientific_path,
        qname_count=qname_count,
        placer_exit_code="0",
        calls=calls,
        best_call=best_call,
        detected=detected,
        extract_elapsed_s=0.0,
        placer_elapsed_s=0.0,
        sample_elapsed_s=0.0,
        joint_metrics=joint_metrics,
    )


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


def evaluate_sample_rows(
    *,
    sampled_rows: Sequence[TruthRow],
    contig_lengths: Dict[str, int],
    args: argparse.Namespace,
    outdir: Path,
    evaluator: Optional[Callable[..., Dict[str, str]]] = None,
) -> tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    sampled_manifest: List[Dict[str, str]] = []
    tasks: List[tuple[int, TruthRow, WindowSpec]] = []
    for index, row in enumerate(sampled_rows):
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
                "truth_length_ins": str(row.length_ins),
                "window_start": str(window.window_start),
                "window_end": str(window.window_end),
                "window_region": window.region,
            }
        )
        tasks.append((index, row, window))

    eval_fn = evaluator
    if eval_fn is None:
        eval_fn = evaluate_existing_sample if args.reuse_existing_samples else run_placer_sample

    ordered_results: List[Optional[Dict[str, str]]] = [None] * len(tasks)
    worker_count = max(1, args.workers)

    if worker_count == 1 or len(tasks) <= 1:
        for index, row, window in tasks:
            ordered_results[index] = eval_fn(
                row=row,
                window=window,
                args=args,
                outdir=outdir,
            )
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=worker_count) as executor:
            future_to_index = {
                executor.submit(
                    eval_fn,
                    row=row,
                    window=window,
                    args=args,
                    outdir=outdir,
                ): index
                for index, row, window in tasks
            }
            for future in concurrent.futures.as_completed(future_to_index):
                ordered_results[future_to_index[future]] = future.result()

    return sampled_manifest, [result for result in ordered_results if result is not None]


def build_summary(
    *,
    seed: Optional[int],
    sampled_rows: Sequence[TruthRow],
    truth_pool_size: int,
    detected_count: int,
    failure_stage_counts: Dict[str, int],
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
        "match_distance_bp": args.match_distance_bp,
        "detected_count": detected_count,
        "undetected_count": len(sampled_rows) - detected_count,
        "detection_rate": detected_count / len(sampled_rows),
        "failure_stage_counts": failure_stage_counts,
        "workers": args.workers,
        "samtools_threads": args.threads,
        "total_elapsed_s": total_elapsed_s,
        "ground_truth": str(args.ground_truth),
        "bam": str(args.bam),
        "ref": str(args.ref),
        "te": str(args.te),
        "placer": str(args.placer),
        "outdir": str(outdir),
    }
    if args.sampled_truth_tsv is not None:
        summary["sampled_truth_tsv"] = str(args.sampled_truth_tsv)
    return summary


def main(argv: Sequence[str]) -> int:
    started = time.time()
    args = parse_args(argv)
    if args.sampled_truth_tsv is not None:
        args.sampled_truth_tsv = require_file(args.sampled_truth_tsv, "sampled truth TSV")
    args.ground_truth = require_file(args.ground_truth, "ground truth table")
    args.bam = require_file(args.bam, "BAM")
    args.ref = require_file(args.ref, "reference FASTA")
    args.te = require_file(args.te, "TE FASTA")
    args.placer = require_executable(args.placer, "PLACER executable")
    if args.match_distance_bp < 0:
        fail("--match-distance-bp must be >= 0")
    if args.min_length_ins is not None and args.min_length_ins < 0:
        fail("--min-length-ins must be >= 0")
    if args.threads <= 0:
        fail("--threads must be > 0")
    if args.workers <= 0:
        fail("--workers must be > 0")
    if args.reuse_existing_samples and args.sampled_truth_tsv is None:
        fail("--reuse-existing-samples requires --sampled-truth-tsv")

    if args.sampled_truth_tsv is None:
        seed = args.seed if args.seed is not None else random.SystemRandom().randrange(1, 2**31)
    else:
        seed = args.seed

    if args.reuse_existing_samples:
        if args.outdir is None:
            fail("--reuse-existing-samples requires --outdir")
        outdir = args.outdir.expanduser().resolve()
        if not outdir.is_dir():
            fail(f"existing output directory not found: {outdir}")
    else:
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

    sampled_manifest, results = evaluate_sample_rows(
        sampled_rows=sampled_rows,
        contig_lengths=contig_lengths,
        args=args,
        outdir=outdir,
    )

    write_tsv(outdir / "sampled_truth.tsv", sampled_manifest)
    write_tsv(outdir / "evaluation.tsv", results)

    detected_count = sum(1 for row in results if row["detected"] == "1")
    failure_stage_counts: Dict[str, int] = {}
    for row in results:
        stage = row["failure_stage"]
        failure_stage_counts[stage] = failure_stage_counts.get(stage, 0) + 1
    total_elapsed_s = time.time() - started
    summary = build_summary(
        seed=seed,
        sampled_rows=sampled_rows,
        truth_pool_size=len(truth_pool),
        detected_count=detected_count,
        failure_stage_counts=failure_stage_counts,
        args=args,
        outdir=outdir,
        total_elapsed_s=total_elapsed_s,
    )
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
