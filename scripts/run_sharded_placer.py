#!/usr/bin/env python3
"""Shard BAM by contig/region, run PLACER per shard, and merge scientific outputs.

Design goals:
- Keep all PLACER modules enabled (no algorithm feature drop).
- Improve scalability on very large BAMs through data sharding.
- Preserve accuracy:
  * `contig` mode is exact (no overlap, no duplicated loci).
  * `region` mode uses overlap padding + core-range filtering at merge time.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import dataclasses
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple


SUMMARY_INT_KEYS = [
    "total_reads",
    "gate1_passed",
    "processed_bins",
    "components",
    "event_consensus_calls",
    "genotype_calls",
    "final_pass_calls",
]

EXPECTED_SCIENTIFIC_HEADER = [
    "chrom",
    "pos",
    "bp_left",
    "bp_right",
    "te",
    "family",
    "subfamily",
    "strand",
    "insert_len",
    "support_reads",
    "alt_struct_reads",
    "ref_span_reads",
    "low_mapq_ref_span_reads",
    "gt",
    "af",
    "gq",
    "best_te_identity",
    "best_te_query_coverage",
    "cross_family_margin",
    "tsd_type",
    "tsd_len",
    "left_flank_align_len",
    "right_flank_align_len",
    "consensus_len",
    "qc",
]


@dataclasses.dataclass(frozen=True)
class ContigInfo:
    name: str
    length: int
    mapped: int
    unmapped: int


@dataclasses.dataclass(frozen=True)
class ShardSpec:
    shard_id: int
    label: str
    chrom: str
    core_start: int
    core_end: int
    fetch_start: int
    fetch_end: int

    @property
    def region(self) -> str:
        if self.fetch_start <= 1 and self.fetch_end <= 0:
            return self.chrom
        return f"{self.chrom}:{self.fetch_start}-{self.fetch_end}"

    @property
    def core_region(self) -> str:
        return f"{self.chrom}:{self.core_start}-{self.core_end}"


@dataclasses.dataclass
class ShardResult:
    spec: ShardSpec
    workdir: Path
    scientific_path: Path
    summary: Dict[str, str]
    n_rows_raw: int
    elapsed_s: float


@dataclasses.dataclass
class ShardProgressState:
    state: str
    stage: str
    submitted_s: float
    started_s: float
    stage_started_s: float
    updated_s: float
    note: str = ""
    elapsed_s: float = 0.0
    rows: int = 0


def format_duration(seconds: float) -> str:
    total = max(0, int(seconds))
    hours, rem = divmod(total, 3600)
    minutes, secs = divmod(rem, 60)
    if hours > 0:
        return f"{hours}h{minutes:02d}m{secs:02d}s"
    if minutes > 0:
        return f"{minutes}m{secs:02d}s"
    return f"{secs}s"


def format_progress_bar(done: int, total: int, width: int = 28) -> str:
    denom = max(1, total)
    filled = min(width, int(width * max(0, done) / denom))
    return "[" + ("#" * filled) + ("-" * (width - filled)) + "]"


class ProgressTracker:
    def __init__(self, shards: Sequence[ShardSpec]) -> None:
        now = time.time()
        self._created_s = now
        self._lock = threading.Lock()
        self._state: Dict[str, ShardProgressState] = {
            spec.label: ShardProgressState(
                state="queued",
                stage="queued",
                submitted_s=now,
                started_s=0.0,
                stage_started_s=now,
                updated_s=now,
                note=spec.region,
            )
            for spec in shards
        }

    def update_stage(self, spec: ShardSpec, stage: str, note: str = "") -> None:
        now = time.time()
        with self._lock:
            rec = self._state[spec.label]
            if rec.started_s <= 0.0:
                rec.started_s = now
            rec.state = "running"
            rec.stage = stage
            rec.stage_started_s = now
            rec.updated_s = now
            if note:
                rec.note = note

    def mark_done(self, spec: ShardSpec, rows: int, elapsed_s: float) -> None:
        now = time.time()
        with self._lock:
            rec = self._state[spec.label]
            if rec.started_s <= 0.0:
                rec.started_s = now
            rec.state = "done"
            rec.stage = "done"
            rec.stage_started_s = now
            rec.updated_s = now
            rec.elapsed_s = elapsed_s
            rec.rows = rows
            rec.note = f"rows={rows}"

    def mark_failed(self, spec: ShardSpec, message: str) -> None:
        now = time.time()
        with self._lock:
            rec = self._state[spec.label]
            if rec.started_s <= 0.0:
                rec.started_s = now
            rec.state = "failed"
            rec.stage = "failed"
            rec.stage_started_s = now
            rec.updated_s = now
            if message:
                rec.note = message

    def counts(self) -> Dict[str, int]:
        with self._lock:
            done = 0
            failed = 0
            running = 0
            queued = 0
            for rec in self._state.values():
                if rec.state == "done":
                    done += 1
                elif rec.state == "failed":
                    failed += 1
                elif rec.state == "running":
                    running += 1
                else:
                    queued += 1
        return {
            "done": done,
            "failed": failed,
            "running": running,
            "queued": queued,
        }

    def render(self, total: int, max_active: int = 4) -> str:
        now = time.time()
        with self._lock:
            done = 0
            failed = 0
            running = 0
            queued = 0
            for rec in self._state.values():
                if rec.state == "done":
                    done += 1
                elif rec.state == "failed":
                    failed += 1
                elif rec.state == "running":
                    running += 1
                else:
                    queued += 1
            counts = {
                "done": done,
                "failed": failed,
                "running": running,
                "queued": queued,
            }
            active = [
                (label, rec)
                for label, rec in self._state.items()
                if rec.state == "running"
            ]
        finished = counts["done"] + counts["failed"]
        bar = format_progress_bar(finished, total)
        line = (
            f"[sharded] progress {bar} {finished}/{max(1, total)} "
            f"done={counts['done']} running={counts['running']} "
            f"queued={counts['queued']} failed={counts['failed']} "
            f"elapsed={format_duration(now - self._created_s)}"
        )
        if not active:
            return line

        active.sort(key=lambda item: item[1].stage_started_s)
        stage_counts: Dict[str, int] = {}
        for _, rec in active:
            stage_counts[rec.stage] = stage_counts.get(rec.stage, 0) + 1
        stage_text = ",".join(f"{name}={count}" for name, count in sorted(stage_counts.items()))

        samples: List[str] = []
        for label, rec in active[:max_active]:
            stage_age = format_duration(now - rec.stage_started_s)
            samples.append(f"{label}:{rec.stage}:{stage_age}")
        return f"{line} | stages={stage_text} | active={' ; '.join(samples)}"


def run_command(
    cmd: Sequence[str],
    *,
    cwd: Optional[Path] = None,
    env: Optional[Dict[str, str]] = None,
    stdout_path: Optional[Path] = None,
    stderr_path: Optional[Path] = None,
    append: bool = False,
) -> None:
    stdout_fh = None
    stderr_fh = None
    try:
        file_mode = "a" if append else "w"
        if stdout_path is not None:
            stdout_path.parent.mkdir(parents=True, exist_ok=True)
            stdout_fh = stdout_path.open(file_mode, encoding="utf-8")
        if stderr_path is not None:
            stderr_path.parent.mkdir(parents=True, exist_ok=True)
            stderr_fh = stderr_path.open(file_mode, encoding="utf-8")
        subprocess.run(
            list(cmd),
            cwd=str(cwd) if cwd else None,
            env=env,
            check=True,
            stdout=stdout_fh if stdout_fh else subprocess.DEVNULL,
            stderr=stderr_fh if stderr_fh else subprocess.DEVNULL,
        )
    finally:
        if stdout_fh:
            stdout_fh.close()
        if stderr_fh:
            stderr_fh.close()


def parse_idxstats_text(text: str) -> List[ContigInfo]:
    contigs: List[ContigInfo] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        name = fields[0]
        if name == "*":
            continue
        try:
            length = int(fields[1])
            mapped = int(fields[2])
            unmapped = int(fields[3])
        except ValueError:
            continue
        if length <= 0:
            continue
        contigs.append(ContigInfo(name=name, length=length, mapped=mapped, unmapped=unmapped))
    return contigs


def load_contigs(samtools: str, bam: Path, min_mapped: int) -> List[ContigInfo]:
    cp = subprocess.run(
        [samtools, "idxstats", str(bam)],
        check=True,
        text=True,
        capture_output=True,
    )
    contigs = parse_idxstats_text(cp.stdout)
    out = [c for c in contigs if c.mapped >= min_mapped]
    return out


def sanitize_label(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s)


def build_shards(
    contigs: Sequence[ContigInfo],
    *,
    mode: str,
    region_size: int,
    overlap_bp: int,
) -> List[ShardSpec]:
    shards: List[ShardSpec] = []
    sid = 0
    for contig in contigs:
        if mode == "contig":
            sid += 1
            label = sanitize_label(f"{sid:04d}_{contig.name}")
            shards.append(
                ShardSpec(
                    shard_id=sid,
                    label=label,
                    chrom=contig.name,
                    core_start=1,
                    core_end=contig.length,
                    fetch_start=1,
                    fetch_end=contig.length,
                )
            )
            continue

        start = 1
        while start <= contig.length:
            core_start = start
            core_end = min(contig.length, start + region_size - 1)
            fetch_start = max(1, core_start - overlap_bp)
            fetch_end = min(contig.length, core_end + overlap_bp)
            sid += 1
            label = sanitize_label(
                f"{sid:04d}_{contig.name}_{core_start}_{core_end}"
            )
            shards.append(
                ShardSpec(
                    shard_id=sid,
                    label=label,
                    chrom=contig.name,
                    core_start=core_start,
                    core_end=core_end,
                    fetch_start=fetch_start,
                    fetch_end=fetch_end,
                )
            )
            start += region_size
    return shards


def parse_scientific(path: Path) -> Tuple[Dict[str, str], List[str], List[Dict[str, str]]]:
    summary: Dict[str, str] = {}
    header: List[str] = []
    rows: List[Dict[str, str]] = []

    in_table = False
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not in_table:
                if line.startswith("#chrom\t"):
                    header = line[1:].split("\t")
                    in_table = True
                    continue
                if not line or line.startswith("#"):
                    continue
                kv = line.split("\t", 1)
                if len(kv) == 2:
                    summary[kv[0]] = kv[1]
                continue

            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < len(header):
                fields = fields + ([""] * (len(header) - len(fields)))
            elif len(fields) > len(header):
                fields = fields[: len(header)]
            rows.append(dict(zip(header, fields)))

    if not header:
        raise RuntimeError(f"failed to parse table header from {path}")
    if header != EXPECTED_SCIENTIFIC_HEADER:
        missing = [col for col in EXPECTED_SCIENTIFIC_HEADER if col not in header]
        unexpected = [col for col in header if col not in EXPECTED_SCIENTIFIC_HEADER]
        problems: List[str] = []
        if missing:
            problems.append(f"missing={','.join(missing)}")
        if unexpected:
            problems.append(f"unexpected={','.join(unexpected)}")
        detail = " ".join(problems) if problems else "column_order_mismatch"
        raise RuntimeError(f"unexpected scientific.txt schema in {path}: {detail}")
    return summary, header, rows


def parse_int(value: str, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def parse_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def same_call_locus(a: Dict[str, str], b: Dict[str, str], dedup_bp: int) -> bool:
    if a.get("chrom", "") != b.get("chrom", ""):
        return False
    return abs(parse_int(a.get("pos", "-1"), -1) - parse_int(b.get("pos", "-1"), -1)) <= dedup_bp


def prefer_new_call(cur: Dict[str, str], incumbent: Dict[str, str]) -> bool:
    cur_support = parse_int(cur.get("support_reads", "0"), 0)
    old_support = parse_int(incumbent.get("support_reads", "0"), 0)
    if cur_support != old_support:
        return cur_support > old_support

    cur_gq = parse_int(cur.get("gq", "0"), 0)
    old_gq = parse_int(incumbent.get("gq", "0"), 0)
    if cur_gq != old_gq:
        return cur_gq > old_gq

    cur_margin = parse_float(cur.get("cross_family_margin", "0"), 0.0)
    old_margin = parse_float(incumbent.get("cross_family_margin", "0"), 0.0)
    if cur_margin != old_margin:
        return cur_margin > old_margin

    cur_te_identity = parse_float(cur.get("best_te_identity", "0"), 0.0)
    old_te_identity = parse_float(incumbent.get("best_te_identity", "0"), 0.0)
    if cur_te_identity != old_te_identity:
        return cur_te_identity > old_te_identity

    cur_te = cur.get("te", "")
    old_te = incumbent.get("te", "")
    if cur_te != old_te:
        if old_te in {"", "UNK", "NA"}:
            return True
        if cur_te in {"", "UNK", "NA"}:
            return False

    cur_consensus_len = parse_int(cur.get("consensus_len", "0"), 0)
    old_consensus_len = parse_int(incumbent.get("consensus_len", "0"), 0)
    if cur_consensus_len != old_consensus_len:
        return cur_consensus_len > old_consensus_len

    cur_query_coverage = parse_float(cur.get("best_te_query_coverage", "0"), 0.0)
    old_query_coverage = parse_float(incumbent.get("best_te_query_coverage", "0"), 0.0)
    if cur_query_coverage != old_query_coverage:
        return cur_query_coverage > old_query_coverage

    return parse_int(cur.get("pos", "0"), 0) < parse_int(incumbent.get("pos", "0"), 0)


def dedup_calls(
    rows: List[Dict[str, str]],
    *,
    chrom_order: Dict[str, int],
    dedup_bp: int,
) -> List[Dict[str, str]]:
    def sort_key(row: Dict[str, str]) -> Tuple[int, int, int, int, str, str]:
        chrom = row.get("chrom", "")
        return (
            chrom_order.get(chrom, 10**9),
            parse_int(row.get("pos", "-1"), -1),
            parse_int(row.get("bp_left", "-1"), -1),
            parse_int(row.get("bp_right", "-1"), -1),
            chrom,
            row.get("te", ""),
        )

    ordered = sorted(rows, key=sort_key)
    out: List[Dict[str, str]] = []
    for row in ordered:
        if out and same_call_locus(row, out[-1], dedup_bp):
            if prefer_new_call(row, out[-1]):
                out[-1] = row
            continue
        out.append(row)
    return out


def sum_summary_key(shards: Sequence[ShardResult], key: str) -> int:
    total = 0
    for shard in shards:
        total += parse_int(shard.summary.get(key, "0"), 0)
    return total


def merge_shard_results(
    shards: Sequence[ShardResult],
    *,
    mode: str,
    dedup_bp: int,
    chrom_order: Dict[str, int],
    region_size: int,
    overlap_bp: int,
    out_path: Path,
) -> None:
    if not shards:
        raise RuntimeError("no shard results to merge")

    base_header: Optional[List[str]] = None
    merged_rows: List[Dict[str, str]] = []

    for shard in shards:
        _, header, rows = parse_scientific(shard.scientific_path)
        if base_header is None:
            base_header = header
        elif header != base_header:
            raise RuntimeError(
                f"inconsistent scientific header in {shard.scientific_path}"
            )

        for row in rows:
            pos = parse_int(row.get("pos", "-1"), -1)
            if pos < 0:
                continue
            if mode == "region":
                if pos < shard.spec.core_start or pos > shard.spec.core_end:
                    continue
            merged_rows.append(row)

    assert base_header is not None
    deduped = dedup_calls(merged_rows, chrom_order=chrom_order, dedup_bp=dedup_bp)

    if mode == "contig":
        summary_values = {k: sum_summary_key(shards, k) for k in SUMMARY_INT_KEYS}
    else:
        summary_values = {k: -1 for k in SUMMARY_INT_KEYS}
    summary_values["final_pass_calls"] = len(deduped)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        out.write("#PLACER sharded merge summary\n")
        for key in [
            "total_reads",
            "gate1_passed",
            "processed_bins",
            "components",
            "event_consensus_calls",
            "genotype_calls",
            "final_pass_calls",
        ]:
            out.write(f"{key}\t{summary_values.get(key, 0)}\n")
        out.write(f"merge_mode\t{mode}\n")
        out.write(f"merge_dedup_distance_bp\t{dedup_bp}\n")
        if mode == "region":
            out.write(f"region_size_bp\t{region_size}\n")
            out.write(f"region_overlap_bp\t{overlap_bp}\n")
        out.write("schema_version\t1.0.0-sharded\n")
        out.write("\n")
        out.write("#" + "\t".join(base_header) + "\n")
        for row in deduped:
            out.write("\t".join(row.get(col, "") for col in base_header) + "\n")


def run_single_shard(
    spec: ShardSpec,
    *,
    bam: Path,
    ref: Path,
    te: Optional[Path],
    placer_bin: Path,
    samtools_bin: str,
    shard_root: Path,
    samtools_threads: int,
    extra_env: Dict[str, str],
    progress_cb: Optional[Callable[[ShardSpec, str, str], None]] = None,
) -> ShardResult:
    t0 = time.time()
    workdir = shard_root / spec.label
    workdir.mkdir(parents=True, exist_ok=True)

    shard_bam = workdir / "input.shard.bam"
    extract_stdout = workdir / "extract.stdout.log"
    extract_stderr = workdir / "extract.stderr.log"
    run_stdout = workdir / "placer.stdout.log"
    run_stderr = workdir / "placer.stderr.log"

    for log_path in (extract_stdout, extract_stderr, run_stdout, run_stderr):
        if log_path.exists():
            log_path.unlink()

    def update(stage: str, note: str = "") -> None:
        if progress_cb is not None:
            progress_cb(spec, stage, note)

    try:
        update("extract", spec.region)
        run_command(
            [
                samtools_bin,
                "view",
                "-@",
                str(max(1, samtools_threads)),
                "-b",
                "-o",
                str(shard_bam),
                str(bam),
                spec.region,
            ],
            stdout_path=extract_stdout,
            stderr_path=extract_stderr,
        )
        update("index", shard_bam.name)
        run_command(
            [
                samtools_bin,
                "index",
                "-@",
                str(max(1, samtools_threads)),
                str(shard_bam),
            ],
            stdout_path=extract_stdout,
            stderr_path=extract_stderr,
            append=True,
        )

        cmd = [str(placer_bin), str(shard_bam), str(ref)]
        if te is not None:
            cmd.append(str(te))
        env = os.environ.copy()
        env.update(extra_env)
        update("placer", spec.core_region)
        run_command(
            cmd,
            cwd=workdir,
            env=env,
            stdout_path=run_stdout,
            stderr_path=run_stderr,
        )

        scientific = workdir / "scientific.txt"
        if not scientific.exists():
            raise RuntimeError(f"scientific.txt not generated for shard {spec.label}")

        summary, _, rows = parse_scientific(scientific)
        elapsed_s = time.time() - t0
        update("postprocess", f"rows={len(rows)}")
        return ShardResult(
            spec=spec,
            workdir=workdir,
            scientific_path=scientific,
            summary=summary,
            n_rows_raw=len(rows),
            elapsed_s=elapsed_s,
        )
    except Exception as exc:
        if progress_cb is not None:
            progress_cb(spec, "failed", str(exc))
        raise


def parse_env_kv(env_args: Sequence[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in env_args:
        if "=" not in item:
            raise ValueError(f"invalid --env (expected KEY=VALUE): {item}")
        k, v = item.split("=", 1)
        k = k.strip()
        if not k:
            raise ValueError(f"invalid --env key: {item}")
        out[k] = v
    return out


def maybe_build_repo_placer(placer_arg: str, build_if_needed: bool) -> Path:
    requested = Path(placer_arg)
    requested_abs = (requested if requested.is_absolute() else Path.cwd() / requested).resolve()
    if not build_if_needed:
        return requested_abs

    repo_root = Path(__file__).resolve().parent.parent
    repo_binary = (repo_root / "build" / "placer").resolve()
    if requested_abs != repo_binary:
        return requested_abs

    build_helper = repo_root / "scripts" / "build_latest_placer.sh"
    if not build_helper.exists():
        return requested_abs

    print("[sharded] checking repo build before shard runs", file=sys.stderr, flush=True)
    cp = subprocess.run(
        [str(build_helper)],
        check=True,
        text=True,
        capture_output=True,
    )
    if cp.stderr:
        print(cp.stderr, file=sys.stderr, end="", flush=True)

    built_bin = cp.stdout.strip()
    if not built_bin:
        raise RuntimeError("build_latest_placer.sh did not return a binary path")
    return Path(built_bin).resolve()


def write_manifest(path: Path, shards: Sequence[ShardResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as out:
        out.write(
            "shard_id\tlabel\tregion\tcore_region\traw_calls\telapsed_s\tworkdir\tscientific\n"
        )
        for s in shards:
            out.write(
                "\t".join(
                    [
                        str(s.spec.shard_id),
                        s.spec.label,
                        s.spec.region,
                        s.spec.core_region,
                        str(s.n_rows_raw),
                        f"{s.elapsed_s:.3f}",
                        str(s.workdir),
                        str(s.scientific_path),
                    ]
                )
                + "\n"
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run PLACER by BAM shards and merge scientific outputs."
    )
    parser.add_argument("--bam", required=True, help="Input BAM (indexed).")
    parser.add_argument("--ref", required=True, help="Reference FASTA.")
    parser.add_argument("--te", default="", help="TE FASTA (optional).")
    parser.add_argument("--placer", default="build/placer", help="PLACER binary path.")
    parser.add_argument("--samtools", default="samtools", help="samtools executable.")
    parser.add_argument(
        "--mode",
        choices=["contig", "region"],
        default="contig",
        help="contig: exact by chromosome; region: fixed windows with overlap.",
    )
    parser.add_argument("--region-size", type=int, default=50_000_000)
    parser.add_argument("--overlap-bp", type=int, default=200_000)
    parser.add_argument("--min-mapped-reads", type=int, default=1)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) // 2))
    parser.add_argument("--samtools-threads", type=int, default=1)
    parser.add_argument("--dedup-bp", type=int, default=50)
    parser.add_argument("--outdir", default="sharded_placer_out")
    parser.add_argument(
        "--build-if-needed",
        dest="build_if_needed",
        action="store_true",
        default=True,
        help="Auto-build build/placer before running shards when using the repo binary path.",
    )
    parser.add_argument(
        "--no-build-if-needed",
        dest="build_if_needed",
        action="store_false",
        help="Skip the pre-run build freshness check.",
    )
    parser.add_argument(
        "--progress-heartbeat-s",
        type=float,
        default=30.0,
        help="Print a progress heartbeat if no shard completes within this many seconds.",
    )
    parser.add_argument(
        "--env",
        action="append",
        default=[],
        help="Extra env for PLACER subprocess, can repeat. Example: --env PLACER_PARALLEL=1",
    )
    parser.add_argument(
        "--keep-shard-bam",
        action="store_true",
        help="Keep per-shard BAM/BAM.bai files after merge.",
    )
    args = parser.parse_args()

    bam = Path(args.bam).resolve()
    ref = Path(args.ref).resolve()
    te = Path(args.te).resolve() if args.te else None
    placer_bin = maybe_build_repo_placer(args.placer, args.build_if_needed)
    outdir = Path(args.outdir).resolve()
    shard_root = outdir / "shards"
    merged_scientific = outdir / "scientific.sharded.txt"
    manifest_path = outdir / "shard_manifest.tsv"

    if not bam.exists():
        print(f"[sharded] BAM not found: {bam}", file=sys.stderr)
        return 2
    if not ref.exists():
        print(f"[sharded] REF not found: {ref}", file=sys.stderr)
        return 2
    if te is not None and not te.exists():
        print(f"[sharded] TE FASTA not found: {te}", file=sys.stderr)
        return 2
    if not placer_bin.exists():
        print(f"[sharded] PLACER binary not found: {placer_bin}", file=sys.stderr)
        return 2
    if shutil.which(args.samtools) is None:
        print(f"[sharded] samtools not found in PATH: {args.samtools}", file=sys.stderr)
        return 2

    # Ensure BAM index exists and is valid for idxstats.
    try:
        subprocess.run(
            [args.samtools, "idxstats", str(bam)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError:
        print(f"[sharded] BAM index missing/invalid, building: {bam}", file=sys.stderr)
        subprocess.run([args.samtools, "index", str(bam)], check=True)

    contigs = load_contigs(args.samtools, bam, max(0, args.min_mapped_reads))
    if not contigs:
        print("[sharded] no contigs pass min-mapped filter", file=sys.stderr)
        return 2

    shards = build_shards(
        contigs,
        mode=args.mode,
        region_size=max(1, args.region_size),
        overlap_bp=max(0, args.overlap_bp),
    )
    if not shards:
        print("[sharded] no shards generated", file=sys.stderr)
        return 2

    extra_env = parse_env_kv(args.env)
    outdir.mkdir(parents=True, exist_ok=True)
    shard_root.mkdir(parents=True, exist_ok=True)
    tracker = ProgressTracker(shards)

    print(
        f"[sharded] mode={args.mode} contigs={len(contigs)} shards={len(shards)} "
        f"workers={max(1, args.workers)}"
        f" heartbeat={max(0.0, args.progress_heartbeat_s):.1f}s",
        flush=True,
    )
    print(tracker.render(len(shards)), flush=True)

    shard_results: List[ShardResult] = []
    failures: List[Tuple[ShardSpec, str]] = []

    def _submit_one(spec: ShardSpec) -> ShardResult:
        return run_single_shard(
            spec,
            bam=bam,
            ref=ref,
            te=te,
            placer_bin=placer_bin,
            samtools_bin=args.samtools,
            shard_root=shard_root,
            samtools_threads=max(1, args.samtools_threads),
            extra_env=extra_env,
            progress_cb=tracker.update_stage,
        )

    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        future_to_spec = {ex.submit(_submit_one, spec): spec for spec in shards}
        total = len(shards)
        pending = set(future_to_spec)
        heartbeat_s = max(0.0, args.progress_heartbeat_s)
        while pending:
            timeout = heartbeat_s if heartbeat_s > 0.0 else None
            done_futs, pending = concurrent.futures.wait(
                pending,
                timeout=timeout,
                return_when=concurrent.futures.FIRST_COMPLETED,
            )
            if not done_futs:
                print(tracker.render(total), flush=True)
                continue
            for fut in done_futs:
                spec = future_to_spec[fut]
                try:
                    result = fut.result()
                    shard_results.append(result)
                    tracker.mark_done(spec, result.n_rows_raw, result.elapsed_s)
                    counts = tracker.counts()
                    finished = counts["done"] + counts["failed"]
                    print(
                        f"[sharded] done {format_progress_bar(finished, total)} "
                        f"{finished}/{total} {spec.label} "
                        f"rows={result.n_rows_raw} time={result.elapsed_s:.2f}s",
                        flush=True,
                    )
                except Exception as exc:  # noqa: BLE001
                    failures.append((spec, str(exc)))
                    tracker.mark_failed(spec, str(exc))
                    counts = tracker.counts()
                    finished = counts["done"] + counts["failed"]
                    print(
                        f"[sharded] FAIL {format_progress_bar(finished, total)} "
                        f"{finished}/{total} {spec.label}: {exc}",
                        file=sys.stderr,
                        flush=True,
                    )

    if failures:
        print(f"[sharded] {len(failures)} shard(s) failed; abort merge", file=sys.stderr)
        for spec, msg in failures:
            print(f"  - {spec.label}: {msg}", file=sys.stderr)
        return 3

    chrom_order = {c.name: i for i, c in enumerate(contigs)}
    shard_results.sort(key=lambda x: x.spec.shard_id)
    merge_shard_results(
        shard_results,
        mode=args.mode,
        dedup_bp=max(0, args.dedup_bp),
        chrom_order=chrom_order,
        region_size=max(1, args.region_size),
        overlap_bp=max(0, args.overlap_bp),
        out_path=merged_scientific,
    )
    write_manifest(manifest_path, shard_results)

    if not args.keep_shard_bam:
        for result in shard_results:
            bam_path = result.workdir / "input.shard.bam"
            bai_path = result.workdir / "input.shard.bam.bai"
            if bam_path.exists():
                bam_path.unlink()
            if bai_path.exists():
                bai_path.unlink()

    print(f"[sharded] merged scientific: {merged_scientific}", flush=True)
    print(f"[sharded] manifest: {manifest_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
