#!/usr/bin/env python3

from __future__ import annotations

import json
import random
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List


ROUND_COUNT = 5
SAMPLE_SIZE = 150
MIN_DETECTION_RATE = 0.90


@dataclass(frozen=True)
class ReleaseRoundSummary:
    round_index: int
    seed: int
    detection_rate: float
    pass_round: bool
    summary_json_path: str


def fail(message: str) -> "NoReturn":
    raise SystemExit(message)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def require_file(path: Path, label: str) -> Path:
    resolved = path.resolve()
    if not resolved.is_file():
        fail(f"{label} not found: {resolved}")
    return resolved


def require_executable(path: Path, label: str) -> Path:
    resolved = require_file(path, label)
    if not (resolved.stat().st_mode & 0o111):
        fail(f"{label} is not executable: {resolved}")
    return resolved


def generate_unique_seeds(round_count: int) -> List[int]:
    rng = random.SystemRandom()
    seeds: List[int] = []
    seen = set()
    while len(seeds) < round_count:
        seed = rng.randrange(1, 2**31)
        if seed in seen:
            continue
        seen.add(seed)
        seeds.append(seed)
    return seeds


def evaluate_release_gate(
    rounds: List[ReleaseRoundSummary],
    *,
    round_count: int,
    min_detection_rate: float,
) -> Dict[str, Any]:
    failed_round_indices = [round_.round_index for round_ in rounds if not round_.pass_round]
    return {
        "required_rounds": round_count,
        "required_detection_rate_strict_gt": min_detection_rate,
        "all_rounds_pass": len(rounds) == round_count and not failed_round_indices,
        "failed_round_indices": failed_round_indices,
        "rounds": [asdict(round_) for round_ in rounds],
    }


def load_detection_rate(summary_path: Path) -> float:
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    detection_rate = payload.get("detection_rate")
    if not isinstance(detection_rate, (int, float)):
        fail(f"invalid detection_rate in {summary_path}")
    return float(detection_rate)


def run_round(
    *,
    round_index: int,
    seed: int,
    root: Path,
    out_root: Path,
    truth_table: Path,
    bam: Path,
    ref_fasta: Path,
    te_fasta: Path,
    placer_bin: Path,
) -> ReleaseRoundSummary:
    round_dir = out_root / f"round_{round_index:02d}"
    command = [
        sys.executable,
        str(root / "scripts" / "random_truth_interval_eval.py"),
        "--ground-truth",
        str(truth_table),
        "--bam",
        str(bam),
        "--ref",
        str(ref_fasta),
        "--te",
        str(te_fasta),
        "--placer",
        str(placer_bin),
        "--sample-size",
        str(SAMPLE_SIZE),
        "--seed",
        str(seed),
        "--threads",
        "1",
        "--workers",
        "1",
        "--outdir",
        str(round_dir),
    ]
    subprocess.run(command, check=True)

    summary_path = round_dir / "summary.json"
    if not summary_path.is_file():
        fail(f"missing summary.json for round {round_index}: {summary_path}")

    detection_rate = load_detection_rate(summary_path)
    return ReleaseRoundSummary(
        round_index=round_index,
        seed=seed,
        detection_rate=detection_rate,
        pass_round=detection_rate > MIN_DETECTION_RATE,
        summary_json_path=str(summary_path),
    )


def main() -> int:
    root = repo_root()
    truth_table = require_file(root / "test_data" / "Yohann_D23.table.txt", "ground truth table")
    bam = require_file(root / "test_data" / "final_merged_Yohann_D23.bam", "BAM")
    require_file(root / "test_data" / "final_merged_Yohann_D23.bam.bai", "BAM index")
    ref_fasta = require_file(
        root / "test_data" / "Astatotilapia_calliptera.fAstCal1.3.dna.toplevel.fa",
        "reference FASTA",
    )
    require_file(
        root / "test_data" / "Astatotilapia_calliptera.fAstCal1.3.dna.toplevel.fa.fai",
        "reference FASTA index",
    )
    te_fasta = require_file(
        root / "test_data" / "MWCichlidTE-3.2.splitted.namefixed_renamed_all_all.fa",
        "TE FASTA",
    )
    require_file(
        root / "test_data" / "MWCichlidTE-3.2.splitted.namefixed_renamed_all_all.fa.fai",
        "TE FASTA index",
    )
    placer_bin = require_executable(root / "build" / "placer", "PLACER executable")

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_root = (root / "placer_out" / "release_readiness" / stamp).resolve()
    out_root.mkdir(parents=True, exist_ok=False)

    rounds: List[ReleaseRoundSummary] = []
    for round_index, seed in enumerate(generate_unique_seeds(ROUND_COUNT), start=1):
        rounds.append(
            run_round(
                round_index=round_index,
                seed=seed,
                root=root,
                out_root=out_root,
                truth_table=truth_table,
                bam=bam,
                ref_fasta=ref_fasta,
                te_fasta=te_fasta,
                placer_bin=placer_bin,
            )
        )

    aggregate = evaluate_release_gate(
        rounds,
        round_count=ROUND_COUNT,
        min_detection_rate=MIN_DETECTION_RATE,
    )
    aggregate["output_root"] = str(out_root)

    aggregate_path = out_root / "release_readiness_summary.json"
    aggregate_path.write_text(
        json.dumps(aggregate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(aggregate, indent=2, sort_keys=True))
    print(f"release_readiness_summary_json\t{aggregate_path}")
    return 0 if aggregate["all_rounds_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
