#!/usr/bin/env python3
"""Merge completed PLACER shard directories into one scientific output."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from run_sharded_placer import merge_existing_shard_directory


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Merge existing PLACER shard scientific.txt files."
    )
    parser.add_argument(
        "shard_root",
        help="Directory containing completed shard work directories.",
    )
    parser.add_argument(
        "--output",
        default="",
        help="Merged scientific output path. Defaults to <shard_root.parent>/scientific.sharded.txt.",
    )
    parser.add_argument("--dedup-bp", type=int, default=50)
    parser.add_argument(
        "--delete-shards",
        action="store_true",
        help="Delete source shard directories after the merged output is written successfully.",
    )
    args = parser.parse_args()

    shard_root = Path(args.shard_root).resolve()
    output = Path(args.output).resolve() if args.output else (
        shard_root.parent / "scientific.sharded.txt"
    )

    try:
        results = merge_existing_shard_directory(
            shard_root,
            out_path=output,
            dedup_bp=max(0, args.dedup_bp),
            delete_source_shards=args.delete_shards,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"[finalize] failed: {exc}", file=sys.stderr)
        return 2

    print(
        f"[finalize] merged {len(results)} shards into {output}",
        flush=True,
    )
    if args.delete_shards:
        print(f"[finalize] deleted source shard directories under {shard_root}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
