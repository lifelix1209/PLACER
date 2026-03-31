#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_HELPER="$REPO_ROOT/scripts/build_latest_placer.sh"

usage() {
    cat <<'EOF'
Usage:
  scripts/run_placer_latest.sh <input.bam> <ref.fa> <te.fa>

Description:
  1) Check whether build/placer is missing or older than source files.
  2) Rebuild automatically when needed.
  3) Run PLACER from the current working directory so job outputs stay local.

Defaults for debug visibility:
  PLACER_PROGRESS_INTERVAL=10000
  PLACER_LOG_STAGE_BINS=1
  PLACER_LOG_STAGE_COMPONENTS=1

Override any of them in the environment before launching if needed.
EOF
}

if [[ $# -lt 3 ]]; then
    usage >&2
    exit 1
fi

if [[ ! -x "$BUILD_HELPER" ]]; then
    echo "[run-latest] build helper not executable: $BUILD_HELPER" >&2
    exit 2
fi

export PLACER_PROGRESS_INTERVAL="${PLACER_PROGRESS_INTERVAL:-10000}"
export PLACER_LOG_STAGE_BINS="${PLACER_LOG_STAGE_BINS:-1}"
export PLACER_LOG_STAGE_COMPONENTS="${PLACER_LOG_STAGE_COMPONENTS:-1}"

PLACER_BIN="$("$BUILD_HELPER")"

echo "[run-latest] binary=$PLACER_BIN" >&2
echo "[run-latest] cwd=$(pwd)" >&2
echo "[run-latest] progress_interval=$PLACER_PROGRESS_INTERVAL log_stage_bins=$PLACER_LOG_STAGE_BINS log_stage_components=$PLACER_LOG_STAGE_COMPONENTS" >&2

exec "$PLACER_BIN" "$@"
