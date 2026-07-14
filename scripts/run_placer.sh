#!/usr/bin/env bash
#SBATCH -J placer_te
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

require_env() {
    local name="$1"
    local value="${!name:-}"
    if [[ -z "$value" ]]; then
        echo "[run_placer] missing required environment variable: $name" >&2
        exit 2
    fi
    printf '%s\n' "$value"
}

OUTPUT_DIR="$(require_env PLACER_OUTPUT_DIR)"
BAM="$(require_env PLACER_BAM)"
REF="$(require_env PLACER_REF)"
TE="$(require_env PLACER_TE)"
PLACER_BIN="${PLACER_BIN:-$OUTPUT_DIR/placer_bin}"

if [[ ! -x "$PLACER_BIN" ]]; then
    echo "[run_placer] PLACER binary is not executable: $PLACER_BIN" >&2
    exit 2
fi

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

SCRATCH_DIR="${SLURM_TMPDIR:-/tmp}/placer_inputs_$$"
mkdir -p "$SCRATCH_DIR"

cp "$BAM" "$SCRATCH_DIR/"
cp "${PLACER_BAM_BAI:-$BAM.bai}" "$SCRATCH_DIR/"
cp "$REF" "$SCRATCH_DIR/"
cp "${PLACER_REF_FAI:-$REF.fai}" "$SCRATCH_DIR/"
cp "$TE" "$SCRATCH_DIR/"
cp "${PLACER_TE_FAI:-$TE.fai}" "$SCRATCH_DIR/"

export PLACER_PARALLEL=1
export PLACER_PARALLEL_WORKERS="${PLACER_PARALLEL_WORKERS:-96}"
export PLACER_PARALLEL_QUEUE_MAX_TASKS="${PLACER_PARALLEL_QUEUE_MAX_TASKS:-384}"
export PLACER_PARALLEL_RESULT_BUFFER_MAX="${PLACER_PARALLEL_RESULT_BUFFER_MAX:-384}"
export PLACER_BAM_THREADS="${PLACER_BAM_THREADS:-2}"
export PLACER_PROGRESS_INTERVAL="${PLACER_PROGRESS_INTERVAL:-100000}"
export PLACER_LOG_PARALLEL_PROGRESS="${PLACER_LOG_PARALLEL_PROGRESS:-1}"
export PLACER_LOG_STAGE_BINS="${PLACER_LOG_STAGE_BINS:-0}"
export PLACER_LOG_STAGE_COMPONENTS="${PLACER_LOG_STAGE_COMPONENTS:-0}"

"$PLACER_BIN" \
    "$SCRATCH_DIR/$(basename "$BAM")" \
    "$SCRATCH_DIR/$(basename "$REF")" \
    "$SCRATCH_DIR/$(basename "$TE")"

echo "Done! Results saved to $OUTPUT_DIR/"
