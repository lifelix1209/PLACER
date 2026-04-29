#!/usr/bin/env bash
# Yohann D23 production submission
# Usage:
#   bash submit_cichlid_yohann_d23.sh <bam> <ref> <te> <out_root> [run_name]

set -euo pipefail

THIS_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export REPO_ROOT="${REPO_ROOT:-$(cd "$THIS_SCRIPT_DIR/.." && pwd)}"
export SCRIPT_DIR="${REPO_ROOT}/scripts"
SLURM_SCRIPT="${PLACER_SLURM_SCRIPT:-${SCRIPT_DIR}/submit_placer_urika_d23.slurm}"

usage() {
    cat >&2 <<'EOF'
Usage:
  bash submit_cichlid_yohann_d23.sh <bam> <ref> <te> <out_root> [run_name]

Or provide the same values via environment variables:
  PLACER_BAM
  PLACER_REF
  PLACER_TE
  PLACER_OUT_ROOT
  PLACER_RUN_NAME (optional)
EOF
}

if [[ $# -gt 0 && $# -ne 4 && $# -ne 5 ]]; then
    usage
    exit 2
fi

if [[ $# -ge 4 ]]; then
    export PLACER_BAM="$1"
    export PLACER_REF="$2"
    export PLACER_TE="$3"
    export PLACER_OUT_ROOT="$4"
fi

if [[ $# -eq 5 ]]; then
    export PLACER_RUN_NAME="$5"
fi

if [[ -z "${PLACER_BAM:-}" || -z "${PLACER_REF:-}" || -z "${PLACER_TE:-}" || -z "${PLACER_OUT_ROOT:-}" ]]; then
    usage
    exit 2
fi

export PLACER_RUN_NAME="${PLACER_RUN_NAME:-yohann_d23}"
export PLACER_SLURM_PARTITION="${PLACER_SLURM_PARTITION:-2204}"
export PLACER_SLURM_CPUS_PER_TASK="${PLACER_SLURM_CPUS_PER_TASK:-96}"
export PLACER_SLURM_MEM="${PLACER_SLURM_MEM:-320G}"
export PLACER_SLURM_TIME="${PLACER_SLURM_TIME:-7-00:00:00}"

default_shard_workers=48
if [[ "$PLACER_SLURM_CPUS_PER_TASK" =~ ^[0-9]+$ && "$PLACER_SLURM_CPUS_PER_TASK" -lt "$default_shard_workers" ]]; then
    default_shard_workers="$PLACER_SLURM_CPUS_PER_TASK"
fi

# Production path: outer exact contig sharding only.
# Keep internal PLACER_PARALLEL disabled until the real-data regression is resolved.
export PLACER_SHARD_WORKERS="${PLACER_SHARD_WORKERS:-$default_shard_workers}"
export PLACER_SAMTOOLS_THREADS="${PLACER_SAMTOOLS_THREADS:-1}"
export PLACER_SHARD_HEARTBEAT_S="${PLACER_SHARD_HEARTBEAT_S:-60}"
export PLACER_BAM_THREADS="${PLACER_BAM_THREADS:-1}"
export PLACER_PROGRESS_INTERVAL="${PLACER_PROGRESS_INTERVAL:-5000}"
export PLACER_LOG_STAGE_BINS="${PLACER_LOG_STAGE_BINS:-1}"
export PLACER_LOG_STAGE_COMPONENTS="${PLACER_LOG_STAGE_COMPONENTS:-1}"

mkdir -p "$PLACER_OUT_ROOT"

echo "[submit] submitting Yohann D23 job..."
echo "[submit] REPO_ROOT: $REPO_ROOT"
echo "[submit] SLURM_SCRIPT: $SLURM_SCRIPT"
echo "[submit] BAM: $PLACER_BAM"
echo "[submit] REF: $PLACER_REF"
echo "[submit] TE: $PLACER_TE"
echo "[submit] OUT: $PLACER_OUT_ROOT"
echo "[submit] SHARD_WORKERS: $PLACER_SHARD_WORKERS"
echo "[submit] SAMTOOLS_THREADS: $PLACER_SAMTOOLS_THREADS"
echo "[submit] PLACER_BAM_THREADS: $PLACER_BAM_THREADS"
echo "[submit] HEARTBEAT_S: $PLACER_SHARD_HEARTBEAT_S"
echo "[submit] PARTITION: $PLACER_SLURM_PARTITION"
echo "[submit] CPUS_PER_TASK: $PLACER_SLURM_CPUS_PER_TASK"
echo "[submit] MEM: $PLACER_SLURM_MEM"
echo "[submit] TIME: $PLACER_SLURM_TIME"

sbatch \
    --partition="$PLACER_SLURM_PARTITION" \
    --cpus-per-task="$PLACER_SLURM_CPUS_PER_TASK" \
    --mem="$PLACER_SLURM_MEM" \
    --time="$PLACER_SLURM_TIME" \
    --export=ALL \
    "$SLURM_SCRIPT"

echo "[submit] job submitted successfully"
