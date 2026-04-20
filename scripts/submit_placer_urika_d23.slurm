#!/usr/bin/env bash
#SBATCH -J placer_batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=3-00:00:00
#SBATCH --output=placer_batch_%j.out
#SBATCH --error=placer_batch_%j.err

set -euo pipefail

# Override cluster-specific settings such as partition/QOS at submit time, e.g.
# PLACER_BAM=/data/input.bam PLACER_REF=/data/ref.fa PLACER_TE=/data/te.fa \
#   PLACER_OUT_ROOT=/data/placer_out PLACER_RUN_NAME=yohann_d23 \
#   sbatch --cpus-per-task=64 --mem=192G --export=ALL,PLACER_BAM,PLACER_REF,PLACER_TE,PLACER_OUT_ROOT,PLACER_RUN_NAME \
#   scripts/submit_placer_urika_d23.slurm

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

require_env() {
    local name="$1"
    local value="${!name:-}"
    if [[ -z "$value" ]]; then
        echo "[slurm] missing required environment variable: $name" >&2
        exit 2
    fi
    printf '%s\n' "$value"
}

BAM="$(require_env PLACER_BAM)"
REF="$(require_env PLACER_REF)"
TE="$(require_env PLACER_TE)"
BAM_BAI="${PLACER_BAM_BAI:-${BAM}.bai}"
REF_FAI="${PLACER_REF_FAI:-${REF}.fai}"
TE_FAI="${PLACER_TE_FAI:-${TE}.fai}"

OUT_ROOT="${PLACER_OUT_ROOT:-$REPO_ROOT/placer_out/slurm_runs}"
RUN_NAME="${PLACER_RUN_NAME:-placer_job${SLURM_JOB_ID:-manual}}"
RUN_DIR="$OUT_ROOT/$RUN_NAME"

mkdir -p "$RUN_DIR"

exec >>"$RUN_DIR/job.stdout.log"
exec 2>>"$RUN_DIR/job.stderr.log"

echo "[slurm] start_time=$(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "[slurm] host=$(hostname)"
echo "[slurm] repo_root=$REPO_ROOT"
echo "[slurm] run_dir=$RUN_DIR"
echo "[slurm] slurm_job_id=${SLURM_JOB_ID:-NA}"
echo "[slurm] slurm_cpus_per_task=${SLURM_CPUS_PER_TASK:-NA}"

for path in \
    "$BAM" \
    "$BAM_BAI" \
    "$REF" \
    "$REF_FAI" \
    "$TE" \
    "$TE_FAI" \
    "$REPO_ROOT/scripts/run_sharded_placer.py" \
    "$REPO_ROOT/scripts/build_latest_placer.sh"; do
    if [[ ! -e "$path" ]]; then
        echo "[slurm] missing required path: $path" >&2
        exit 2
    fi
done

TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-1}"
PLACER_SHARD_WORKERS="${PLACER_SHARD_WORKERS:-$TOTAL_CPUS}"
PLACER_SAMTOOLS_THREADS="${PLACER_SAMTOOLS_THREADS:-1}"
PLACER_SHARD_HEARTBEAT_S="${PLACER_SHARD_HEARTBEAT_S:-60}"

export PLACER_BUILD_JOBS="${PLACER_BUILD_JOBS:-$TOTAL_CPUS}"
export PLACER_PROGRESS_INTERVAL="${PLACER_PROGRESS_INTERVAL:-5000}"
export PLACER_LOG_STAGE_BINS="${PLACER_LOG_STAGE_BINS:-1}"
export PLACER_LOG_STAGE_COMPONENTS="${PLACER_LOG_STAGE_COMPONENTS:-1}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "[slurm] bam=$BAM"
echo "[slurm] ref=$REF"
echo "[slurm] te=$TE"
echo "[slurm] PLACER_BUILD_JOBS=$PLACER_BUILD_JOBS"
echo "[slurm] PLACER_PARALLEL=disabled"
echo "[slurm] PLACER_SHARD_WORKERS=$PLACER_SHARD_WORKERS"
echo "[slurm] PLACER_SAMTOOLS_THREADS=$PLACER_SAMTOOLS_THREADS"
echo "[slurm] PLACER_SHARD_HEARTBEAT_S=$PLACER_SHARD_HEARTBEAT_S"
echo "[slurm] PLACER_PROGRESS_INTERVAL=$PLACER_PROGRESS_INTERVAL"
echo "[slurm] PLACER_LOG_STAGE_BINS=$PLACER_LOG_STAGE_BINS"
echo "[slurm] PLACER_LOG_STAGE_COMPONENTS=$PLACER_LOG_STAGE_COMPONENTS"

if command -v git >/dev/null 2>&1 && git -C "$REPO_ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "[slurm] git_commit=$(git -C "$REPO_ROOT" rev-parse --short HEAD)"
fi

cd "$RUN_DIR"
python3 "$REPO_ROOT/scripts/run_sharded_placer.py" \
    --bam "$BAM" \
    --ref "$REF" \
    --te "$TE" \
    --placer "$REPO_ROOT/build/placer" \
    --mode contig \
    --workers "$PLACER_SHARD_WORKERS" \
    --samtools-threads "$PLACER_SAMTOOLS_THREADS" \
    --outdir "$RUN_DIR" \
    --progress-heartbeat-s "$PLACER_SHARD_HEARTBEAT_S" \
    --resume \
    --env "PLACER_PARALLEL=0" \
    --env "PLACER_BAM_THREADS=1" \
    --env "PLACER_PROGRESS_INTERVAL=$PLACER_PROGRESS_INTERVAL" \
    --env "PLACER_LOG_STAGE_BINS=$PLACER_LOG_STAGE_BINS" \
    --env "PLACER_LOG_STAGE_COMPONENTS=$PLACER_LOG_STAGE_COMPONENTS"

echo "[slurm] finished_time=$(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "[slurm] scientific_sharded_txt=$RUN_DIR/scientific.sharded.txt"
echo "[slurm] shard_manifest_tsv=$RUN_DIR/shard_manifest.tsv"
