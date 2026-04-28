#!/usr/bin/env bash
#SBATCH -J placer_batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=320G
#SBATCH --time=7-00:00:00
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
    "$REPO_ROOT/scripts/build_latest_placer.sh"; do
    if [[ ! -e "$path" ]]; then
        echo "[slurm] missing required path: $path" >&2
        exit 2
    fi
done

TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-1}"
DEFAULT_TASKGRAPH_WORKERS=48
if [[ "$TOTAL_CPUS" =~ ^[0-9]+$ && "$TOTAL_CPUS" -lt "$DEFAULT_TASKGRAPH_WORKERS" ]]; then
    DEFAULT_TASKGRAPH_WORKERS="$TOTAL_CPUS"
fi
PLACER_TASKGRAPH_WORKERS="${PLACER_TASKGRAPH_WORKERS:-$DEFAULT_TASKGRAPH_WORKERS}"
PLACER_TASKGRAPH_QUEUE_MAX_TASKS="${PLACER_TASKGRAPH_QUEUE_MAX_TASKS:-$((PLACER_TASKGRAPH_WORKERS * 4))}"
PLACER_BAM_THREADS="${PLACER_BAM_THREADS:-2}"

export PLACER_BUILD_JOBS="${PLACER_BUILD_JOBS:-$TOTAL_CPUS}"
export PLACER_PROGRESS_INTERVAL="${PLACER_PROGRESS_INTERVAL:-5000}"
export PLACER_LOG_STAGE_BINS="${PLACER_LOG_STAGE_BINS:-1}"
export PLACER_LOG_STAGE_COMPONENTS="${PLACER_LOG_STAGE_COMPONENTS:-1}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

if [[ -n "${HTSLIB_ROOT:-}" ]]; then
    export PATH="$HTSLIB_ROOT/bin:$PATH"
    export LD_LIBRARY_PATH="$HTSLIB_ROOT/lib:${LD_LIBRARY_PATH:-}"
elif [[ -n "${CONDA_PREFIX:-}" ]]; then
    export PATH="$CONDA_PREFIX/bin:$PATH"
    export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
fi

echo "[slurm] bam=$BAM"
echo "[slurm] ref=$REF"
echo "[slurm] te=$TE"
echo "[slurm] htslib_root=${HTSLIB_ROOT:-NA}"
echo "[slurm] htslib_include_dir=${HTSLIB_INCLUDE_DIR:-NA}"
echo "[slurm] htslib_library=${HTSLIB_LIBRARY:-NA}"
echo "[slurm] PLACER_BUILD_JOBS=$PLACER_BUILD_JOBS"
echo "[slurm] PLACER_TASKGRAPH_WORKERS=$PLACER_TASKGRAPH_WORKERS"
echo "[slurm] PLACER_TASKGRAPH_QUEUE_MAX_TASKS=$PLACER_TASKGRAPH_QUEUE_MAX_TASKS"
echo "[slurm] PLACER_BAM_THREADS=$PLACER_BAM_THREADS"
echo "[slurm] PLACER_PROGRESS_INTERVAL=$PLACER_PROGRESS_INTERVAL"
echo "[slurm] PLACER_LOG_STAGE_BINS=$PLACER_LOG_STAGE_BINS"
echo "[slurm] PLACER_LOG_STAGE_COMPONENTS=$PLACER_LOG_STAGE_COMPONENTS"

if command -v git >/dev/null 2>&1 && git -C "$REPO_ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "[slurm] git_commit=$(git -C "$REPO_ROOT" rev-parse --short HEAD)"
fi

cd "$RUN_DIR"
PLACER_PARALLEL=1 \
PLACER_PARALLEL_WORKERS="$PLACER_TASKGRAPH_WORKERS" \
PLACER_PARALLEL_QUEUE_MAX_TASKS="$PLACER_TASKGRAPH_QUEUE_MAX_TASKS" \
PLACER_BAM_THREADS="$PLACER_BAM_THREADS" \
PLACER_PROGRESS_INTERVAL="$PLACER_PROGRESS_INTERVAL" \
PLACER_LOG_STAGE_BINS="$PLACER_LOG_STAGE_BINS" \
PLACER_LOG_STAGE_COMPONENTS="$PLACER_LOG_STAGE_COMPONENTS" \
PLACER_LOG_PARALLEL_PROGRESS=1 \
"$REPO_ROOT/build/placer" "$BAM" "$REF" "$TE"

echo "[slurm] finished_time=$(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "[slurm] scientific_txt=$RUN_DIR/scientific.txt"
