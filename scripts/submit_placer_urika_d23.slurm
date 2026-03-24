#!/usr/bin/env bash
#SBATCH -J placer_urika_d23
#SBATCH --partition=2204
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=3-00:00:00
#SBATCH --output=placer_urika_d23_%j.out
#SBATCH --error=placer_urika_d23_%j.err

set -euo pipefail

# Default sizing targets partition 2204. Override at submit time if needed, e.g.
# sbatch --partition=2004 --cpus-per-task=48 --mem=256G scripts/submit_placer_urika_d23.slurm

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

BAM="/mnt/home1/miska/hl725/scratch/projects/ciclid/urika/data/final_merged_Yohann_D23.bam"
BAM_BAI="/mnt/home1/miska/hl725/scratch/projects/ciclid/urika/data/final_merged_Yohann_D23.bam.bai"
REF="/mnt/home1/miska/hl725/scratch/projects/ciclid/urika/references/Astatotilapia_calliptera.fAstCal1.3.dna.toplevel.fa"
REF_FAI="/mnt/home1/miska/hl725/scratch/projects/ciclid/urika/references/Astatotilapia_calliptera.fAstCal1.3.dna.toplevel.fa.fai"
TE="/mnt/home1/miska/hl725/scratch/projects/ciclid/urika/references/MWCichlidTE-3.2.splitted.namefixed_renamed_all_all.fa"
TE_FAI="/mnt/home1/miska/hl725/scratch/projects/ciclid/urika/references/MWCichlidTE-3.2.splitted.namefixed_renamed_all_all.fa.fai"

OUT_ROOT="${PLACER_OUT_ROOT:-/mnt/home1/miska/hl725/scratch/projects/ciclid/urika/placer_runs}"
RUN_NAME="final_merged_Yohann_D23_job${SLURM_JOB_ID:-manual}"
RUN_DIR="$OUT_ROOT/$RUN_NAME"

mkdir -p "$RUN_DIR"

exec > >(tee -a "$RUN_DIR/job.stdout.log")
exec 2> >(tee -a "$RUN_DIR/job.stderr.log" >&2)

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
    "$REPO_ROOT/scripts/run_placer_latest.sh" \
    "$REPO_ROOT/scripts/build_latest_placer.sh"; do
    if [[ ! -e "$path" ]]; then
        echo "[slurm] missing required path: $path" >&2
        exit 2
    fi
done

TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-1}"
if (( TOTAL_CPUS >= 8 )); then
    export PLACER_BAM_THREADS=4
    export PLACER_PARALLEL_WORKERS=$((TOTAL_CPUS - 4))
elif (( TOTAL_CPUS >= 4 )); then
    export PLACER_BAM_THREADS=2
    export PLACER_PARALLEL_WORKERS=$((TOTAL_CPUS - 2))
else
    export PLACER_BAM_THREADS=1
    export PLACER_PARALLEL_WORKERS=1
fi
if (( PLACER_PARALLEL_WORKERS < 1 )); then
    export PLACER_PARALLEL_WORKERS=1
fi

export PLACER_BUILD_JOBS="${PLACER_BUILD_JOBS:-$TOTAL_CPUS}"
export PLACER_PARALLEL="${PLACER_PARALLEL:-1}"
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
echo "[slurm] PLACER_PARALLEL=$PLACER_PARALLEL"
echo "[slurm] PLACER_BAM_THREADS=$PLACER_BAM_THREADS"
echo "[slurm] PLACER_PARALLEL_WORKERS=$PLACER_PARALLEL_WORKERS"
echo "[slurm] PLACER_PROGRESS_INTERVAL=$PLACER_PROGRESS_INTERVAL"
echo "[slurm] PLACER_LOG_STAGE_BINS=$PLACER_LOG_STAGE_BINS"
echo "[slurm] PLACER_LOG_STAGE_COMPONENTS=$PLACER_LOG_STAGE_COMPONENTS"

if command -v git >/dev/null 2>&1 && git -C "$REPO_ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "[slurm] git_commit=$(git -C "$REPO_ROOT" rev-parse --short HEAD)"
fi

cd "$RUN_DIR"
"$REPO_ROOT/scripts/run_placer_latest.sh" "$BAM" "$REF" "$TE"

echo "[slurm] finished_time=$(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "[slurm] scientific_txt=$RUN_DIR/scientific.txt"
