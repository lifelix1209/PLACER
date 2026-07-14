#!/usr/bin/env bash
#SBATCH -J placer_sim_truth
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Build a simulation-based TE-insertion truth set and benchmark PLACER against it,
# end to end, on a compute node.
#
#   1. plant known TE insertions into the reference        (simulate_te_insertions.py)
#   2. simulate ONT reads from the modified genome         (badread, one job per contig)
#   3. align reads back to the UNTOUCHED reference         (minimap2 map-ont -> sorted BAM)
#   4. call insertions                                     (PLACER binary)
#   5. score calls vs the planted truth                    (benchmark_sim_truth.py)
#
# The truth is exact (breakpoint, family, strand, length), so this measures real
# precision AND recall -- the controlled complement to the assembly-based GraffiTE
# truth (assemble_d2_flye.sh -> graffite_truthset.sh), which is precision-strong but
# recall-incomplete at D2's short read N50.
#
# COVERAGE SWEEP: SIM_QUANTITY may be a comma-separated list of badread coverages
# (e.g. SIM_QUANTITY=10x,25x,40x). The insertions are planted once (shared truth)
# and steps 2-5 are repeated per coverage into cov_<tag>/ subdirectories, with an
# aggregated benchmark_summary.tsv across coverages. This is what stresses the
# depth-sensitive parts of the decision layer -- beta-binomial genotyping and its
# per-sample overdispersion estimate, and the worst-case-lFDR / conformal-null
# calibration -- which a single fixed coverage does not exercise.
#
# Configure with environment variables (submit like:
#   sbatch --export=ALL,PLACER_BIN=/path/placer,SIM_QUANTITY=10x,25x,40x scripts/simulate_te_truthset.slurm )
set -euo pipefail

# ---- config (override via --export / env) -----------------------------------
REPO="${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
REF="${REF:-$REPO/test_data/fAstCal671.primary.fa}"
TE="${TE:-$REPO/test_data/MWCichlidTE-3.2.splitted.namefixed_renamed_all_all.fa}"
PLACER_BIN="${PLACER_BIN:-$REPO/build/placer}"
OUT="${OUT:-$REPO/analysis/sim_truthset}"

SIM_CONTIGS="${SIM_CONTIGS:-}"          # comma list, empty = all contigs in REF
SIM_NUM="${SIM_NUM:-150}"               # insertions per contig
SIM_QUANTITY="${SIM_QUANTITY:-25x}"     # badread coverage; comma list = sweep (e.g. 10x,25x,40x)
SIM_LENGTH="${SIM_LENGTH:-7000,5000}"   # badread read length mean,stdev (D2-like)
SIM_SEED="${SIM_SEED:-7}"
CHUNK_MB="${CHUNK_MB:-5}"               # simulate reads from ~this many Mb per badread job
MATCH_RADIUS="${MATCH_RADIUS:-100}"     # breakpoint match tolerance bp
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"

BADREAD="${BADREAD:-badread}"           # badread executable (pip install git+https://github.com/rrwick/Badread)

# ---- preflight --------------------------------------------------------------
for tool in "$BADREAD" minimap2 samtools python3; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not found on PATH" >&2; exit 1; }
done
for f in "$REF" "$TE" "$PLACER_BIN"; do
  [[ -e "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done
[[ -x "$PLACER_BIN" ]] || { echo "ERROR: PLACER binary not executable: $PLACER_BIN" >&2; exit 1; }

mkdir -p "$OUT"
MOD_FA="$OUT/modified.fa"
TRUTH="$OUT/truth.tsv"
IFS=',' read -ra QUANTITIES <<< "$SIM_QUANTITY"
echo "[config] REF=$REF"
echo "[config] TE=$TE"
echo "[config] contigs=${SIM_CONTIGS:-ALL} num=$SIM_NUM coverages=${QUANTITIES[*]} length=$SIM_LENGTH threads=$THREADS"
echo "[config] out=$OUT"

# ---- 1. plant insertions (once; shared across the coverage sweep) -----------
echo "[1/5] planting insertions"
python3 "$REPO/scripts/simulate_te_insertions.py" \
  --ref "$REF" --te "$TE" \
  ${SIM_CONTIGS:+--contigs "$SIM_CONTIGS"} \
  --num "$SIM_NUM" --seed "$SIM_SEED" \
  --out-fa "$MOD_FA" --out-truth "$TRUTH"

# ---- split the modified genome into ~CHUNK_MB pieces (once) ------------------
# badread is single-threaded, so we split the modified genome into ~CHUNK_MB pieces
# and simulate them independently. Chunk record names are arbitrary: minimap2 aligns
# by sequence content, so reads still land at their true reference coordinates (the
# only cost is a handful of reads lost at chunk seams). This saturates all cores and
# avoids the largest chromosome becoming the wall-time bottleneck. The chunks depend
# only on the planted genome, so they are shared across coverages.
SPLIT="$OUT/chunks"
rm -rf "$SPLIT"; mkdir -p "$SPLIT"
python3 - "$MOD_FA" "$SPLIT" "$CHUNK_MB" <<'PY'
import sys, pathlib
mod, outdir, chunk_mb = sys.argv[1], pathlib.Path(sys.argv[2]), int(sys.argv[3])
step = chunk_mb * 1_000_000
idx = 0
def emit(seq):
    global idx
    for i in range(0, len(seq), step):
        piece = seq[i:i + step]
        with open(outdir / f"chunk{idx:05d}.fa", "w") as f:
            f.write(f">chunk{idx:05d}\n")
            for j in range(0, len(piece), 60):
                f.write(piece[j:j + 60] + "\n")
        idx += 1
buf = []
for line in open(mod):
    if line.startswith(">"):
        if buf: emit("".join(buf)); buf = []
    else:
        buf.append(line.strip())
if buf: emit("".join(buf))
print(f"  {idx} chunks")
PY

sim_one() {
  local fa="$1" seed="$2" qty="$3" len="$4" badread="$5"
  "$badread" simulate --seed "$seed" --reference "$fa" --quantity "$qty" --length "$len" \
    2>"${fa%.fa}.badread.log" | gzip -c > "${fa%.fa}.fq.gz"
}
export -f sim_one

# ---- steps 2-5 per coverage -------------------------------------------------
SUMMARY="$OUT/benchmark_summary.tsv"
printf "coverage\tprecision\trecall\tf1\ttp\tfp\tfn\temitted\n" > "$SUMMARY"

# Parse benchmark_sim_truth.py's human-readable report into a summary row.
summary_row() {
  local qty="$1" report="$2"
  awk -v qty="$qty" '
    /^precision:/            {p=$2}
    /^recall:/               {r=$2}
    /^f1:/                   {f=$2}
    /true positives/         {tp=$NF}
    /false positives:/       {fp=$NF}
    /false negatives/        {fn=$NF}
    /PLACER calls emitted:/  {em=$NF}
    END {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", qty, p, r, f, tp, fp, fn, em}
  ' "$report" >> "$SUMMARY"
}

for QTY in "${QUANTITIES[@]}"; do
  QTAG="${QTY//[^0-9a-zA-Z]/}"
  COV_OUT="$OUT/cov_${QTAG}"
  mkdir -p "$COV_OUT"
  READS="$COV_OUT/reads.fq.gz"
  BAM="$COV_OUT/sim.sorted.bam"
  RUN="$COV_OUT/placer_run"; mkdir -p "$RUN"
  echo "=================================================================="
  echo "[coverage $QTY]"

  # ---- 2. simulate reads (fixed-size chunks, run in parallel) ---------------
  echo "[2/5] simulating ONT reads with badread ($QTY, length $SIM_LENGTH, ~${CHUNK_MB}Mb chunks)"
  # distinct seed per chunk (offset by coverage index) keeps runs reproducible yet independent
  i=0
  for fa in "$SPLIT"/*.fa; do echo "$fa $((SIM_SEED + i))"; i=$((i + 1)); done \
    | xargs -P "$THREADS" -L1 bash -c 'sim_one "$0" "$1" "'"$QTY"'" "'"$SIM_LENGTH"'" "'"$BADREAD"'"'
  cat "$SPLIT"/*.fq.gz > "$READS"
  echo "  reads -> $READS"

  # ---- 3. align back to the untouched reference -----------------------------
  echo "[3/5] aligning to reference (minimap2 map-ont)"
  minimap2 -ax map-ont -t "$THREADS" "$REF" "$READS" \
    | samtools sort -@ "$THREADS" -o "$BAM" -
  samtools index -@ "$THREADS" "$BAM"

  # ---- 4. run PLACER --------------------------------------------------------
  echo "[4/5] running PLACER"
  ( cd "$RUN"
    export PLACER_PARALLEL=1
    export PLACER_PARALLEL_WORKERS="${PLACER_PARALLEL_WORKERS:-$THREADS}"
    export PLACER_BAM_THREADS="${PLACER_BAM_THREADS:-4}"
    "$PLACER_BIN" "$BAM" "$REF" "$TE" )

  # ---- 5. benchmark ---------------------------------------------------------
  echo "[5/5] benchmarking vs planted truth (radius ${MATCH_RADIUS}bp)"
  python3 "$REPO/scripts/benchmark_sim_truth.py" \
    --scientific "$RUN/scientific.txt" \
    --truth "$TRUTH" \
    --radius "$MATCH_RADIUS" \
    --out-matched "$COV_OUT/matched.tsv" | tee "$COV_OUT/benchmark.txt"
  summary_row "$QTY" "$COV_OUT/benchmark.txt"
done

echo "=================================================================="
echo "coverage sweep summary ($SUMMARY):"
column -t "$SUMMARY"
echo "DONE. Per-coverage results in $OUT/cov_<tag>/ (benchmark.txt, matched.tsv, placer_run/scientific.txt); truth in $TRUTH."
