#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT_DIR/build}"
PLACER_BIN="$BUILD_DIR/placer"

BAM="$ROOT_DIR/test_data/sim_te_benchmark/sim_te_benchmark.bam"
REF="$ROOT_DIR/test_data/sim_te_benchmark/ref.fa"
TE="$ROOT_DIR/test_data/sim_te_benchmark/te_tldr.fa"
TRUTH_STRONG="$ROOT_DIR/test_data/sim_te_benchmark/truth_positive_strong.tsv"
TRUTH_ALL="$ROOT_DIR/test_data/sim_te_benchmark/truth_positive_all.tsv"
TRUTH_NEG="$ROOT_DIR/test_data/sim_te_benchmark/truth_negative_fp.tsv"

if [[ ! -x "$PLACER_BIN" ]]; then
    echo "[tune] placer binary not found: $PLACER_BIN" >&2
    echo "[tune] build first: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j" >&2
    exit 1
fi

TMP_DIR="${TMPDIR:-/tmp}/placer_tune_$$"
mkdir -p "$TMP_DIR"
trap 'rm -rf "$TMP_DIR"' EXIT

declare -a POA_MIN_READS_VALUES=(2 3 4 5)
declare -a EVIDENCE_ALPHA_VALUES=(0.08 0.10 0.12)
declare -a BREAKPOINT_MAD_MAX_VALUES=(80 120 160)

RESULTS_TSV="$ROOT_DIR/tuning_results.tsv"
echo -e "poa_min_reads\tevidence_alpha\tbreakpoint_mad_max\tcalls\tstrong_hits\tstrong_miss\tall_hits\tall_miss\tneg_fp_hits\tneg_clean" > "$RESULTS_TSV"

best_poa=0
best_alpha=""
best_mad=0
best_calls=-1
best_strong_hits=-1
best_all_hits=-1
best_neg_fp=999999

total_runs=$(( ${#POA_MIN_READS_VALUES[@]} * ${#EVIDENCE_ALPHA_VALUES[@]} * ${#BREAKPOINT_MAD_MAX_VALUES[@]} ))
run_id=0

score_positive_truth() {
    local calls_file="$1"
    local truth_file="$2"
    awk 'BEGIN{FS="\t"} \
         FNR==NR{n++; c[n]=$1; p[n]=$2+0; next} \
         FNR==1{next} \
         {hit=0; qchr=$1; qpos=$2+0; for(i=1;i<=n;i++){if(c[i]==qchr && (p[i]-qpos<=200) && (qpos-p[i]<=200)){hit=1; break}} if(hit)h++; else m++;} \
         END{print h+0, m+0}' \
         "$calls_file" "$truth_file"
}

score_negative_truth() {
    local calls_file="$1"
    local truth_file="$2"
    awk 'BEGIN{FS="\t"} \
         FNR==NR{n++; c[n]=$1; p[n]=$2+0; next} \
         FNR==1{next} \
         {fp=0; qchr=$1; qpos=$2+0; for(i=1;i<=n;i++){if(c[i]==qchr && (p[i]-qpos<=200) && (qpos-p[i]<=200)){fp=1; break}} if(fp)fph++; else clean++;} \
         END{print fph+0, clean+0}' \
         "$calls_file" "$truth_file"
}

for poa_min_reads in "${POA_MIN_READS_VALUES[@]}"; do
    for evidence_alpha in "${EVIDENCE_ALPHA_VALUES[@]}"; do
        for breakpoint_mad_max in "${BREAKPOINT_MAD_MAX_VALUES[@]}"; do
            run_id=$((run_id + 1))
            echo "[tune] run ${run_id}/${total_runs}: poa_min_reads=${poa_min_reads}, evidence_alpha=${evidence_alpha}, breakpoint_mad_max=${breakpoint_mad_max}"

            (
                cd "$ROOT_DIR"
                PLACER_ASSEMBLY_POA_MIN_READS="$poa_min_reads" \
                PLACER_EVIDENCE_MIN_SUPPORT_ALPHA="$evidence_alpha" \
                PLACER_EVIDENCE_BREAKPOINT_MAD_MAX="$breakpoint_mad_max" \
                "$PLACER_BIN" "$BAM" "$REF" "$TE" >/dev/null 2>&1
            )

            calls_file="$TMP_DIR/calls.tsv"
            awk 'BEGIN{FS="\t"; OFS="\t"} $1 ~ /^chr/ {te=$6; sub(/:.*/, "", te); print $1,$3,te}' \
                "$ROOT_DIR/scientific.txt" > "$calls_file"

            calls_count=$(wc -l < "$calls_file" | tr -d ' ')
            read -r strong_hits strong_miss < <(score_positive_truth "$calls_file" "$TRUTH_STRONG")
            read -r all_hits all_miss < <(score_positive_truth "$calls_file" "$TRUTH_ALL")
            read -r neg_fp_hits neg_clean < <(score_negative_truth "$calls_file" "$TRUTH_NEG")

            echo -e "${poa_min_reads}\t${evidence_alpha}\t${breakpoint_mad_max}\t${calls_count}\t${strong_hits}\t${strong_miss}\t${all_hits}\t${all_miss}\t${neg_fp_hits}\t${neg_clean}" \
                >> "$RESULTS_TSV"

            better=0
            if (( all_hits > best_all_hits )); then
                better=1
            elif (( all_hits == best_all_hits )); then
                if (( neg_fp_hits < best_neg_fp )); then
                    better=1
                elif (( neg_fp_hits == best_neg_fp )); then
                    if (( strong_hits > best_strong_hits )); then
                        better=1
                    elif (( strong_hits == best_strong_hits && calls_count > best_calls )); then
                        better=1
                    fi
                fi
            fi

            if (( better == 1 )); then
                best_poa=$poa_min_reads
                best_alpha=$evidence_alpha
                best_mad=$breakpoint_mad_max
                best_calls=$calls_count
                best_strong_hits=$strong_hits
                best_all_hits=$all_hits
                best_neg_fp=$neg_fp_hits
            fi
        done
    done
done

echo "[tune] done"
echo "[tune] best: poa_min_reads=${best_poa}, evidence_alpha=${best_alpha}, breakpoint_mad_max=${best_mad}"
echo "[tune] best metrics: calls=${best_calls}, strong_hits=${best_strong_hits}, all_hits=${best_all_hits}, neg_fp_hits=${best_neg_fp}"
echo "[tune] results: $RESULTS_TSV"
