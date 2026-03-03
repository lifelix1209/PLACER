#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT_DIR/build}"
PLACER_BIN="${PLACER_BIN:-$BUILD_DIR/placer}"

DATA_DIR="${DATA_DIR:-$ROOT_DIR/test_data/sim_te_benchmark}"
BAM="${BAM:-$DATA_DIR/sim_te_benchmark.bam}"
REF="${REF:-$DATA_DIR/ref.fa}"
TE="${TE:-$DATA_DIR/te_tldr.fa}"
TRUTH_EVENTS="${TRUTH_EVENTS:-$DATA_DIR/truth_events.tsv}"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/benchmark_stress_out}"
MATCH_WINDOW_BP="${MATCH_WINDOW_BP:-200}"

mkdir -p "$OUT_DIR"

if [[ ! -x "$PLACER_BIN" ]]; then
    echo "[stress] placer binary not found: $PLACER_BIN" >&2
    echo "[stress] build first: cmake --build build -j" >&2
    exit 1
fi
if [[ ! -f "$TRUTH_EVENTS" ]]; then
    echo "[stress] truth events file not found: $TRUTH_EVENTS" >&2
    exit 1
fi

echo "[stress] running PLACER benchmark callset generation..."
(
    cd "$ROOT_DIR"
    "$PLACER_BIN" "$BAM" "$REF" "$TE" >"$OUT_DIR/placer.stdout" 2>"$OUT_DIR/placer.stderr"
    cp scientific.txt "$OUT_DIR/scientific.txt"
)

CALLS_TSV="$OUT_DIR/calls.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} $1 ~ /^chr/ {print $1, $3}' "$OUT_DIR/scientific.txt" > "$CALLS_TSV"

subset_truth() {
    local name="$1"
    local awk_expr="$2"
    local out="$OUT_DIR/truth_${name}.tsv"
    awk -F'\t' -v OFS='\t' "$awk_expr" "$TRUTH_EVENTS" > "$out"
    echo "$out"
}

score_positive_subset() {
    local calls="$1"
    local truth="$2"
    local window_bp="$3"
    awk -v w="$window_bp" 'BEGIN{FS=OFS="\t"}
         FNR==NR{n++; c[n]=$1; p[n]=$2+0; next}
         {hit=0; for(i=1;i<=n;i++){if(c[i]==$1 && (p[i]-$2<=w) && ($2-p[i]<=w)){hit=1; break}}
          if(hit)h++; else m++; t++}
         END{print t+0, h+0, m+0}' "$calls" "$truth"
}

score_negative_subset() {
    local calls="$1"
    local truth="$2"
    local window_bp="$3"
    awk -v w="$window_bp" 'BEGIN{FS=OFS="\t"}
         FNR==NR{n++; c[n]=$1; p[n]=$2+0; next}
         {fp=0; for(i=1;i<=n;i++){if(c[i]==$1 && (p[i]-$2<=w) && ($2-p[i]<=w)){fp=1; break}}
          if(fp)f++; else clean++; t++}
         END{print t+0, f+0, clean+0}' "$calls" "$truth"
}

POS_SPLIT="$(subset_truth "pos_split_indel" '
    NR>1 && $5=="POS" && (tolower($8) ~ /split/ || $1 ~ /SPLIT/) {print $2, $3, $1}
')"
POS_SHORT="$(subset_truth "pos_short_ins" '
    NR>1 && $5=="POS" && (tolower($8) ~ /short|300bp/ || $1 ~ /SHORT_INS/) {print $2, $3, $1}
')"
POS_LOW_SUPPORT="$(subset_truth "pos_low_support_or_low_kmer" '
    NR>1 && $5=="POS" && (tolower($8) ~ /low_support|low_kmer/ || $1 ~ /LOW_SUPPORT|LOW_KMER/) {print $2, $3, $1}
')"
NEG_LOW_COMPLEX="$(subset_truth "neg_low_complexity" '
    NR>1 && $5=="NEG" && (tolower($8) ~ /polya|low_complex/ || $1 ~ /LOW_COMPLEX|POLYA/) {print $2, $3, $1}
')"

RESULTS="$OUT_DIR/stress_subset_metrics.tsv"
echo -e "subset\ttruth_total\thit_or_fp\tmiss_or_clean\trecall_or_fp_rate\twindow_bp" > "$RESULTS"

append_pos_metric() {
    local subset="$1"
    local truth_file="$2"
    read -r total hit miss < <(score_positive_subset "$CALLS_TSV" "$truth_file" "$MATCH_WINDOW_BP")
    local recall="0.0000"
    if [[ "$total" -gt 0 ]]; then
        recall="$(awk -v h="$hit" -v t="$total" 'BEGIN{printf "%.4f", h/t}')"
    fi
    echo -e "${subset}\t${total}\t${hit}\t${miss}\t${recall}\t${MATCH_WINDOW_BP}" >> "$RESULTS"
}

append_neg_metric() {
    local subset="$1"
    local truth_file="$2"
    read -r total fp clean < <(score_negative_subset "$CALLS_TSV" "$truth_file" "$MATCH_WINDOW_BP")
    local fp_rate="0.0000"
    if [[ "$total" -gt 0 ]]; then
        fp_rate="$(awk -v f="$fp" -v t="$total" 'BEGIN{printf "%.4f", f/t}')"
    fi
    echo -e "${subset}\t${total}\t${fp}\t${clean}\t${fp_rate}\t${MATCH_WINDOW_BP}" >> "$RESULTS"
}

append_pos_metric "split_indel_supported_pos" "$POS_SPLIT"
append_pos_metric "short_insertion_pos" "$POS_SHORT"
append_pos_metric "low_support_or_low_kmer_pos" "$POS_LOW_SUPPORT"
append_neg_metric "low_complexity_neg" "$NEG_LOW_COMPLEX"

echo "[stress] done"
echo "[stress] scientific: $OUT_DIR/scientific.txt"
echo "[stress] calls:      $CALLS_TSV"
echo "[stress] metrics:    $RESULTS"
