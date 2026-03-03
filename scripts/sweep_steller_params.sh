#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PLACER_BIN="${PLACER_BIN:-$ROOT_DIR/build/placer}"

DEFAULT_DATA_ROOT="$ROOT_DIR/../sTELLeR_supplementary"
DATA_ROOT="${DATA_ROOT:-$DEFAULT_DATA_ROOT}"
BAM="${BAM:-$DATA_ROOT/testdata.bam}"
REF="${REF:-$DATA_ROOT/sTELLeR_supplementary/refgenome/HG38_chr22.fa}"
TE_FASTA="${TE_FASTA:-$DATA_ROOT/fasta/TEfastasequences.fa}"
TRUTH="${TRUTH:-$DATA_ROOT/testbamTRUTH.txt}"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/benchmark_steller_sweep}"
MATCH_WINDOW_BP="${MATCH_WINDOW_BP:-200}"
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-0}"
PROGRESS_EVERY="${PROGRESS_EVERY:-25}"

mkdir -p "$OUT_DIR"

if [[ ! -x "$PLACER_BIN" ]]; then
    echo "[steller-sweep] placer binary not found: $PLACER_BIN" >&2
    exit 1
fi
for f in "$BAM" "$REF" "$TE_FASTA" "$TRUTH"; do
    if [[ ! -f "$f" ]]; then
        echo "[steller-sweep] required file missing: $f" >&2
        exit 1
    fi
done

RESULTS="$OUT_DIR/sweep_results.tsv"
echo -e "run_id\temit_low_conf\tlow_conf_min\tlow_conf_max_tier\tevidence_alpha\tbreakpoint_mad_max\tassembly_min_len\tassembly_min_identity\tte_vote_min\tte_identity_min\ttruth_total\tcall_total\ttruth_hits_pos\ttruth_hits_family\tmatched_calls_pos\tfp_calls\tprecision_pos\trecall_pos\tf1_pos\twindow_bp" > "$RESULTS"

score_calls() {
    local calls_tsv="$1"
    local truth_tsv="$2"
    local window_bp="$3"
    awk -v w="$window_bp" 'BEGIN{
            FS=OFS="\t";
        }
        function norm_chr(c, t) {
            t = c;
            sub(/^chr/, "", t);
            sub(/^CHR/, "", t);
            return toupper(t);
        }
        function norm_family(f, u) {
            u = toupper(f);
            if (index(u, "ALU") == 1) return "ALU";
            if (index(u, "L1") == 1) return "L1";
            if (index(u, "SVA") == 1) return "SVA";
            if (index(u, "HERV") == 1) return "HERV";
            return u;
        }
        FNR==NR{
            n_truth++;
            truth_chr[n_truth] = norm_chr($1);
            truth_pos[n_truth] = $2 + 0;
            truth_family[n_truth] = norm_family($3);
            next;
        }
        {
            n_calls++;
            c_chr = norm_chr($1);
            c_pos = $2 + 0;
            c_family = norm_family($3);
            best_i = 0;
            best_d = 1e18;
            for (i = 1; i <= n_truth; i++) {
                if (truth_chr[i] != c_chr) continue;
                d = truth_pos[i] - c_pos;
                if (d < 0) d = -d;
                if (d <= w && d < best_d) {
                    best_d = d;
                    best_i = i;
                }
            }
            if (best_i > 0) {
                matched_calls++;
                truth_hit[best_i] = 1;
                if (c_family == truth_family[best_i]) {
                    truth_family_hit[best_i] = 1;
                }
            }
        }
        END{
            truth_hits = 0;
            truth_family_hits = 0;
            for (i = 1; i <= n_truth; i++) {
                if (truth_hit[i]) truth_hits++;
                if (truth_family_hit[i]) truth_family_hits++;
            }
            fp_calls = n_calls - matched_calls;
            precision = (n_calls > 0) ? matched_calls / n_calls : 0.0;
            recall = (n_truth > 0) ? truth_hits / n_truth : 0.0;
            f1 = ((precision + recall) > 0) ? (2.0 * precision * recall / (precision + recall)) : 0.0;
            printf("%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n",
                n_truth, n_calls, truth_hits, truth_family_hits, fp_calls, precision, recall, f1);
        }' "$truth_tsv" "$calls_tsv"
}

extract_calls() {
    local scientific="$1"
    local calls_out="$2"
    awk -F'\t' 'BEGIN{OFS="\t"} $1 !~ /^#/ && NF >= 40 && $2 ~ /^-?[0-9]+$/ {print $1,$3,$6}' \
        "$scientific" > "$calls_out"
}

emit_low_conf_values=(1)
low_conf_min_values=(1 2)
low_conf_tier_values=(2 3)
evidence_alpha_values=(0.005 0.01 0.02 0.04)
mad_max_values=(120 200)
assembly_min_len_values=(20 40 60)
assembly_min_identity_values=(0.30 0.40 0.50)
te_vote_min_values=(0.30 0.40)
te_identity_min_values=(0.20 0.30)

total_runs=$(( ${#emit_low_conf_values[@]} * ${#low_conf_min_values[@]} * ${#low_conf_tier_values[@]} * ${#evidence_alpha_values[@]} * ${#mad_max_values[@]} * ${#assembly_min_len_values[@]} * ${#assembly_min_identity_values[@]} * ${#te_vote_min_values[@]} * ${#te_identity_min_values[@]} ))
run_id=0

for emit_low_conf in "${emit_low_conf_values[@]}"; do
for low_conf_min in "${low_conf_min_values[@]}"; do
for low_conf_tier in "${low_conf_tier_values[@]}"; do
for evidence_alpha in "${evidence_alpha_values[@]}"; do
for mad_max in "${mad_max_values[@]}"; do
for asm_min_len in "${assembly_min_len_values[@]}"; do
for asm_min_ident in "${assembly_min_identity_values[@]}"; do
for te_vote_min in "${te_vote_min_values[@]}"; do
for te_ident_min in "${te_identity_min_values[@]}"; do
    run_id=$((run_id + 1))
    if (( run_id == 1 || run_id % PROGRESS_EVERY == 0 || run_id == total_runs )); then
        echo "[steller-sweep] run ${run_id}/${total_runs} (alpha=${evidence_alpha}, mad=${mad_max}, asm_len=${asm_min_len}, asm_id=${asm_min_ident}, vote=${te_vote_min}, te_id=${te_ident_min})"
    fi

    scientific_file="$OUT_DIR/scientific.run_${run_id}.txt"
    calls_file="$OUT_DIR/calls.run_${run_id}.tsv"

    (
        cd "$ROOT_DIR"
        PLACER_EMIT_LOW_CONFIDENCE_CALLS="$emit_low_conf" \
        PLACER_LOW_CONF_MIN_SUPPORT_READS="$low_conf_min" \
        PLACER_LOW_CONF_MAX_TIER="$low_conf_tier" \
        PLACER_EVIDENCE_MIN_SUPPORT_ALPHA="$evidence_alpha" \
        PLACER_EVIDENCE_BREAKPOINT_MAD_MAX="$mad_max" \
        PLACER_ASSEMBLY_MIN_FRAGMENT_LEN="$asm_min_len" \
        PLACER_ASSEMBLY_MIN_CONSENSUS_LEN="$asm_min_len" \
        PLACER_ASSEMBLY_MIN_IDENTITY_EST="$asm_min_ident" \
        PLACER_TE_VOTE_FRACTION_MIN="$te_vote_min" \
        PLACER_TE_MEDIAN_IDENTITY_MIN="$te_ident_min" \
        "$PLACER_BIN" "$BAM" "$REF" "$TE_FASTA" >/dev/null 2>/dev/null
        cp scientific.txt "$scientific_file"
    )

    extract_calls "$scientific_file" "$calls_file"
    read -r truth_total call_total truth_hits_pos truth_hits_family fp_calls precision_pos recall_pos f1_pos < <(
        score_calls "$calls_file" "$TRUTH" "$MATCH_WINDOW_BP"
    )
    matched_calls_pos=$(( call_total - fp_calls ))

    echo -e "${run_id}\t${emit_low_conf}\t${low_conf_min}\t${low_conf_tier}\t${evidence_alpha}\t${mad_max}\t${asm_min_len}\t${asm_min_ident}\t${te_vote_min}\t${te_ident_min}\t${truth_total}\t${call_total}\t${truth_hits_pos}\t${truth_hits_family}\t${matched_calls_pos}\t${fp_calls}\t${precision_pos}\t${recall_pos}\t${f1_pos}\t${MATCH_WINDOW_BP}" >> "$RESULTS"
    if [[ "$KEEP_INTERMEDIATE" != "1" ]]; then
        rm -f "$scientific_file" "$calls_file"
    fi
done
done
done
done
done
done
done
done
done

BEST="$OUT_DIR/sweep_top20.tsv"
SORTED_TMP="$OUT_DIR/.sweep_sorted.tmp.tsv"
tail -n +2 "$RESULTS" | sort -t $'\t' \
    -k13,13nr \
    -k14,14nr \
    -k16,16n \
    -k19,19nr \
    -k12,12n > "$SORTED_TMP"
{
    head -n 1 "$RESULTS"
    awk 'NR<=20 {print}' "$SORTED_TMP"
} > "$BEST"
rm -f "$SORTED_TMP"

echo "[steller-sweep] done"
echo "[steller-sweep] all results: $RESULTS"
echo "[steller-sweep] top20: $BEST"
sed -n '1,25p' "$BEST"
