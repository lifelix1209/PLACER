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
OUT_DIR="${OUT_DIR:-$ROOT_DIR/benchmark_steller_out}"
MATCH_WINDOW_BP="${MATCH_WINDOW_BP:-200}"

mkdir -p "$OUT_DIR"

if [[ ! -x "$PLACER_BIN" ]]; then
    echo "[steller-benchmark] placer binary not found: $PLACER_BIN" >&2
    exit 1
fi
for f in "$BAM" "$REF" "$TE_FASTA" "$TRUTH"; do
    if [[ ! -f "$f" ]]; then
        echo "[steller-benchmark] required file missing: $f" >&2
        exit 1
    fi
done

extract_calls() {
    local scientific="$1"
    local out="$2"
    awk -F'\t' 'BEGIN{OFS="\t"} $1 !~ /^#/ && NF >= 40 && $2 ~ /^-?[0-9]+$/ {print $1,$3,$6}' \
        "$scientific" > "$out"
}

score_mode() {
    local calls="$1"
    local truth="$2"
    local mode="$3"
    local window_bp="$4"
    local out="$5"
    awk -v w="$window_bp" -v mode="$mode" 'BEGIN{
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
                    matched_calls_family++;
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
            precision = (n_calls > 0) ? matched_calls / n_calls : 0.0;
            recall = (n_truth > 0) ? truth_hits / n_truth : 0.0;
            family_recall = (n_truth > 0) ? truth_family_hits / n_truth : 0.0;
            family_precision = (n_calls > 0) ? matched_calls_family / n_calls : 0.0;

            printf("%s\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n",
                mode, n_truth, n_calls, truth_hits, truth_family_hits,
                matched_calls, precision, recall, family_precision, family_recall, w);
        }' "$truth" "$calls" >> "$out"
}

run_mode() {
    local mode="$1"
    local scientific="$OUT_DIR/scientific.${mode}.txt"
    local calls="$OUT_DIR/calls.${mode}.tsv"
    local stderr_file="$OUT_DIR/placer.${mode}.stderr"

    case "$mode" in
        default)
            (
                cd "$ROOT_DIR"
                "$PLACER_BIN" "$BAM" "$REF" "$TE_FASTA" > /dev/null 2> "$stderr_file"
                cp scientific.txt "$scientific"
            )
            ;;
        relaxed)
            (
                cd "$ROOT_DIR"
                PLACER_EMIT_LOW_CONFIDENCE_CALLS=1 \
                PLACER_LOW_CONF_MIN_SUPPORT_READS=1 \
                PLACER_LOW_CONF_MAX_TIER=3 \
                PLACER_EVIDENCE_MIN_SUPPORT_ALPHA=0.02 \
                PLACER_EVIDENCE_BREAKPOINT_MAD_MAX=200 \
                PLACER_ASSEMBLY_MIN_FRAGMENT_LEN=40 \
                PLACER_ASSEMBLY_MIN_CONSENSUS_LEN=40 \
                PLACER_ASSEMBLY_MIN_IDENTITY_EST=0.40 \
                PLACER_PURE_SOFTCLIP_MIN_READS=2 \
                PLACER_PURE_SOFTCLIP_MIN_FRAGMENTS=2 \
                PLACER_PURE_SOFTCLIP_MIN_IDENTITY=0.20 \
                "$PLACER_BIN" "$BAM" "$REF" "$TE_FASTA" > /dev/null 2> "$stderr_file"
                cp scientific.txt "$scientific"
            )
            ;;
        *)
            echo "[steller-benchmark] unknown mode: $mode" >&2
            exit 1
            ;;
    esac

    extract_calls "$scientific" "$calls"
}

RESULTS="$OUT_DIR/metrics.tsv"
echo -e "mode\ttruth_total\tcall_total\ttruth_hits_pos\ttruth_hits_family\tmatched_calls_pos\tprecision_pos\trecall_pos\tprecision_family\trecall_family\twindow_bp" > "$RESULTS"

run_mode default
score_mode "$OUT_DIR/calls.default.tsv" "$TRUTH" "default" "$MATCH_WINDOW_BP" "$RESULTS"

run_mode relaxed
score_mode "$OUT_DIR/calls.relaxed.tsv" "$TRUTH" "relaxed" "$MATCH_WINDOW_BP" "$RESULTS"

echo "[steller-benchmark] done"
echo "[steller-benchmark] results: $RESULTS"
cat "$RESULTS"
