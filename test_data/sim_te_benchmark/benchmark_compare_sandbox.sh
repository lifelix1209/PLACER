#!/usr/bin/env bash
set -euo pipefail

# Reproducible TLDR vs PLACER benchmark in sandboxed terminals.
# Key constraints handled here:
# 1) `conda run` may fail with NoWritableEnvsDirError in restricted sessions.
# 2) MAFFT may fail writing to /dev/stderr in sandboxed terminals.
#    We solve this by patching TLDR's MAFFT invocation to use `--progress <file>`
#    instead of relying on MAFFT's default `/dev/stderr` progress sink.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${ROOT_DIR}/test_data/sim_te_benchmark"

TLDR_ENV_PREFIX="${TLDR_ENV_PREFIX:-/opt/anaconda3/envs/tldr-dev-test}"
TLDR_PY="${TLDR_PY:-${TLDR_ENV_PREFIX}/bin/python}"
TLDR_ENTRY="${TLDR_ENTRY:-/Users/hanzhang/Desktop/File/softwares/opensource/tldr_optimized/tldr/tldr}"
PLACER_BIN="${PLACER_BIN:-${ROOT_DIR}/build/placer}"
OUT_DIR="${OUT_DIR:-/tmp/placer_tldr_compare}"
WORKERS="${WORKERS:-8}"

mkdir -p "${OUT_DIR}"
mkdir -p "${OUT_DIR}/bin"

export PATH="${TLDR_ENV_PREFIX}/bin:${PATH}"

echo "[Preflight] Binary checks"
for t in exonerate minimap2 samtools mafft; do
    if ! command -v "${t}" >/dev/null 2>&1; then
        echo "[ERROR] missing dependency: ${t}" >&2
        exit 3
    fi
    echo "  - ${t}: $(command -v "${t}")"
done
echo "  - python: ${TLDR_PY}"

if [[ ! -x "${TLDR_PY}" ]]; then
    echo "[ERROR] python missing: ${TLDR_PY}" >&2
    exit 4
fi
if [[ ! -f "${TLDR_ENTRY}" ]]; then
    echo "[ERROR] tldr entry missing: ${TLDR_ENTRY}" >&2
    exit 5
fi
if [[ ! -x "${PLACER_BIN}" ]]; then
    echo "[ERROR] placer binary missing: ${PLACER_BIN}" >&2
    echo "        Build first: cmake --build build -j" >&2
    exit 6
fi

echo "[Preflight] MAFFT sandbox sanity test"
TEST_FA="${OUT_DIR}/mafft_test.fa"
TEST_OUT="${OUT_DIR}/mafft_test.out"
cat > "${TEST_FA}" <<'EOF'
>a
ACGT
>b
ACGT
EOF
if ! mafft --randomseed 1 --progress "${OUT_DIR}/mafft_test.progress.log" "${TEST_FA}" > "${TEST_OUT}" 2>"${OUT_DIR}/mafft_test.err"; then
    echo "[ERROR] mafft failed with explicit --progress log" >&2
    exit 7
fi
if [[ ! -s "${TEST_OUT}" ]]; then
    echo "[ERROR] mafft produced empty output" >&2
    exit 8
fi

PATCHED_TLDR="${OUT_DIR}/tldr_patched.py"
python3 - "${TLDR_ENTRY}" "${PATCHED_TLDR}" <<'PY'
import sys
from pathlib import Path

src_path = Path(sys.argv[1])
dst_path = Path(sys.argv[2])
text = src_path.read_text()
needle = "    args = ['mafft', '--randomseed', '1', in_fn]\n"
replacement = (
    "    progress_log = os.path.abspath(out_fn + '.progress.log')\n"
    "    args = ['mafft', '--randomseed', '1', '--progress', progress_log, in_fn]\n"
)
if needle not in text:
    raise SystemExit("Failed to patch TLDR: MAFFT invocation pattern not found")
dst_path.write_text(text.replace(needle, replacement, 1))
PY
chmod +x "${PATCHED_TLDR}"

TLDR_OUT="${OUT_DIR}/tldr"
PLACER_STDOUT="${OUT_DIR}/placer.stdout"
PLACER_STDERR="${OUT_DIR}/placer.stderr"
PLACER_SCI="${OUT_DIR}/placer.scientific.txt"
TLDR_TIME="${OUT_DIR}/tldr.time.txt"
PLACER_TIME="${OUT_DIR}/placer.time.txt"

echo "[Config] workers=${WORKERS}"

echo "[Run] TLDR"
/usr/bin/time -p -o "${TLDR_TIME}" \
    "${TLDR_PY}" "${PATCHED_TLDR}" \
    -b "${DATA_DIR}/sim_te_benchmark.bam" \
    -e "${DATA_DIR}/te_tldr.fa" \
    -r "${DATA_DIR}/ref.fa" \
    -p "${WORKERS}" \
    -o "${TLDR_OUT}" \
    > "${TLDR_OUT}.stdout" 2> "${TLDR_OUT}.stderr"

echo "[Run] PLACER"
PLACER_PARALLEL=1 \
PLACER_PARALLEL_WORKERS="${WORKERS}" \
PLACER_BAM_THREADS="${WORKERS}" \
/usr/bin/time -p -o "${PLACER_TIME}" \
    "${PLACER_BIN}" \
    "${DATA_DIR}/sim_te_benchmark.bam" \
    "${DATA_DIR}/ref.fa" \
    "${DATA_DIR}/te_tldr.fa" \
    > "${PLACER_STDOUT}" 2> "${PLACER_STDERR}"
cp "${ROOT_DIR}/scientific.txt" "${PLACER_SCI}"

awk -F'\t' 'NR>1 && $1!="" {pos=int(($3+$4)/2); print $2"\t"pos}' "${TLDR_OUT}.table.txt" > "${OUT_DIR}/tldr.calls.all.tsv"
awk -F'\t' 'NR>1 && $1!="" && $26=="PASS" {pos=int(($3+$4)/2); print $2"\t"pos}' "${TLDR_OUT}.table.txt" > "${OUT_DIR}/tldr.calls.pass.tsv"
awk -F'\t' '$1 ~ /^chr/ {print $1"\t"$3}' "${PLACER_SCI}" > "${OUT_DIR}/placer.calls.tsv"

calc_metrics() {
    local calls_file="$1"
    local label="$2"
    local calls strong all_pos neg tp fp precision recall f1
    calls="$(wc -l < "${calls_file}" | tr -d ' ')"
    local truth_total
    truth_total="$(awk -F'\t' 'NR>1 {n++} END{print n+0}' "${DATA_DIR}/truth_positive_all.tsv")"
    strong="$(awk -F'\t' 'FNR==NR {c[$1":"$2]=1; next} NR==1 {next} {hit=0; for (k in c){split(k,a,":"); if (a[1]==$1 && (a[2]-$2<=200 && $2-a[2]<=200)) {hit=1; break}} if(hit) s++} END{print s+0}' "${calls_file}" "${DATA_DIR}/truth_positive_strong.tsv")"
    all_pos="$(awk -F'\t' 'FNR==NR {c[$1":"$2]=1; next} NR==1 {next} {hit=0; for (k in c){split(k,a,":"); if (a[1]==$1 && (a[2]-$2<=200 && $2-a[2]<=200)) {hit=1; break}} if(hit) s++} END{print s+0}' "${calls_file}" "${DATA_DIR}/truth_positive_all.tsv")"
    neg="$(awk -F'\t' 'FNR==NR {c[$1":"$2]=1; next} NR==1 {next} {hit=0; for (k in c){split(k,a,":"); if (a[1]==$1 && (a[2]-$2<=200 && $2-a[2]<=200)) {hit=1; break}} if(hit) s++} END{print s+0}' "${calls_file}" "${DATA_DIR}/truth_negative_fp.tsv")"
    tp="$(awk -F'\t' 'FNR==NR {if(NR>1) t[++n]=$1":"$2; next} {hit=0; split($0,a,"\t"); for(i=1;i<=n;i++){split(t[i],b,":"); if (b[1]==a[1] && (b[2]-a[2]<=200 && a[2]-b[2]<=200)) {hit=1; break}} if(hit) s++} END{print s+0}' "${DATA_DIR}/truth_positive_all.tsv" "${calls_file}")"
    fp=$((calls - tp))
    precision="$(awk -v tp="${tp}" -v calls="${calls}" 'BEGIN{if(calls>0) printf "%.4f", tp/calls; else printf "0.0000"}')"
    recall="$(awk -v hit="${all_pos}" -v total="${truth_total}" 'BEGIN{if(total>0) printf "%.4f", hit/total; else printf "0.0000"}')"
    f1="$(awk -v p="${precision}" -v r="${recall}" 'BEGIN{if((p+r)>0) printf "%.4f", 2*p*r/(p+r); else printf "0.0000"}')"
    echo "${label}: calls=${calls} strong_hit=${strong} all_hit=${all_pos} neg_hit=${neg} tp_calls=${tp} fp_calls=${fp} precision=${precision} recall=${recall} f1=${f1}"
}

echo "[Summary] Metrics"
calc_metrics "${OUT_DIR}/tldr.calls.all.tsv" "TLDR_ALL"
calc_metrics "${OUT_DIR}/tldr.calls.pass.tsv" "TLDR_PASS"
calc_metrics "${OUT_DIR}/placer.calls.tsv" "PLACER"

echo "[Summary] Runtime (seconds)"
awk '$1=="real"{print "  - TLDR: "$2}' "${TLDR_TIME}"
awk '$1=="real"{print "  - PLACER: "$2}' "${PLACER_TIME}"

echo "[Summary] TLDR filter distribution"
awk -F'\t' 'NR>1 && $1!="" {c[$26]++} END{for(k in c) print "  - "k": "c[k]}' "${TLDR_OUT}.table.txt" | sort

echo "[Done] Outputs in ${OUT_DIR}"
