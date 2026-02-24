#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  scripts/bootstrap_te_reclassify.sh \
    --bam <input.bam> \
    --ref <ref.fa> \
    --base-te <base_te.fa> \
    --pass1-fasta <pass1_bootstrap_consensus.fasta> \
    [--placer ./build/placer] \
    [--outdir bootstrap_pass2] \
    [--min-len 80] \
    [-- <extra placer args>]

Description:
  1) Deduplicate pass-1 bootstrap consensus sequences by exact sequence.
  2) Merge deduplicated consensus with base TE library.
  3) Run PLACER pass-2 with merged TE library in an isolated output folder.

Notes:
  - Pass-2 run forces PLACER_BOOTSTRAP_EXPORT=0 to avoid recursive export.
  - Additional placer args after '--' are forwarded to placer.
EOF
}

abspath() {
    local p="$1"
    if [[ "$p" = /* ]]; then
        printf '%s\n' "$p"
    else
        printf '%s\n' "$(cd "$(dirname "$p")" && pwd)/$(basename "$p")"
    fi
}

PLACER_BIN="./build/placer"
OUTDIR="bootstrap_pass2"
MIN_LEN=80
BAM=""
REF=""
BASE_TE=""
PASS1_FASTA=""
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam)
            BAM="${2:-}"; shift 2 ;;
        --ref)
            REF="${2:-}"; shift 2 ;;
        --base-te)
            BASE_TE="${2:-}"; shift 2 ;;
        --pass1-fasta)
            PASS1_FASTA="${2:-}"; shift 2 ;;
        --placer)
            PLACER_BIN="${2:-}"; shift 2 ;;
        --outdir)
            OUTDIR="${2:-}"; shift 2 ;;
        --min-len)
            MIN_LEN="${2:-}"; shift 2 ;;
        --help|-h)
            usage; exit 0 ;;
        --)
            shift
            EXTRA_ARGS+=("$@")
            break ;;
        *)
            EXTRA_ARGS+=("$1")
            shift ;;
    esac
done

if [[ -z "$BAM" || -z "$REF" || -z "$BASE_TE" || -z "$PASS1_FASTA" ]]; then
    usage
    exit 1
fi

if [[ ! -x "$PLACER_BIN" ]]; then
    echo "[bootstrap] placer binary not executable: $PLACER_BIN" >&2
    exit 2
fi
if [[ ! -f "$BAM" ]]; then
    echo "[bootstrap] BAM not found: $BAM" >&2
    exit 2
fi
if [[ ! -f "$REF" ]]; then
    echo "[bootstrap] reference not found: $REF" >&2
    exit 2
fi
if [[ ! -f "$BASE_TE" ]]; then
    echo "[bootstrap] base TE library not found: $BASE_TE" >&2
    exit 2
fi
if [[ ! -f "$PASS1_FASTA" ]]; then
    echo "[bootstrap] pass1 fasta not found: $PASS1_FASTA" >&2
    exit 2
fi

if ! [[ "$MIN_LEN" =~ ^[0-9]+$ ]]; then
    echo "[bootstrap] --min-len must be an integer, got: $MIN_LEN" >&2
    exit 2
fi

mkdir -p "$OUTDIR"

PLACER_BIN_ABS="$(abspath "$PLACER_BIN")"
BAM_ABS="$(abspath "$BAM")"
REF_ABS="$(abspath "$REF")"
BASE_TE_ABS="$(abspath "$BASE_TE")"
PASS1_FASTA_ABS="$(abspath "$PASS1_FASTA")"
OUTDIR_ABS="$(abspath "$OUTDIR")"

DEDUP_FASTA="$OUTDIR_ABS/bootstrap_incremental_dedup.fasta"
MERGED_TE="$OUTDIR_ABS/te_bootstrap_merged.fasta"
PASS2_DIR="$OUTDIR_ABS/pass2_run"
PASS2_SCI="$OUTDIR_ABS/scientific_pass2.txt"
REPORT_TXT="$OUTDIR_ABS/bootstrap_report.txt"

awk -v min_len="$MIN_LEN" '
BEGIN { RS=">"; ORS="" }
NR == 1 { next }
{
    seq = $0;
    sub(/^[^\n]*\n/, "", seq);
    gsub(/[ \t\r\n]/, "", seq);
    seq = toupper(seq);
    if (length(seq) < min_len) {
        next;
    }
    if (!(seq in seen)) {
        seen[seq] = 1;
        ++n;
        printf(">BOOTSTRAP_SEQ_%d\n%s\n", n, seq);
    }
}
' "$PASS1_FASTA_ABS" > "$DEDUP_FASTA"

BASE_COUNT="$(grep -c '^>' "$BASE_TE_ABS" || true)"
DEDUP_COUNT="$(grep -c '^>' "$DEDUP_FASTA" || true)"

cat "$BASE_TE_ABS" "$DEDUP_FASTA" > "$MERGED_TE"
MERGED_COUNT="$(grep -c '^>' "$MERGED_TE" || true)"

mkdir -p "$PASS2_DIR"
(
    cd "$PASS2_DIR"
    PLACER_BOOTSTRAP_EXPORT=0 \
    "$PLACER_BIN_ABS" "$BAM_ABS" "$REF_ABS" "$MERGED_TE" "${EXTRA_ARGS[@]}"
)

if [[ -f "$PASS2_DIR/scientific.txt" ]]; then
    cp "$PASS2_DIR/scientific.txt" "$PASS2_SCI"
fi

{
    echo "bootstrap_min_len=$MIN_LEN"
    echo "base_library=$BASE_TE_ABS"
    echo "base_entries=$BASE_COUNT"
    echo "pass1_fasta=$PASS1_FASTA_ABS"
    echo "dedup_entries=$DEDUP_COUNT"
    echo "merged_library=$MERGED_TE"
    echo "merged_entries=$MERGED_COUNT"
    echo "pass2_scientific=$PASS2_SCI"
} > "$REPORT_TXT"

echo "[bootstrap] done"
echo "[bootstrap] dedup_entries=$DEDUP_COUNT merged_entries=$MERGED_COUNT"
echo "[bootstrap] pass2 scientific: $PASS2_SCI"
echo "[bootstrap] report: $REPORT_TXT"
