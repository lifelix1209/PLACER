#!/usr/bin/env bash
# Step 1/4 of the assembly-based TE truth set.
# Assemble the D2 ONT reads with Flye so GraffiTE can derive a caller-independent
# truth set of non-reference TE insertions.
#
# D2 is ONT PromethION (~25x, read N50 ~6.3 kb). At this read length the assembly
# will be fragmented in long/repetitive TE regions, so treat the resulting truth
# as high-precision but incomplete for recall (pair it with a simulation truth).
#
# Usage: assemble_d2_flye.sh <input.bam> <out_dir> [threads] [read_type]
#   read_type: --nano-hq (R10 / recent basecalls, default) or --nano-raw (older R9)
set -euo pipefail

BAM="${1:?path to D2 BAM required}"
OUT="${2:?output directory required}"
THREADS="${3:-16}"
READ_TYPE="${4:---nano-hq}"

for tool in samtools flye; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "ERROR: '$tool' not found on PATH. Install e.g. 'mamba install -c bioconda $tool'." >&2
    exit 1
  }
done
[[ -f "$BAM" ]] || { echo "ERROR: input BAM not found: $BAM" >&2; exit 1; }

mkdir -p "$OUT"
READS="$OUT/d2.reads.fq.gz"

if [[ ! -s "$READS" ]]; then
  echo "[1/2] Extracting reads -> $READS (dropping secondary/supplementary; unmapped kept)"
  samtools fastq -@ "$THREADS" -F 0x900 -n "$BAM" | gzip -c > "$READS"
else
  echo "[1/2] Reusing existing $READS"
fi

echo "[2/2] Running Flye ($READ_TYPE, $THREADS threads) -> $OUT/flye"
flye "$READ_TYPE" "$READS" --threads "$THREADS" --out-dir "$OUT/flye"

echo "DONE. Assembly: $OUT/flye/assembly.fasta"
