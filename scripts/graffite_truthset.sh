#!/usr/bin/env bash
# Step 2/4 of the assembly-based TE truth set.
# Run GraffiTE in assembly mode to call non-reference TE insertion polymorphisms
# from the D2 assembly vs the reference. Genotyping is disabled: we only need the
# assembly-derived pME calls as truth, not graph genotyping.
#
# Usage: graffite_truthset.sh <assembly.fasta> <reference.fa> <te_library.fa> <out_dir> [sample] [main.nf]
set -euo pipefail

ASM="${1:?assembly FASTA required}"
REF="${2:?reference FASTA required}"
LIB="${3:?TE library FASTA required}"
OUT="${4:?output directory required}"
SAMPLE="${5:-D2}"
MAIN_NF="${6:-.research/repos/GraffiTE/main.nf}"

command -v nextflow >/dev/null 2>&1 || {
  echo "ERROR: nextflow not found. Install: 'curl -s https://get.nextflow.io | bash' (needs Java 11+)." >&2
  exit 1
}
command -v singularity >/dev/null 2>&1 || command -v docker >/dev/null 2>&1 || {
  echo "ERROR: GraffiTE needs singularity or docker for its containers; neither found." >&2
  exit 1
}
for f in "$ASM" "$REF" "$LIB" "$MAIN_NF"; do
  [[ -f "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

mkdir -p "$OUT"
CSV="$OUT/assemblies.csv"
ASM_ABS="$(cd "$(dirname "$ASM")" && pwd)/$(basename "$ASM")"
printf 'path,sample\n%s,%s\n' "$ASM_ABS" "$SAMPLE" > "$CSV"
echo "Wrote $CSV"

nextflow run "$MAIN_NF" \
  --assemblies "$CSV" \
  --reference "$REF" \
  --TE_library "$LIB" \
  --genotype false \
  --out "$OUT/graffite"

echo "DONE. GraffiTE output under $OUT/graffite (truth VCF: pangenome.vcf / *.vcf)."
