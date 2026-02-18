# sim_te_benchmark Simulation Documentation

## 1. Purpose

`test_data/sim_te_benchmark` is a deterministic synthetic benchmark for TE insertion detection.
It focuses on three evaluation axes:

1. baseline detectability (`BASELINE`)
2. false-positive suppression (`F1_FP`)
3. expected false-negative stress (`F2_FN`)

The dataset is template-driven (not a generic read-error simulator), so each event is auditable.

## 2. TE Sequence Provenance (Curated Real TE Consensus)

This benchmark no longer uses repeated toy seeds to build TE inserts.
It now uses curated Dfam consensus sequences:

- `L1HS_5end` (`DF000000226`, v4) + `L1HS_3end` (`DF000000225`, v4), concatenated as `L1:L1HS`
- `AluYa5` (`DF000000053`, v4), used as `Alu:AluYa5`
- `SVA_F` (`DF000001072`, v4), used as `SVA:SVA_F`

Rationale:

- Dfam is a curated TE family database and RepeatMasker standard source.
- Human non-reference polymorphic insertions are dominated by LINE-1, Alu, and SVA families.
- Recent L1HS in Dfam are represented by curated 5' and 3' models, so we concatenate those two curated parts for a long template.
- To keep TLDR remapping stable in this synthetic benchmark, long L1/SVA inserts use higher-complexity consensus segments with short poly(A) tails.

## 3. Current Design (Long ONT-style Inserts)

This version models longer ONT-like insertion signatures:

- insertion/clip length range: **300 bp to 2000 bp**
- aligned flank design: **600 bp left + 600 bp right** for insertion templates (`>500bp` per side)
- long insertion CIGARs such as:
  - `600M2000I600M` (L1)
  - `600M2000I600M` (L1, secondary locus)
  - `600M1600I600M` (L1)
  - `600M900I600M` (non-TE random negative)
  - `600M310I600M` (Alu)
  - `600M300I600M` (Alu)
- long soft-clips / split support:
  - `600M350S` (polyA challenge)
  - `600M500S` + `500S600M` (split + supplementary)

## 4. Files and Roles

Directory: `test_data/sim_te_benchmark`

- `ref.fa` / `ref.fa.fai`: single contig reference (`chrTEST`)
- `te_consensus_dfam.fa`: source consensus records with Dfam accession/version in headers
- `te_tldr.fa`: benchmark TE library emitted by generator
  - `L1:L1HS` length 3038 (L1HS_5end + L1HS_3end)
  - `Alu:AluYa5` length 311
  - `SVA:SVA_F` length 1375
- `truth_events.tsv`: master truth table (8 events)
- `truth_positive_strong.tsv`: baseline positives
- `truth_positive_all.tsv`: all positives
- `truth_negative_fp.tsv`: negatives
- `read_manifest.tsv`: full synthetic alignment manifest
- `sim_te_benchmark.bam` / `.bai`: final sorted/indexed BAM
- `generate_sim_te_benchmark.py`: deterministic generator

## 5. Event Logic

Ground-truth events are defined in `truth_events.tsv`.
Each event has fixed position, label, challenge type, support count, and a specific read template.

- `TP_L1_INS_STRONG`: POS, `600M2000I600M`, support=12
- `TP_L1_INS_SECONDARY`: POS, `600M2000I600M`, support=12
- `TP_ALU_INS_STRONG`: POS, `600M310I600M`, support=10
- `TP_L1_SPLIT_CLIP`: POS, `600M500S` + supplementary `500S600M`, support=6 read names
- `F1_FP_RANDOM_INS`: NEG, `600M900I600M`, non-TE random inserted sequence, support=10
- `F1_FP_POLYA_SOFTCLIP`: NEG, `600M350S`, polyA soft-clip tail, support=8
- `F2_FN_LOW_SUPPORT_L1`: POS, `600M1600I600M`, low support (2 reads)
- `F2_FN_SHORT_INS_ALU`: POS, `600M300I600M`, lower-end long insertion, support=7

Background reads:

- 140 reads named `BG_0000..BG_0139`
- CIGAR `400M`, mapQ 60
- exact reference slices

## 6. Generation Pipeline

The generator script:

1. loads `ref.fa`
2. loads curated Dfam consensus templates from `te_consensus_dfam.fa`
3. constructs event insert sequences:
   - L1/SVA: high-complexity consensus segment + short poly(A)
   - Alu: near full-length consensus + poly(A)
   - negatives: random DNA or polyA-only softclip
4. emits background + event reads into SAM records
5. writes truth tables and `read_manifest.tsv`
6. converts SAM to BAM and runs `samtools sort/index`

BAM header therefore shows:

- `samtools view -b ... sim_te_benchmark.sam`
- `samtools sort ... sim_te_benchmark.bam`

## 7. Reproducibility

Regenerate all artifacts with:

```bash
python3 test_data/sim_te_benchmark/generate_sim_te_benchmark.py
```

Default seed is fixed (`20260218`) for deterministic outputs.

## 8. Quick Integrity Checks

Recommended checks after regeneration:

1. `read_manifest.tsv` row count and event support counts
2. TE FASTA lengths (`3038/311/1375`)
3. CIGAR distribution includes the long-insert templates above
4. BAM header contains the expected SAM->BAM->sort lineage

## 9. References

1. Dfam API version endpoint (release metadata): https://dfam.org/api/version
2. Dfam family records used:
   - L1HS_5end: https://dfam.org/api/families/DF000000226
   - L1HS_3end: https://dfam.org/api/families/DF000000225
   - AluYa5: https://dfam.org/api/families/DF000000053
   - SVA_F: https://dfam.org/api/families/DF000001072
3. Dfam as RepeatMasker standard library source:
   - https://www.repeatmasker.org/dev/RepeatModeler/
   - L1HS split-model note in repeat browser: https://www.repeatbrowser.org/?dfam=DF000000226
4. Human polymorphic non-reference TE classes (L1/Alu/SVA):
   - Ewing AD et al. *Cell Genomics* (2023): https://pmc.ncbi.nlm.nih.gov/articles/PMC9910659/
5. Hallmarks of LINE1-mediated insertions (polyA/TSD/5' truncation):
   - Raiz J et al. *Nucleic Acids Research* (2012): https://pmc.ncbi.nlm.nih.gov/articles/PMC3439880/
   - Occurrence frequencies in long-read validated insertions: *Mobile DNA* (2025): https://pubmed.ncbi.nlm.nih.gov/40187185/
