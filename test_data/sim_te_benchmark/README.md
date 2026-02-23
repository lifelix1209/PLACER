# Synthetic TE Benchmark Dataset

## Files
- `sim_te_benchmark.bam` / `sim_te_benchmark.bam.bai`: synthetic benchmark BAM
- `ref.fa` / `ref.fa.fai`: reference FASTA used by BAM
- `te_tldr.fa`: TE library used to generate TE-like insertions
- `te_consensus_dfam.fa`: curated TE consensus templates (Dfam accessions)
- `truth_events.tsv`: full event table
- `truth_positive_strong.tsv`: baseline positive truth set
- `truth_positive_all.tsv`: all positive events (including F2/FN challenge)
- `truth_negative_fp.tsv`: negative events for F1/FP challenge
- `read_manifest.tsv`: all synthetic SAM records
- `generate_sim_te_benchmark.py`: deterministic dataset generator

## TE sequence source

This dataset now uses **real curated TE consensus sequences** instead of repeated toy seeds.

- `L1HS_5end` (`DF000000226`, v4) + `L1HS_3end` (`DF000000225`, v4)
- `AluYa5` (`DF000000053`, v4)
- `SVA_F` (`DF000001072`, v4)

All are from Dfam API family records (release 3.9, March 2025).

## Event design

| event_id | truth_label | challenge_type | scenario |
|---|---|---|---|
| TP_L1_INS_STRONG | POS | BASELINE | long_insertion_I_2000bp |
| TP_ALU_INS_STRONG | POS | BASELINE | long_insertion_I_310bp |
| TP_L1_SPLIT_CLIP | POS | BASELINE | split_plus_long_softclip_500bp |
| F1_FP_RANDOM_INS | NEG | F1_FP | non_te_long_insertion_900bp |
| F1_FP_POLYA_SOFTCLIP | NEG | F1_FP | polyA_softclip_only_350bp |
| F2_FN_LOW_SUPPORT_L1 | POS | F2_FN | te_insertion_low_support_1600bp |
| F2_FN_SHORT_INS_ALU | POS | F2_FN | te_insertion_low_end_300bp |
| TP_L1_INS_SECONDARY | POS | BASELINE | long_insertion_I_2000bp_locus2 |

Insertion and clip sizes now span **300 bp to 2000 bp** to better reflect long ONT-like events.
Insertion templates use **600 bp left + 600 bp right** aligned flanks (`>500 bp` per side).

Insertion construction is TE-specific:

- L1/SVA: high-complexity consensus segments + short poly(A) tail
- Alu: near full-length AluYa5 + poly(A) extension
- negatives keep non-TE random insertion / polyA-only softclip challenges

## Regeneration

Rebuild all benchmark artifacts (manifest/truth/TE/BAM) with:

```bash
python3 test_data/sim_te_benchmark/generate_sim_te_benchmark.py
```

## Benchmark matching recommendation

- Compare by `(chrom, pos)` with `±200bp` slack.
- Use `truth_positive_strong.tsv` for baseline recall.
- Use `truth_negative_fp.tsv` to audit false positives.
- Use `truth_positive_all.tsv` for stress recall including low-support events.
