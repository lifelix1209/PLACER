# Synthetic TE Benchmark Dataset

## Files
- `sim_te_benchmark.bam` / `sim_te_benchmark.bam.bai`: synthetic benchmark BAM
- `ref.fa` / `ref.fa.fai`: reference FASTA used by BAM
- `te.fa`: TE library used to generate TE-like insertions
- `truth_events.tsv`: full event table
- `truth_positive_strong.tsv`: baseline positive truth set
- `truth_positive_all.tsv`: all positive events (including F2/FN challenge)
- `truth_negative_fp.tsv`: negative events for F1/FP challenge
- `read_manifest.tsv`: all synthetic SAM records

## Event design

| event_id | truth_label | challenge_type | scenario |
|---|---|---|---|
| TP_L1_INS_STRONG | POS | BASELINE | large_insertion_I |
| TP_ALU_INS_STRONG | POS | BASELINE | large_insertion_I |
| TP_L1_SPLIT_CLIP | POS | BASELINE | split_plus_softclip |
| F1_FP_RANDOM_INS | NEG | F1_FP | non_te_large_insertion |
| F1_FP_POLYA_SOFTCLIP | NEG | F1_FP | polyA_softclip_only |
| F2_FN_LOW_SUPPORT_L1 | POS | F2_FN | te_insertion_low_support |
| F2_FN_SHORT_INS_ALU | POS | F2_FN | te_insertion_short_I |
| F2_FN_LOWMAPQ_ERVL | POS | F2_FN | te_insertion_low_mapq |

## Benchmark matching recommendation

- Compare by `(chrom, pos)` with `±200bp` slack.
- Use `truth_positive_strong.tsv` for baseline recall.
- Use `truth_negative_fp.tsv` to audit false positives.
- Use `truth_positive_all.tsv` for stress recall including low-support events.
