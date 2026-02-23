# Benchmark Notes

## Real-Data Observation

In ONT real-data runs, TE classification was over-filtered before assembly.

Observed summary:

- total extracted fragments: `9962`
- fragments with any k-mer match: `538` (`5.4%`)
- fragments with `kmer_support >= 0.50`: `0`
- fragments with `kmer_support >= 0.30`: `5`

Root cause:

- `kmer_support` in current quick classifier is `best_hit_kmers / total_kmers`, not true alignment identity.
- Old hard thresholds (`te_vote_fraction_min=0.60`, `te_median_identity_min=0.50`) were too strict for noisy ONT fragments.
- Candidates were dropped before assembly, so assembly could not rescue low-identity TE calls.

## Current Design (Implemented)

Two-stage TE decision:

- Stage A (pre-assembly): weak gate
  - `te_name` non-empty
  - `fragment_count >= max(1, te_min_fragments_for_vote / 2)`
- Stage B (post-assembly): strong gate
  - `PASS_CLASSIC` if standard vote thresholds pass
  - `PASS_ASM_RESCUE` if rescue thresholds pass and assembly identity is high enough

Pure soft-clip components are evaluated in Stage B with separate read/fragment/identity constraints.

## Updated Default Parameters

- `te_kmer_size = 13`
- `te_vote_fraction_min = 0.40`
- `te_median_identity_min = 0.30`
- `te_min_fragments_for_vote = 2`
- `te_rescue_vote_fraction_min = 0.25`
- `te_rescue_median_identity_min = 0.20`
- `te_pure_softclip_min_reads = 6`
- `te_pure_softclip_min_fragments = 6`
- `te_pure_softclip_min_identity = 0.35`

## Useful Runtime Overrides

```bash
PLACER_TE_KMER_SIZE=13 \
PLACER_TE_VOTE_FRACTION_MIN=0.40 \
PLACER_TE_MEDIAN_IDENTITY_MIN=0.30 \
PLACER_TE_RESCUE_VOTE_FRACTION_MIN=0.25 \
PLACER_TE_RESCUE_MEDIAN_IDENTITY_MIN=0.20 \
PLACER_PURE_SOFTCLIP_MIN_READS=6 \
PLACER_PURE_SOFTCLIP_MIN_FRAGMENTS=6 \
PLACER_PURE_SOFTCLIP_MIN_IDENTITY=0.35 \
./build/placer <input.bam> <ref.fa> <te.fa>
```

## Sanity Check on `sim_te_benchmark`

With current defaults and tuning script:

- best params in search grid remain:
  - `PLACER_ASSEMBLY_POA_MIN_READS=2`
  - `PLACER_EVIDENCE_MIN_SUPPORT_ALPHA=0.08`
  - `PLACER_EVIDENCE_BREAKPOINT_MAD_MAX=80`
- best metrics:
  - `calls=6`
  - `strong_hits=3`
  - `all_hits=6`
  - `neg_fp_hits=0`
