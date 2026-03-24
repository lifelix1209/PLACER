# Component And Proposal Layer

This document describes the current proposal layer only.
It does not describe final TE acceptance.

Relevant code:

- `src/pipeline/pipeline.cpp`
- `src/component/insert_fragment_module.cpp`
- `src/component/te_quick_classifier.cpp`
- `include/pipeline.h`

## 1. Role In The Current Pipeline

The component stage exists to propose loci and collect local sequence material.

It is responsible for:

- grouping local evidence into candidate loci
- extracting insertion-related read segments
- producing junction diagnostics used downstream
- accelerating TE library search with fragment-level classification

It is not responsible for:

- final TE family assignment
- final locus acceptance
- uncertain/non-TE labeling
- confidence calibration

Final acceptance happens later on the event-level path in
`Pipeline::process_bin_records(...)`.

## 2. Component Construction

`LinearBinComponentModule::build(...)` works on Gate1-passed reads.

Current output object: `ComponentCall`

Main fields:

- `chrom`, `tid`
- `bin_start`, `bin_end`
- `anchor_pos`
- `read_indices`
- `soft_clip_read_indices`
- `split_sa_read_indices`
- `insertion_read_indices`
- `breakpoint_candidates`

Current behavior:

- evidence points are collected from soft clips, long insertions, and SA hints
- local density windows are built from those points
- reads are assigned to the strongest compatible window
- each surviving window becomes one proposal `ComponentCall`

`ComponentCall` is a candidate locus container, not a call record.

## 3. Fragment Extraction

`CigarInsertionFragmentModule` and `SplitSAFragmentModule` extract
`InsertionFragment` objects from each proposal component.

Fragment sources:

- `kClipRefLeft`
- `kClipRefRight`
- `kCigarInsertion`
- `kSplitSa`

Important fragment fields:

- `fragment_id`
- `read_id`
- `source`
- `sequence`
- `anchor_len`
- `ref_side`
- `ref_junc_pos`
- `nm`
- `split_sa_reliable`

Current downstream uses:

- event breakpoint inference
- event consensus input construction
- TE library shortlist acceleration
- denovo parent-pool evidence scanning

## 4. Fragment-Level TE Search

`TEKmerQuickClassifierModule::classify(...)` performs fragment-level TE lookup.

Current semantics:

- multi-k unique-kmer index
- optional low-complexity soft-clip suppression
- optional low-support rescue by direct identity estimation on shortlist entries
- optional TSV export through `ins_fragment_hits_tsv_path`

Returned object: `FragmentTEHit`

Main fields:

- `te_name`
- `fragment_len`
- `aligned_len_est`
- `kmer_support`
- `coverage`
- `multik_support`
- `rescue_used`
- `hit_kmers`
- `total_kmers`

This stage is a search accelerator only.

It can propose likely TE templates, but it does not decide the final family for
the call. Final family assignment comes from
`TEKmerQuickClassifierModule::align_insert_sequence(...)` on the segmented
`insert_seq`.

## 5. Boundary To Final Path

The handoff from proposal layer to final path is:

1. candidate `ComponentCall`
2. local read recollection
3. `EventReadEvidence`
4. `EventConsensus`
5. `EventSegmentation`
6. insert-sequence TE alignment
7. boundary check
8. final TE emission

Any logic that stops before `EventReadEvidence` or before segmented
`insert_seq` alignment is upstream-only.
