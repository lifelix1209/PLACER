# PLACER Current Algorithm

This document describes the algorithm that is on the current terminal path.
It is intentionally limited to the code that still participates in final TE
acceptance.

Code baseline:

- `src/main.cpp`
- `src/pipeline/pipeline.cpp`
- `src/pipeline/decision_policy.cpp`
- `src/gate1/gate1_module.cpp`
- `src/component/insert_fragment_module.cpp`
- `src/component/te_quick_classifier.cpp`
- `src/component/tsd_detector.cpp`
- `include/pipeline.h`

## 1. Scope

Current PLACER emits only final PASS non-reference TE insertion calls.

`scientific.txt` contains:

- summary rows:
  - `total_reads`
  - `gate1_passed`
  - `processed_bins`
  - `components`
  - `event_consensus_calls`
  - `genotype_calls`
  - `final_pass_calls`
  - `schema_version`
- final call rows:
  - `chrom`, `pos`, `bp_left`, `bp_right`
  - `te`, `family`, `subfamily`, `strand`, `insert_len`
  - `support_reads`, `alt_struct_reads`, `ref_span_reads`, `low_mapq_ref_span_reads`
  - `gt`, `af`, `gq`
  - `best_te_identity`, `best_te_query_coverage`, `cross_family_margin`
  - `tsd_type`, `tsd_len`
  - `left_flank_align_len`, `right_flank_align_len`, `consensus_len`
  - `qc`

The final file does not contain candidate calls, uncertain states, generic
insertion states, or post-hoc confidence labels.

## 2. Terminal Path

The current terminal path is:

1. Gate1 proposal filtering
2. Component proposal construction
3. Local fragment extraction
4. Event-level local recollection
5. Event existence genotyping
6. Event consensus construction
7. Event consensus segmentation
8. Insert-sequence TE alignment
9. Boundary consistency check
10. Final TE acceptance and emission

Only steps 4-10 decide whether a final call exists.

## 3. Stage Details

### 3.1 Gate1

`SignalFirstGate1Module` removes obviously irrelevant reads and keeps reads with
structural insertion evidence or sufficient mapping quality.

This stage is proposal-only. It cannot emit final calls.

Current semantics:

- long `I` CIGAR segments are direct proposal signals, even when read-level
  `MAPQ` is low
- soft-clip-specific flank and global `NM` vetoes apply to clip-only proposal
  reads, not to reads that already carry a direct long-insertion signal

### 3.2 Component Proposal

`LinearBinComponentModule::build(...)` converts local evidence points into
candidate windows and candidate loci.

Current semantics:

- read-to-window assignment prefers breakpoint-specific evidence over proxy
  evidence
- long `I` CIGAR evidence outranks soft clips, and soft clips outrank SA-only
  hints
- within the same evidence class, stronger local support still wins

Each `ComponentCall` is only a proposal object:

- `chrom`, `tid`
- `bin_start`, `bin_end`
- `anchor_pos`
- `read_indices`
- per-class read indices
- `breakpoint_candidates`

This stage cannot accept or reject a TE event.

### 3.3 Local Fragment Extraction

`CigarInsertionFragmentModule` and `SplitSAFragmentModule` extract:

- reference-left soft clips
- reference-right soft clips
- long `I` CIGAR segments
- split/SA-derived junction fragments

These fragments are used for:

- candidate-side breakpoint hints
- event string construction
- TE library shortlist acceleration

They are not final TE evidence by themselves.

### 3.4 Event-Level Local Recollection

`collect_event_read_evidence(...)` re-fetches reads around the proposed locus and
recounts only event-level evidence.

The current evidence object is `EventReadEvidence`:

- `bp_left`, `bp_right`
- `alt_split_reads`
- `alt_indel_reads`
- `alt_left_clip_reads`
- `alt_right_clip_reads`
- `alt_struct_reads`
- `ref_span_reads`
- `low_mapq_ref_span_reads`
- `support_qnames`
- `ref_span_qnames`

Current semantics:

- breakpoint resolution still defines the narrow ref-span recount interval
- alt-support recount may use the enclosing component breakpoint envelope so
  nearby clip support is not dropped just because proxy breakpoints are less
  precise than the winning direct-indel breakpoint

This replaces the old generic evidence summary.

Breakpoint recount is local to `ComponentCall::anchor_pos`.
If recollection contains multiple breakpoint clusters, this stage compares
breakpoint hypotheses built from split/indel/clip evidence and resolves the
event breakpoint before recomputing `alt_struct_reads` and `ref_span_reads`.
Candidate breakpoint hypotheses are ranked by breakpoint-bearing evidence
strength first, then by anchor locality. Split/indel-derived hypotheses
therefore outrank clip-only clusters when the latter are more abundant but less
specific.

### 3.5 Event Existence

`genotype_call(...)` invokes `genotype_event_from_alt_vs_ref(...)` on
`alt_struct_reads` vs `ref_span_reads`.

The final path continues only when:

- `GT != 0/0`
- `GT != ./.`
- `GQ >= 20`

If the best genotype is reference or uninformative, the locus is rejected.

### 3.6 Event Consensus

`build_event_consensus(...)` builds an event consensus from event strings that
cover the event with full context when available, and otherwise falls back to
one-sided clip-derived event strings for rescue.

Implementation:

- full-context event strings (`CIGAR` insertion / reliable split) are preferred
  and clip-only strings do not mix into that POA once full-context strings
  exist
- if no full-context event strings exist, one-sided clip strings enter POA as a
  rescue path
- input is event strings, not insert-only fragments
- abPOA builds the consensus

Typical rejection reasons:

- `NO_EVENT_STRING_READS`
- `INSUFFICIENT_EVENT_READS`
- `EMPTY_EVENT_CONSENSUS`

### 3.7 Event Segmentation

`segment_event_consensus(...)` decomposes the event consensus into:

- `left_flank_seq`
- `insert_seq`
- `right_flank_seq`

It also records:

- `left_ref_start`, `left_ref_end`
- `right_ref_start`, `right_ref_end`
- `left_flank_align_len`, `right_flank_align_len`
- `left_flank_identity`, `right_flank_identity`

The final path continues only when tripartite segmentation succeeds.

Current semantics:

- left/right flank candidates are scored independently
- near-best repeated placements are carried forward
- final disambiguation happens at the left/right pair level, not by single-side
  early veto

### 3.8 Insert-Sequence TE Alignment

`align_insert_seq_to_te(...)` calls
`TEKmerQuickClassifierModule::align_insert_sequence(...)`.

Current semantics:

- k-mer index still builds the TE shortlist from the full `insert_seq`
- final family ambiguity is resolved primarily by base-level `insert_seq`
  alignment, with shortlist support used only as a bounded secondary prior
- final subfamily comes from the best-aligned template inside the winning family
- raw fragment votes do not decide final family

Current hard conditions in this stage:

- `insert_seq_len >= 80`
- `best_te_identity >= 0.80`
- `best_te_query_coverage >= 0.80`
- `cross_family_margin >= 0.10`

Typical rejection reasons:

- `INSERT_SEQ_TOO_SHORT`
- `NO_TE_ALIGNMENT_SHORTLIST`
- `NO_TE_ALIGNMENT_MATCH`
- `TE_ALIGNMENT_LOW_IDENTITY`
- `TE_ALIGNMENT_LOW_QUERY_COVERAGE`
- `TE_ALIGNMENT_CROSS_FAMILY_AMBIGUOUS`

### 3.9 Boundary Consistency

`check_boundary_consistency(...)` validates the reference-side connection implied
by the segmented event.

Allowed boundary classes on the current path:

- `TSD`
- `BLUNT`
- `SMALL_DEL`

If the left/right connection cannot be explained as one coherent insertion
event, the locus is rejected.

### 3.10 Final Emission

`emit_final_te_call(...)` materializes the final row only after all hard checks
have passed:

- event existence
- event consensus
- event segmentation
- insert-sequence TE alignment
- boundary consistency

`finalize_final_calls(...)` then deduplicates by locus and keeps the strongest
call by:

1. `support_reads`
2. `gq`
3. `cross_family_margin`
4. `best_te_identity`
5. resolved TE name preference over `NA`/`UNK`
6. `event_consensus_len`
7. `best_te_query_coverage`

## 4. Proposal/Internal Code That Still Exists

The repository still contains proposal-layer or internal-only code that is not
part of final acceptance:

- `component_module_.build()`
- `ins_fragment_module_.extract()`
- fragment-level low-k rescue inside `TEKmerQuickClassifierModule::classify(...)`

Their current role is limited to candidate generation, debug export, or search
acceleration.

They must not be interpreted as final TE calling logic.

## 5. Removed Historical Logic

The following historical semantics are no longer on the current path:

- breakpoint-first generic insertion acceptance
- generic evidence hard filtering
- component-level TE vote final acceptance
- anchor-locked TE consensus final acceptance
- placeability scoring and tiering
- post-hoc confidence calibration
- open-set `TE_CERTAIN / TE_UNCERTAIN / NON_TE`
- rescue-path final acceptance labels

## 6. Current Regression Coverage

Key tests covering the current path:

- `tests/test_event_read_evidence.cpp`
- `tests/test_event_segmentation.cpp`
- `tests/test_insert_te_alignment.cpp`
- `tests/test_terminal_path_without_legacy_assembly.cpp`
- `tests/test_ref_span_dominant_rejection.cpp`
- `tests/test_run_sharded_placer.py`
