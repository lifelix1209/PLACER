# PLACER

Current release: `0.0.4`

PLACER detects non-reference transposable element (TE) insertions from long-read BAM alignments (ONT/PacBio). The current pipeline includes evidence scoring, abPOA-only local assembly, and likelihood-based genotyping.

## Clone

abPOA is required and tracked as a git submodule.

```bash
git clone --recursive <your-repo-url>
cd PLACER
```

If you already cloned without submodules:

```bash
git submodule update --init --recursive
```

## Build

Cross-platform user-space install, configure, compile, and test:

```bash
./build.sh
```

On Windows PowerShell:

```powershell
.\build.ps1
```

Or call the shared Python entrypoint directly on any platform:

```bash
python scripts/build.py
```

By default this uses an existing abPOA checkout or initializes the submodule
when needed, creates `build/.venv`, installs `requirements.txt` there,
configures CMake into `build/`, builds `placer`, and runs CTest. The script
does not require administrator privileges and does not install system packages.

Useful options:

```bash
./build.sh --no-tests
./build.sh --no-venv
./build.sh --skip-py-deps
./build.sh --build-dir build_release --build-type Release --jobs 16
```

You can also build and immediately run PLACER:

```bash
./build.sh <input.bam> <ref.fa> <te.fa>
```

On Windows, use the same arguments with `.\build.ps1`.

Manual build:

```bash
python3 -m pip install -r requirements.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

If htslib discovery is ambiguous, you can override paths with
`HTSLIB_ROOT`, `HTSLIB_INCLUDE_DIR`, or `HTSLIB_LIBRARY`.

## Test

```bash
(cd build && ctest --output-on-failure)
```

Current CTest covers tracked engineering regression tests only. It does not
decide whether the algorithm is release-ready on the real dataset.

## Algorithm Overview

PLACER calls non-reference TE insertions from long-read alignments by keeping
the decision process local to each candidate locus:

1. **Read scan and breakpoint seeds.** The BAM is streamed once. Soft clips,
   split-read anchors, and insertion/deletion CIGAR signals seed candidate
   breakpoint hypotheses in genomic bins.
2. **Local component assembly.** Nearby breakpoint hypotheses are grouped into
   components. For each feasible component, PLACER gathers the supporting reads
   and builds a local insertion consensus with abPOA.
3. **Breakpoint closure.** The consensus is aligned back to the reference to
   classify the event as a closed insertion, one-sided insertion, unplaced
   insert, or rejected pre-segmentation case. Boundary/TSD evidence is scored
   separately from TE evidence.
4. **TE sequence evidence.** The inferred insert sequence is compared with the
   TE library using quick family prefiltering followed by exact alignment.
   Calls may be resolved to a subfamily, resolved to family only, or retained as
   unknown TE evidence when the sequence is TE-like but not confidently named.
5. **Genotype and joint explanation.** Alternative structural reads and
   reference spanning reads are genotyped with a likelihood model. PLACER scores
   reference, non-TE insertion, unknown TE insertion, and resolved TE insertion
   hypotheses for ranking, then uses a single event-explanation comparison to
   decide whether TE, non-TE insertion, reference, or artifact best explains the
   evidence.
6. **Final output tiers.** TE calls with closed breakpoints are emitted as
   `PASS_TE_CLOSED`. TE calls that are supported but not breakpoint-closed are
   emitted as `PASS_TE_IMPRECISE`, so downstream users can keep closed and
   imprecise calls separate.
7. **Reference poly-N veto.** Final TE call candidates whose final breakpoint
   coordinate or breakpoint bounds fall on `N` bases in the indexed reference
   are dropped before final de-duplication. This prevents ambiguous reference
   assembly gaps from being reported as accepted TE insertions.

The current policy deliberately avoids retethering a missed truth locus to the
nearest unrelated PLACER call. One-sided rescues are constrained to the local
component's own evidence, and evidence-tier calls remain distinguishable from
final high-confidence calls.

## Repository Hygiene

This repository is intended to publish source code only.

- Local outputs such as `build/`, `placer_out/`, shard outputs, bootstrap outputs, and assistant notes are ignored by `.gitignore`.
- `test_data/` is reserved for small reproducible fixtures. Large local references and TE libraries should stay outside the tracked set and are ignored when kept under their current local filenames.
- Keep real BAM/VCF/IGV screenshots and validation exports outside the tracked tree.
- Generated QC plots and TSVs, including `test_data/poly_n_distribution/`, are local analysis outputs and are ignored.
- Before pushing to GitHub, check `git status --short` and avoid staging large binaries, generated plots, local BAMs, machine-local absolute paths, or validation exports.

## Run

```bash
./build/placer <input.bam> <ref.fa> <te.fa>
```

Limit the initial BAM scan to an indexed region while keeping event-level
local fetches on the same original BAM:

```bash
./build/placer --region chr1:1-50000000 <input.bam> <ref.fa> <te.fa>
./build/placer --region chr1 <input.bam> <ref.fa> <te.fa>
./build/placer --final-fdr-q 0.05 <input.bam> <ref.fa> <te.fa>
./build/placer --min-final-raw-cigar-insert-len-bp 50 <input.bam> <ref.fa> <te.fa>
```

`--min-final-raw-cigar-insert-len-bp` removes final calls without a local raw
CIGAR `I` operation longer than the threshold. The default is `50`; set it to
`0` to disable this final raw-CIGAR insertion filter.

Default output:

- `scientific.txt`

Optional debug outputs (off by default to reduce I/O on large runs):

```bash
PLACER_INS_FRAGMENTS_FASTA_PATH=ins_fragments.fasta \
PLACER_INS_FRAGMENT_HITS_TSV_PATH=ins_fragment_hits.tsv \
./build/placer <input.bam> <ref.fa> <te.fa>
```

## SLURM Run

Use `scripts/run_placer.sh` as the only supported cluster submission
entrypoint. It is environment-variable driven so the repository does not carry
machine-local paths.

```bash
PLACER_OUTPUT_DIR=/data/placer_out/yohann_d23 \
PLACER_BAM=/data/input.bam \
PLACER_REF=/data/ref.fa \
PLACER_TE=/data/te.fa \
PLACER_BIN=/data/placer_out/yohann_d23/placer_bin \
sbatch -p <partition> scripts/run_placer.sh
```

Required environment variables:

- `PLACER_OUTPUT_DIR`
- `PLACER_BAM`
- `PLACER_REF`
- `PLACER_TE`

Optional overrides:

- `PLACER_BIN`
- `PLACER_BAM_BAI`
- `PLACER_REF_FAI`
- `PLACER_TE_FAI`
- `PLACER_PARALLEL_WORKERS`
- `PLACER_PARALLEL_QUEUE_MAX_TASKS`
- `PLACER_PARALLEL_RESULT_BUFFER_MAX`
- `PLACER_BAM_THREADS`
- `PLACER_FINAL_FDR_Q`
- `PLACER_MIN_FINAL_RAW_CIGAR_INSERT_LEN_BP`

`scientific.txt` now contains only final PASS TE insertion calls. The row schema is:

- `chrom`, `pos`, `bp_left`, `bp_right`
- `te`, `family`, `subfamily`, `strand`, `insert_len`
- `support_reads`, `alt_struct_reads`, `raw_cigar_insert_reads`,
  `max_raw_cigar_insert_len`, `ref_span_reads`, `low_mapq_ref_span_reads`
- `gt`, `af`, `gq`
- `best_te_identity`, `best_te_query_coverage`, `cross_family_margin`
- `te_sequence_model_label`, `te_sequence_model_score`
- `te_sequence_model_gc`, `te_sequence_model_entropy`
- `te_sequence_model_tandem_fraction`
- `te_sequence_model_low_complexity_fraction`
- `te_sequence_model_jsd_k5`, `te_sequence_model_jsd_k6`
- `te_sequence_model_k9_containment`
- `te_annotation_confidence`, `te_annotation_class`, `te_annotation_order`
- `te_annotation_intervals`
- `te_annotation_residual_fraction`, `te_annotation_masked_fraction`
- `family_status` (`COMMITTED` or `ABSTAINED`; independent of the family name,
  so a library family named `Unknown` remains a valid committed label)
- `tsd_type`, `tsd_len`
- `left_flank_align_len`, `right_flank_align_len`, `consensus_len`
- optional `insert_seq` when `PLACER_DEBUG_OUTPUT_INSERT_SEQ=1`
- `qc`
- posterior/local-FDR diagnostics: `te_posterior`, `non_te_posterior`,
  `artifact_posterior`, `posterior_qc`, `latent_mechanism`, `lfdr`,
  `worst_case_lfdr`, `lfdr_qc`
- robust mechanistic diagnostics:
  `mechanistic_lower_log_bf_te_vs_artifact`,
  `mechanistic_lower_log_bf_te_vs_non_te`,
  `mechanistic_ref_conflict_signal`, `mechanistic_ambiguity_width`,
  `mechanistic_blocks`
- final conformal-null selector diagnostics:
  `conformal_null_p`, `conformal_by_threshold`,
  `conformal_dominated_nulls`, `conformal_null_count`, `conformal_qc`

`evidence_ledger.tsv` contains candidate-level diagnostics for evaluated and
ledger-only hypotheses, including rejected/evidence-only rows. It is the
precision-first audit table; `scientific.txt` remains final PASS-only output.
Ledger rows include both narrow breakpoint coordinates (`bp_left`, `bp_right`)
and coverage-tier coordinates (`coverage_left`, `coverage_right`). Broad
owner-bin coverage is only assigned to candidate rows whose TE explanation is
`PASS_TE*`; non-TE, artifact, ambiguous, and pre-expensive rows keep narrow
breakpoint coverage. This targets TLDR-like TE evidence coverage without
over-covering PLACER-only artifact regions. The ledger also reports
split/indel/left-clip/right-clip counts, full/partial consensus input counts,
left/right anchor counts, consensus length, posterior probabilities, and
`support_qname_count` so low-allele-fraction and partial-anchor rescues can be
audited at the event-community level without writing full read-name lists to
TSV output.

`coverage_candidates.tsv` is a compact high-coverage candidate-region table
derived from the evidence ledger. It reports `coverage_left` and
`coverage_right` for downstream review or recall-oriented benchmarking.

Final TE calls are selected by a sample-local conformal-null FDR procedure.
All non-final evidence ledger rows with a mechanistic certificate are used as
internal null controls. A final candidate receives a dominance conformal
p-value: the add-one-smoothed fraction of null controls that are at least as
TE-like in every positive evidence coordinate and no more reference-supported
than the candidate. PLACER also computes the add-one-smoothed upper-tail
p-value for raw `alt_struct_reads` and uses the maximum of the dominance and
structural-support p-values. Candidate p-values are selected with the
Benjamini-Yekutieli step-up rule at `PLACER_FINAL_FDR_Q` or `--final-fdr-q`.
The default target is `0.10`. After FDR selection, selected candidates that are
Pareto-dominated by another selected candidate in the same run are removed.

Low-allele-fraction structural insertions are not promoted solely because a
single candidate has a favorable alt/ref count. The finalization stage builds a
local event graph from ledger rows that share breakpoint proximity and
supporting reads. A low-AF focal row must have a same-event structural neighbor
connected by read-support overlap before it receives the
`EVENT_COMMUNITY_STABLE` certificate used by the e-value fallback. This keeps
the recall-oriented rescue tied to event-community stability rather than a
single empirical threshold.

When no single ledger row is individually promotable, PLACER can also form a
community-level structural candidate from a connected read-overlap component of
weak event rows. The aggregate counts unique supporting read names as alt
support and uses a robust reference-span summary instead of summing overlapping
reference evidence. These calls remain `UNKNOWN` family/subfamily and are marked
`EVENT_COMMUNITY_AGGREGATED`; they are still filtered by the existing
e-value/FDR fallbacks. The reported breakpoint interval is widened to the local
event-resolution radius so the output represents the uncertainty of the
community, not an over-precise single row.

A single low-AF row can also receive `EVENT_INTERVAL_STABLE` only when it already
has a non-degenerate breakpoint interval with boundary evidence. This represents
within-row spatial stability and keeps point-like isolated low-AF rows rejected.

`null_controls.tsv` contains internal empirical-null diagnostics derived from
the same non-final evidence ledger rows. Rows with `row_type=null_control` form
the sample-local null tail. Rows with `row_type=final_call_tail_query` report
the final call's conservative mechanistic score, legacy empirical upper-tail
probability, and conformal selector diagnostics.

The `robust_mechanistic_*`, posterior, and legacy LFDR columns are retained for
debugging and ablation only. They no longer decide final TE emission.

## Large BAM Performance Run

For one large BAM, the preferred performance path is one PLACER process with
internal dynamic bin scheduling. The process streams the BAM once,
materializes exact `(tid, bin_index)` owner-bin tasks, and lets worker threads
dynamically consume a bounded cost-aware queue:

```bash
PLACER_PARALLEL=1 \
PLACER_PARALLEL_WORKERS=96 \
PLACER_PARALLEL_QUEUE_MAX_TASKS=384 \
PLACER_PARALLEL_RESULT_BUFFER_MAX=384 \
PLACER_BAM_THREADS=2 \
PLACER_LOG_PARALLEL_PROGRESS=1 \
./build/placer <input.bam> <ref.fa> <te.fa>
```

Short bins finishing early do not strand CPU resources: their workers
immediately take more ready bin tasks. Long genomic regions receive more CPU
time because they generate more bin tasks. A single unusually heavy bin may
still create a tail; if profiling shows that pattern, split that stage further
inside PLACER rather than returning to external process sharding.

Workers may finish out of order, but final results are reduced in
deterministic owner-bin order.

Useful environment variables:

```bash
PLACER_PARALLEL=1 \
PLACER_PARALLEL_WORKERS=8 \
PLACER_PARALLEL_QUEUE_MAX_TASKS=16 \
PLACER_PARALLEL_RESULT_BUFFER_MAX=16 \
PLACER_LOG_PARALLEL_PROGRESS=1 \
./build/placer <input.bam> <ref.fa> <te.fa>
```

Interpretation:

- `[Pipeline] processed=...` reports BAM scan progress.
- `[Pipeline][parallel] ...` reports raw bins discovered, snapshot-window size, executor backlog, running tasks, completed tasks, and reducer-committed bins.
- If scan progress rises while reducer-committed bins stall, the hot path is exact task execution rather than BAM reading.

## Validation Utilities

Two local validation helpers are kept under `scripts/`. They are not cluster
submission entrypoints and write only ignored local outputs under `placer_out/`
unless `--outdir` is provided.

Randomly sample PASS truth loci, build per-locus subset BAMs, run PLACER, and
write `evaluation.tsv` plus `summary.json`:

```bash
python3 scripts/random_truth_interval_eval.py \
  --ground-truth <truth.tsv> \
  --bam <input.bam> \
  --ref <ref.fa> \
  --te <te.fa> \
  --placer ./build/placer \
  --sample-size 10 \
  --seed 20260330 \
  --window-flank-bp 2000 \
  --match-distance-bp 100
```

Randomly sample PASS truth loci and only create read-complete local BAMs for
IGV review:

```bash
python3 scripts/rebuild_igv_validation_bams.py \
  --ground-truth <truth.tsv> \
  --bam <input.bam> \
  --ref <ref.fa> \
  --sample-size 10 \
  --seed 20260330 \
  --window-flank-bp 2000
```

The BAM helper writes `sampled_truth.tsv`, `igv_bams.tsv`, and one
`Sxx/subset.bam` plus index per sampled truth locus.

Summarize reference `N` runs before a production run:

```bash
python3 scripts/poly_n_distribution.py <ref.fa>
```

This writes run-level, genome-wide, and per-contig summaries plus PNG plots
under `placer_out/poly_n_distribution/` by default.

## Tuning

Current pipeline semantics:

- Candidate discovery is only a proposal layer.
- Final output requires event-level structural evidence, event consensus, tripartite segmentation, insert-sequence TE alignment, boundary consistency, and non-reference genotype.
- `scientific.txt` contains the final call table.

Useful environment overrides for the current main path:

```bash
PLACER_TE_KMER_SIZE=13 \
PLACER_TE_KMER_SIZES=9,11,13 \
PLACER_TSD_MIN_LEN=3 \
PLACER_TSD_MAX_LEN=50 \
PLACER_TSD_FLANK_WINDOW=150 \
PLACER_TSD_BG_P_MAX=0.05 \
PLACER_GENOTYPE_MIN_DEPTH=3 \
PLACER_GENOTYPE_ERROR_RATE=0.02 \
PLACER_EVENT_CONSENSUS_POA_MIN_READS=2 \
PLACER_EVENT_CONSENSUS_POA_MAX_READS=48 \
./build/placer <input.bam> <ref.fa> <te.fa>
```

## Notes

- On Apple Silicon, abPOA builds through the top-level CMake SIMDE definitions.
