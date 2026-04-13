# PLACER

Current release: `0.0.3`

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

```bash
python3 -m pip install -r requirements.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Test

```bash
(cd build && ctest --output-on-failure)
```

Current CTest covers the C++ pipeline units plus Python regression tests for shard merging and random-truth evaluation.

## Repository Hygiene

This repository is intended to publish source code only.

- Local outputs such as `build/`, `placer_out/`, shard outputs, bootstrap outputs, and assistant notes are ignored by `.gitignore`.
- `test_data/` is reserved for small reproducible fixtures. Large local references and TE libraries should stay outside the tracked set and are ignored when kept under their current local filenames.
- Keep real BAM/VCF/IGV screenshots and validation exports outside the tracked tree.
- Before pushing to GitHub, run:

```bash
./scripts/repo_audit.sh
```

The audit script fails if tracked files include large binaries, data-like outputs, or machine-local absolute paths.

## Run

```bash
./build/placer <input.bam> <ref.fa> <te.fa>
```

Default output:

- `scientific.txt`

Optional debug outputs (off by default to reduce I/O on large runs):

```bash
PLACER_INS_FRAGMENTS_FASTA_PATH=ins_fragments.fasta \
PLACER_INS_FRAGMENT_HITS_TSV_PATH=ins_fragment_hits.tsv \
./build/placer <input.bam> <ref.fa> <te.fa>
```

## SLURM Run

Use the repository submission script as the only supported cluster entrypoint.

```bash
PLACER_BAM=/data/input.bam
PLACER_REF=/data/ref.fa
PLACER_TE=/data/te.fa
PLACER_OUT_ROOT=/data/placer_out
PLACER_RUN_NAME=yohann_d23

sbatch \
  --cpus-per-task 64 \
  --mem 192G \
  --export=ALL,PLACER_BAM,PLACER_REF,PLACER_TE,PLACER_OUT_ROOT,PLACER_RUN_NAME \
  scripts/submit_placer_urika_d23.slurm
```

Required environment variables:

- `PLACER_BAM`
- `PLACER_REF`
- `PLACER_TE`

Optional overrides:

- `PLACER_OUT_ROOT`
- `PLACER_RUN_NAME`
- `PLACER_BAM_BAI`
- `PLACER_REF_FAI`
- `PLACER_TE_FAI`

The submission script writes run logs into:

- `<PLACER_OUT_ROOT>/<PLACER_RUN_NAME>/job.stdout.log`
- `<PLACER_OUT_ROOT>/<PLACER_RUN_NAME>/job.stderr.log`

Expected log chain:

- `[slurm] ...`
- `[run-latest] ...`
- `[PLACER] run started`

If `[PLACER] run started` never appears in the run logs, the job did not reach the current repository binary entrypoint.

`scientific.txt` now contains only final PASS TE insertion calls. The row schema is:

- `chrom`, `pos`, `bp_left`, `bp_right`
- `te`, `family`, `subfamily`, `strand`, `insert_len`
- `support_reads`, `alt_struct_reads`, `ref_span_reads`, `low_mapq_ref_span_reads`
- `gt`, `af`, `gq`
- `best_te_identity`, `best_te_query_coverage`, `cross_family_margin`
- `tsd_type`, `tsd_len`
- `left_flank_align_len`, `right_flank_align_len`, `consensus_len`
- `qc`

## Sharded Run (Large BAM)

For very large BAMs on a single node, use exact contig sharding as the primary path:

```bash
python3 scripts/run_sharded_placer.py \
  --bam <input.bam> \
  --ref <ref.fa> \
  --te <te.fa> \
  --placer ./build/placer \
  --mode contig \
  --workers 32 \
  --resume \
  --outdir sharded_placer_out
```

Outputs:

- `sharded_placer_out/scientific.sharded.txt` (canonical merged exact callset)
- `sharded_placer_out/shard_manifest.tsv` (live shard state + final per-shard timing)
- `sharded_placer_out/shards/<label>/scientific.txt` (per-contig shard result)

Operational notes:

- `--mode contig` is the exact production path for large BAMs.
- Restarting the same `--outdir` with `--resume` reuses completed shard results.
- Large-BAM SLURM runs should launch the sharded orchestrator directly rather than `run_placer_latest.sh`.

Accuracy notes:

- `--mode contig` is exact (no overlap artifacts).
- `--mode region` uses overlap + core-range merge filtering for boundary safety.

Memory/backpressure tuning for internal parallel mode:

```bash
PLACER_PARALLEL=1 \
PLACER_PARALLEL_WORKERS=8 \
PLACER_PARALLEL_QUEUE_MAX_TASKS=64 \
./build/placer <input.bam> <ref.fa> <te.fa>
```

## Random Truth Evaluation

Run a deterministic local cohort evaluation against a truth table:

```bash
python3 scripts/random_truth_interval_eval.py \
  --ground-truth <truth.tsv> \
  --bam <input.bam> \
  --ref <ref.fa> \
  --te <te.fa> \
  --placer ./build/placer \
  --sample-size 10 \
  --seed 20260330 \
  --threads 1 \
  --workers 1 \
  --outdir placer_out/random_truth_eval_serial
```

Replay the same sampled cohort in parallel:

```bash
python3 scripts/random_truth_interval_eval.py \
  --ground-truth <truth.tsv> \
  --bam <input.bam> \
  --ref <ref.fa> \
  --te <te.fa> \
  --placer ./build/placer \
  --sampled-truth-tsv placer_out/random_truth_eval_serial/sampled_truth.tsv \
  --seed 20260330 \
  --threads 1 \
  --workers 4 \
  --outdir placer_out/random_truth_eval_parallel
```

`evaluation.tsv` records per-sample detection calls, timing fields, and parsed joint-decision diagnostics from `run.log`.

## Tuning

Current pipeline semantics:

- Candidate discovery is only a proposal layer.
- Final output requires event-level structural evidence, event consensus, tripartite segmentation, insert-sequence TE alignment, boundary consistency, and non-reference genotype.
- `scientific.txt` and `scientific.sharded.txt` share the same call columns.

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
