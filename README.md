# PLACER

Current release: `0.0.1`

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

Current CTest includes `test_decision_policy` (runtime decision matrix checks).

## Repository Hygiene

This repository is intended to publish source code only.

- Local outputs such as `build/`, `placer_out/`, shard outputs, bootstrap outputs, and assistant notes are ignored by `.gitignore`.
- Keep real BAM/VCF/IGV screenshots and validation exports outside the tracked tree.
- Before pushing to GitHub, run:

```bash
./scripts/repo_audit.sh
```

The audit script fails if tracked files include large binaries, data-like outputs, or machine-local absolute paths.

## Run

```bash
./build/placer <input.bam> <ref.fa> [te.fa]
```

Default output:

- `scientific.txt`

Optional debug outputs (off by default to reduce I/O on large runs):

```bash
PLACER_INS_FRAGMENTS_FASTA_PATH=ins_fragments.fasta \
PLACER_INS_FRAGMENT_HITS_TSV_PATH=ins_fragment_hits.tsv \
./build/placer <input.bam> <ref.fa> [te.fa]
```

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

For very large BAMs, run by shard and merge final calls without disabling modules:

```bash
python3 scripts/run_sharded_placer.py \
  --bam <input.bam> \
  --ref <ref.fa> \
  --te <te.fa> \
  --placer ./build/placer \
  --mode contig \
  --workers 8 \
  --outdir sharded_placer_out
```

Outputs:

- `sharded_placer_out/scientific.sharded.txt` (merged callset)
- `sharded_placer_out/shard_manifest.tsv` (per-shard runtime/row stats)

Accuracy notes:

- `--mode contig` is recommended and exact (no overlap artifacts).
- `--mode region` uses overlap + core-range merge filtering for boundary safety.

Memory/backpressure tuning for internal parallel mode:

```bash
PLACER_PARALLEL=1 \
PLACER_PARALLEL_WORKERS=8 \
PLACER_PARALLEL_QUEUE_MAX_TASKS=64 \
./build/placer <input.bam> <ref.fa> <te.fa>
```

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
