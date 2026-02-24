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
ctest --test-dir build --output-on-failure
```

Current CTest includes `test_decision_policy` (runtime decision matrix checks).

## Run

```bash
./build/placer <input.bam> <ref.fa> [te.fa]
```

Default outputs are written to repository root:

- `scientific.txt`
- `ins_fragments.fasta`
- `ins_fragment_hits.tsv`

`scientific.txt` includes a dedicated `insertion_qc` field (separate from `te_qc`).

## Real-Data Tuning (ONT)

Current default TE-classification settings are tuned for noisy long-read data:

- `te_kmer_size = 13`
- `te_vote_fraction_min = 0.40`
- `te_median_identity_min = 0.30`
- `te_rescue_vote_fraction_min = 0.25`
- `te_rescue_median_identity_min = 0.20`

Pipeline behavior:

- Stage A (pre-assembly): weak TE gate to avoid early false negatives.
- Stage B (post-assembly): final TE decision combines vote support and assembly identity.

Useful environment overrides:

```bash
PLACER_TE_KMER_SIZE=13 \
PLACER_TE_VOTE_FRACTION_MIN=0.40 \
PLACER_TE_MEDIAN_IDENTITY_MIN=0.30 \
PLACER_TE_RESCUE_VOTE_FRACTION_MIN=0.25 \
PLACER_TE_RESCUE_MEDIAN_IDENTITY_MIN=0.20 \
./build/placer <input.bam> <ref.fa> <te.fa>
```

## Bootstrap Reclassification (Pass-1 / Pass-2)

When TE reference quality is low, you can keep uncertain calls and bootstrap an incremental TE library.

Pass-1 export controls:

```bash
PLACER_BOOTSTRAP_EXPORT=1 \
PLACER_BOOTSTRAP_EXPORT_INCLUDE_NON_TE=1 \
PLACER_BOOTSTRAP_MIN_CONSENSUS_LEN=80 \
PLACER_BOOTSTRAP_FASTA_PATH=pass1_bootstrap_consensus.fasta \
PLACER_BOOTSTRAP_TSV_PATH=pass1_bootstrap_calls.tsv \
./build/placer <input.bam> <ref.fa> <te.fa>
```

Then run pass-2 with merged base+bootstrap library:

```bash
scripts/bootstrap_te_reclassify.sh \
  --bam <input.bam> \
  --ref <ref.fa> \
  --base-te <te.fa> \
  --pass1-fasta pass1_bootstrap_consensus.fasta \
  --placer ./build/placer \
  --outdir bootstrap_pass2 \
  --min-len 80
```

The script generates:

- `bootstrap_pass2/te_bootstrap_merged.fasta`
- `bootstrap_pass2/scientific_pass2.txt`
- `bootstrap_pass2/bootstrap_report.txt`

Optional low-confidence acceptance (default off):

```bash
PLACER_EMIT_LOW_CONFIDENCE_CALLS=1 \
PLACER_LOW_CONF_MIN_SUPPORT_READS=2 \
PLACER_LOW_CONF_MAX_TIER=2 \
./build/placer <input.bam> <ref.fa> <te.fa>
```

Posterior calibration knobs (optional):

```bash
PLACER_TE_CONF_CERTAIN_MIN=0.85 \
PLACER_TE_CONF_UNCERTAIN_MIN=0.35 \
PLACER_TE_CONF_BIAS=-3.0 \
PLACER_TE_CONF_W_TOP1=2.4 \
PLACER_TE_CONF_W_MARGIN=2.0 \
PLACER_TE_CONF_W_ASM_IDENTITY=1.8 \
PLACER_TE_CONF_W_SUPPORT=0.8 \
./build/placer <input.bam> <ref.fa> <te.fa>
```

`scientific.txt` now reports calibrated TE confidence as `te_conf_prob`.

## Notes

- On Apple Silicon, abPOA builds through the top-level CMake SIMDE definitions.
- Tuned benchmark helper: `scripts/tune_sim_benchmark.sh`.
