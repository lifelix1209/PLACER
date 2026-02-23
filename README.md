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

`0.0.1` currently has no registered CTest cases, so `ctest` may report `No tests were found`.

## Run

```bash
./build/placer <input.bam> <ref.fa> [te.fa]
```

Default outputs are written to repository root:

- `scientific.txt`
- `ins_fragments.fasta`
- `ins_fragment_hits.tsv`

## Notes

- On Apple Silicon, abPOA builds through the top-level CMake SIMDE definitions.
- Tuned benchmark helper: `scripts/tune_sim_benchmark.sh`.
