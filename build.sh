#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"
BUILD_HELPER="$REPO_ROOT/scripts/build_latest_placer.sh"

if [[ ! -x "$BUILD_HELPER" ]]; then
    echo "[build] build helper not executable: $BUILD_HELPER" >&2
    exit 2
fi

echo "=== Building PLACER ==="
PLACER_BIN="$("$BUILD_HELPER")"
echo "=== Build complete ==="

if [[ "${PLACER_RUN_TESTS:-1}" != "0" ]]; then
    echo "Running tests..."
    ctest --test-dir "$REPO_ROOT/build" --output-on-failure
    echo "=== All tests passed ==="
else
    echo "Skipping tests because PLACER_RUN_TESTS=0"
fi

if [[ $# -gt 0 && $# -lt 3 ]]; then
    echo "[build] missing run arguments" >&2
    echo "Usage: ./build.sh <input.bam> <ref.fa> <te.fa>" >&2
    exit 2
fi

if [[ $# -ge 3 ]]; then
    echo ""
    echo "=== Running PLACER ==="
    exec "$PLACER_BIN" "$@"
fi

echo ""
echo "=== Run skipped ==="
echo "Provide BAM/REF paths to run PLACER:"
echo "  ./build.sh <input.bam> <ref.fa> <te.fa>"
echo ""
echo "For an explicit conda htslib build, run for example:"
echo "  HTSLIB_ROOT=/path/to/conda ./build.sh"
echo "  HTSLIB_INCLUDE_DIR=/path/to/conda/include HTSLIB_LIBRARY=/path/to/conda/lib/libhts.so ./build.sh"
