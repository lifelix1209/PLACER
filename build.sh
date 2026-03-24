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

if [[ $# -ge 2 ]]; then
    echo ""
    echo "=== Running PLACER ==="
    exec "$PLACER_BIN" "$@"
fi

echo ""
echo "=== Run skipped ==="
echo "Provide BAM/REF paths to run PLACER:"
echo "  ./build.sh <input.bam> <ref.fa> [te.fa]"
