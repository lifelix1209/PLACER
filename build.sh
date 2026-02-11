#!/bin/bash
# Build script for PLACER Phase 1

set -e

BUILD_DIR=build

echo "=== Building PLACER Phase 1 ==="

if command -v nproc >/dev/null 2>&1; then
    JOBS=$(nproc)
elif command -v sysctl >/dev/null 2>&1; then
    JOBS=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
else
    JOBS=4
fi

# Create build directory
mkdir -p $BUILD_DIR
cd $BUILD_DIR

# Configure with CMake
echo "Configuring with CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
echo "Building..."
cmake --build . -j "$JOBS"

echo "=== Build complete ==="

# Run tests
echo "Running tests..."
ctest --output-on-failure

echo "=== All tests passed ==="

# Optional demo run (user-provided input)
if [ $# -ge 2 ]; then
    echo ""
    echo "=== Running demo ==="
    ./placer "$1" "$2" "${3:-}"
else
    echo ""
    echo "=== Demo skipped ==="
    echo "Provide BAM/REF paths to run demo:"
    echo "  ./build.sh <input.bam> <ref.fa> [te.fa]"
fi
