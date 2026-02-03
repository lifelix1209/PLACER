#!/bin/bash
# Build script for PLACER Phase 1

set -e

BUILD_DIR=build

echo "=== Building PLACER Phase 1 ==="

# Create build directory
mkdir -p $BUILD_DIR
cd $BUILD_DIR

# Configure with CMake
echo "Configuring with CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
echo "Building..."
make -j$(nproc)

echo "=== Build complete ==="

# Run tests
echo "Running tests..."
ctest --output-on-failure

echo "=== All tests passed ==="

# Run demo
echo ""
echo "=== Running demo on test data ==="
./placer /mnt/home1/miska/hl725/projects/tldr_optimized/test/test.bam
