#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="${PLACER_BUILD_DIR:-$REPO_ROOT/build}"
BUILD_TYPE="${PLACER_CMAKE_BUILD_TYPE:-Release}"
BINARY_PATH="$BUILD_DIR/placer"

detect_jobs() {
    if [[ -n "${PLACER_BUILD_JOBS:-}" ]]; then
        printf '%s\n' "$PLACER_BUILD_JOBS"
        return
    fi
    if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
        printf '%s\n' "$SLURM_CPUS_PER_TASK"
        return
    fi
    if command -v nproc >/dev/null 2>&1; then
        nproc
        return
    fi
    if command -v sysctl >/dev/null 2>&1; then
        sysctl -n hw.ncpu 2>/dev/null || echo 4
        return
    fi
    echo 4
}

find_newer_source() {
    local binary="$1"
    local path=""
    local newer=""
    local -a watch_list=(
        "$REPO_ROOT/CMakeLists.txt"
        "$REPO_ROOT/src"
        "$REPO_ROOT/include"
        "$REPO_ROOT/third_party/abPOA/CMakeLists.txt"
        "$REPO_ROOT/third_party/abPOA/src"
        "$REPO_ROOT/third_party/abPOA/include"
    )

    for path in "${watch_list[@]}"; do
        if [[ ! -e "$path" ]]; then
            continue
        fi
        if [[ -f "$path" ]]; then
            if [[ "$path" -nt "$binary" ]]; then
                printf '%s\n' "$path"
                return 0
            fi
            continue
        fi
        newer="$(find "$path" \
            -type f \
            \( -name '*.cpp' -o -name '*.c' -o -name '*.h' -o -name 'CMakeLists.txt' \) \
            -newer "$binary" \
            -print \
            -quit 2>/dev/null || true)"
        if [[ -n "$newer" ]]; then
            printf '%s\n' "$newer"
            return 0
        fi
    done

    return 1
}

needs_configure=0
needs_build=0
reason=""

if [[ ! -f "$BUILD_DIR/CMakeCache.txt" ]]; then
    needs_configure=1
    needs_build=1
    reason="missing CMake cache"
fi

if [[ ! -x "$BINARY_PATH" ]]; then
    needs_build=1
    if [[ -z "$reason" ]]; then
        reason="missing placer binary"
    fi
fi

if [[ -x "$BINARY_PATH" ]]; then
    if newer_source="$(find_newer_source "$BINARY_PATH")"; then
        needs_build=1
        if [[ -z "$reason" ]]; then
            reason="source newer than binary: $newer_source"
        fi
    fi
fi

mkdir -p "$BUILD_DIR"

if [[ "$needs_configure" -eq 1 ]]; then
    echo "[build-latest] configuring build dir: $BUILD_DIR" >&2
    cmake -S "$REPO_ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
fi

if [[ "$needs_build" -eq 1 ]]; then
    echo "[build-latest] building latest placer ($reason)" >&2
    cmake --build "$BUILD_DIR" -j "$(detect_jobs)"
else
    echo "[build-latest] build is up to date: $BINARY_PATH" >&2
fi

if [[ ! -x "$BINARY_PATH" ]]; then
    echo "[build-latest] placer binary not found after build: $BINARY_PATH" >&2
    exit 2
fi

printf '%s\n' "$BINARY_PATH"
