#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="${PLACER_BUILD_DIR:-$REPO_ROOT/build}"
BUILD_TYPE="${PLACER_CMAKE_BUILD_TYPE:-Release}"
BINARY_PATH="$BUILD_DIR/placer"

ensure_abpoa_submodule() {
    local abpoa_cmake="$REPO_ROOT/third_party/abPOA/CMakeLists.txt"
    local abpoa_header="$REPO_ROOT/third_party/abPOA/include/abpoa.h"

    if [[ -f "$abpoa_cmake" && -f "$abpoa_header" ]]; then
        return
    fi

    if ! command -v git >/dev/null 2>&1; then
        echo "[build-latest] abPOA submodule is missing and git is not available" >&2
        exit 2
    fi

    echo "[build-latest] initializing abPOA submodule" >&2
    git -C "$REPO_ROOT" submodule update --init --recursive

    if [[ ! -f "$abpoa_cmake" || ! -f "$abpoa_header" ]]; then
        echo "[build-latest] abPOA submodule is still incomplete after initialization" >&2
        exit 2
    fi
}

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

cmake_configure_args=(
    -S "$REPO_ROOT"
    -B "$BUILD_DIR"
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
)

explicit_htslib_config=0
if [[ -n "${HTSLIB_ROOT:-}" ]]; then
    cmake_configure_args+=("-DHTSLIB_ROOT=$HTSLIB_ROOT")
    explicit_htslib_config=1
fi
if [[ -n "${HTSLIB_INCLUDE_DIR:-}" ]]; then
    cmake_configure_args+=("-DHTSLIB_INCLUDE_DIR=$HTSLIB_INCLUDE_DIR")
    explicit_htslib_config=1
fi
if [[ -n "${HTSLIB_LIBRARY:-}" ]]; then
    cmake_configure_args+=("-DHTSLIB_LIBRARY=$HTSLIB_LIBRARY")
    explicit_htslib_config=1
fi

ensure_abpoa_submodule

needs_configure=0
needs_build=0
reason=""

if [[ ! -f "$BUILD_DIR/CMakeCache.txt" ]]; then
    needs_configure=1
    needs_build=1
    reason="missing CMake cache"
fi

if [[ "${PLACER_FORCE_CONFIGURE:-0}" == "1" ]]; then
    needs_configure=1
    needs_build=1
    reason="PLACER_FORCE_CONFIGURE=1"
fi

if [[ "$explicit_htslib_config" -eq 1 ]]; then
    needs_configure=1
    needs_build=1
    if [[ -z "$reason" ]]; then
        reason="explicit htslib configuration"
    fi
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
    cmake "${cmake_configure_args[@]}"
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
