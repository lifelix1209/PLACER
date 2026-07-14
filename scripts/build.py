#!/usr/bin/env python3
"""Cross-platform user-space build entrypoint for PLACER."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]


def default_jobs() -> str:
    return str(os.cpu_count() or 4)


def env_flag(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() not in {"0", "false", "no", "off"}


def resolve_path(value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (REPO_ROOT / path).resolve()


def run_command(command: list[str], *, env: dict[str, str] | None = None) -> None:
    printable = " ".join(command)
    print(f"+ {printable}", flush=True)
    subprocess.run(command, cwd=REPO_ROOT, env=env, check=True)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "One-click install/build for PLACER. This cross-platform script "
            "runs entirely in user space. No administrator privileges are required."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--build-dir",
        default=os.environ.get("PLACER_BUILD_DIR", str(REPO_ROOT / "build")),
        help="CMake build directory. Default: build/",
    )
    parser.add_argument(
        "--build-type",
        default=os.environ.get("PLACER_CMAKE_BUILD_TYPE", "Release"),
        help="CMake build type/configuration. Default: Release",
    )
    parser.add_argument(
        "--jobs",
        default=os.environ.get("PLACER_BUILD_JOBS", default_jobs()),
        help="Parallel build jobs. Default: detected CPU count",
    )
    parser.add_argument(
        "--no-tests",
        action="store_true",
        default=not env_flag("PLACER_RUN_TESTS", True),
        help="Build only; skip CTest.",
    )
    parser.add_argument(
        "--no-venv",
        action="store_true",
        default=not env_flag("PLACER_USE_VENV", True),
        help="Install Python deps into the active Python env.",
    )
    parser.add_argument(
        "--skip-py-deps",
        action="store_true",
        default=env_flag("PLACER_SKIP_PY_DEPS", False),
        help="Do not install Python dependencies.",
    )
    parser.add_argument(
        "--venv-dir",
        default=os.environ.get("PLACER_VENV_DIR"),
        help="Python virtual environment path. Default: <build-dir>/.venv",
    )
    parser.add_argument(
        "run_args",
        nargs="*",
        metavar="PLACER_ARG",
        help="Optional PLACER run arguments: <input.bam> <ref.fa> <te.fa> [extra args...]",
    )
    parser.epilog = """Environment overrides:
  PLACER_BUILD_DIR
  PLACER_CMAKE_BUILD_TYPE
  PLACER_BUILD_JOBS
  PLACER_RUN_TESTS=0
  PLACER_USE_VENV=0
  PLACER_SKIP_PY_DEPS=1
  PLACER_VENV_DIR
  HTSLIB_ROOT, HTSLIB_INCLUDE_DIR, HTSLIB_LIBRARY

By default this script uses an existing abPOA checkout or initializes
submodules if needed, creates build/.venv, installs requirements.txt there,
configures CMake into build/, compiles build/placer, and runs CTest.
"""
    return parser.parse_args(argv)


def prepare_source_tree() -> None:
    print("=== Preparing PLACER source tree ===", flush=True)
    abpoa_cmake = REPO_ROOT / "third_party" / "abPOA" / "CMakeLists.txt"
    if abpoa_cmake.is_file():
        print("abPOA submodule is already present: third_party/abPOA", flush=True)
        return

    if (REPO_ROOT / ".git").exists():
        print("Initializing git submodules...", flush=True)
        run_command(["git", "-C", str(REPO_ROOT), "submodule", "update", "--init", "--recursive"])

    if not abpoa_cmake.is_file():
        raise SystemExit(
            "[build] missing required abPOA checkout at third_party/abPOA\n"
            "[build] clone with --recursive, or run git submodule update --init --recursive "
            "in a git checkout."
        )


def install_python_deps(args: argparse.Namespace, build_dir: Path) -> str:
    python_bin = sys.executable
    if args.skip_py_deps:
        print(
            "Skipping Python dependency install because PLACER_SKIP_PY_DEPS=1 "
            "or --skip-py-deps was set",
            flush=True,
        )
        return python_bin

    print("=== Installing Python dependencies ===", flush=True)
    if args.no_venv:
        print(f"Using active Python environment: {python_bin}", flush=True)
        run_command([python_bin, "-m", "pip", "install", "-r", str(REPO_ROOT / "requirements.txt")])
        return python_bin

    venv_dir = resolve_path(args.venv_dir) if args.venv_dir else build_dir / ".venv"
    run_command([python_bin, "-m", "venv", str(venv_dir)])
    python_bin = str(venv_dir / ("Scripts/python.exe" if os.name == "nt" else "bin/python"))
    print(f"Using Python virtual environment: {venv_dir}", flush=True)
    run_command([python_bin, "-m", "pip", "install", "-r", str(REPO_ROOT / "requirements.txt")])
    return python_bin


def build_placer(args: argparse.Namespace, build_dir: Path, python_bin: str) -> Path:
    print("=== Building PLACER ===", flush=True)
    run_command(
        [
            "cmake",
            "-S",
            str(REPO_ROOT),
            "-B",
            str(build_dir),
            f"-DCMAKE_BUILD_TYPE={args.build_type}",
            f"-DPLACER_PYTHON_EXECUTABLE={python_bin}",
        ]
    )
    run_command(
        [
            "cmake",
            "--build",
            str(build_dir),
            "--config",
            args.build_type,
            "--parallel",
            str(args.jobs),
        ]
    )
    print("=== Build complete ===", flush=True)
    suffix = ".exe" if os.name == "nt" else ""
    config_bin = build_dir / args.build_type / f"placer{suffix}"
    return config_bin if config_bin.exists() else build_dir / f"placer{suffix}"


def run_tests(args: argparse.Namespace, build_dir: Path) -> None:
    if args.no_tests:
        print("Skipping tests because PLACER_RUN_TESTS=0 or --no-tests was set", flush=True)
        return
    print("Running tests...", flush=True)
    run_command(["ctest", "--test-dir", str(build_dir), "--build-config", args.build_type, "--output-on-failure"])
    print("=== All tests passed ===", flush=True)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if 0 < len(args.run_args) < 3:
        raise SystemExit("[build] missing run arguments\nUsage: build.py [options] <input.bam> <ref.fa> <te.fa>")

    build_dir = resolve_path(args.build_dir)
    prepare_source_tree()
    python_bin = install_python_deps(args, build_dir)
    placer_bin = build_placer(args, build_dir, python_bin)
    run_tests(args, build_dir)

    if len(args.run_args) >= 3:
        print("\n=== Running PLACER ===", flush=True)
        run_command([str(placer_bin), *args.run_args])
        return 0

    print(
        "\n=== Run skipped ===\n"
        "Provide BAM/REF paths to run PLACER:\n"
        "  python scripts/build.py <input.bam> <ref.fa> <te.fa>\n"
        "\n"
        "You can override htslib discovery explicitly with:\n"
        "  HTSLIB_ROOT=/path/to/conda python scripts/build.py\n"
        "  HTSLIB_INCLUDE_DIR=/path/to/conda/include "
        "HTSLIB_LIBRARY=/path/to/conda/lib/libhts.so python scripts/build.py",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
