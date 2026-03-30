#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

status=0

echo "[repo_audit] root: $ROOT_DIR"

if [[ -n "$(git status --short)" ]]; then
    echo "[repo_audit] warning: working tree is not clean"
    git status --short
fi

python3 - <<'PY' || status=$?
import os
import re
import subprocess
import sys
from pathlib import Path

root = Path.cwd()
tracked = subprocess.check_output(["git", "ls-files", "-z"], text=False).split(b"\0")
tracked_paths = [Path(p.decode()) for p in tracked if p]

suspicious_exts = {
    ".bam", ".bai", ".cram", ".crai", ".sam",
    ".vcf", ".bcf", ".bed",
    ".png", ".jpg", ".jpeg", ".pdf", ".tif", ".tiff",
}
suspicious_paths = []
large_paths = []
for path in tracked_paths:
    if path.parts[:1] == ("placer_out",):
        suspicious_paths.append(str(path))
        continue
    if path.suffix.lower() in suspicious_exts:
        suspicious_paths.append(str(path))
    try:
        size = (root / path).stat().st_size
    except FileNotFoundError:
        continue
    if size > 5 * 1024 * 1024:
        large_paths.append((str(path), size))

text_like_suffixes = {
    ".md", ".txt", ".sh", ".py", ".cpp", ".c", ".h", ".hpp",
    ".cmake", ".yml", ".yaml", ".toml", ".json", ".slurm",
}
text_like_names = {"CMakeLists.txt", ".gitignore"}
abs_pattern = re.compile(r"(/Users/|/mnt/|/home/|[A-Za-z]:\\\\)")
absolute_path_hits = []
for path in tracked_paths:
    if path.name not in text_like_names and path.suffix.lower() not in text_like_suffixes:
        continue
    try:
        text = (root / path).read_text(encoding="utf-8")
    except (UnicodeDecodeError, FileNotFoundError):
        continue
    for lineno, line in enumerate(text.splitlines(), start=1):
        if path == Path("scripts/repo_audit.sh") and line.strip().startswith("abs_pattern = re.compile("):
            continue
        if abs_pattern.search(line):
            absolute_path_hits.append(f"{path}:{lineno}:{line.strip()}")

failed = False
if suspicious_paths:
    failed = True
    print("[repo_audit] tracked data-like files that should stay out of the public repo:")
    for item in suspicious_paths:
        print(f"  - {item}")

if large_paths:
    failed = True
    print("[repo_audit] tracked files larger than 5 MiB:")
    for path, size in large_paths:
        print(f"  - {path} ({size} bytes)")

if absolute_path_hits:
    failed = True
    print("[repo_audit] tracked files containing absolute local paths:")
    for item in absolute_path_hits:
        print(f"  - {item}")

if failed:
    sys.exit(1)

print("[repo_audit] passed")
PY

exit "$status"
