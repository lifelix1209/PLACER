#!/usr/bin/env python3

import os
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class ClusterScriptTest(unittest.TestCase):
    def test_build_latest_forwards_explicit_htslib_paths_to_cmake(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fake_bin = tmp / "bin"
            fake_bin.mkdir()
            calls_path = tmp / "cmake.calls"
            fake_cmake = fake_bin / "cmake"
            fake_cmake.write_text(
                textwrap.dedent(
                    """\
                    #!/usr/bin/env bash
                    set -euo pipefail
                    printf '%s\\0' "$@" >> "$CALLS_PATH"
                    printf '\\n' >> "$CALLS_PATH"
                    if [[ "${1:-}" == "-S" ]]; then
                        build_dir=""
                        while [[ $# -gt 0 ]]; do
                            if [[ "$1" == "-B" ]]; then
                                shift
                                build_dir="$1"
                                break
                            fi
                            shift
                        done
                        mkdir -p "$build_dir"
                        : > "$build_dir/CMakeCache.txt"
                        exit 0
                    fi
                    if [[ "${1:-}" == "--build" ]]; then
                        build_dir="$2"
                        mkdir -p "$build_dir"
                        printf '#!/usr/bin/env bash\\nexit 0\\n' > "$build_dir/placer"
                        chmod +x "$build_dir/placer"
                        exit 0
                    fi
                    exit 2
                    """
                ),
                encoding="utf-8",
            )
            fake_cmake.chmod(0o755)

            env = os.environ.copy()
            env.update(
                {
                    "PATH": f"{fake_bin}{os.pathsep}{env.get('PATH', '')}",
                    "CALLS_PATH": str(calls_path),
                    "PLACER_BUILD_DIR": str(tmp / "build"),
                    "PLACER_BUILD_JOBS": "5",
                    "HTSLIB_ROOT": "/opt/conda",
                    "HTSLIB_INCLUDE_DIR": "/opt/conda/include",
                    "HTSLIB_LIBRARY": "/opt/conda/lib/libhts.so",
                }
            )

            cp = subprocess.run(
                [str(REPO_ROOT / "scripts" / "build_latest_placer.sh")],
                cwd=REPO_ROOT,
                env=env,
                text=True,
                capture_output=True,
                check=True,
            )

            self.assertEqual(cp.stdout.strip(), str(tmp / "build" / "placer"))
            calls = calls_path.read_text(encoding="utf-8").replace("\0", " ")
            self.assertIn("-DHTSLIB_ROOT=/opt/conda", calls)
            self.assertIn("-DHTSLIB_INCLUDE_DIR=/opt/conda/include", calls)
            self.assertIn("-DHTSLIB_LIBRARY=/opt/conda/lib/libhts.so", calls)

    def test_cichlid_submit_wrapper_uses_repo_root_script_path(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fake_bin = tmp / "bin"
            fake_bin.mkdir()
            sbatch_args = tmp / "sbatch.args"
            fake_sbatch = fake_bin / "sbatch"
            fake_sbatch.write_text(
                textwrap.dedent(
                    """\
                    #!/usr/bin/env bash
                    set -euo pipefail
                    printf '%s\\n' "$@" > "$SBATCH_ARGS"
                    echo "Submitted batch job 123"
                    """
                ),
                encoding="utf-8",
            )
            fake_sbatch.chmod(0o755)

            out_root = tmp / "out"
            env = os.environ.copy()
            env.update(
                {
                    "PATH": f"{fake_bin}{os.pathsep}{env.get('PATH', '')}",
                    "SBATCH_ARGS": str(sbatch_args),
                }
            )

            cp = subprocess.run(
                [
                    "bash",
                    str(REPO_ROOT / "scripts" / "submit_cichlid_yohann_d23.sh"),
                    "/data/input.bam",
                    "/data/ref.fa",
                    "/data/te.fa",
                    str(out_root),
                    "run1",
                ],
                cwd=REPO_ROOT,
                env=env,
                text=True,
                capture_output=True,
                check=True,
            )

            args = sbatch_args.read_text(encoding="utf-8").splitlines()
            self.assertIn(str(REPO_ROOT / "scripts" / "submit_placer_urika_d23.slurm"), args)
            self.assertNotIn("scripts/scripts", "\n".join(args))
            self.assertIn("[submit] TASKGRAPH_WORKERS:", cp.stdout)
            self.assertIn("[submit] TASKGRAPH_QUEUE_MAX_TASKS:", cp.stdout)
            self.assertNotIn("SHARD_WORKERS", cp.stdout)


if __name__ == "__main__":
    unittest.main()
