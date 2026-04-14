#!/usr/bin/env python3

import os
import shutil
import stat
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SUBMIT_SCRIPT = REPO_ROOT / "scripts" / "submit_placer_urika_d23.slurm"
RUN_SCRIPT = REPO_ROOT / "scripts" / "run_placer_latest.sh"
SHARDED_SCRIPT = REPO_ROOT / "scripts" / "run_sharded_placer.py"


def make_executable(path: Path) -> None:
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def write_script(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    make_executable(path)


class SubmitPlacerSlurmTest(unittest.TestCase):
    def make_fake_repo(self, root: Path) -> Path:
        repo = root / "repo"
        (repo / "scripts").mkdir(parents=True, exist_ok=True)
        shutil.copy2(SUBMIT_SCRIPT, repo / "scripts" / "submit_placer_urika_d23.slurm")
        shutil.copy2(RUN_SCRIPT, repo / "scripts" / "run_placer_latest.sh")
        shutil.copy2(SHARDED_SCRIPT, repo / "scripts" / "run_sharded_placer.py")
        make_executable(repo / "scripts" / "submit_placer_urika_d23.slurm")
        make_executable(repo / "scripts" / "run_placer_latest.sh")
        make_executable(repo / "scripts" / "run_sharded_placer.py")
        return repo

    def install_fake_build_chain(self, repo: Path) -> Path:
        stamp_path = repo / "build_helper.invoked"
        build_dir = repo / "build"
        placer_bin = build_dir / "placer"

        write_script(
            repo / "scripts" / "build_latest_placer.sh",
            f"""#!/usr/bin/env bash
set -euo pipefail
mkdir -p "{build_dir}"
printf '1\\n' > "{stamp_path}"
printf '%s\\n' "{placer_bin}"
""",
        )

        write_script(
            placer_bin,
            """#!/usr/bin/env bash
set -euo pipefail
echo "[PLACER] run started" >&2
echo "[PLACER] pipeline finished" >&2
printf 'ok\\n' > scientific.txt
echo "[PLACER] wrote scientific.txt path=$(pwd)/scientific.txt" >&2
""",
        )

        return stamp_path

    def install_fake_sharded_runner(self, repo: Path) -> None:
        write_script(
            repo / "scripts" / "run_sharded_placer.py",
            """#!/usr/bin/env python3
import json
import os
import pathlib
import sys

outdir = pathlib.Path(sys.argv[sys.argv.index("--outdir") + 1])
outdir.mkdir(parents=True, exist_ok=True)
(outdir / "scientific.sharded.txt").write_text("ok\\n", encoding="utf-8")
(outdir / "shard_manifest.tsv").write_text("label\\tchrom\\n", encoding="utf-8")
(outdir / "runner.argv.json").write_text(json.dumps(sys.argv[1:]), encoding="utf-8")
(outdir / "runner.env.json").write_text(
    json.dumps(
        {
            "PLACER_PARALLEL": os.environ.get("PLACER_PARALLEL", ""),
            "PLACER_PARALLEL_WORKERS": os.environ.get("PLACER_PARALLEL_WORKERS", ""),
            "PLACER_BAM_THREADS": os.environ.get("PLACER_BAM_THREADS", ""),
        }
    ),
    encoding="utf-8",
)
print("[sharded] mode=contig contigs=2 shards=2 workers=8 heartbeat=30.0s")
print("[sharded] merged scientific: " + str(outdir / "scientific.sharded.txt"))
print("[sharded] manifest: " + str(outdir / "shard_manifest.tsv"))
""",
        )

    def write_fake_inputs(self, repo: Path) -> dict[str, str]:
        data_dir = repo / "fake_data"
        data_dir.mkdir(parents=True, exist_ok=True)

        bam = data_dir / "input.bam"
        ref = data_dir / "ref.fa"
        te = data_dir / "te.fa"

        bam.write_text("", encoding="utf-8")
        (data_dir / "input.bam.bai").write_text("", encoding="utf-8")
        ref.write_text(">chr1\nACGT\n", encoding="utf-8")
        (data_dir / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
        te.write_text(">TE1\nACGT\n", encoding="utf-8")
        (data_dir / "te.fa.fai").write_text("TE1\t4\t5\t4\t5\n", encoding="utf-8")

        return {
            "PLACER_BAM": str(bam),
            "PLACER_REF": str(ref),
            "PLACER_TE": str(te),
        }

    def test_run_placer_latest_requires_te_argument(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            repo = self.make_fake_repo(Path(tmpdir))
            stamp_path = self.install_fake_build_chain(repo)

            completed = subprocess.run(
                ["bash", str(repo / "scripts" / "run_placer_latest.sh"), "input.bam", "ref.fa"],
                cwd=repo,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 1)
            self.assertIn("Usage:", completed.stderr)
            self.assertFalse(stamp_path.exists())

    def test_submit_script_runs_without_slurm_mempernode_and_emits_log_chain(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            repo = self.make_fake_repo(Path(tmpdir))
            self.install_fake_build_chain(repo)
            self.install_fake_sharded_runner(repo)
            env = os.environ.copy()
            env.update(self.write_fake_inputs(repo))
            env["PLACER_OUT_ROOT"] = str(repo / "placer_out")
            env["SLURM_CPUS_PER_TASK"] = "64"
            env["SLURM_JOB_ID"] = "3240474"
            env.pop("SLURM_MEMPERNODE", None)

            completed = subprocess.run(
                ["bash", str(repo / "scripts" / "submit_placer_urika_d23.slurm")],
                cwd=repo,
                capture_output=True,
                text=True,
                env=env,
                check=False,
            )

            run_dir = repo / "placer_out" / "placer_job3240474"
            stdout_log = run_dir / "job.stdout.log"
            stderr_log = run_dir / "job.stderr.log"

            self.assertEqual(completed.returncode, 0, msg=completed.stdout + completed.stderr)
            self.assertTrue(stdout_log.is_file())
            self.assertTrue(stderr_log.is_file())
            stdout_text = stdout_log.read_text(encoding="utf-8")
            argv_text = (run_dir / "runner.argv.json").read_text(encoding="utf-8")
            env_text = (run_dir / "runner.env.json").read_text(encoding="utf-8")
            self.assertIn("[slurm] slurm_cpus_per_task=64", stdout_text)
            self.assertIn("[sharded] mode=contig", stdout_text)
            self.assertIn("[slurm] PLACER_PARALLEL=disabled", stdout_text)
            self.assertTrue((run_dir / "scientific.sharded.txt").is_file())
            self.assertTrue((run_dir / "shard_manifest.tsv").is_file())
            self.assertNotIn("PLACER_PARALLEL=1", argv_text)
            self.assertIn('"PLACER_PARALLEL": ""', env_text)
            self.assertIn('"PLACER_PARALLEL_WORKERS": ""', env_text)
            self.assertIn(
                f"[slurm] scientific_sharded_txt={run_dir / 'scientific.sharded.txt'}",
                stdout_text,
            )
            self.assertIn(
                f"[slurm] shard_manifest_tsv={run_dir / 'shard_manifest.tsv'}",
                stdout_text,
            )


if __name__ == "__main__":
    unittest.main()
