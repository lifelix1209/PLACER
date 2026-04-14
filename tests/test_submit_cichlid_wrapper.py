#!/usr/bin/env python3

import json
import os
import stat
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WRAPPER_SCRIPT = REPO_ROOT / "submit_cichlid_yohann_d23.sh"


def make_executable(path: Path) -> None:
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def write_script(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    make_executable(path)


class SubmitCichlidWrapperTest(unittest.TestCase):
    def test_wrapper_accepts_cli_paths_and_submits_tracked_slurm_script(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            repo = Path(tmpdir) / "repo"
            repo.mkdir(parents=True, exist_ok=True)
            wrapper = repo / "submit_cichlid_yohann_d23.sh"
            wrapper.write_text(WRAPPER_SCRIPT.read_text(encoding="utf-8"), encoding="utf-8")
            make_executable(wrapper)

            slurm_script = repo / "scripts" / "submit_placer_urika_d23.slurm"
            write_script(
                slurm_script,
                """#!/usr/bin/env bash
set -euo pipefail
""",
            )

            capture_json = repo / "sbatch.capture.json"
            fake_bin = repo / "fake_bin"
            fake_bin.mkdir(parents=True, exist_ok=True)
            write_script(
                fake_bin / "sbatch",
                f"""#!/usr/bin/env python3
import json
import os
import sys
from pathlib import Path

Path({str(capture_json)!r}).write_text(
    json.dumps(
        {{
            "argv": sys.argv[1:],
            "REPO_ROOT": os.environ.get("REPO_ROOT", ""),
            "PLACER_BAM": os.environ.get("PLACER_BAM", ""),
            "PLACER_REF": os.environ.get("PLACER_REF", ""),
            "PLACER_TE": os.environ.get("PLACER_TE", ""),
            "PLACER_OUT_ROOT": os.environ.get("PLACER_OUT_ROOT", ""),
            "PLACER_RUN_NAME": os.environ.get("PLACER_RUN_NAME", ""),
            "PLACER_SHARD_WORKERS": os.environ.get("PLACER_SHARD_WORKERS", ""),
        }}
    ),
    encoding="utf-8",
)
""",
            )

            data_dir = repo / "inputs"
            data_dir.mkdir(parents=True, exist_ok=True)
            bam = data_dir / "sample.bam"
            ref = data_dir / "ref.fa"
            te = data_dir / "te.fa"
            out_root = data_dir / "placer_out"
            for path in (bam, ref, te):
                path.write_text("", encoding="utf-8")

            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env['PATH']}"

            completed = subprocess.run(
                [
                    "bash",
                    str(wrapper),
                    str(bam),
                    str(ref),
                    str(te),
                    str(out_root),
                    "yohann_d23",
                ],
                cwd=repo,
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, msg=completed.stdout + completed.stderr)
            payload = json.loads(capture_json.read_text(encoding="utf-8"))
            self.assertEqual(payload["REPO_ROOT"], str(repo))
            self.assertEqual(payload["PLACER_BAM"], str(bam))
            self.assertEqual(payload["PLACER_REF"], str(ref))
            self.assertEqual(payload["PLACER_TE"], str(te))
            self.assertEqual(payload["PLACER_OUT_ROOT"], str(out_root))
            self.assertEqual(payload["PLACER_RUN_NAME"], "yohann_d23")
            self.assertEqual(payload["PLACER_SHARD_WORKERS"], "8")
            self.assertEqual(
                payload["argv"][-1],
                str(slurm_script),
            )


if __name__ == "__main__":
    unittest.main()
