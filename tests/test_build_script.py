from pathlib import Path
import os
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class BuildScriptTests(unittest.TestCase):
    def test_python_entrypoint_documents_cross_platform_user_space_build(self):
        result = subprocess.run(
            ["python3", str(ROOT / "scripts" / "build.py"), "--help"],
            check=True,
            cwd=ROOT,
            text=True,
            capture_output=True,
        )

        self.assertIn("One-click install/build", result.stdout)
        self.assertIn("cross-platform", result.stdout)
        self.assertIn("No administrator privileges", result.stdout)
        self.assertIn("--no-venv", result.stdout)
        self.assertIn("--no-tests", result.stdout)
        self.assertIn("PLACER_SKIP_PY_DEPS", result.stdout)
        self.assertIn("build/.venv", result.stdout)

    def test_shell_and_powershell_wrappers_delegate_to_python_entrypoint(self):
        shell_text = (ROOT / "build.sh").read_text()
        powershell_text = (ROOT / "build.ps1").read_text()

        self.assertIn("scripts/build.py", shell_text)
        self.assertIn("scripts/build.py", powershell_text)

    def test_build_helpers_avoid_sudo_or_system_package_installers(self):
        checked_paths = [
            ROOT / "build.sh",
            ROOT / "build.ps1",
            ROOT / "scripts" / "build.py",
            ROOT / "README.md",
        ]

        combined = "\n".join(path.read_text().lower() for path in checked_paths)
        forbidden_tokens = ["sudo", "apt-get", "apt ", "brew install", "yum ", "dnf ", "pacman "]
        for token in forbidden_tokens:
            self.assertNotIn(token, combined)

    def test_existing_abpoa_checkout_skips_git_submodule_update(self):
        self.assertTrue((ROOT / "third_party" / "abPOA" / "CMakeLists.txt").exists())

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fake_bin = tmp_path / "bin"
            fake_bin.mkdir()
            build_dir = tmp_path / "build"

            git_script = fake_bin / "git"
            git_script.write_text(
                "#!/usr/bin/env bash\n"
                "echo 'git should not be called when abPOA is present' >&2\n"
                "exit 13\n",
            )
            cmake_script = fake_bin / "cmake"
            cmake_script.write_text(
                "#!/usr/bin/env bash\n"
                "printf '%s\\n' \"$@\" >> \"$PLACER_FAKE_CMAKE_LOG\"\n"
                "exit 0\n",
            )
            git_script.chmod(0o755)
            cmake_script.chmod(0o755)
            cmake_log = tmp_path / "cmake.log"
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env['PATH']}"
            env["PLACER_FAKE_CMAKE_LOG"] = str(cmake_log)

            result = subprocess.run(
                [
                    "python3",
                    str(ROOT / "scripts" / "build.py"),
                    "--skip-py-deps",
                    "--no-tests",
                    "--jobs",
                    "1",
                    "--build-dir",
                    str(build_dir),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                env=env,
            )
            cmake_output = cmake_log.read_text()

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("abPOA submodule is already present", result.stdout)
        self.assertIn("--parallel", cmake_output)


if __name__ == "__main__":
    unittest.main()
