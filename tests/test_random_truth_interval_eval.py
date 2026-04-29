#!/usr/bin/env python3

import importlib.util
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_random_truth_eval_module():
    module_path = REPO_ROOT / "scripts" / "random_truth_interval_eval.py"
    spec = importlib.util.spec_from_file_location("random_truth_interval_eval", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class RandomTruthIntervalEvalTest(unittest.TestCase):
    def test_parse_args_accepts_batch_placer_run(self):
        module = load_random_truth_eval_module()

        args = module.parse_args(
            [
                "--ground-truth",
                "/tmp/truth.tsv",
                "--bam",
                "/tmp/input.bam",
                "--ref",
                "/tmp/ref.fa",
                "--te",
                "/tmp/te.fa",
                "--batch-placer-run",
            ]
        )

        self.assertTrue(args.batch_placer_run)

    def test_parse_batch_placer_stdout(self):
        module = load_random_truth_eval_module()

        statuses = module.parse_batch_placer_stdout(
            "sample_id\texit_code\telapsed_s\n"
            "S01\t0\t1.25\n"
            "S02\t2\t0.5\n"
        )

        self.assertEqual(statuses["S01"].exit_code, 0)
        self.assertEqual(statuses["S01"].elapsed_s, 1.25)
        self.assertEqual(statuses["S02"].exit_code, 2)
        self.assertEqual(statuses["S02"].elapsed_s, 0.5)


if __name__ == "__main__":
    unittest.main()
