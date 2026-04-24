#!/usr/bin/env python3

import importlib.util
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "scripts" / "check_release_readiness.py"

SPEC = importlib.util.spec_from_file_location("check_release_readiness", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Failed to load module spec from {MODULE_PATH}")
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class CheckReleaseReadinessTest(unittest.TestCase):
    def test_generate_unique_seeds_returns_requested_unique_count(self):
        seeds = MODULE.generate_unique_seeds(5)
        self.assertEqual(len(seeds), 5)
        self.assertEqual(len(set(seeds)), 5)
        self.assertTrue(all(isinstance(seed, int) and seed > 0 for seed in seeds))

    def test_evaluate_release_gate_requires_all_rounds_to_clear_strict_threshold(self):
        rounds = [
            MODULE.ReleaseRoundSummary(1, 11, 0.91, True, "/tmp/r1/summary.json"),
            MODULE.ReleaseRoundSummary(2, 22, 0.94, True, "/tmp/r2/summary.json"),
            MODULE.ReleaseRoundSummary(3, 33, 0.97, True, "/tmp/r3/summary.json"),
            MODULE.ReleaseRoundSummary(4, 44, 0.95, True, "/tmp/r4/summary.json"),
            MODULE.ReleaseRoundSummary(5, 55, 0.90, False, "/tmp/r5/summary.json"),
        ]

        verdict = MODULE.evaluate_release_gate(
            rounds,
            round_count=5,
            min_detection_rate=0.90,
        )

        self.assertFalse(verdict["all_rounds_pass"])
        self.assertEqual(verdict["required_rounds"], 5)
        self.assertEqual(verdict["required_detection_rate_strict_gt"], 0.90)
        self.assertEqual(verdict["failed_round_indices"], [5])

    def test_evaluate_release_gate_passes_only_when_every_round_is_strictly_above_threshold(self):
        rounds = [
            MODULE.ReleaseRoundSummary(1, 11, 0.91, True, "/tmp/r1/summary.json"),
            MODULE.ReleaseRoundSummary(2, 22, 0.94, True, "/tmp/r2/summary.json"),
            MODULE.ReleaseRoundSummary(3, 33, 0.97, True, "/tmp/r3/summary.json"),
            MODULE.ReleaseRoundSummary(4, 44, 0.95, True, "/tmp/r4/summary.json"),
            MODULE.ReleaseRoundSummary(5, 55, 0.93, True, "/tmp/r5/summary.json"),
        ]

        verdict = MODULE.evaluate_release_gate(
            rounds,
            round_count=5,
            min_detection_rate=0.90,
        )

        self.assertTrue(verdict["all_rounds_pass"])
        self.assertEqual(verdict["failed_round_indices"], [])


if __name__ == "__main__":
    unittest.main()
