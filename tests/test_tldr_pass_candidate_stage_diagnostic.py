#!/usr/bin/env python3

import importlib.util
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_module():
    module_path = REPO_ROOT / "scripts" / "tldr_pass_candidate_stage_diagnostic.py"
    spec = importlib.util.spec_from_file_location("tldr_pass_candidate_stage_diagnostic", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class TldrPassCandidateStageDiagnosticTest(unittest.TestCase):
    def test_candidate_class_prefers_te_and_pre_expensive_labels(self):
        module = load_module()

        self.assertEqual(
            module.candidate_class(
                {
                    "final_qc": "PASS_TE_CLOSED|PASS_INSERT_TE_ALIGNMENT",
                    "candidate_retention_reason": "",
                }
            ),
            "pass_te_candidate",
        )
        self.assertEqual(
            module.candidate_class(
                {
                    "final_qc": "NOT_EVALUATED_PRE_EXPENSIVE_STAGE",
                    "candidate_retention_reason": "LEDGER_ONLY_PRE_EXPENSIVE_STAGE",
                }
            ),
            "pre_expensive_stage_candidate",
        )
        self.assertEqual(
            module.candidate_class(
                {
                    "final_qc": "REFERENCE_OR_ARTIFACT",
                    "candidate_retention_reason": "",
                }
            ),
            "reference_or_artifact",
        )

    def test_summarize_pass_rows_counts_candidate_class_presence_once_per_tldr(self):
        module = load_module()
        pass_rows = [
            module.TldrPassRow(index=0, chrom="chr1", midpoint=100.0),
            module.TldrPassRow(index=1, chrom="chr1", midpoint=1000.0),
            module.TldrPassRow(index=2, chrom="chr2", midpoint=500.0),
        ]
        candidate_hits = {
            0: {"pass_te_candidate": 3, "pre_expensive_stage_candidate": 2},
            1: {"pass_te_candidate": 10},
        }
        final_matched = {0}

        rows = module.summarize_pass_rows(pass_rows, candidate_hits, final_matched)

        by_metric = {row["metric"]: row for row in rows}
        self.assertEqual(by_metric["total_tldr_pass"]["all"], 3)
        self.assertEqual(by_metric["final_matched_tldr_pass"]["all"], 1)
        self.assertEqual(by_metric["coverage_candidate_any_class"]["missed_final"], 1)
        self.assertEqual(by_metric["has_pass_te_candidate"]["all"], 2)
        self.assertEqual(by_metric["has_pass_te_candidate"]["missed_final"], 1)
        self.assertEqual(by_metric["has_pre_expensive_stage_candidate"]["final_matched"], 1)

    def test_combination_rows_summarize_class_sets_by_final_status(self):
        module = load_module()
        pass_rows = [
            module.TldrPassRow(index=0, chrom="chr1", midpoint=100.0),
            module.TldrPassRow(index=1, chrom="chr1", midpoint=200.0),
            module.TldrPassRow(index=2, chrom="chr1", midpoint=300.0),
        ]
        candidate_hits = {
            0: {"pass_te_candidate": 1, "reference_or_artifact": 1},
            1: {"reference_or_artifact": 2},
        }
        rows = module.combination_rows(pass_rows, candidate_hits, final_matched={0})

        by_combo = {(row["status"], row["class_combination"]): row for row in rows}
        self.assertEqual(by_combo[("final_matched", "pass_te_candidate+reference_or_artifact")]["count"], 1)
        self.assertEqual(by_combo[("missed_final", "reference_or_artifact")]["count"], 1)
        self.assertEqual(by_combo[("missed_final", "no_coverage_candidate")]["count"], 1)


if __name__ == "__main__":
    unittest.main()
