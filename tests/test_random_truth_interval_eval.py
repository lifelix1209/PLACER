#!/usr/bin/env python3

import importlib.util
import sys
import tempfile
import unittest
from argparse import Namespace
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
        self.assertTrue(module.use_pipeline_placer_runner(args))

    def test_default_eval_uses_pipeline_placer_runner(self):
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
            ]
        )

        self.assertTrue(module.use_pipeline_placer_runner(args))

    def test_evaluate_sample_rows_runs_pipeline_placer_once_by_default(self):
        module = load_random_truth_eval_module()

        row = module.TruthRow(
            sample_id="S01",
            uuid="u1",
            chrom="chr1",
            start=100,
            end=105,
            strand="+",
            family="L1",
            subfamily="L1-1",
            filter_value="PASS",
            raw={},
            length_ins=100,
        )
        args = Namespace(
            reuse_existing_samples=False,
            window_flank_bp=10,
            workers=1,
        )
        calls = []
        prepared = {
            "S01": module.PreparedSubset(
                sample_id="S01",
                sample_dir=Path("/tmp/out/S01"),
                log_path=Path("/tmp/out/S01/run.log"),
                read_names_path=Path("/tmp/out/S01/window_qnames.txt"),
                subset_bam=Path("/tmp/out/S01/subset.bam"),
                qname_count=2,
                extract_elapsed_s=0.1,
            )
        }

        original_prepare = module.prepare_sample_subset_bams
        original_run_pipeline = module.run_pipeline_placer_samples
        original_evaluate = module.evaluate_pipeline_sample
        try:
            def fake_prepare_sample_subset_bams(**kwargs):
                calls.append("prepare")
                return prepared

            def fake_run_pipeline_placer_samples(**kwargs):
                calls.append(("pipeline", sorted(kwargs["prepared_subsets"])))
                return {
                    "S01": module.PipelinePlacerStatus(
                        sample_id="S01",
                        exit_code=0,
                        elapsed_s=0.2,
                    )
                }

            def fake_evaluate_pipeline_sample(**kwargs):
                calls.append(("evaluate", kwargs["row"].sample_id))
                return {"sample_id": kwargs["row"].sample_id}

            module.prepare_sample_subset_bams = fake_prepare_sample_subset_bams
            module.run_pipeline_placer_samples = fake_run_pipeline_placer_samples
            module.evaluate_pipeline_sample = fake_evaluate_pipeline_sample

            _, results = module.evaluate_sample_rows(
                sampled_rows=[row],
                contig_lengths={"chr1": 1000},
                args=args,
                outdir=Path("/tmp/out"),
            )
        finally:
            module.prepare_sample_subset_bams = original_prepare
            module.run_pipeline_placer_samples = original_run_pipeline
            module.evaluate_pipeline_sample = original_evaluate

        self.assertEqual(calls, ["prepare", ("pipeline", ["S01"]), ("evaluate", "S01")])
        self.assertEqual(results, [{"sample_id": "S01"}])

    def test_parse_args_defaults_match_distance_to_100bp(self):
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
            ]
        )

        self.assertEqual(args.match_distance_bp, 100)

    def test_parse_args_accepts_placer_raw_cigar_insert_threshold(self):
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
                "--placer-min-final-raw-cigar-insert-len-bp",
                "80",
            ]
        )

        self.assertEqual(args.placer_min_final_raw_cigar_insert_len_bp, 80)

    def test_build_placer_env_sets_raw_cigar_insert_threshold(self):
        module = load_random_truth_eval_module()

        args = Namespace(placer_min_final_raw_cigar_insert_len_bp=80)

        env = module.build_placer_env(args)

        self.assertEqual(env["PLACER_MIN_FINAL_RAW_CIGAR_INSERT_LEN_BP"], "80")

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

    def test_load_sampled_truth_rows_accepts_igv_sampled_calls_tsv(self):
        module = load_random_truth_eval_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            sampled_calls = Path(tmpdir) / "sampled_calls.tsv"
            sampled_calls.write_text(
                "chrom\twindow_start0\twindow_end\tname\tscore\tstrand\tlabel\t"
                "subset\tsource\tcandidate_start\tcandidate_end\tsource_id\t"
                "family\tsubfamily\tinsert_len\ttldr_filter\n"
                "chr1\t900\t1200\tT_ONLY_01|chr1:1000-1004|TLDR\t0\t+\t"
                "T_ONLY_01\tTLDR_ONLY\tTLDR\t1000\t1004\tu1\tL1\tL1-1\t591\tPASS\n",
                encoding="utf-8",
            )

            rows = module.load_sampled_truth_rows(sampled_calls)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].sample_id, "T_ONLY_01")
        self.assertEqual(rows[0].uuid, "u1")
        self.assertEqual(rows[0].chrom, "chr1")
        self.assertEqual(rows[0].start, 1000)
        self.assertEqual(rows[0].end, 1004)
        self.assertEqual(rows[0].filter_value, "PASS")
        self.assertEqual(rows[0].length_ins, 591)
        self.assertEqual(rows[0].raw["subset"], "TLDR_ONLY")

    def test_load_truth_rows_preserves_length_ins_without_min_filter(self):
        module = load_random_truth_eval_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            truth = Path(tmpdir) / "truth.tsv"
            truth.write_text(
                "UUID\tChrom\tStart\tEnd\tStrand\tFamily\tSubfamily\t"
                "LengthIns\tFilter\n"
                "u1\tchr1\t100\t101\t+\tL1\tL1-1\t591\tPASS\n",
                encoding="utf-8",
            )

            rows = module.load_truth_rows(truth, "PASS")

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].length_ins, 591)

    def test_extract_region_qnames_keeps_supplementary_and_secondary_window_evidence(self):
        module = load_random_truth_eval_module()

        region_stdout = (
            "primary\t0\tchr1\t100\t60\t100M\t*\t0\t0\tACGT\tIIII\n"
            "secondary\t256\tchr1\t120\t0\t100M\t*\t0\t0\tACGT\tIIII\n"
            "supplementary\t2048\tchr1\t130\t60\t100M\t*\t0\t0\tACGT\tIIII\n"
            "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"
            "primary\t0\tchr1\t140\t60\t100M\t*\t0\t0\tACGT\tIIII\n"
        )

        qnames = module.extract_region_qnames_from_region_view(
            region_stdout,
            "chr1:100-200",
        )

        self.assertEqual(qnames, ["primary", "secondary", "supplementary"])

    def test_evidence_ledger_interval_prefers_coverage_envelope(self):
        module = load_random_truth_eval_module()

        row = {
            "bp_left": "100",
            "bp_right": "110",
            "coverage_left": "50",
            "coverage_right": "250",
        }

        self.assertEqual(module.evidence_ledger_interval(row), (50, 250))

    def test_build_result_row_captures_robust_fields_and_same_event_ledger(self):
        module = load_random_truth_eval_module()
        truth = module.TruthRow(
            sample_id="S01",
            uuid="u1",
            chrom="chr1",
            start=1000,
            end=1004,
            strand="+",
            family="L1",
            subfamily="L1-1",
            filter_value="PASS",
            raw={},
            length_ins=591,
        )
        window = module.WindowSpec(
            chrom="chr1",
            truth_start=1000,
            truth_end=1004,
            window_start=900,
            window_end=1200,
        )
        best_call = {
            "chrom": "chr1",
            "interval_distance_bp": "0",
            "pos": "1002",
            "bp_left": "1000",
            "bp_right": "1004",
            "te": "L1-1",
            "family": "L1",
            "subfamily": "L1-1",
            "gt": "0/1",
            "af": "0.50",
            "gq": "60",
            "qc": "PASS_TE_CLOSED",
            "robust_mechanistic_qc": "PASS_TE_LFDR",
            "robust_mechanistic_worst_case_lfdr": "0.04",
            "mechanistic_lower_log_bf_te_vs_artifact": "7.0",
            "mechanistic_lower_log_bf_te_vs_non_te": "6.5",
            "conformal_null_p": "0.01",
            "conformal_by_threshold": "0.02",
            "conformal_qc": "PASS_CONFORMAL_FDR",
        }
        evidence_ledger_rows = [
            {
                "chrom": "chr1",
                "bp_left": "998",
                "bp_right": "1005",
                "final_qc": "TE_LFDR_HIGH",
                "robust_mechanistic_qc": "TE_LFDR_HIGH",
                "robust_mechanistic_worst_case_lfdr": "0.42",
            },
            {
                "chrom": "chr1",
                "bp_left": "5000",
                "bp_right": "5010",
                "final_qc": "TE_LFDR_HIGH",
                "robust_mechanistic_qc": "TE_LFDR_HIGH",
            },
        ]

        remote_call = dict(best_call)
        remote_call.update(
            {
                "chrom": "chr2",
                "pos": "5000",
                "bp_left": "5000",
                "bp_right": "5000",
                "robust_mechanistic_worst_case_lfdr": "0.20",
            }
        )
        same_window_extra_call = dict(best_call)
        same_window_extra_call.update(
            {
                "pos": "1150",
                "bp_left": "1150",
                "bp_right": "1150",
                "robust_mechanistic_worst_case_lfdr": "0.08",
            }
        )

        result = module.build_result_row(
            row=truth,
            window=window,
            sample_dir=Path("/tmp/S01"),
            read_names_path=Path("/tmp/S01/window_qnames.txt"),
            subset_bam=Path("/tmp/S01/subset.bam"),
            scientific_path=Path("/tmp/S01/scientific.txt"),
            qname_count=4,
            placer_exit_code="0",
            calls=[best_call, remote_call, same_window_extra_call],
            best_call=best_call,
            detected=True,
            extract_elapsed_s=0.1,
            placer_elapsed_s=0.2,
            sample_elapsed_s=0.3,
            joint_metrics={},
            scientific_summary={"event_consensus_calls": 1, "genotype_calls": 1, "final_pass_calls": 1},
            evidence_ledger_path=Path("/tmp/S01/evidence_ledger.tsv"),
            evidence_ledger_rows=evidence_ledger_rows,
            null_controls_path=Path("/tmp/S01/null_controls.tsv"),
            null_control_rows=[
                {
                    "row_type": "null_control",
                    "empirical_null_upper_tail_p": "0.80",
                },
                {
                    "row_type": "final_call_tail_query",
                    "empirical_null_upper_tail_p": "0.05",
                },
            ],
        )

        self.assertEqual(result["placer_call_count"], "3")
        self.assertEqual(result["window_placer_call_count"], "2")
        self.assertEqual(result["nearest_call_robust_mechanistic_qc"], "PASS_TE_LFDR")
        self.assertEqual(result["nearest_call_robust_mechanistic_worst_case_lfdr"], "0.04")
        self.assertEqual(result["nearest_call_conformal_qc"], "PASS_CONFORMAL_FDR")
        self.assertEqual(result["nearest_call_conformal_null_p"], "0.01")
        self.assertEqual(result["same_event_evidence_ledger_count"], "1")
        self.assertEqual(result["same_event_evidence_only_count"], "1")
        self.assertEqual(result["same_event_robust_lfdr_high_count"], "1")
        self.assertEqual(result["evidence_ledger_row_count"], "2")
        self.assertEqual(result["final_call_robust_worst_case_lfdr_count"], "3")
        self.assertAlmostEqual(
            float(result["final_call_robust_worst_case_lfdr_max"]),
            0.20,
        )
        self.assertEqual(result["final_call_conformal_null_p_count"], "3")
        self.assertEqual(result["final_call_conformal_null_p_max"], "0.01")
        self.assertEqual(result["null_control_row_count"], "1")
        self.assertEqual(result["final_call_empirical_null_tail_p_count"], "1")
        self.assertEqual(result["final_call_empirical_null_tail_p_max"], "0.05")

    def test_build_summary_reports_final_and_evidence_metrics(self):
        module = load_random_truth_eval_module()
        args = Namespace(
            truth_filter="PASS",
            min_length_ins=None,
            window_flank_bp=2000,
            match_distance_bp=500,
            workers=4,
            threads=2,
            batch_placer_run=False,
            reuse_existing_samples=False,
            sampled_truth_tsv=None,
            ground_truth=Path("/tmp/truth.tsv"),
            bam=Path("/tmp/input.bam"),
            ref=Path("/tmp/ref.fa"),
            te=Path("/tmp/te.fa"),
            placer=Path("/tmp/placer"),
        )
        args._pipeline_placer_statuses = {
            "S01": module.PipelinePlacerStatus("S01", 0, 1.25),
            "S02": module.PipelinePlacerStatus("S02", 0, 0.75),
        }
        rows = [
            module.TruthRow("S01", "u1", "chr1", 10, 10, "+", "L1", "L1-1", "PASS", {}, 100),
            module.TruthRow("S02", "u2", "chr1", 20, 20, "+", "L1", "L1-2", "PASS", {}, 100),
        ]
        results = [
            {
                "detected": "1",
                "placer_call_count": "2",
                "window_placer_call_count": "1",
                "failure_stage": "DETECTED",
                "same_event_evidence_ledger_count": "1",
                "same_event_evidence_only_count": "0",
                "same_event_robust_lfdr_high_count": "0",
                "evidence_ledger_row_count": "3",
                "final_call_robust_worst_case_lfdr_count": "2",
                "final_call_robust_worst_case_lfdr_sum": "0.12",
                "final_call_robust_worst_case_lfdr_max": "0.08",
                "null_control_row_count": "3",
                "final_call_empirical_null_tail_p_count": "2",
                "final_call_empirical_null_tail_p_sum": "0.09",
                "final_call_empirical_null_tail_p_max": "0.05",
                "final_call_conformal_null_p_count": "2",
                "final_call_conformal_null_p_sum": "0.03",
                "final_call_conformal_null_p_max": "0.02",
            },
            {
                "detected": "0",
                "placer_call_count": "0",
                "window_placer_call_count": "0",
                "failure_stage": "GENOTYPED_BUT_NO_FINAL_CALL",
                "same_event_evidence_ledger_count": "2",
                "same_event_evidence_only_count": "2",
                "same_event_robust_lfdr_high_count": "2",
                "evidence_ledger_row_count": "4",
                "final_call_robust_worst_case_lfdr_count": "0",
                "final_call_robust_worst_case_lfdr_sum": "0",
                "final_call_robust_worst_case_lfdr_max": "NA",
                "null_control_row_count": "4",
                "final_call_empirical_null_tail_p_count": "0",
                "final_call_empirical_null_tail_p_sum": "0",
                "final_call_empirical_null_tail_p_max": "NA",
                "final_call_conformal_null_p_count": "0",
                "final_call_conformal_null_p_sum": "0",
                "final_call_conformal_null_p_max": "NA",
            },
        ]

        summary = module.build_summary(
            seed=7,
            sampled_rows=rows,
            truth_pool_size=20,
            detected_count=1,
            failure_stage_counts={"DETECTED": 1, "GENOTYPED_BUT_NO_FINAL_CALL": 1},
            args=args,
            outdir=Path("/tmp/out"),
            total_elapsed_s=1.5,
            results=results,
        )

        self.assertEqual(summary["final_pass_count"], 1)
        self.assertEqual(summary["all_final_pass_count"], 2)
        self.assertEqual(summary["same_event_final_count"], 1)
        self.assertAlmostEqual(summary["truth_coverage_rate"], 0.5)
        self.assertAlmostEqual(summary["truth_fraction_of_placer_final_calls"], 1.0)
        self.assertFalse(summary["passes_truth_coverage_95pct"])
        self.assertTrue(summary["passes_truth_fraction_80pct"])
        self.assertEqual(summary["same_event_evidence_ledger_count"], 3)
        self.assertEqual(summary["same_event_evidence_only_count"], 2)
        self.assertEqual(summary["same_event_robust_lfdr_high_count"], 2)
        self.assertEqual(summary["evidence_ledger_row_count"], 7)
        self.assertAlmostEqual(summary["final_call_robust_worst_case_lfdr_mean"], 0.06)
        self.assertAlmostEqual(summary["final_call_robust_worst_case_lfdr_max"], 0.08)
        self.assertEqual(summary["truth_subset_metrics"]["NA"]["final_pass_count"], 1)
        self.assertEqual(summary["truth_subset_metrics"]["NA"]["all_final_pass_count"], 2)
        self.assertEqual(summary["null_control_row_count"], 7)
        self.assertAlmostEqual(summary["final_call_empirical_null_tail_p_mean"], 0.045)
        self.assertAlmostEqual(summary["final_call_empirical_null_tail_p_max"], 0.05)
        self.assertAlmostEqual(summary["final_call_conformal_null_p_mean"], 0.015)
        self.assertAlmostEqual(summary["final_call_conformal_null_p_max"], 0.02)
        self.assertEqual(summary["placer_run_mode"], "pipeline")
        self.assertEqual(
            summary["placer_pipeline_manifest"],
            "/tmp/out/placer_pipeline_manifest.tsv",
        )
        self.assertEqual(summary["placer_pipeline_log"], "/tmp/out/placer_pipeline.log")
        self.assertAlmostEqual(summary["placer_pipeline_elapsed_s"], 2.0)
        self.assertNotIn("batch_manifest", summary)
        self.assertNotIn("batch_placer_elapsed_s", summary)


if __name__ == "__main__":
    unittest.main()
