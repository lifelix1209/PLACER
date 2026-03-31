#!/usr/bin/env python3

import csv
import importlib.util
import sys
import tempfile
import threading
import time
import types
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "scripts" / "random_truth_interval_eval.py"

SPEC = importlib.util.spec_from_file_location("random_truth_interval_eval", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Failed to load module spec from {MODULE_PATH}")
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def make_truth_row(sample_id: str, start: int, end: int, *, length_ins: int = 0) -> object:
    return MODULE.TruthRow(
        sample_id=sample_id,
        uuid=f"uuid_{sample_id}",
        chrom="chr1",
        start=start,
        end=end,
        strand="+",
        family="Fam",
        subfamily="Subfam",
        filter_value="PASS",
        raw={},
        length_ins=length_ins,
    )


def write_truth_table(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "UUID",
        "Chrom",
        "Start",
        "End",
        "Strand",
        "Family",
        "Subfamily",
        "LengthIns",
        "Filter",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class RandomTruthIntervalEvalTest(unittest.TestCase):
    def test_load_truth_rows_applies_min_length_ins_filter(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            truth_path = Path(tmpdir) / "truth.tsv"
            write_truth_table(
                truth_path,
                [
                    {
                        "UUID": "uuid_short",
                        "Chrom": "chr1",
                        "Start": "100",
                        "End": "120",
                        "Strand": "+",
                        "Family": "Fam",
                        "Subfamily": "Sub",
                        "LengthIns": "800",
                        "Filter": "PASS",
                    },
                    {
                        "UUID": "uuid_long",
                        "Chrom": "chr1",
                        "Start": "200",
                        "End": "220",
                        "Strand": "-",
                        "Family": "Fam",
                        "Subfamily": "Sub",
                        "LengthIns": "1400",
                        "Filter": "PASS",
                    },
                    {
                        "UUID": "uuid_filtered",
                        "Chrom": "chr1",
                        "Start": "300",
                        "End": "320",
                        "Strand": "+",
                        "Family": "Fam",
                        "Subfamily": "Sub",
                        "LengthIns": "2000",
                        "Filter": "LowQual",
                    },
                ],
            )

            rows = MODULE.load_truth_rows(
                truth_path,
                truth_filter="PASS",
                min_length_ins=1000,
            )

        self.assertEqual([row.uuid for row in rows], ["uuid_long"])
        self.assertEqual(rows[0].length_ins, 1400)

    def test_evaluate_sample_rows_manifest_includes_truth_length_ins(self):
        rows = [make_truth_row("S09", 900, 930, length_ins=1500)]
        args = types.SimpleNamespace(
            window_flank_bp=25,
            reuse_existing_samples=False,
            workers=1,
        )
        contig_lengths = {"chr1": 2000}

        def evaluator(*, row, window, args, outdir):
            return {
                "sample_id": row.sample_id,
                "truth_length_ins": str(row.length_ins),
                "window_region": window.region,
            }

        with tempfile.TemporaryDirectory() as tmpdir:
            manifest, results = MODULE.evaluate_sample_rows(
                sampled_rows=rows,
                contig_lengths=contig_lengths,
                args=args,
                outdir=Path(tmpdir),
                evaluator=evaluator,
            )

        self.assertEqual(manifest[0]["truth_length_ins"], "1500")
        self.assertEqual(results[0]["truth_length_ins"], "1500")

    def test_evaluate_sample_rows_parallel_preserves_order_and_uses_multiple_workers(self):
        rows = [
            make_truth_row("S01", 100, 110),
            make_truth_row("S02", 200, 210),
            make_truth_row("S03", 300, 310),
        ]
        args = types.SimpleNamespace(
            window_flank_bp=25,
            reuse_existing_samples=False,
            workers=3,
        )
        contig_lengths = {"chr1": 1000}
        delays = {"S01": 0.20, "S02": 0.20, "S03": 0.05}

        lock = threading.Lock()
        active_workers = 0
        max_active_workers = 0

        def evaluator(*, row, window, args, outdir):
            nonlocal active_workers, max_active_workers
            with lock:
                active_workers += 1
                max_active_workers = max(max_active_workers, active_workers)
            time.sleep(delays[row.sample_id])
            with lock:
                active_workers -= 1
            return {
                "sample_id": row.sample_id,
                "window_region": window.region,
            }

        with tempfile.TemporaryDirectory() as tmpdir:
            manifest, results = MODULE.evaluate_sample_rows(
                sampled_rows=rows,
                contig_lengths=contig_lengths,
                args=args,
                outdir=Path(tmpdir),
                evaluator=evaluator,
            )

        self.assertEqual([row["sample_id"] for row in manifest], ["S01", "S02", "S03"])
        self.assertEqual([row["sample_id"] for row in results], ["S01", "S02", "S03"])
        self.assertGreaterEqual(max_active_workers, 2)

    def test_parse_joint_component_metrics_picks_last_joint_line(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "run.log"
            log_path.write_text(
                "[Pipeline][joint] best_kind=REFERENCE runner_kind=NON_TE_INSERTION "
                "final_qc=REJECT_EVENT_EXISTENCE emit_te=0 emit_unknown=0\n"
                "[Pipeline][joint] best_kind=TE_RESOLVED runner_kind=NON_TE_INSERTION "
                "final_qc=PASS_FINAL_TE_CALL emit_te=1 emit_unknown=0\n",
                encoding="utf-8",
            )

            metrics = MODULE.parse_joint_component_metrics(log_path)

        self.assertEqual(metrics["joint_best_kind"], "TE_RESOLVED")
        self.assertEqual(metrics["joint_runner_up_kind"], "NON_TE_INSERTION")
        self.assertEqual(metrics["joint_final_qc"], "PASS_FINAL_TE_CALL")
        self.assertEqual(metrics["joint_emit_te"], "1")
        self.assertEqual(metrics["joint_emit_unknown"], "0")

    def test_parse_joint_component_metrics_extracts_stage_fields(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "run.log"
            log_path.write_text(
                "[Pipeline][joint] consensus_qc=PASS_EVENT_CONSENSUS consensus_len=1700 "
                "seg_qc=NO_RIGHT_FLANK_MATCH seg_insert_len=0 "
                "te_qc=NO_TE_ALIGNMENT te_pass=0 "
                "best_kind=NON_TE_INSERTION runner_kind=TE_UNKNOWN "
                "final_qc=NON_TE_INSERTION emit_te=0 emit_unknown=0\n",
                encoding="utf-8",
            )

            metrics = MODULE.parse_joint_component_metrics(log_path)

        self.assertEqual(metrics["consensus_qc"], "PASS_EVENT_CONSENSUS")
        self.assertEqual(metrics["consensus_len"], "1700")
        self.assertEqual(metrics["seg_qc"], "NO_RIGHT_FLANK_MATCH")
        self.assertEqual(metrics["seg_insert_len"], "0")
        self.assertEqual(metrics["te_qc"], "NO_TE_ALIGNMENT")
        self.assertEqual(metrics["te_pass"], "0")

    def test_classify_failure_stage_returns_event_but_no_consensus(self):
        stage = MODULE.classify_failure_stage(
            detected=False,
            joint_metrics={
                "joint_best_kind": "NA",
                "joint_runner_up_kind": "NA",
                "joint_final_qc": "REJECT_EVENT_EXISTENCE",
                "joint_emit_te": "0",
                "joint_emit_unknown": "0",
                "consensus_qc": "INSUFFICIENT_EVENT_READS",
                "consensus_len": "NA",
                "seg_qc": "NA",
                "seg_insert_len": "NA",
                "te_qc": "NA",
                "te_pass": "NA",
            },
        )

        self.assertEqual(stage, "EVENT_BUT_NO_CONSENSUS")

    def test_classify_failure_stage_returns_consensus_but_no_segmentation(self):
        stage = MODULE.classify_failure_stage(
            detected=False,
            joint_metrics={
                "joint_best_kind": "NON_TE_INSERTION",
                "joint_runner_up_kind": "REFERENCE",
                "joint_final_qc": "NON_TE_INSERTION",
                "joint_emit_te": "0",
                "joint_emit_unknown": "0",
                "consensus_qc": "PASS_EVENT_CONSENSUS",
                "consensus_len": "1850",
                "seg_qc": "NO_TRIPARTITE_EVENT_SEGMENTATION",
                "seg_insert_len": "NA",
                "te_qc": "NA",
                "te_pass": "NA",
            },
        )

        self.assertEqual(stage, "CONSENSUS_BUT_NO_SEGMENTATION")

    def test_classify_failure_stage_returns_segmented_but_no_te_alignment(self):
        stage = MODULE.classify_failure_stage(
            detected=False,
            joint_metrics={
                "joint_best_kind": "NON_TE_INSERTION",
                "joint_runner_up_kind": "REFERENCE",
                "joint_final_qc": "NON_TE_INSERTION",
                "joint_emit_te": "0",
                "joint_emit_unknown": "0",
                "consensus_qc": "PASS_EVENT_CONSENSUS",
                "consensus_len": "1900",
                "seg_qc": "PASS_EVENT_SEGMENTATION",
                "seg_insert_len": "1250",
                "te_qc": "TE_ALIGNMENT_LOW_QUERY_COVERAGE",
                "te_pass": "0",
            },
        )

        self.assertEqual(stage, "SEGMENTED_BUT_NO_TE_ALIGNMENT")

    def test_classify_failure_stage_returns_te_evidence_but_non_te_wins(self):
        stage = MODULE.classify_failure_stage(
            detected=False,
            joint_metrics={
                "joint_best_kind": "NON_TE_INSERTION",
                "joint_runner_up_kind": "TE_UNKNOWN",
                "joint_final_qc": "NON_TE_INSERTION",
                "joint_emit_te": "0",
                "joint_emit_unknown": "0",
                "consensus_qc": "PASS_EVENT_CONSENSUS",
                "consensus_len": "1600",
                "seg_qc": "PASS_EVENT_SEGMENTATION",
                "seg_insert_len": "1100",
                "te_qc": "PASS_INSERT_TE_ALIGNMENT_UNKNOWN",
                "te_pass": "1",
            },
        )

        self.assertEqual(stage, "TE_EVIDENCE_BUT_NON_TE_WINS")

    def test_build_result_row_includes_joint_component_metrics(self):
        row = make_truth_row("S07", 500, 520)
        window = MODULE.WindowSpec(
            chrom="chr1",
            truth_start=500,
            truth_end=520,
            window_start=450,
            window_end=570,
        )

        result = MODULE.build_result_row(
            row=row,
            window=window,
            sample_dir=Path("/tmp/S07"),
            read_names_path=Path("/tmp/S07/window_qnames.txt"),
            subset_bam=Path("/tmp/S07/subset.bam"),
            scientific_path=Path("/tmp/S07/scientific.txt"),
            qname_count=8,
            placer_exit_code="0",
            calls=[],
            best_call=None,
            detected=False,
            extract_elapsed_s=0.0,
            placer_elapsed_s=0.0,
            sample_elapsed_s=0.0,
            joint_metrics={
                "joint_best_kind": "TE_UNKNOWN",
                "joint_runner_up_kind": "NON_TE_INSERTION",
                "joint_final_qc": "PASS_FINAL_TE_CALL_UNKNOWN",
                "joint_emit_te": "1",
                "joint_emit_unknown": "1",
            },
        )

        self.assertEqual(result["joint_best_kind"], "TE_UNKNOWN")
        self.assertEqual(result["joint_runner_up_kind"], "NON_TE_INSERTION")
        self.assertEqual(result["joint_final_qc"], "PASS_FINAL_TE_CALL_UNKNOWN")
        self.assertEqual(result["joint_emit_te"], "1")
        self.assertEqual(result["joint_emit_unknown"], "1")

    def test_build_result_row_records_timing_fields(self):
        row = make_truth_row("S08", 700, 710)
        window = MODULE.WindowSpec(
            chrom="chr1",
            truth_start=700,
            truth_end=710,
            window_start=650,
            window_end=760,
        )

        result = MODULE.build_result_row(
            row=row,
            window=window,
            sample_dir=Path("/tmp/S08"),
            read_names_path=Path("/tmp/S08/window_qnames.txt"),
            subset_bam=Path("/tmp/S08/subset.bam"),
            scientific_path=Path("/tmp/S08/scientific.txt"),
            qname_count=6,
            placer_exit_code="0",
            calls=[],
            best_call=None,
            detected=False,
            extract_elapsed_s=0.12,
            placer_elapsed_s=1.34,
            sample_elapsed_s=1.46,
            joint_metrics={
                "joint_best_kind": "NA",
                "joint_runner_up_kind": "NA",
                "joint_final_qc": "NA",
                "joint_emit_te": "NA",
                "joint_emit_unknown": "NA",
            },
        )

        self.assertEqual(result["extract_elapsed_s"], "0.120000")
        self.assertEqual(result["placer_elapsed_s"], "1.340000")
        self.assertEqual(result["sample_elapsed_s"], "1.460000")

    def test_build_summary_records_parallel_runtime_metadata(self):
        args = types.SimpleNamespace(
            ground_truth=Path("/tmp/truth.tsv"),
            bam=Path("/tmp/input.bam"),
            ref=Path("/tmp/ref.fa"),
            te=Path("/tmp/te.fa"),
            placer=Path("/tmp/placer"),
            truth_filter="PASS",
            min_length_ins=None,
            window_flank_bp=2000,
            match_distance_bp=1000,
            threads=1,
            workers=4,
            sampled_truth_tsv=None,
        )

        summary = MODULE.build_summary(
            seed=20260330,
            sampled_rows=[make_truth_row("S01", 100, 110)],
            truth_pool_size=200,
            detected_count=1,
            failure_stage_counts={"DETECTED": 1},
            args=args,
            outdir=Path("/tmp/out"),
            total_elapsed_s=12.5,
        )

        self.assertEqual(summary["workers"], 4)
        self.assertEqual(summary["samtools_threads"], 1)
        self.assertEqual(summary["total_elapsed_s"], 12.5)

    def test_build_summary_records_long_te_metadata(self):
        args = types.SimpleNamespace(
            ground_truth=Path("/tmp/truth.tsv"),
            bam=Path("/tmp/input.bam"),
            ref=Path("/tmp/ref.fa"),
            te=Path("/tmp/te.fa"),
            placer=Path("/tmp/placer"),
            truth_filter="PASS",
            min_length_ins=1000,
            window_flank_bp=2000,
            match_distance_bp=1000,
            threads=1,
            workers=1,
            sampled_truth_tsv=None,
        )

        summary = MODULE.build_summary(
            seed=20260330,
            sampled_rows=[make_truth_row("S01", 100, 110, length_ins=1500)],
            truth_pool_size=42,
            detected_count=6,
            failure_stage_counts={
                "DETECTED": 6,
                "CONSENSUS_BUT_NO_SEGMENTATION": 4,
            },
            args=args,
            outdir=Path("/tmp/out"),
            total_elapsed_s=10.0,
        )

        self.assertEqual(summary["min_length_ins"], 1000)
        self.assertEqual(summary["failure_stage_counts"]["DETECTED"], 6)
        self.assertEqual(summary["failure_stage_counts"]["CONSENSUS_BUT_NO_SEGMENTATION"], 4)


if __name__ == "__main__":
    unittest.main()
