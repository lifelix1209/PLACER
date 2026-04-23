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


def write_alignment(
    out_handle,
    *,
    qname: str,
    flag: int,
    pos0: int,
    mapq: int = 60,
) -> None:
    import pysam

    seg = pysam.AlignedSegment()
    seg.query_name = qname
    seg.flag = flag
    seg.reference_id = 0
    seg.reference_start = pos0
    seg.mapping_quality = mapq
    seg.cigar = ((0, 4),)
    seg.query_sequence = "ACGT"
    seg.query_qualities = pysam.qualitystring_to_array("IIII")
    out_handle.write(seg)


class RandomTruthIntervalEvalTest(unittest.TestCase):
    def test_resolve_samtools_threads_uses_all_detected_cpus_by_default(self):
        original_cpu_count = MODULE.os.cpu_count
        MODULE.os.cpu_count = lambda: 12
        try:
            self.assertEqual(MODULE.resolve_samtools_threads(None), 12)
        finally:
            MODULE.os.cpu_count = original_cpu_count

    def test_resolve_samtools_threads_preserves_explicit_user_value(self):
        self.assertEqual(MODULE.resolve_samtools_threads(3), 3)

    def test_extract_primary_qnames_from_region_view_keeps_unique_primary_names(self):
        region_stdout = "\n".join(
            [
                "readA\t0\tchr1\t101\t60\t50M\t*\t0\t0\tACGT\t*",
                "readA\t2048\tchr1\t130\t60\t50M\t*\t0\t0\tACGT\t*",
                "readB\t256\tchr1\t140\t60\t50M\t*\t0\t0\tACGT\t*",
                "readC\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*",
                "readD\t0\tchr1\t155\t60\t50M\t*\t0\t0\tACGT\t*",
            ]
        )

        qnames = MODULE.extract_primary_qnames_from_region_view(
            region_stdout,
            "chr1:100-200",
        )

        self.assertEqual(qnames, ["readA", "readD"])

    def test_build_qname_membership_maps_preserves_order_and_shared_qnames(self):
        sample_to_qnames, qname_to_sample_ids, union_qnames = MODULE.build_qname_membership_maps(
            {
                "S01": ["readA", "shared", "readB"],
                "S02": ["shared", "readC"],
            }
        )

        self.assertEqual(sample_to_qnames["S01"], ["readA", "shared", "readB"])
        self.assertEqual(sample_to_qnames["S02"], ["shared", "readC"])
        self.assertEqual(qname_to_sample_ids["shared"], ["S01", "S02"])
        self.assertEqual(union_qnames, ["readA", "shared", "readB", "readC"])

    def test_split_union_subset_bam_duplicates_shared_qnames(self):
        import pysam

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            union_bam = tmp / "union_subset.bam"
            s01_bam = tmp / "S01.subset.bam"
            s02_bam = tmp / "S02.subset.bam"

            header = {
                "HD": {"VN": "1.6"},
                "SQ": [{"SN": "chr1", "LN": 1000}],
            }
            with pysam.AlignmentFile(union_bam, "wb", header=header) as out_handle:
                write_alignment(out_handle, qname="shared", flag=0, pos0=100)
                write_alignment(out_handle, qname="shared", flag=2048, pos0=500)
                write_alignment(out_handle, qname="only_s1", flag=0, pos0=140)
                write_alignment(out_handle, qname="only_s2", flag=0, pos0=220)

            record_counts = MODULE.split_union_subset_bam(
                union_subset_bam=union_bam,
                qname_to_sample_ids={
                    "shared": ["S01", "S02"],
                    "only_s1": ["S01"],
                    "only_s2": ["S02"],
                },
                sample_to_subset_bam={
                    "S01": s01_bam,
                    "S02": s02_bam,
                },
            )

            self.assertEqual(record_counts, {"S01": 3, "S02": 3})
            with pysam.AlignmentFile(s01_bam, "rb") as handle:
                self.assertEqual(
                    [record.query_name for record in handle],
                    ["shared", "shared", "only_s1"],
                )
            with pysam.AlignmentFile(s02_bam, "rb") as handle:
                self.assertEqual(
                    [record.query_name for record in handle],
                    ["shared", "shared", "only_s2"],
                )

    def test_evaluate_sample_rows_prepares_subsets_once_before_running_samples(self):
        rows = [
            make_truth_row("S01", 100, 110),
            make_truth_row("S02", 200, 210),
        ]
        args = types.SimpleNamespace(
            window_flank_bp=25,
            reuse_existing_samples=False,
            workers=1,
        )
        contig_lengths = {"chr1": 1000}
        trace: list[str] = []

        original_prepare = getattr(MODULE, "prepare_sample_subset_bams", None)
        original_run = MODULE.run_placer_sample

        def fake_prepare(*, sampled_rows, windows_by_sample_id, args, outdir):
            trace.append("prepare")
            return {
                "S01": {"qname_count": 2, "extract_elapsed_s": 0.1},
                "S02": {"qname_count": 3, "extract_elapsed_s": 0.2},
            }

        def fake_run(*, row, window, args, outdir):
            trace.append(f"run:{row.sample_id}")
            return {
                "sample_id": row.sample_id,
                "window_region": window.region,
            }

        MODULE.prepare_sample_subset_bams = fake_prepare
        MODULE.run_placer_sample = fake_run
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                manifest, results = MODULE.evaluate_sample_rows(
                    sampled_rows=rows,
                    contig_lengths=contig_lengths,
                    args=args,
                    outdir=Path(tmpdir),
                )
        finally:
            if original_prepare is None:
                delattr(MODULE, "prepare_sample_subset_bams")
            else:
                MODULE.prepare_sample_subset_bams = original_prepare
            MODULE.run_placer_sample = original_run

        self.assertEqual(trace, ["prepare", "run:S01", "run:S02"])
        self.assertEqual([row["sample_id"] for row in manifest], ["S01", "S02"])
        self.assertEqual([row["sample_id"] for row in results], ["S01", "S02"])

    def test_run_placer_sample_uses_prepared_subset_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "S01"
            sample_dir.mkdir()
            subset_bam = sample_dir / "subset.bam"
            subset_bam.write_bytes(b"BAM")
            (sample_dir / "window_qnames.txt").write_text("readA\nreadB\n", encoding="utf-8")

            row = make_truth_row("S01", 100, 110)
            window = MODULE.WindowSpec(
                chrom="chr1",
                truth_start=100,
                truth_end=110,
                window_start=80,
                window_end=130,
            )
            args = types.SimpleNamespace(
                placer=Path("/tmp/fake_placer"),
                ref=Path("/tmp/ref.fa"),
                te=Path("/tmp/te.fa"),
                match_distance_bp=1000,
            )

            original_subprocess_run = MODULE.subprocess.run
            original_parse_scientific_txt = MODULE.parse_scientific_txt
            original_choose_best_call = MODULE.choose_best_call
            original_parse_joint_component_metrics = MODULE.parse_joint_component_metrics

            calls = []

            def fake_run(cmd, cwd=None, check=False, text=True, capture_output=True):
                calls.append((cmd, cwd))
                return types.SimpleNamespace(returncode=0, stdout="", stderr="")

            MODULE.subprocess.run = fake_run
            MODULE.parse_scientific_txt = lambda path: []
            MODULE.choose_best_call = lambda row, calls: None
            MODULE.parse_joint_component_metrics = lambda path: {
                "joint_best_kind": "NA",
                "joint_runner_up_kind": "NA",
                "joint_final_qc": "NA",
                "joint_emit_te": "NA",
                "joint_emit_unknown": "NA",
            }
            try:
                result = MODULE.run_placer_sample(
                    row=row,
                    window=window,
                    args=args,
                    outdir=Path(tmpdir),
                )
            finally:
                MODULE.subprocess.run = original_subprocess_run
                MODULE.parse_scientific_txt = original_parse_scientific_txt
                MODULE.choose_best_call = original_choose_best_call
                MODULE.parse_joint_component_metrics = original_parse_joint_component_metrics

            self.assertEqual(calls[0][0][1], str(subset_bam))
            self.assertEqual(result["window_qname_count"], "2")

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

    def test_parse_scientific_summary_extracts_pipeline_counters(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            scientific_path = Path(tmpdir) / "scientific.txt"
            scientific_path.write_text(
                "#PLACER streaming pipeline summary\n"
                "total_reads\t72\n"
                "gate1_passed\t53\n"
                "processed_bins\t1\n"
                "components\t7\n"
                "event_consensus_calls\t5\n"
                "genotype_calls\t5\n"
                "final_pass_calls\t1\n"
                "schema_version\t1.0.0\n"
                "\n"
                "#chrom\tpos\tbp_left\tbp_right\tte\tfamily\tsubfamily\tstrand\tinsert_len\t"
                "support_reads\talt_struct_reads\tref_span_reads\tlow_mapq_ref_span_reads\t"
                "gt\taf\tgq\tbest_te_identity\tbest_te_query_coverage\tcross_family_margin\t"
                "tsd_type\ttsd_len\tleft_flank_align_len\tright_flank_align_len\tconsensus_len\tqc\n",
                encoding="utf-8",
            )

            summary = MODULE.parse_scientific_summary(scientific_path)

        self.assertEqual(summary["event_consensus_calls"], 5)
        self.assertEqual(summary["genotype_calls"], 5)
        self.assertEqual(summary["final_pass_calls"], 1)

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
            scientific_summary=None,
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
            scientific_summary=None,
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
            scientific_summary=None,
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
            scientific_summary=None,
        )

        self.assertEqual(stage, "TE_EVIDENCE_BUT_NON_TE_WINS")

    def test_classify_failure_stage_uses_scientific_summary_when_joint_metrics_missing(self):
        stage = MODULE.classify_failure_stage(
            detected=False,
            joint_metrics={
                "joint_best_kind": "NA",
                "joint_runner_up_kind": "NA",
                "joint_final_qc": "NA",
                "joint_emit_te": "NA",
                "joint_emit_unknown": "NA",
                "consensus_qc": "NA",
                "consensus_len": "NA",
                "seg_qc": "NA",
                "seg_insert_len": "NA",
                "te_qc": "NA",
                "te_pass": "NA",
            },
            scientific_summary={
                "event_consensus_calls": 5,
                "genotype_calls": 5,
                "final_pass_calls": 0,
            },
            best_call=None,
        )

        self.assertEqual(stage, "GENOTYPED_BUT_NO_FINAL_CALL")

    def test_classify_failure_stage_reports_final_call_elsewhere(self):
        stage = MODULE.classify_failure_stage(
            detected=False,
            joint_metrics={
                "joint_best_kind": "NA",
                "joint_runner_up_kind": "NA",
                "joint_final_qc": "NA",
                "joint_emit_te": "NA",
                "joint_emit_unknown": "NA",
                "consensus_qc": "NA",
                "consensus_len": "NA",
                "seg_qc": "NA",
                "seg_insert_len": "NA",
                "te_qc": "NA",
                "te_pass": "NA",
            },
            scientific_summary={
                "event_consensus_calls": 7,
                "genotype_calls": 7,
                "final_pass_calls": 2,
            },
            best_call={
                "interval_distance_bp": "2401",
            },
        )

        self.assertEqual(stage, "FINAL_CALL_ELSEWHERE")

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
