#!/usr/bin/env python3

import importlib.util
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    module_path = REPO_ROOT / relative_path
    scripts_dir = str(REPO_ROOT / "scripts")
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    spec = importlib.util.spec_from_file_location(name, module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class RebuildIgvValidationBamsTest(unittest.TestCase):
    def test_parse_args_defaults_to_ten_pass_truth_loci_with_2kb_flanks(self):
        module = load_module(
            "rebuild_igv_validation_bams",
            "scripts/rebuild_igv_validation_bams.py",
        )

        args = module.parse_args(
            [
                "--ground-truth",
                "/tmp/truth.tsv",
                "--bam",
                "/tmp/input.bam",
                "--ref",
                "/tmp/ref.fa",
            ]
        )

        self.assertEqual(args.sample_size, 10)
        self.assertEqual(args.truth_filter, "PASS")
        self.assertEqual(args.window_flank_bp, 2000)

    def test_build_sample_manifest_clips_windows_to_reference_bounds(self):
        random_eval = load_module(
            "random_truth_interval_eval",
            "scripts/random_truth_interval_eval.py",
        )
        module = load_module(
            "rebuild_igv_validation_bams",
            "scripts/rebuild_igv_validation_bams.py",
        )
        row = random_eval.TruthRow(
            sample_id="S01",
            uuid="truth-1",
            chrom="chr1",
            start=100,
            end=120,
            strand="+",
            family="LINE",
            subfamily="L1",
            filter_value="PASS",
            raw={},
            length_ins=500,
        )

        manifest, windows = module.build_sample_manifest(
            [row],
            {"chr1": 1000},
            flank_bp=200,
        )

        self.assertEqual(manifest[0]["window_region"], "chr1:1-320")
        self.assertEqual(windows["S01"].window_start, 1)
        self.assertEqual(windows["S01"].window_end, 320)

    def test_build_bam_manifest_records_igv_paths(self):
        random_eval = load_module(
            "random_truth_interval_eval",
            "scripts/random_truth_interval_eval.py",
        )
        module = load_module(
            "rebuild_igv_validation_bams",
            "scripts/rebuild_igv_validation_bams.py",
        )
        row = random_eval.TruthRow(
            sample_id="S01",
            uuid="truth-1",
            chrom="chr1",
            start=100,
            end=120,
            strand="+",
            family="LINE",
            subfamily="L1",
            filter_value="PASS",
            raw={},
            length_ins=500,
        )
        window = random_eval.WindowSpec(
            chrom="chr1",
            truth_start=100,
            truth_end=120,
            window_start=1,
            window_end=320,
        )
        subset = random_eval.PreparedSubset(
            sample_id="S01",
            sample_dir=Path("/tmp/out/S01"),
            log_path=Path("/tmp/out/S01/run.log"),
            read_names_path=Path("/tmp/out/S01/window_qnames.txt"),
            subset_bam=Path("/tmp/out/S01/subset.bam"),
            qname_count=12,
            extract_elapsed_s=0.5,
        )

        manifest = module.build_bam_manifest(
            [row],
            {"S01": window},
            {"S01": subset},
        )

        self.assertEqual(manifest[0]["subset_bam"], "/tmp/out/S01/subset.bam")
        self.assertEqual(manifest[0]["subset_bai"], "/tmp/out/S01/subset.bam.bai")
        self.assertEqual(manifest[0]["qname_count"], "12")


if __name__ == "__main__":
    unittest.main()
