#!/usr/bin/env python3

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "scripts" / "run_sharded_placer.py"

SPEC = importlib.util.spec_from_file_location("run_sharded_placer", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Failed to load module spec from {MODULE_PATH}")
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def make_summary(**overrides):
    summary = {
        "total_reads": "100",
        "gate1_passed": "40",
        "processed_bins": "10",
        "components": "3",
        "event_consensus_calls": "2",
        "genotype_calls": "2",
        "final_pass_calls": "1",
        "schema_version": "1.0.0",
    }
    summary.update(overrides)
    return summary


def make_row(**overrides):
    row = {
        "chrom": "chr1",
        "pos": "100",
        "bp_left": "95",
        "bp_right": "105",
        "te": "GypsyA",
        "family": "Gypsy",
        "subfamily": "GypsyA",
        "strand": "NA",
        "insert_len": "120",
        "support_reads": "5",
        "alt_struct_reads": "5",
        "ref_span_reads": "1",
        "low_mapq_ref_span_reads": "0",
        "gt": "0/1",
        "af": "0.83",
        "gq": "28",
        "best_te_identity": "0.91",
        "best_te_query_coverage": "0.84",
        "cross_family_margin": "0.20",
        "tsd_type": "TSD",
        "tsd_len": "5",
        "left_flank_align_len": "65",
        "right_flank_align_len": "67",
        "consensus_len": "210",
        "qc": "PASS",
    }
    row.update(overrides)
    return row


def write_scientific(path: Path, summary, rows, header=None):
    columns = header or MODULE.EXPECTED_SCIENTIFIC_HEADER
    with path.open("w", encoding="utf-8") as handle:
        handle.write("#PLACER streaming pipeline summary\n")
        for key, value in summary.items():
            handle.write(f"{key}\t{value}\n")
        handle.write("\n")
        handle.write("#" + "\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(col, "")) for col in columns) + "\n")


class RunShardedPlacerTest(unittest.TestCase):
    def test_parse_scientific_rejects_unexpected_schema(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "unexpected.scientific.txt"
            malformed_header = [
                col for col in MODULE.EXPECTED_SCIENTIFIC_HEADER
                if col != "cross_family_margin"
            ]
            malformed_header.insert(18, "unexpected_margin")
            write_scientific(path, make_summary(), [make_row()], header=malformed_header)

            with self.assertRaisesRegex(RuntimeError, "unexpected scientific.txt schema"):
                MODULE.parse_scientific(path)

    def test_merge_shard_results_dedups_by_locus_using_new_metrics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            shard_one_path = root / "shard1.scientific.txt"
            shard_two_path = root / "shard2.scientific.txt"
            merged_path = root / "scientific.sharded.txt"

            row_one = make_row(
                pos="100",
                te="GypsyA",
                family="Gypsy",
                subfamily="GypsyA",
                cross_family_margin="0.20",
                best_te_identity="0.91",
                best_te_query_coverage="0.84",
                consensus_len="210",
            )
            row_two = make_row(
                pos="103",
                te="L1HS",
                family="L1",
                subfamily="L1HS",
                cross_family_margin="0.37",
                best_te_identity="0.89",
                best_te_query_coverage="0.88",
                consensus_len="190",
            )

            write_scientific(shard_one_path, make_summary(total_reads="120", event_consensus_calls="3"), [row_one])
            write_scientific(shard_two_path, make_summary(total_reads="80", event_consensus_calls="4"), [row_two])

            shard_one = MODULE.ShardResult(
                spec=MODULE.ShardSpec(
                    shard_id=1,
                    label="0001_chr1",
                    chrom="chr1",
                    core_start=1,
                    core_end=1000,
                    fetch_start=1,
                    fetch_end=1000,
                ),
                workdir=root,
                scientific_path=shard_one_path,
                summary=make_summary(total_reads="120", event_consensus_calls="3"),
                n_rows_raw=1,
                elapsed_s=1.0,
            )
            shard_two = MODULE.ShardResult(
                spec=MODULE.ShardSpec(
                    shard_id=2,
                    label="0002_chr1",
                    chrom="chr1",
                    core_start=1,
                    core_end=1000,
                    fetch_start=1,
                    fetch_end=1000,
                ),
                workdir=root,
                scientific_path=shard_two_path,
                summary=make_summary(total_reads="80", event_consensus_calls="4"),
                n_rows_raw=1,
                elapsed_s=1.0,
            )

            MODULE.merge_shard_results(
                [shard_one, shard_two],
                mode="contig",
                dedup_bp=10,
                chrom_order={"chr1": 0},
                region_size=1000,
                overlap_bp=100,
                out_path=merged_path,
            )

            summary, header, rows = MODULE.parse_scientific(merged_path)
            self.assertEqual(header, MODULE.EXPECTED_SCIENTIFIC_HEADER)
            self.assertEqual(summary["event_consensus_calls"], "7")
            self.assertEqual(summary["final_pass_calls"], "1")
            self.assertEqual(summary["schema_version"], "1.0.0-sharded")
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["te"], "L1HS")
            self.assertEqual(rows[0]["family"], "L1")
            self.assertEqual(rows[0]["cross_family_margin"], "0.37")


if __name__ == "__main__":
    unittest.main()
