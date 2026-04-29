#!/usr/bin/env python3

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_module():
    module_path = REPO_ROOT / "scripts" / "blast_insert_sequence.py"
    spec = importlib.util.spec_from_file_location("blast_insert_sequence", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class BlastInsertSequenceTest(unittest.TestCase):
    def test_normalize_sequence_keeps_iupac_bases_and_uppercases(self):
        module = load_module()

        seq = module.normalize_sequence(" acgtNNry-\n\t")

        self.assertEqual(seq, "ACGTNNRY")

    def test_load_sequence_file_accepts_fasta(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / "insert.fa"
            fasta.write_text(">insert\nacgt\nnn\n", encoding="utf-8")

            seq = module.load_sequence_from_file(fasta)

        self.assertEqual(seq, "ACGTNN")

    def test_parse_blast_rows_computes_query_coverage(self):
        module = load_module()

        rows = module.parse_blast_tabular(
            "query1\tTE_A\t95.000\t90\t2\t1\t1\t90\t10\t99\t1e-20\t180\t100\t500\n"
            "query1\tTE_B\t80.000\t50\t10\t0\t20\t69\t1\t50\t1e-5\t75\t100\t300\n"
        )

        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0].subject_id, "TE_A")
        self.assertEqual(rows[0].query_coverage_pct, 90.0)
        self.assertEqual(rows[1].query_coverage_pct, 50.0)


if __name__ == "__main__":
    unittest.main()
