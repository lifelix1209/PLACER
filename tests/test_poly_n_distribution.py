from pathlib import Path
import csv
import contextlib
import gzip
import io
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from poly_n_distribution import bin_poly_n_lengths, main, scan_poly_n_runs


class PolyNDistributionTests(unittest.TestCase):
    def test_scan_poly_n_runs_merges_runs_across_fasta_lines(self):
        fasta = Path(self._testMethodName + ".fa")
        self.addCleanup(fasta.unlink, missing_ok=True)
        fasta.write_text(">chr1\nAAANN\nNNTTA\n>chr2 description\nnCGNN\n", encoding="ascii")

        runs = list(scan_poly_n_runs(fasta))

        self.assertEqual(
            [(r.contig, r.start, r.end, r.length) for r in runs],
            [
                ("chr1", 3, 7, 4),
                ("chr2", 0, 1, 1),
                ("chr2", 3, 5, 2),
            ],
        )

    def test_scan_poly_n_runs_accepts_gzipped_fasta(self):
        fasta = Path(self._testMethodName + ".fa.gz")
        self.addCleanup(fasta.unlink, missing_ok=True)
        with gzip.open(fasta, "wt", encoding="ascii") as handle:
            handle.write(">chrGz\nACNN\nNNTA\n")

        runs = list(scan_poly_n_runs(fasta))

        self.assertEqual(
            [(r.contig, r.start, r.end, r.length) for r in runs],
            [("chrGz", 2, 6, 4)],
        )

    def test_bin_poly_n_lengths_uses_ranges_not_exact_lengths(self):
        summary = bin_poly_n_lengths([1, 10, 11, 100, 101, 1000, 1001])

        self.assertEqual([row.label for row in summary], ["1-10 bp", "11-100 bp", "101-1000 bp", ">1000 bp"])
        self.assertEqual([row.run_count for row in summary], [2, 2, 2, 1])
        self.assertEqual([row.total_bp for row in summary], [11, 111, 1101, 1001])

    def test_main_writes_poly_n_summary_split_by_contig(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / "ref.fa"
            fasta.write_text(
                f">chr1\nNNNAANNNNNNNNNNN\n>chr2\nA{'N' * 101}A{'N' * 1001}\n>chr3\nACGTACGT\n",
                encoding="ascii",
            )

            with contextlib.redirect_stdout(io.StringIO()):
                status = main([str(fasta), "--outdir", str(tmp), "--prefix", "ref"])

            self.assertEqual(status, 0)
            by_contig_tsv = tmp / "ref.poly_n_by_contig_summary.tsv"
            by_contig_png = tmp / "ref.poly_n_by_contig.png"
            self.assertTrue(by_contig_tsv.is_file())
            self.assertTrue(by_contig_png.is_file())
            self.assertEqual(by_contig_png.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
            with by_contig_tsv.open(newline="", encoding="ascii") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(
                [(row["contig"], row["bin"], row["run_count"], row["total_N_bp"]) for row in rows],
                [
                    ("chr1", "1-10 bp", "1", "3"),
                    ("chr1", "11-100 bp", "1", "11"),
                    ("chr1", "101-1000 bp", "0", "0"),
                    ("chr1", ">1000 bp", "0", "0"),
                    ("chr2", "1-10 bp", "0", "0"),
                    ("chr2", "11-100 bp", "0", "0"),
                    ("chr2", "101-1000 bp", "1", "101"),
                    ("chr2", ">1000 bp", "1", "1001"),
                    ("chr3", "1-10 bp", "0", "0"),
                    ("chr3", "11-100 bp", "0", "0"),
                    ("chr3", "101-1000 bp", "0", "0"),
                    ("chr3", ">1000 bp", "0", "0"),
                ],
            )


if __name__ == "__main__":
    unittest.main()
