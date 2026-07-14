#!/usr/bin/env python3
"""Summarize poly-N runs in a FASTA reference genome."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence, TextIO


DEFAULT_BINS = (
    (1, 10, "1-10 bp"),
    (11, 100, "11-100 bp"),
    (101, 1000, "101-1000 bp"),
    (1001, None, ">1000 bp"),
)

BIN_COLORS = ("#4477AA", "#EE7733", "#228833", "#CC6677")


@dataclass(frozen=True)
class PolyNRun:
    contig: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class BinSummary:
    label: str
    min_bp: int
    max_bp: int | None
    run_count: int
    total_bp: int


@contextmanager
def open_fasta_text(fasta_path: Path) -> Iterator[TextIO]:
    if fasta_path.suffix == ".gz":
        with gzip.open(fasta_path, "rt", encoding="ascii") as handle:
            yield handle
    else:
        with fasta_path.open("rt", encoding="ascii") as handle:
            yield handle


def scan_fasta_contigs(fasta_path: Path) -> list[str]:
    """Return contig identifiers in FASTA order."""
    contigs: list[str] = []
    with open_fasta_text(fasta_path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">"):
                contigs.append(line[1:].split()[0])
    return contigs


def scan_poly_n_runs(fasta_path: Path) -> Iterator[PolyNRun]:
    """Yield consecutive N/n runs as 0-based half-open coordinates."""
    contig: str | None = None
    pos = 0
    run_start: int | None = None

    with open_fasta_text(fasta_path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if contig is not None and run_start is not None:
                    yield PolyNRun(contig, run_start, pos)
                contig = line[1:].split()[0]
                pos = 0
                run_start = None
                continue

            if contig is None:
                raise ValueError(f"Sequence data found before a FASTA header in {fasta_path}")

            for base in line:
                if base == "N" or base == "n":
                    if run_start is None:
                        run_start = pos
                elif run_start is not None:
                    yield PolyNRun(contig, run_start, pos)
                    run_start = None
                pos += 1

    if contig is not None and run_start is not None:
        yield PolyNRun(contig, run_start, pos)


def bin_poly_n_lengths(lengths: Iterable[int]) -> list[BinSummary]:
    counts = [0] * len(DEFAULT_BINS)
    totals = [0] * len(DEFAULT_BINS)

    for length in lengths:
        if length <= 0:
            continue
        for index, (min_bp, max_bp, _label) in enumerate(DEFAULT_BINS):
            if length >= min_bp and (max_bp is None or length <= max_bp):
                counts[index] += 1
                totals[index] += length
                break

    return [
        BinSummary(label, min_bp, max_bp, counts[index], totals[index])
        for index, (min_bp, max_bp, label) in enumerate(DEFAULT_BINS)
    ]


def summarize_runs_by_contig(contigs: Sequence[str], runs: Sequence[PolyNRun]) -> dict[str, list[BinSummary]]:
    lengths_by_contig: dict[str, list[int]] = {contig: [] for contig in contigs}
    for run in runs:
        lengths_by_contig.setdefault(run.contig, []).append(run.length)

    return {contig: bin_poly_n_lengths(lengths) for contig, lengths in lengths_by_contig.items()}


def write_runs_tsv(runs: Sequence[PolyNRun], output_path: Path) -> None:
    with output_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["contig", "start_0based", "end_0based_exclusive", "start_1based", "end_1based", "length_bp"])
        for run in runs:
            writer.writerow([run.contig, run.start, run.end, run.start + 1, run.end, run.length])


def write_summary_tsv(summary: Sequence[BinSummary], output_path: Path) -> None:
    with output_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["bin", "min_bp", "max_bp", "run_count", "total_N_bp"])
        for row in summary:
            writer.writerow([row.label, row.min_bp, "" if row.max_bp is None else row.max_bp, row.run_count, row.total_bp])


def write_contig_summary_tsv(summary_by_contig: Mapping[str, Sequence[BinSummary]], output_path: Path) -> None:
    with output_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["contig", "bin", "min_bp", "max_bp", "run_count", "total_N_bp"])
        for contig, summary in summary_by_contig.items():
            for row in summary:
                writer.writerow(
                    [
                        contig,
                        row.label,
                        row.min_bp,
                        "" if row.max_bp is None else row.max_bp,
                        row.run_count,
                        row.total_bp,
                    ]
                )


def configure_matplotlib_cache(output_path: Path) -> None:
    # Keep Matplotlib cache inside the workspace when the user home directory is not writable.
    os.environ.setdefault("MPLCONFIGDIR", str(output_path.parent / ".mplconfig"))
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("XDG_CACHE_HOME", str(output_path.parent / ".cache"))


def plot_summary(summary: Sequence[BinSummary], output_path: Path) -> None:
    configure_matplotlib_cache(output_path)
    import matplotlib.pyplot as plt

    labels = [row.label for row in summary]
    counts = [row.run_count for row in summary]
    total_bp = [row.total_bp for row in summary]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    bar_color = "#4477AA"
    edge_color = "#222222"

    axes[0].bar(labels, counts, color=bar_color, edgecolor=edge_color)
    axes[0].set_title("Poly-N run count")
    axes[0].set_ylabel("Number of runs")

    axes[1].bar(labels, total_bp, color="#66AA55", edgecolor=edge_color)
    axes[1].set_title("Total N bases")
    axes[1].set_ylabel("Total bp")

    for ax in axes:
        ax.set_xlabel("Run length bin")
        ax.tick_params(axis="x", rotation=25)
        ax.grid(axis="y", alpha=0.25)
        for container in ax.containers:
            ax.bar_label(container, fmt="%d", padding=3, fontsize=8)

    fig.suptitle("Reference genome poly-N distribution")
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_contig_summary(
    summary_by_contig: Mapping[str, Sequence[BinSummary]], output_path: Path, max_contigs: int = 30
) -> None:
    configure_matplotlib_cache(output_path)
    import matplotlib.pyplot as plt

    contig_rows = [
        (contig, sum(row.total_bp for row in summary), summary)
        for contig, summary in summary_by_contig.items()
        if sum(row.total_bp for row in summary) > 0
    ]
    contig_rows.sort(key=lambda item: item[1], reverse=True)
    contig_rows = contig_rows[:max_contigs]

    if not contig_rows:
        fig, ax = plt.subplots(figsize=(8, 2.5), layout="constrained")
        ax.axis("off")
        ax.text(0.5, 0.5, "No poly-N runs found", ha="center", va="center", fontsize=12)
        fig.savefig(output_path, dpi=300)
        plt.close(fig)
        return

    contigs = [contig for contig, _total, _summary in contig_rows]
    y_positions = list(range(len(contigs)))
    fig_height = max(4.0, 0.28 * len(contigs) + 1.4)
    fig, ax = plt.subplots(figsize=(10, fig_height), layout="constrained")

    left = [0] * len(contigs)
    bin_labels = [label for _min_bp, _max_bp, label in DEFAULT_BINS]
    for bin_index, label in enumerate(bin_labels):
        values = [summary[bin_index].total_bp for _contig, _total, summary in contig_rows]
        ax.barh(
            y_positions,
            values,
            left=left,
            color=BIN_COLORS[bin_index % len(BIN_COLORS)],
            edgecolor="#222222",
            linewidth=0.4,
            label=label,
        )
        left = [offset + value for offset, value in zip(left, values)]

    ax.set_yticks(y_positions, labels=contigs)
    ax.invert_yaxis()
    ax.set_xlabel("Total N bp")
    ax.set_ylabel("Contig")
    ax.set_title(f"Top {len(contigs)} contigs by total N bases")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(title="Run length bin", loc="lower right")

    for y_position, total_bp in zip(y_positions, left):
        ax.text(total_bp, y_position, f" {total_bp}", va="center", fontsize=8)

    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path, help="Input reference FASTA.")
    parser.add_argument("--outdir", type=Path, default=Path("placer_out/poly_n_distribution"), help="Output directory.")
    parser.add_argument("--prefix", default=None, help="Output filename prefix. Default: FASTA stem.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    fasta = args.fasta.expanduser().resolve()
    if not fasta.is_file():
        raise SystemExit(f"FASTA not found: {fasta}")

    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or fasta.name

    contigs = scan_fasta_contigs(fasta)
    runs = list(scan_poly_n_runs(fasta))
    summary = bin_poly_n_lengths(run.length for run in runs)
    summary_by_contig = summarize_runs_by_contig(contigs, runs)

    runs_tsv = outdir / f"{prefix}.poly_n_runs.tsv"
    summary_tsv = outdir / f"{prefix}.poly_n_summary.tsv"
    contig_summary_tsv = outdir / f"{prefix}.poly_n_by_contig_summary.tsv"
    plot_png = outdir / f"{prefix}.poly_n_histogram.png"
    contig_plot_png = outdir / f"{prefix}.poly_n_by_contig.png"

    write_runs_tsv(runs, runs_tsv)
    write_summary_tsv(summary, summary_tsv)
    write_contig_summary_tsv(summary_by_contig, contig_summary_tsv)
    plot_summary(summary, plot_png)
    plot_contig_summary(summary_by_contig, contig_plot_png)

    print(f"FASTA: {fasta}")
    print(f"Poly-N runs: {len(runs)}")
    print(f"Runs TSV: {runs_tsv}")
    print(f"Summary TSV: {summary_tsv}")
    print(f"Per-contig summary TSV: {contig_summary_tsv}")
    print(f"Histogram PNG: {plot_png}")
    print(f"Per-contig PNG: {contig_plot_png}")
    print()
    print("bin\trun_count\ttotal_N_bp")
    for row in summary:
        print(f"{row.label}\t{row.run_count}\t{row.total_bp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
