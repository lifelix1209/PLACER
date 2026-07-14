#!/usr/bin/env python3
"""Publication figure for the TE-caller benchmark (benchmark_callers.py output).

Two panels sharing one colour-per-caller mapping:
  A. F1 per caller (horizontal bars, sorted)         -- headline magnitude
  B. precision-recall scatter with iso-F1 contours   -- the precision/recall trade-off

Design: colourblind-safe categorical palette assigned to callers in a fixed order;
identity is carried by direct labels (never colour alone); marker shape separates
TE-aware callers (circle) from general SV callers (square). Static, print-ready
(300 dpi PNG + vector PDF).
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# colourblind-safe categorical palette (validated, light surface), fixed order
PALETTE = ["#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7", "#e34948", "#e87ba4", "#eb6834"]
INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e6e5e2"
SURFACE = "#fcfcfb"

# which callers classify TE family (circle) vs general SV callers (square)
TE_AWARE = {"PLACER", "placer_py", "TLDR", "sTELLeR"}


def load(path: Path) -> list[dict]:
    rows = []
    with path.open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            r["precision"] = float(r["precision"])
            r["recall"] = float(r["recall"])
            r["f1"] = float(r["f1"])
            rows.append(r)
    return rows


def f1_contours(ax):
    r = np.linspace(0.001, 1, 300)
    p = np.linspace(0.001, 1, 300)
    R, P = np.meshgrid(r, p)
    F = 2 * P * R / (P + R)
    cs = ax.contour(R, P, F, levels=[0.2, 0.4, 0.6, 0.8, 0.9], colors=GRID, linewidths=1, zorder=1)
    ax.clabel(cs, fmt=lambda v: f"F1={v:.1f}", fontsize=7, colors=MUTED)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True, type=Path, help="benchmark_callers.py TSV")
    ap.add_argument("--out", required=True, type=Path, help="output PNG (a .pdf sibling is also written)")
    ap.add_argument("--radius", type=int, default=100)
    ap.add_argument("--title", default="Long-read TE-insertion caller benchmark (simulated chr8, 25x ONT)")
    args = ap.parse_args()

    rows = load(args.results)
    # fixed colour per caller by input order (colour follows the entity)
    colour = {r["caller"]: PALETTE[i % len(PALETTE)] for i, r in enumerate(rows)}

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    })
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={"width_ratios": [1, 1.15]})
    fig.suptitle(args.title, fontsize=12, fontweight="bold", x=0.5, y=0.99)

    # --- Panel A: F1 bars, sorted ascending so best sits on top ---
    srt = sorted(rows, key=lambda r: r["f1"])
    y = np.arange(len(srt))
    axA.barh(y, [r["f1"] for r in srt], color=[colour[r["caller"]] for r in srt],
             height=0.62, zorder=3)
    for yi, r in zip(y, srt):
        axA.text(r["f1"] + 0.012, yi, f"{r['f1']:.2f}", va="center", ha="left",
                 fontsize=9, color=INK)
    axA.set_yticks(y)
    axA.set_yticklabels([r["caller"] for r in srt], fontsize=9.5, color=INK)
    axA.set_xlim(0, 1.0)
    axA.set_xlabel("F1 score", fontsize=10)
    axA.set_title("A   F1 by caller", loc="left", fontsize=11, fontweight="bold")
    axA.grid(axis="x", color=GRID, linewidth=1, zorder=0)
    axA.set_axisbelow(True)
    for s in ("top", "right", "left"):
        axA.spines[s].set_visible(False)
    axA.tick_params(length=0)

    # --- Panel B: precision-recall scatter with F1 contours ---
    def draw_points(ax, label_which, offsets):
        for r in rows:
            marker = "o" if r["caller"] in TE_AWARE else "s"
            ax.scatter(r["recall"], r["precision"], s=130, marker=marker,
                       color=colour[r["caller"]], edgecolor=SURFACE, linewidth=1.5, zorder=5)
        for r in rows:
            if r["caller"] not in label_which:
                continue
            dx, dy = offsets.get(r["caller"], (8, 6))
            ax.annotate(r["caller"], (r["recall"], r["precision"]),
                        textcoords="offset points", xytext=(dx, dy), fontsize=9,
                        color=INK, fontweight="medium", zorder=6)

    f1_contours(axB)
    # on the main axes label the precision~1.0 spread; the top-right high-recall
    # cluster is labelled in the zoom inset so labels never overprint each other.
    draw_points(axB, {"TLDR", "SVIM", "sTELLeR", "placer_py"},
                {"TLDR": (8, 5), "SVIM": (8, 5), "sTELLeR": (-24, -18), "placer_py": (6, 8)})
    axB.set_xlim(0, 1.02)
    axB.set_ylim(0, 1.05)
    axB.set_xlabel("Recall", fontsize=10)
    axB.set_ylabel("Precision", fontsize=10)
    axB.set_title("B   Precision vs recall", loc="left", fontsize=11, fontweight="bold")
    for s in ("top", "right"):
        axB.spines[s].set_visible(False)
    axB.tick_params(length=0)

    # zoom inset for the three clustered high performers
    axi = axB.inset_axes([0.10, 0.30, 0.44, 0.44])
    f1_contours(axi)
    draw_points(axi, {"PLACER", "cuteSV", "Sniffles2"},
                {"PLACER": (7, 5), "cuteSV": (7, -14), "Sniffles2": (-64, 5)})
    axi.set_xlim(0.86, 1.005)
    axi.set_ylim(0.93, 1.008)
    axi.set_facecolor(SURFACE)
    axi.tick_params(length=0, labelsize=7)
    axi.set_xticks([0.9, 1.0]); axi.set_yticks([0.95, 1.0])
    for s in axi.spines.values():
        s.set_edgecolor(GRID)
    axB.indicate_inset_zoom(axi, edgecolor=MUTED, alpha=0.5, linewidth=1)

    # shared legend for marker shape
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=MUTED, markersize=9,
               label="TE-aware caller (classifies family)"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=MUTED, markersize=9,
               label="general SV caller (insertion only)"),
    ]
    axB.legend(handles=handles, loc="lower left", frameon=False, fontsize=8.5,
               handletextpad=0.4, borderpad=0.2)

    fig.text(0.5, 0.005,
             f"Truth = 120 planted TE insertions; a call is a true positive within {args.radius} bp of an unclaimed truth breakpoint.",
             ha="center", fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    fig.savefig(args.out, dpi=300, facecolor=SURFACE)
    fig.savefig(args.out.with_suffix(".pdf"), facecolor=SURFACE)
    print(f"wrote {args.out} and {args.out.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
