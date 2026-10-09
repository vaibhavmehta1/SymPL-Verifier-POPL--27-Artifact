#!/usr/bin/env python3
"""
plot_scaling.py
===============

Plot the verifier's dimension-scaling results: wall-clock verification time vs.
phase-space dimension D, in two side-by-side panels --

    left   coupled   (N-DOF harmonic chain; the dimensions interact)
    right  regular   (N independent harmonic oscillators)

One colour per integrator, shared log-log axes.  Reads the CSV emitted by
`benchmark.py --plot` (columns method,variant,N,D,time_s,status) and writes a
PNG or PDF (format inferred from the output extension).

Run with a python whose matplotlib works (NOT `sage -python`):

    python3 plot_scaling.py scaling.csv scaling.pdf [budget_seconds]
"""

import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter


ORDER = ["Symplectic Euler", "Leapfrog",
         "Implicit midpoint", "Gauss-Legendre 2 (GL-2)"]

# Methods to omit from the scaling plot (substring match, case-insensitive).
SKIP = ["stoermer", "stormer", "verlet"]
PANELS = [("Coupled Oscillators", "Coupled Oscillators" ),
          ("Independent Oscillators", "Independent Oscillators")]


def load(csv_path):
    series = {}
    for row in csv.DictReader(open(csv_path)):
        variant = row.get("variant", "coupled")
        series.setdefault((row["method"], variant), []).append(
            (int(row["D"]), float(row["time_s"]), row["status"]))
    for key in series:
        series[key].sort(key=lambda r: r[0])
    return series


def main(argv):
    if len(argv) < 3:
        print("usage: python3 plot_scaling.py <csv> <out.png|out.pdf> [budget_s]")
        return 2
    csv_path, out_path = argv[1], argv[2]
    budget = float(argv[3]) if len(argv) > 3 else 120.0

    series = load(csv_path)
    if not series:
        print("no rows in", csv_path)
        return 1

    methods = sorted({m for m, _ in series
                      if not any(s in m.lower() for s in SKIP)},
                     key=lambda m: ORDER.index(m) if m in ORDER else 99)
    present = {v for _, v in series}
    panels = [(v, title) for v, title in PANELS if v in present]
    if not panels:
        panels = [(next(iter(present)), "Scaling")]
    all_N = sorted({n for pts in series.values() for n, _, _ in pts})

    cmap = plt.get_cmap("tab10")
    color_of = {m: cmap(i % 10) for i, m in enumerate(methods)}

    plt.rcParams.update({
        "font.size": 23, "axes.titlesize": 28, "axes.labelsize": 25,
        "xtick.labelsize": 21, "ytick.labelsize": 21, "legend.fontsize": 20,
    })
    fig, axes = plt.subplots(1, len(panels), figsize=(7.0 * len(panels), 5.8),
                             sharex=True, sharey=True, squeeze=False)
    axes = list(axes[0])

    for ax, (variant, title) in zip(axes, panels):
        for method in methods:
            pts = series.get((method, variant))
            if not pts:
                continue
            color = color_of[method]
            xs = [n for n, _, _ in pts]
            ys = [t for _, t, _ in pts]
            ax.plot(xs, ys, "-o", color=color, lw=2.3, markersize=7,
                    markeredgecolor="white", markeredgewidth=0.8, zorder=2)

        ax.axhline(budget, ls="--", color="0.45", lw=1.4, zorder=1)
        # ax.text(all_N[0], budget * 1.06, "2-min budget",
        #         ha="left", va="bottom", fontsize=15, color="0.35")
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xticks(all_N)
        ax.get_xaxis().set_major_formatter(ScalarFormatter())
        ax.set_xlabel("Dimension")
        ax.set_title(title)
        ax.grid(False)
        ax.margins(x=0.06)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    axes[0].set_ylabel("Time (s)")

    handles = [Line2D([], [], color=color_of[m], lw=2.8, marker="o",
                      markersize=8, markeredgecolor="white", label=m)
               for m in methods]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=20,
               framealpha=0.95, bbox_to_anchor=(0.5, -0.01))

    # fig.suptitle("Time vs. Phase-space dimension",
    #              y=0.98, fontsize=28)
    fig.tight_layout(rect=(0, 0.14, 1, 0.93))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print("wrote", out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
