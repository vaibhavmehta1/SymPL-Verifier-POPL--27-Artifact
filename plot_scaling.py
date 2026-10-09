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

Run with a Python interpreter that has a working Matplotlib installation:

    python3 plot_scaling.py scaling.csv scaling.pdf [budget_seconds]
"""

import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter


ORDER = ["Symplectic Euler", "Leapfrog", "Stormer-Verlet",
         "Implicit midpoint", "Gauss-Legendre 2 (GL-2)"]

PANELS = [("coupled", "Coupled Oscillators"),
          ("regular", "Independent Oscillators")]
VARIANT_ALIASES = {title.lower(): key for key, title in PANELS}
METHOD_ALIASES = {name.lower(): name for name in ORDER}
METHOD_ALIASES.update({"leapfrog (kdk)": "Leapfrog",
                       "stoermer-verlet (dkd)": "Stormer-Verlet"})


def load(csv_path):
    series = {}
    with open(csv_path, newline="") as source:
        for row in csv.DictReader(source):
            variant = row.get("variant", "coupled").strip().lower()
            variant = VARIANT_ALIASES.get(variant, variant)
            method = row["method"].strip()
            method = METHOD_ALIASES.get(method.lower(), method)
            status = row["status"].strip().upper()
            if status == "SKIP":
                continue
            series.setdefault((method, variant), []).append(
                (int(row["D"]), float(row["time_s"]), status))
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

    methods = sorted({m for m, _ in series},
                     key=lambda m: (ORDER.index(m) if m in ORDER else 99, m))
    present = {v for _, v in series}
    panels = [(v, title) for v, title in PANELS if v in present]
    panels += [(v, v) for v in sorted(present - {v for v, _ in PANELS})]
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
            # Failed or timed-out checks are not successful verification times.
            ys = [t if status == "PASS" and t > 0 else float("nan")
                  for _, t, status in pts]
            ax.plot(xs, ys, "-o", color=color, lw=2.3, markersize=7,
                    markeredgecolor="white", markeredgewidth=0.8, zorder=2)
            for status, marker in (("TIMEOUT", "^"), ("FAIL", "x"), ("ERROR", "x")):
                marked = [(n, t) for n, t, s in pts if s == status and t > 0]
                if marked:
                    ax.plot([n for n, _ in marked], [t for _, t in marked],
                            linestyle="none", marker=marker, color=color,
                            markersize=9, zorder=3)

        ax.axhline(budget, ls="--", color="0.45", lw=1.4, zorder=1)
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
    statuses = {s for pts in series.values() for _, _, s in pts}
    for label, marker, present_status in (
            ("Timeout", "^", "TIMEOUT" in statuses),
            ("Failed/error", "x", bool(statuses & {"FAIL", "ERROR"}))):
        if present_status:
            handles.append(Line2D([], [], color="0.4", linestyle="none",
                                  marker=marker, markersize=9, label=label))
    handles.append(Line2D([], [], color="0.45", ls="--", label="Time budget"))
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=14,
               framealpha=0.95, bbox_to_anchor=(0.5, -0.01))

    legend_space = 0.07 * ((len(handles) + 1) // 2)
    fig.tight_layout(rect=(0, legend_space, 1, 1))
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
