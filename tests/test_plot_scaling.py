"""Plotting regressions; run with a Python that has Matplotlib installed."""

import csv
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import matplotlib.pyplot as plt

import benchmark
import plot_scaling


class ScalingPlotTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.addCleanup(plt.close, "all")
        self.root = Path(self.tmp.name)

    def write_csv(self, rows):
        path = self.root / "input.csv"
        with path.open("w", newline="") as output:
            writer = csv.writer(output)
            writer.writerow(["method", "variant", "N", "D", "time_s", "status"])
            writer.writerows(rows)
        return path

    def render(self, path):
        output = self.root / "nested" / "plot.png"
        with patch.object(plot_scaling.plt, "close"):
            self.assertEqual(plot_scaling.main(["plot_scaling.py", str(path),
                                                str(output), "30"]), 0)
        self.assertTrue(output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
        return plt.gcf()

    def test_current_runner_data_has_both_panels_and_every_method(self):
        # Exercise the real CSV writer and plotting subprocess, not a duplicate
        # of the runner's metadata/serialization conventions.
        benches = [b for b in benchmark.registry()
                   if b.family == "scaling" and b.meta["N"] == 2]
        output = self.root / "runner.pdf"
        benchmark.write_scaling_plot([(b, "PASS", 1.0) for b in benches],
                                     str(output), 30)
        self.assertTrue(output.read_bytes().startswith(b"%PDF"))
        fig = self.render(output.with_suffix(".csv"))
        self.assertEqual([ax.get_title() for ax in fig.axes],
                         ["Coupled Oscillators", "Independent Oscillators"])
        for ax in fig.axes:
            curves = [line for line in ax.lines if line.get_marker() == "o"]
            self.assertEqual(len(curves), 5)
            self.assertTrue(all(list(line.get_ydata()) == [1.0] for line in curves))

    def test_legacy_csv_keeps_both_panels_and_normalizes_methods(self):
        path = self.write_csv([
            ["Leapfrog ", "Coupled Oscillators", 2, 4, 1, "PASS"],
            ["Leapfrog", "Independent Oscillators", 2, 4, 2, "PASS"],
        ])
        fig = self.render(path)
        self.assertEqual(len(fig.axes), 2)
        self.assertEqual([t.get_text() for t in fig.legends[0].texts],
                         ["Leapfrog", "Time budget"])

    def test_unsuccessful_runs_are_not_drawn_as_passes(self):
        path = self.write_csv([
            ["symplectic Euler", "coupled", 2, 4, 1, "PASS"],
            ["symplectic Euler", "coupled", 4, 8, 30, "TIMEOUT"],
            ["symplectic Euler", "coupled", 8, 16, 2, "FAIL"],
            ["symplectic Euler", "coupled", 16, 32, 0, "ERROR"],
            ["symplectic Euler", "coupled", 32, 64, 0, "SKIP"],
        ])
        fig = self.render(path)
        lines = fig.axes[0].lines
        passes = next(line for line in lines if line.get_marker() == "o")
        self.assertEqual(list(passes.get_xdata()), [4, 8, 16, 32])
        self.assertEqual(passes.get_ydata()[0], 1)
        self.assertTrue(all(math.isnan(y) for y in passes.get_ydata()[1:]))
        timeout = next(line for line in lines if line.get_marker() == "^")
        self.assertEqual(list(timeout.get_xdata()), [8])
        self.assertEqual(list(timeout.get_ydata()), [30])
        failed = next(line for line in lines if line.get_marker() == "x")
        self.assertEqual(list(failed.get_xdata()), [16])

    def test_unknown_variants_are_not_silently_dropped(self):
        path = self.write_csv([
            ["custom method", "coupled", 2, 4, 1, "PASS"],
            ["custom method", "custom variant", 2, 4, 1, "PASS"],
        ])
        self.assertEqual([ax.get_title() for ax in self.render(path).axes],
                         ["Coupled Oscillators", "custom variant"])


if __name__ == "__main__":
    unittest.main()
