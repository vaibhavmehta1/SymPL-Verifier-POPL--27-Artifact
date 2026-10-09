#!/usr/bin/env python3
"""
benchmark.py
============

Heavy-duty evaluation harness for the symplectic verifier.

It runs the verifier over a curated benchmark suite, organised into three
families, TIMES each run on the wall clock, and -- crucially -- KILLS any run
that exceeds a hard per-benchmark budget (default 120 s).  The kill is done at
the process-GROUP level (setsid + killpg): the verifier offloads heavy algebra
to a Maxima/ECL child process, and an in-process `signal.alarm` cannot interrupt
a blocking C call inside that child (see the project memory).  Only tearing down
the whole process group reliably reclaims a wedged CAS computation.

Benchmark families
------------------
  TEXTBOOK   classic integrators (symplectic Euler, leapfrog, Stoermer-Verlet,
             implicit midpoint, GL-2, and Tao's explicit doubling scheme), each
             on increasingly complex Hamiltonians:
                 H1 harmonic oscillator  (1 DOF, linear)
                 H2 pendulum             (1 DOF, transcendental)
                 H3 Henon-Heiles         (2 DOF, coupled polynomial)
                 H4 double pendulum      (2 DOF, NONSEPARABLE, chaotic)
             The kick-drift methods stop at H3 (separable); H4 is reached only
             by implicit midpoint (a deliberate CAS timeout) and by Tao, whose
             doubling makes the nonseparable flow explicit and so verifies.
  CONFORMAL  conformal symplectic methods that contract omega by e^{-gam h}:
                 conformal Stoermer-Verlet (dissipative leapfrog)
                 conformal relativistic    (Franca-Jordan-Vidal RGD, 1st order)
  BESPOKE    the hand-built / domain maps:
                 Forest-Ruth, Tao, Boris, multipole kick, RESPA, exact bend.

Each benchmark is verified in a fresh subprocess:
  * a .spec file is checked by `sage -python parse.py <spec>`
    (exit 0 = VERIFIED, exit 1 = not verified);
  * the exact-bend map is symbolically intractable to simplify to zero, so it is
    certified numerically by `sage -python examples/exact_bend_certificate.py`
    (exit 0 = CERTIFIED) -- a CAS-performance limit, not a correctness one.

Usage
-----
    sage -python benchmark.py                  # run everything
    python3 benchmark.py                        # same (spawns `sage` per item)
    python3 benchmark.py --timeout 60           # tighter budget
    python3 benchmark.py --only textbook        # one family
    python3 benchmark.py --only gl2 henon       # substring filters (AND)
    python3 benchmark.py --list                 # list benchmarks, do not run
"""

import argparse
import glob
import os
import signal
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(ROOT, "examples")
BENCH = os.path.join(EX, "benchmarks")

DEFAULT_TIMEOUT = 120          # hard per-benchmark budget, seconds
KILL_GRACE = 5                 # SIGTERM -> SIGKILL grace, seconds
DEFAULT_PLOT = os.path.join(EX, "benchmarks", "scaling", "scaling.pdf")

SAGE = os.environ.get("SAGE", "sage")


class Bench:
    """One benchmark: a name, its family, and the command that verifies it.

    `expect` is the verdict a correct run must produce: "VERIFIED" for the
    structure-preserving maps, "CERTIFIED" for the numeric exact-bend check.
    """

    def __init__(self, family, name, argv, expect="VERIFIED", note="", meta=None):
        self.family = family
        self.name = name
        self.argv = argv
        self.expect = expect
        self.note = note
        self.meta = meta          # scaling rows carry {"method","N","D"}


def _spec(family, name, relpath, note=""):
    return Bench(family, name, [SAGE, "-python", "parse.py",
                                os.path.join(EX, relpath)], note=note)


def _tb(name, fname):
    return _spec("textbook", name, os.path.join("benchmarks", "textbook", fname))


def _cf(name, fname, note=""):
    return _spec("conformal", name,
                 os.path.join("benchmarks", "conformal", fname), note=note)


def registry():
    """The ordered benchmark suite."""
    b = []

    # -- TEXTBOOK: 5 methods x {harmonic, pendulum, henon-heiles} ----------
    methods = [
        ("symplectic_euler", "symplectic Euler"),
        ("leapfrog",         "leapfrog (KDK)"),
        ("stormer_verlet",   "Stoermer-Verlet (DKD)"),
        ("implicit_midpoint", "implicit midpoint"),
        ("gl2",              "Gauss-Legendre 2 (GL-2)"),
    ]
    hams = [
        ("harmonic",      "H1 harmonic"),
        ("pendulum",      "H2 pendulum"),
        ("henon_heiles",  "H3 Henon-Heiles"),
    ]
    for mkey, mlabel in methods:
        for hkey, hlabel in hams:
            b.append(_tb("{} / {}".format(mlabel, hlabel),
                         "{}_{}.spec".format(mkey, hkey)))
    # H4: double pendulum -- nonseparable & chaotic, so the separable kick-drift
    # methods do not apply.  Two methods reach it:
    #   * implicit midpoint -- its 4-D rational-trig stage solve is CAS-
    #     intractable to simplify to zero, so it is EXPECTED to exhaust the
    #     budget (a deliberate kill demo);
    #   * Tao's doubling -- makes the same nonseparable H EXPLICITLY integrable
    #     (no stage solve), so its map is genuinely symplectic where midpoint
    #     times out.  Its five composed rational-trig flows still overwhelm
    #     Maxima's simplify-to-zero (same CAS wall as exact_bend), so it is
    #     certified NUMERICALLY -- a CAS-performance limit, not a correctness one.
    b.append(_tb("implicit midpoint / H4 double pendulum",
                 "implicit_midpoint_double_pendulum.spec"))
    b.append(Bench("textbook", "Tao (doubling) / H4 double pendulum",
                   [SAGE, "-python",
                    os.path.join(EX, "tao_double_pendulum_certificate.py")],
                   expect="CERTIFIED",
                   note="explicit but CAS-heavy; certified numerically"))

    # -- TAO on H1..H3: Tao's explicit doubling scheme also handles the classic
    # separable Hamiltonians, verified on the same doubled phase space.
    for hkey, hlabel in hams:
        b.append(_tb("Tao (doubling) / {}".format(hlabel),
                     "tao_{}.spec".format(hkey)))

    # -- SCALING (textbook methods on an N-DOF coupled harmonic chain) -----
    # Auto-discovered from examples/benchmarks/scaling/ (regenerate with
    # gen_scaling.py).  Files are named "<method>_chain<NN>.spec".
    morder = ["symplectic_euler", "leapfrog", "stormer_verlet",
              "implicit_midpoint", "gl2"]
    mlabels = dict(zip(morder,
                       ["symplectic Euler", "leapfrog (KDK)",
                        "Stoermer-Verlet (DKD)", "implicit midpoint",
                        "Gauss-Legendre 2 (GL-2)"]))

    vorder = {"coupled": 0, "regular": 1}

    def _parse_scale(path):
        # "<variant>_<method>_chain<NN>.spec"
        stem = os.path.basename(path)[:-5]
        left, n = stem.rsplit("_chain", 1)
        variant, method = left.split("_", 1)
        return variant, method, int(n)

    def _scale_key(path):
        variant, method, N = _parse_scale(path)
        return (morder.index(method) if method in morder else 99,
                vorder.get(variant, 9), N)

    for path in sorted(glob.glob(os.path.join(BENCH, "scaling", "*.spec")),
                       key=_scale_key):
        variant, method, N = _parse_scale(path)
        base = mlabels.get(method, method)
        label = "{} [{}] / N={} (D={})".format(base, variant, N, 2 * N)
        b.append(Bench("scaling", label, [SAGE, "-python", "parse.py", path],
                       meta={"method": base, "variant": variant,
                             "N": N, "D": 2 * N}))

    # -- CONFORMAL ---------------------------------------------------------
    b.append(_cf("conformal Stoermer-Verlet", "conformal_sv.spec",
                 note="contracts omega by e^{-gam h}"))
    b.append(_cf("conformal relativistic (RGD)", "conformal_relativistic.spec",
                 note="contracts omega by e^{-gam h}"))

    # -- BESPOKE -----------------------------------------------------------
    b.append(_spec("bespoke", "Forest-Ruth (4th order)", "forest_ruth.spec"))
    b.append(_spec("bespoke", "Tao (nonseparable)", "tao.spec"))
    b.append(_spec("bespoke", "Boris pusher (volume form)", "boris.spec"))
    b.append(_spec("bespoke", "multipole kick", "multipole_kick.spec"))
    b.append(_spec("bespoke", "RESPA (multiple time step)", "respa.spec"))
    b.append(Bench("bespoke", "exact bend (numeric certificate)",
                   [SAGE, "-python",
                    os.path.join(EX, "exact_bend_certificate.py")],
                   expect="CERTIFIED",
                   note="symbolic simplify intractable; certified numerically"))
    return b


def run_one(bench, timeout):
    """Run one benchmark in its own process group; kill the group on timeout.

    Returns (status, elapsed_seconds, returncode_or_None).  status is one of
    PASS, FAIL, TIMEOUT, ERROR.
    """
    start = time.monotonic()
    try:
        proc = subprocess.Popen(
            bench.argv, cwd=ROOT,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            start_new_session=True,           # new process group => killable as a unit
        )
    except FileNotFoundError as exc:
        return "ERROR", 0.0, None, "cannot launch: {}".format(exc)

    pgid = os.getpgid(proc.pid)
    try:
        out, _ = proc.communicate(timeout=timeout)
        elapsed = time.monotonic() - start
        text = out.decode("utf-8", "replace") if out else ""
        if proc.returncode == 0:
            return "PASS", elapsed, proc.returncode, text
        return "FAIL", elapsed, proc.returncode, text
    except subprocess.TimeoutExpired:
        # Tear down the whole group: SIGTERM, brief grace, then SIGKILL.
        _killpg(pgid, signal.SIGTERM)
        try:
            proc.communicate(timeout=KILL_GRACE)
        except subprocess.TimeoutExpired:
            _killpg(pgid, signal.SIGKILL)
            try:
                proc.communicate(timeout=KILL_GRACE)
            except subprocess.TimeoutExpired:
                pass
        elapsed = time.monotonic() - start
        return "TIMEOUT", elapsed, None, ""


def _killpg(pgid, sig):
    try:
        os.killpg(pgid, sig)
    except ProcessLookupError:
        pass


def _verdict(bench, status):
    """Map a raw run status to PASS/FAIL given what the benchmark should yield."""
    if status == "PASS":
        return "PASS"
    return status


# ---------------------------------------------------------------------------
#  Pretty-printing
# ---------------------------------------------------------------------------

_TTY = sys.stdout.isatty()
_COLOR = {"PASS": "\033[32m", "FAIL": "\033[31m",
          "TIMEOUT": "\033[33m", "ERROR": "\033[35m"}
_RESET = "\033[0m"
_DIM = "\033[2m"


def _paint(plain_padded, key):
    """Wrap an already-padded string in the status colour (TTY only)."""
    if _TTY and key in _COLOR:
        return _COLOR[key] + plain_padded + _RESET
    return plain_padded


def render_table(headers, rows, aligns):
    """Render a box-drawn table.

    `rows` is a list of plain-string tuples; `aligns` gives 'l'/'r'/'c' per
    column.  Any cell whose text equals a known status (PASS/FAIL/...) is
    coloured.  Width maths is done on the plain text, so colour never disturbs
    alignment.
    """
    widths = [len(h) for h in headers]
    for r in rows:
        for i, cell in enumerate(r):
            widths[i] = max(widths[i], len(cell))

    def rule(left, mid, right):
        return left + mid.join("─" * (w + 2) for w in widths) + right

    def pad(cell, i):
        w = widths[i]
        if aligns[i] == "r":
            return cell.rjust(w)
        if aligns[i] == "c":
            return cell.center(w)
        return cell.ljust(w)

    def row(cells, header=False):
        parts = []
        for i, cell in enumerate(cells):
            s = pad(cell, i)
            key = cell.strip()
            if not header:
                s = _paint(s, key)
            parts.append(" " + s + " ")
        return "│" + "│".join(parts) + "│"

    out = [rule("┌", "┬", "┐"),
           row(headers, header=True),
           rule("├", "┼", "┤")]
    out.extend(row(r) for r in rows)
    out.append(rule("└", "┴", "┘"))
    return "\n".join(out)


def _find_plot_python():
    """Locate a python interpreter whose matplotlib actually imports.

    Try an explicit $PLOT_PYTHON override, the running interpreter (including
    virtual environments), Conda environments, then system Python. Each
    candidate is probed with a quick `import matplotlib`.
    """
    cands = []
    if os.environ.get("PLOT_PYTHON"):
        cands.append(os.environ["PLOT_PYTHON"])
    cands.append(sys.executable)
    for base in (os.environ.get("CONDA_PREFIX"),
                 os.path.expanduser("~/anaconda3"),
                 os.path.expanduser("~/miniconda3")):
        if base and os.path.isdir(base):
            cands.append(os.path.join(base, "bin", "python3"))
            cands.extend(sorted(glob.glob(
                os.path.join(base, "envs", "*", "bin", "python3"))))
    cands += ["python3", "/usr/bin/python3"]

    env = {k: v for k, v in os.environ.items()
           if k not in ("PYTHONPATH", "PYTHONHOME")}
    seen = set()
    for c in cands:
        if not c or c in seen:
            continue
        seen.add(c)
        try:
            r = subprocess.run([c, "-c", "import matplotlib"], env=env,
                               stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)
            if r.returncode == 0:
                return c, env
        except Exception:
            continue
    return None, env


def write_scaling_plot(results, plot_path, budget):
    """Emit a CSV of the scaling-family results and render a PDF or PNG.

    Plotting uses an interpreter with a working Matplotlib installation. The
    CSV is always written, so the plot can be regenerated by hand at any time.
    """
    rows = [(b.meta["method"], b.meta.get("variant", "coupled"),
             b.meta["N"], b.meta["D"], t, v)
            for b, v, t in results
            if b.family == "scaling" and b.meta and v != "SKIP"]
    if not rows:
        print("\n  --plot: no scaling benchmarks in this run; "
              "nothing to plot (try --only scaling --plot).")
        return

    out_dir = os.path.dirname(os.path.abspath(plot_path))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    csv_path = os.path.splitext(plot_path)[0] + ".csv"
    with open(csv_path, "w") as fh:
        fh.write("method,variant,N,D,time_s,status\n")
        for method, variant, N, D, t, v in rows:
            fh.write("{},{},{},{},{:.3f},{}\n".format(method, variant, N, D, t, v))
    print("\n  wrote scaling data : {}".format(csv_path))

    # Render with an interpreter that has a working Matplotlib installation.
    py, env = _find_plot_python()
    if py is None:
        print("  no python with a working matplotlib found; CSV kept at {}\n"
              "  (set PLOT_PYTHON, or run: <py> plot_scaling.py {} {} {:.0f})"
              .format(csv_path, csv_path, plot_path, budget))
        return
    proc = subprocess.run([py, os.path.join(ROOT, "plot_scaling.py"),
                           csv_path, plot_path, str(budget)],
                          cwd=ROOT, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if proc.returncode == 0:
        print("  wrote scaling plot : {}  (via {})".format(plot_path, py))
    else:
        print("  could not render plot ({} exit {}); CSV kept at {}"
              .format(py, proc.returncode, csv_path))
        tail = proc.stdout.decode("utf-8", "replace").rstrip().splitlines()[-4:]
        for line in tail:
            print("        | " + line)


def main(argv):
    ap = argparse.ArgumentParser(description="Symplectic-verifier benchmark suite")
    ap.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT,
                    help="hard per-benchmark budget in seconds (default %(default)s)")
    ap.add_argument("--only", nargs="+", default=None, metavar="FILTER",
                    help="run only benchmarks whose family or name contains ALL "
                         "of these substrings (case-insensitive)")
    ap.add_argument("--list", action="store_true",
                    help="list the selected benchmarks and exit")
    ap.add_argument("--verbose", action="store_true",
                    help="echo each benchmark's captured output")
    ap.add_argument("--plot", nargs="?", const=DEFAULT_PLOT, default=None,
                    metavar="PNG",
                    help="after running, write a scaling plot (time vs N) and a "
                         "CSV of the scaling results; with no argument defaults "
                         "to " + DEFAULT_PLOT)
    args = ap.parse_args(argv[1:])

    benches = registry()
    if args.only:
        needles = [s.lower() for s in args.only]
        benches = [b for b in benches
                   if all(n in (b.family + " " + b.name).lower() for n in needles)]
    if not benches:
        print("no benchmarks match filter {!r}".format(args.only))
        return 2

    if args.list:
        for b in benches:
            print("  [{:9s}] {}".format(b.family, b.name))
        print("\n{} benchmark(s).".format(len(benches)))
        return 0

    print("=" * 78)
    print("SYMPLECTIC VERIFIER -- BENCHMARK SUITE")
    print("  benchmarks : {}".format(len(benches)))
    print("  budget     : {:.0f} s per benchmark (process-group kill)".format(args.timeout))
    print("=" * 78)

    results = []
    cur_family = None
    walled = {}            # (method, variant) -> N at which a scaling sweep first timed out
    for b in benches:
        if b.family != cur_family:
            cur_family = b.family
            print("\n  {} ---".format(cur_family.upper()))
        sys.stdout.write("    {:46s} ... ".format(b.name[:46]))
        sys.stdout.flush()

        # once a scaling sweep has timed out, skip its larger dimensions
        skey = (b.meta["method"], b.meta["variant"]) if (
            b.family == "scaling" and b.meta) else None
        if skey is not None and skey in walled and b.meta["N"] > walled[skey]:
            print("{}  (past wall, N>{})".format(
                _paint("{:7s}".format("SKIP"), "SKIP"), walled[skey]))
            results.append((b, "SKIP", 0.0))
            continue

        status, elapsed, rc, text = run_one(b, args.timeout)
        verdict = _verdict(b, status)
        print("{}  {:7.2f}s".format(_paint("{:7s}".format(verdict), verdict), elapsed))
        if skey is not None and verdict == "TIMEOUT":
            walled[skey] = b.meta["N"]
        if args.verbose and text:
            for line in text.rstrip().splitlines():
                print("        | " + line)
        elif verdict in ("FAIL", "ERROR") and text:
            for line in text.rstrip().splitlines()[-3:]:
                print("        | " + line)
        results.append((b, verdict, elapsed))

    # -- results table -----------------------------------------------------
    rows = []
    prev_fam = None
    for b, v, t in results:
        fam = b.family if b.family != prev_fam else ""
        prev_fam = b.family
        tcell = "-" if v == "SKIP" else "{:.2f}s".format(t)
        rows.append((fam, b.name, v, tcell))
    print("\nRESULTS")
    print(render_table(["Family", "Benchmark", "Status", "Time"],
                        rows, ["l", "l", "c", "r"]))

    # -- summary table -----------------------------------------------------
    n_pass = sum(1 for _, v, _ in results if v == "PASS")
    n_fail = sum(1 for _, v, _ in results if v == "FAIL")
    n_to = sum(1 for _, v, _ in results if v == "TIMEOUT")
    n_err = sum(1 for _, v, _ in results if v == "ERROR")
    n_skip = sum(1 for _, v, _ in results if v == "SKIP")
    total_t = sum(t for _, _, t in results)
    by_fam = {}                          # family -> [n_run, n_pass, time]
    for b, v, t in results:
        d = by_fam.setdefault(b.family, [0, 0, 0.0])
        if v != "SKIP":
            d[0] += 1
        d[1] += 1 if v == "PASS" else 0
        d[2] += t
    srows = []
    for fam in ("textbook", "scaling", "conformal", "bespoke"):
        if fam in by_fam:
            n, p, t = by_fam[fam]
            srows.append((fam, "{}/{}".format(p, n), "{:.2f}s".format(t)))
    n_run = n_pass + n_fail + n_to + n_err
    srows.append(("TOTAL", "{}/{}".format(n_pass, n_run),
                  "{:.2f}s".format(total_t)))
    skip_note = ", {} skipped (past wall)".format(n_skip) if n_skip else ""
    print("\nSUMMARY  ({} pass, {} fail, {} timeout, {} error{})"
          .format(n_pass, n_fail, n_to, n_err, skip_note))
    print(render_table(["Family", "Passed", "Time"], srows, ["l", "r", "r"]))

    slow = sorted(results, key=lambda r: -r[2])[:3]
    print("\n  slowest: " + ", ".join(
        "{} ({:.1f}s, {})".format(b.name, t, v) for b, v, t in slow))

    if args.plot is not None:
        write_scaling_plot(results, args.plot, args.timeout)

    return 0 if (n_fail == 0 and n_err == 0) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
