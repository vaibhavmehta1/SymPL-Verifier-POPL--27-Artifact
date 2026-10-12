# Symplectic verifier artifact

Checks whether a program's pullback preserves a differential form: `C* B = A`.
Using `pre = post = dq ^ dp` checks canonical symplecticity.

## Installation

The verifier needs SageMath; plotting additionally needs Matplotlib. Choose one
installation below. There is no project build step. On Windows, use WSL2; the
benchmark runner requires POSIX process-group signals.

**Without Conda — Ubuntu 22.04 or Debian 12:**

```bash
sudo apt-get update
sudo apt-get install -y sagemath python3-matplotlib python3-venv
```

These releases provide SageMath 9.5
([Ubuntu package](https://packages.ubuntu.com/jammy/sagemath),
[Debian package](https://packages.debian.org/bookworm/sagemath)). On Ubuntu,
if the package is unavailable, enable Universe with
`sudo add-apt-repository universe`, then repeat the commands. Package availability
on other OS releases differs; use the
[Sage installation guide](https://doc.sagemath.org/html/en/installation/index.html).

**With Conda — Linux or macOS:**

If Conda is missing, install [Miniforge](https://github.com/conda-forge/miniforge#install)
first. In Bash, using the default installation location:

```bash
curl -fL -o Miniforge3.sh \
  "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3.sh
source "$HOME/miniforge3/etc/profile.d/conda.sh"
```

Then create the environment, following the
[Sage Conda installation route](https://doc.sagemath.org/html/en/installation/conda.html):

```bash
conda create -n symplectic-artifact --override-channels -c conda-forge sage matplotlib
conda activate symplectic-artifact
```

Reactivate it in new terminals. Do not use `pip install sage` as a substitute.

**Check either installation**, from this source directory:

```bash
sage --version
sage -python -c 'from sage.all import SR; print(SR(1))'
python3 -c 'import matplotlib; print(matplotlib.__version__)'
sage -python parse.py examples/benchmarks/textbook/symplectic_euler_pendulum.spec
```

Expect `1`, a Matplotlib version, and `RESULT : VERIFIED` from the last three
commands. If Sage cannot write its user directory, run
`export DOT_SAGE="$PWD/.sage-runtime"` first.

If system Matplotlib has a NumPy import/ABI error, install plotting dependencies
in a separate virtual environment; Sage can remain installed through apt:

```bash
python3 -m venv .venv-plot
.venv-plot/bin/python -m pip install matplotlib
export PLOT_PYTHON="$PWD/.venv-plot/bin/python"
```

Use `.venv-plot/bin/python` for standalone plotting and plotting tests below.
`PLOT_PYTHON` selects it for `benchmark.py --plot`; `SAGE=/path/to/sage` can
similarly select the verifier executable.

## Run

```bash
sage -python symplectic_verifier.py       # core regression tests
sage -python parse.py                    # parser demos
sage -python -m unittest discover -s tests -p test_implicit_solve.py -v
python3 benchmark.py --only textbook harmonic --timeout 30
python3 benchmark.py --list
python3 benchmark.py --timeout 120 > benchmark-results.txt 2>&1
```

The first two finish with `ALL TESTS PASSED` and
`ALL DEMOS BEHAVED AS EXPECTED`. Deliberately negative cases print `FAILED`.
The harmonic selection has six cases. The full registry has 109 cases and can
run for many minutes. Filters match family/display names, case-insensitively;
multiple filters use AND, e.g. `--only textbook gl-2 harmonic`.

Run GL-2 double pendulum with
`python3 benchmark.py --only textbook gl-2 double --timeout 120`.
Its input is `examples/benchmarks/textbook/gl2_double_pendulum.spec`;
symbolic verification may time out.

For a single `.spec`, `parse.py` exits 0 on verification, 1 on a residual;
syntax/runtime errors also exit nonzero. The benchmark runner reports `PASS`,
`FAIL`, `TIMEOUT`, `ERROR`, or `SKIP`. **Timeouts and skips alone still produce
exit 0**: inspect its summary. Larger dimensions are skipped after the first
scaling timeout for each method/variant.

## Reproduce the evaluation (Table 1 and Figure 10)

Section 6 of [the supplied paper](<popl27-paper75 (1).pdf>) has one results
table (Table 1, p. 23) and one plot (Figure 10, p. 24). After installation, run
the following **from the repository root** to measure all registered cases and
create the plot. If plotting uses a separate environment, set `PLOT_PYTHON` as
shown in Installation first.

```bash
mkdir -p results/evaluation

# Save the machine and software versions used for these measurements.
{
  date -u
  uname -a
  command -v lscpu >/dev/null && lscpu
  sage --version
  sage -python -c 'import sys; print(sys.version)'
  "${PLOT_PYTHON:-python3}" -c 'import matplotlib; print(matplotlib.__version__)'
} > results/evaluation/environment.txt 2>&1

# Table 1 cases, Figure 10 cases, and the artifact's additional benchmarks.
python3 benchmark.py --timeout 120 \
  --plot results/evaluation/all-scaling.pdf > results/evaluation/run.log 2>&1

# Select Figure 10's four methods from the fresh scaling measurements.
"${PLOT_PYTHON:-python3}" plot_scaling.py \
  results/evaluation/all-scaling.csv results/evaluation/figure10.pdf 120 --paper
```

This schedules 109 cases. In a second terminal, use
`tail -f results/evaluation/run.log` to monitor progress. Run on an otherwise
idle machine; allow up to roughly four hours at the 120-second budget, though
successful cases and skipping larger dimensions after timeouts usually shorten
the run substantially.

| Paper result | Where to find the fresh result |
| --- | --- |
| Table 1: general-purpose integrators | `run.log`, `RESULTS` rows in family `textbook`; omit the extra `Stoermer-Verlet (DKD)` rows |
| Table 1: conformal integrators | `run.log`, family `conformal` |
| Table 1: bespoke integrators | `run.log`, family `bespoke`; omit the extra `Tao (nonseparable)` row |
| Figure 10: both oscillator panels | `figure10.pdf`; underlying measurements are in `all-scaling.csv` |

Table 1 calls `Forest-Ruth` **Ruth–Forest**, `multipole kick` **SciBMad1**,
`exact bend` **SciBMad2**, and `RESPA (multiple time step)` **r-RESPA**. Its
Euler/leapfrog double-pendulum cells are **N/A**, so no cases run for those cells.
GL-2 double pendulum, Tao double pendulum, and exact bend all use the standard
symbolic verifier. Keep `TIMEOUT` as a timeout when assembling the table;
it is not a successful verification. Check the log for `FAIL` and `ERROR` too.

**Limits:** the runner emits a text results table, not the paper's typeset
Table 1 or its LoC column. The paper counts an assignment to `n` variables as
`n` lines; `wc -l` on a `.spec` does not reproduce that metric. The original
LoC-counting script and hand-simplified exact-bend proof are not included.
These commands rerun the available inputs; they do not establish identical
timings, source-line counts, or full paper fidelity (see below).

To rerun only a Table 1 family, without repeating scaling:

```bash
python3 benchmark.py --only textbook --timeout 120 > results/evaluation/textbook.log 2>&1
python3 benchmark.py --only conformal --timeout 120 > results/evaluation/conformal.log 2>&1
python3 benchmark.py --only bespoke --timeout 120 > results/evaluation/bespoke.log 2>&1
```

### Figure 10 only

These instructions target Section 6, Figure 10 (pp. 24–25) of
[the supplied paper](<popl27-paper75 (1).pdf>). Run from the project root after
installation. Use `PLOT_PYTHON` from the installation section if plotting needs
its own environment.

**Replot the supplied measurements** (quick; no new verification):

```bash
"${PLOT_PYTHON:-python3}" plot_scaling.py \
  examples/benchmarks/scaling/scaling.csv results/paper-reference/figure10.pdf 120 --paper
```

The saved CSV agrees with the paper's reported timeout thresholds. `--paper`
selects its four curves and colour order; timeout points use triangles rather
than the paper's connected circular markers. This recreates the data view,
not a pixel-identical figure or an independent rerun.

**Rerun every Figure 10 scaling series** (fresh measurements):

```bash
mkdir -p results/paper-rerun
python3 benchmark.py --paper-scaling --list  # 64 scheduled cases

# Record the environment alongside the results.
{
  date -u
  uname -a
  lscpu
  sage --version
  sage -python -c 'import sys; print(sys.version)'
} > results/paper-rerun/environment.txt

python3 benchmark.py --paper-scaling --timeout 120 \
  --plot results/paper-rerun/figure10.pdf > results/paper-rerun/run.log 2>&1

# Optional PNG, or retry rendering after fixing Matplotlib.
"${PLOT_PYTHON:-python3}" plot_scaling.py \
  results/paper-rerun/figure10.csv results/paper-rerun/figure10.png 120 --paper
```

The runner emits the PDF, CSV (`method,variant,N,D,time_s,status`), and log.
Read the log summary: timeouts/skips do not cause a nonzero exit. Run one sweep
at a time on an otherwise idle machine. Timings include Sage startup, parsing,
verification, and output. Exact timings and timeout thresholds depend on the
machine and Sage/CAS versions; the PDF does not specify a complete environment.

The experiment is:

| Setting | Figure 10 reproduction |
| --- | --- |
| Methods | Symplectic Euler, leapfrog (KDK), implicit midpoint, GL-2; no Tao or extra DKD Störmer–Verlet curve |
| Oscillators `N` | 2, 4, 8, 16, 32, 64, 128, 256 |
| Physical dimension `D = 2N` | 4, 8, 16, 32, 64, 128, 256, 512; excludes implicit-stage scratch coordinates |
| Independent (`regular`) | `H = sum_i (p_i**2 + q_i**2)/2` |
| Coupled (`coupled`) | Same `H` plus `(k/2)*sum_i (q_(i+1)-q_i)**2`; open chain, symbolic `k` and `h` |
| Budget | 120 seconds per case; larger dimensions skipped after the first timeout per method/variant |

The 64 cases are 4 methods × 2 systems × 8 sizes. The saved paper-method data
contains 48 measured rows (44 passes, 4 timeouts); the remaining 16 were skipped.
Reported first timeout dimensions are:

| Method | Coupled | Independent |
| --- | --- | --- |
| Implicit midpoint | 16 | 256 |
| GL-2 | 4 | 64 |

The explicit curves reach `D=512`. A run can take many minutes; the nominal
64-case budget is 128 minutes before timeout cleanup/plotting, although early
skipping usually makes it much shorter. Inputs are already included; optionally
run `python3 gen_scaling.py` to regenerate them (replaces all scaling `.spec`
files, including the extra DKD method). `--only scaling` runs all 80 cases;
`--paper-scaling` selects only the 64 Figure 10 cases.

For a quick pipeline check, add `--only 'N=2 ('` to the paper command and use a
separate output directory. Plotting regressions run with
`"${PLOT_PYTHON:-python3}" -m unittest discover -s tests -v`.

## Fidelity to the paper

**The scaling inputs match the described experiment, but the verifier is not a
faithful implementation of the entire paper.** The two audited implicit-solver
symbol bugs are fixed: each unresolved solve occurrence gets fresh symbols,
including repeated commands and loop iterations. Caller-supplied `Solve(phi=...)`
symbols are now rejected; omit that argument. The paper's indexed-loop semantics
and Section 5 SymplErr bounds remain unimplemented. The registry includes GL-2
double pendulum and uses symbolic verification for all cases, but a fresh full
Table 1 run has not been validated.

See [PAPER_AUDIT.md](PAPER_AUDIT.md) for the paper-to-code comparison,
counterexamples, and limits of the checks. The earlier audit's
`tests/paper_fidelity.py` driver is not present in this checkout. The included
implicit-solve regressions run with the Sage unittest command above. Check the
GL-2 double-pendulum equations separately with
`sage -python -m unittest discover -s tests -p test_gl2_double_pendulum.py -v`.
Unresolved solution values remain conservative symbolic placeholders, so some
valid identities requiring relations between solution values may not verify.
Their state dependence is not available to `Form.d()` on an unsimplified
returned pullback. The scaling inputs' `N=2` update maps passed independent
exact-arithmetic reference checks.

## Write and extend specifications

Save this as `my_integrator.spec`, then run
`sage -python parse.py my_integrator.spec`:

```text
coords q p
params h
pre  dq ^ dp
post dq ^ dp
program {
    p := p - h*sin(q)
    q := q + h*p
}
```

Assignments execute in order; tuple assignments evaluate their right-hand sides
simultaneously. Declare all state/scratch variables in `coords` and constants in
`params`. Use `**` for powers, `^` for wedge products, and `d(expr)` for exterior
derivatives. Statements need newlines or semicolons; `#` starts a comment.

Implicit assignment `(Q, P) := solve (a, b) -> (b/2, -sin(a)/2)` means
`(Q_new, P_new) = (Q_old, P_old) + h*(b/2, -sin(a)/2)`, with `a, b` standing for
new values. It requires parameter `h` and declared targets `Q, P`; choose unknown
names distinct from coordinates/parameters. See the implicit-midpoint examples.
`while^N` and `for^T` repeat a body a literal number of times; neither automatically
increments a symbolic index. Indexed assertions use the Python `verify_for` API.

- **New integrator/system:** add a `.spec` under `examples/`; existing assignments,
  solves, and loops need no core changes. Use textbook cases as templates.
- **New benchmark:** add it to `benchmark.py:registry()` using `_tb`, `_cf`, or
  `_spec`. Only scaling files are automatically discovered. Negative tests need
  a driver that asserts rejection: `Bench.expect` does not invert exit status.
- **New scalar function:** register its Sage callable in `parse.py:FUNCTIONS`.
- **New command:** implement its pullback in `symplectic_verifier.py:Verifier.pull`
  and add parser syntax in `parse.py`. Add positive and negative regression cases
  to the built-in suites.
- **Change scaling dimensions/methods:** edit `gen_scaling.py` and run
  `python3 gen_scaling.py`. This replaces all scaling `.spec` files. The RESPA
  generators similarly overwrite their inputs; `gen_respa_ff.py` verifies first.

For Python integration, `parse.parse_spec(text)` returns
`(ctx, pre, program, post)` and `Verifier(ctx).verify(pre, program, post)` returns
`(ok, pullback, residual)`. Run drivers with `sage -python`; use exact Sage
constants such as `QQ(1)/3`, since plain Python `1/3` is a float. See the core
self-tests for `Context`, `Form`, `Assign`, `Solve`, and `verify_for` examples.

## Package the artifact

| Include | Purpose |
| --- | --- |
| `symplectic_verifier.py` | Essential verifier library; requires SageMath |
| `parse.py` and selected `.spec` files | Essential for the specification CLI |
| `examples/` and `benchmark.py` | Evaluation inputs and benchmark runner |
| `gen_scaling.py`, `gen_respa_md.py`, `gen_respa_ff.py` | Input regeneration |
| `plot_scaling.py`, `tests/`, `README.md`, `PAPER_AUDIT.md` | Plotting, regression/audit checks, instructions |

Preserve the directory layout. From the project root:

```bash
tar --exclude='__pycache__' --exclude='*.pyc' \
  -czf symplectic-verifier-artifact.tar.gz \
  README.md PAPER_AUDIT.md symplectic_verifier.py parse.py benchmark.py \
  gen_scaling.py gen_respa_md.py gen_respa_ff.py plot_scaling.py tests examples
```

Extract into an empty directory with
`tar -xzf symplectic-verifier-artifact.tar.gz -C /path/to/directory`.
Papers and saved results are optional. Exclude caches, virtual environments,
`.sage-runtime`, editor settings, LaTeX intermediates, and the unrelated `test.py`.
For release, add a license, hardware/software versions, and fresh result logs.
The archive does not contain SageMath: provide a prepared VM/container for offline
use. With Conda, also save `conda env export > environment.yml`.

## Interpretation and validation

`VERIFIED` means Sage's simplifier reduced the symbolic residual to zero. This
is not a separately checked proof certificate, and a residual may reflect a
simplification limit. Smoothness/domain assumptions and implicit-solution
existence are not checked; implicit differentiation assumes `I - h*J_xi` is
invertible and the step size is state-independent.

Exact bend and Tao double pendulum use the standard symbolic verifier on their
`.spec` files. They may report `TIMEOUT` at the chosen budget; this is an accepted
benchmark outcome, not a verification success. There is no numerical fallback.

Validation: core tests and five plotting regressions pass; regenerated inputs
match all 80 included files. The paper audit's ten exact-arithmetic reference
maps pass; the two implicit-solve checks now pass and the indexed-loop check
still fails as documented above. An eight-case
paper-mode run (`N=2`, 120-second budget) produced seven passes and one coupled
GL-2 timeout, with both panels rendered. See `results/paper-smoke/` and
`results/paper-audit/`. These checks used SageMath 9.5/Python 3.10.12 and a separate
Matplotlib 3.10.9 environment on Ubuntu 22.04. The full dimension sweep and fresh
OS/Conda installations were not rerun during this audit.
