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
python3 benchmark.py --only textbook harmonic --timeout 30
python3 benchmark.py --list
python3 benchmark.py --timeout 120 > benchmark-results.txt 2>&1
```

The first two finish with `ALL TESTS PASSED` and
`ALL DEMOS BEHAVED AS EXPECTED`. Deliberately negative cases print `FAILED`.
The harmonic selection has six cases. The full registry has 108 cases and can
run for many minutes. Filters match family/display names, case-insensitively;
multiple filters use AND, e.g. `--only textbook gl-2 harmonic`.

For a single `.spec`, `parse.py` exits 0 on verification, 1 on a residual;
syntax/runtime errors also exit nonzero. The benchmark runner reports `PASS`,
`FAIL`, `TIMEOUT`, `ERROR`, or `SKIP`. **Timeouts and skips alone still produce
exit 0**: inspect its summary. Larger dimensions are skipped after the first
scaling timeout for each method/variant.

## Plot scaling results

```bash
# Small end-to-end run: both variants, all five methods, N=2.
python3 benchmark.py --only scaling 'N=2 (' --timeout 30 --plot results/smoke.pdf

# Full sweep; writes results/scaling.csv alongside the PDF.
python3 benchmark.py --only scaling --timeout 120 --plot results/scaling.pdf

# Replot an existing CSV (PDF or PNG).
python3 plot_scaling.py results/scaling.csv results/scaling.png 120
python3 -m unittest discover -s tests -v
```

The plotter accepts current `coupled`/`regular` variant names and the older
`Coupled Oscillators`/`Independent Oscillators` labels in the bundled CSV.
It plots all methods, with separate panels per variant. Circles/lines represent
passes, triangles mark timeouts, crosses mark positive-duration failures/errors,
and the dashed line shows the budget. Skipped rows are omitted. The CSV remains
available if rendering fails; consult it for full statuses and timings.

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
| `examples/` and `benchmark.py` | Full evaluation inputs, numerical checks, runner |
| `gen_scaling.py`, `gen_respa_md.py`, `gen_respa_ff.py` | Input regeneration |
| `plot_scaling.py`, `tests/`, `README.md` | Plotting, plotting regressions, instructions |

Preserve the directory layout. From the project root:

```bash
tar --exclude='__pycache__' --exclude='*.pyc' \
  -czf symplectic-verifier-artifact.tar.gz \
  README.md symplectic_verifier.py parse.py benchmark.py \
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

The exact-bend and Tao double-pendulum scripts print `CERTIFIED SYMPLECTIC` but
perform numerical sampling (200-bit arithmetic, tolerance `1e-40`), not symbolic
proofs. Report these separately from symbolic verification.

Validated on Ubuntu 22.04 with SageMath 9.5/Python 3.10.12: core tests, parser demos,
six harmonic benchmarks, and the exact-bend numerical check. Four plotting
regressions passed with Matplotlib 3.10.9 in a separate virtual environment.
A fresh 10-case scaling run produced PDF/PNG plots with both panels: nine passes
and one coupled GL-2 timeout at 30 seconds, displayed as a timeout marker.
The bundled legacy CSV also rendered correctly. The full benchmark suite and
fresh OS/Conda installations have not been tested.
