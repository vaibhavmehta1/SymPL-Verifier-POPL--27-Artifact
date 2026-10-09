# Symplectic verifier artifact

This project checks differential-form preservation for programs describing
integration steps. Given a precondition `A`, program `C`, and postcondition `B`,
it computes the pullback `C* B` and checks whether it equals `A`. Setting both
assertions to `dq ^ dp` checks preservation of the canonical symplectic form;
other forms express volume preservation or conformal symplecticity.

## Installation (start here)

The verifier requires **SageMath**, including its symbolic algebra dependencies.
There is no project compilation step or pip package to install. Matplotlib is
optional for verification but required for plots; the setup below includes it.
Use Linux for the evaluation artifact. macOS can also use this Conda setup;
on Windows, run it inside an Ubuntu **WSL2** terminal because the benchmark
runner requires POSIX process-group signals. See the
[Windows/WSL setup in the Sage installation guide](https://doc.sagemath.org/html/en/installation/index.html).

### 1. Prepare a terminal

Use Bash for the commands below. On macOS, enter `bash` in Terminal first.
You need an Internet connection to download the installer and dependencies.
On Ubuntu/Debian (including Ubuntu under WSL2), install the download and archive
utilities if they are missing:

```bash
sudo apt-get update
sudo apt-get install -y curl ca-certificates tar bzip2
```

On macOS, `curl` and `tar` are normally already available. Administrator access
is needed only for the system-package commands above; install Conda and the
project environment as your normal user.

### 2. Install Conda using Miniforge

If `conda --version` already works, skip the installer and proceed to step 3.
Otherwise, download the installer matching your OS and CPU architecture:

```bash
curl -fL -o Miniforge3.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3.sh
```

Accept the license and use the default installation directory
`$HOME/miniforge3`. Then enable Conda in the current Bash session:

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda --version
```

If you chose another installation directory, substitute it in the `source`
command. To enable `conda activate` in future Bash terminals, run `conda init bash`
and reopen the terminal. For other platforms or installer issues, see the
[official Miniforge instructions](https://github.com/conda-forge/miniforge#install).

### 3. Create and activate the SageMath environment

```bash
conda create -n symplectic-artifact --override-channels -c conda-forge sage matplotlib
conda activate symplectic-artifact
```

Confirm the package plan when prompted. This installs SageMath and its algebra
backends, Python, and Matplotlib together. This follows the
[official SageMath Conda installation route](https://doc.sagemath.org/html/en/installation/conda.html).
Do not substitute `pip install sage`: that is not a complete SageMath installation.

The command resolves available package versions; it is not a version lock.
The existing verifier checks were run with SageMath 9.5 and Python 3.10.12;
a newly resolved environment must pass the checks below before release.

### 4. Open or extract the artifact

If you already have this source directory, change into it:

```bash
cd /path/to/POPLArtifact
```

Replace `/path/to/POPLArtifact` with your actual checkout location. Alternatively,
if you received the archive produced by this README, run the following from the
directory containing the archive:

```bash
mkdir -p symplectic-verifier-artifact
tar -xzf symplectic-verifier-artifact.tar.gz -C symplectic-verifier-artifact
cd symplectic-verifier-artifact
```

Run all subsequent project commands from this directory. It should contain
`README.md`, `symplectic_verifier.py`, `parse.py`, and `examples/`.

### 5. Check the installation and run the verifier

With the `symplectic-artifact` environment active:

```bash
sage --version
sage -python -c 'from sage.all import SR; print(SR(1))'
python3 -c 'import matplotlib; print(matplotlib.__version__)'
sage -python parse.py examples/benchmarks/textbook/symplectic_euler_pendulum.spec
```

The import checks should print `1` and a Matplotlib version, respectively.
The specification should finish with `RESULT : VERIFIED` and exit 0. Then run
both built-in regression suites:

```bash
sage -python symplectic_verifier.py
sage -python parse.py
```

Expect `ALL TESTS PASSED` and `ALL DEMOS BEHAVED AS EXPECTED`. Some individual
negative examples intentionally print `FAILED`; the final suite messages are
the success criteria. To check the benchmark runner as well:

```bash
python3 benchmark.py --only textbook harmonic --timeout 30
```

This selects six cases. Inspect the summary for six passes; a slow machine may
need a larger timeout. A full evaluation is described below.

### 6. Resume work or preserve the installed environment

In each new terminal, activate the environment and enter the source directory
before running commands:

```bash
conda activate symplectic-artifact
cd /path/to/POPLArtifact
```

After validating your installation, capture its dependencies from the project
root for the artifact release:

```bash
conda env export -n symplectic-artifact > environment.yml
conda list --explicit > conda-explicit.txt
```

Include those generated files in the release. An evaluator can restore the
exported environment using `conda env create -f environment.yml`, or use
`conda create -n symplectic-artifact --file conda-explicit.txt` on the same
platform. These are alternative installation routes to step 3 and require that
the environment name is not already in use. The explicit package list is
platform-specific. For offline evaluation, supply a prepared VM or container;
a source archive and dependency lists alone do not contain the dependencies.

### Installation troubleshooting

- **`conda` is missing or activation fails:** source the installed
  `etc/profile.d/conda.sh` as in step 2, then activate the environment again.
- **`sage` is missing or `sage.all` cannot be imported:** activate the environment
  and use `sage -python` for verifier scripts. Ordinary system Python may not
  contain SageMath. If needed, the benchmark runner accepts
  `SAGE=/absolute/path/to/sage` in its environment.
- **Sage cannot write its user directory:** from the project root, run
  `export DOT_SAGE="$PWD/.sage-runtime"` before starting Sage. This selects a
  writable runtime directory; exclude it from the packaged artifact.
- **Matplotlib reports a NumPy ABI/import error:** check that the Conda
  environment is active and repeat the Matplotlib import check above. Avoid
  mixing system Python packages into this environment. Set
  `export PLOT_PYTHON="$(command -v python3)"` after the import check succeeds to
  make the benchmark plotter use that interpreter. Plotting is optional for
  verifier runs; the known plot-label limitation is documented below.

## What to package

The **minimal verifier library is `symplectic_verifier.py` plus SageMath**.
For the command-line interface, include **`parse.py` as well**, alongside the
specifications you want to check. There is no compilation step for this project.

| Files | Purpose | Needed when |
| --- | --- | --- |
| `symplectic_verifier.py` | Forms, command representation, pullbacks, verification, built-in tests | Always |
| `parse.py` | Specification parser, CLI, and parser demos | Using `.spec` input |
| `examples/**/*.spec` and `examples/*.spec` | Example programs and benchmark inputs | Demonstrating or evaluating the verifier |
| `benchmark.py` | Timed subprocess runner with per-case timeouts | Reproducing the evaluation |
| `examples/exact_bend_certificate.py` | Numerical sampling check of `examples/exact_bend.spec` | Running the registered exact-bend benchmark |
| `examples/tao_double_pendulum_certificate.py` | Numerical sampling check of `examples/benchmarks/textbook/tao_double_pendulum.spec` | Running the registered Tao double-pendulum benchmark |
| `gen_scaling.py` | Regenerates the scaling specifications | Changing or regenerating the dimension sweep |
| `gen_respa_md.py`, `gen_respa_ff.py` | Generate the molecular-dynamics RESPA examples | Regenerating those examples |
| `plot_scaling.py` | Renders scaling CSV data with Matplotlib | Producing plots |
| `README.md` | Setup, usage, interpretation, and extension instructions | Distributing the artifact |
| `main.tex`, `main.pdf`, `boris.tex`, `boris.pdf` | Accompanying mathematical exposition | Optional documentation; not runtime dependencies |
| Existing `.csv` and plot `.pdf` files | Saved evaluation outputs | Optional reference results; not inputs to verification |

For a research artifact, include the core, parser, **all of `examples/`**, the
runner, generators, plotter, and this README. Keeping the directory layout is
important: the runner and numerical scripts resolve inputs relative to it.

Exclude `__pycache__/`, `*.pyc`, editor/agent configuration directories, and
LaTeX build intermediates (`*.aux`, `*.log`, `*.fls`, `*.fdb_latexmk`, `*.out`,
`*.synctex.gz`). `test.py` is an unrelated Sage vector-space experiment;
the actual verifier tests are inside `symplectic_verifier.py`.

From the project root, this command creates a source artifact containing the
evaluation tools and existing reference results:

```bash
tar --exclude='__pycache__' --exclude='*.pyc' \
    -czf symplectic-verifier-artifact.tar.gz \
    README.md symplectic_verifier.py parse.py benchmark.py \
    gen_scaling.py gen_respa_md.py gen_respa_ff.py plot_scaling.py examples
```

This archive does not bundle SageMath. For an offline evaluation, also provide
a prepared VM or container with SageMath and the optional plotting dependencies.
Before release, choose an appropriate license (none is supplied here), record
the exact environment and hardware used, and include fresh evaluation logs.
Existing saved timings do not establish timings on a new machine.

## Quick start and expected results

```bash
# Core algebra and verifier regression tests, including negative examples.
sage -python symplectic_verifier.py

# Parser demos, also including deliberately rejected examples.
sage -python parse.py

# Verify one specification.
sage -python parse.py examples/benchmarks/textbook/symplectic_euler_pendulum.spec

# Small benchmark selection with a per-case time limit.
python3 benchmark.py --only textbook harmonic --timeout 30
```

The first two commands finish with `ALL TESTS PASSED` and
`ALL DEMOS BEHAVED AS EXPECTED`, respectively. A printed `FAILED` inside a
negative test is expected: forward Euler is deliberately rejected, and Boris
preserves the tested volume form but not the tested canonical two-form.
The single pendulum specification should print `RESULT : VERIFIED` and exit 0.
The harmonic selection contains six cases expected to pass with sufficient time.

For a single specification, `parse.py` exits 0 when verified and 1 when the
computed residual is nonzero. Invalid input or runtime errors also produce a
nonzero exit, with a traceback. The CLI accepts one input filename; it has no
`--help` flag. A reported residual is indexed by a wedge basis and gives the
corresponding coefficient of `C* B - A`.

## Writing a specification

Save this as `my_integrator.spec` and run
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

This is symplectic Euler for the pendulum. Separate assignments execute in
order. In contrast, `(q, p) := (q + h*p, p - h*sin(q))` evaluates both right-hand
sides in the old state; substituting that simultaneous assignment gives the
non-symplectic forward Euler method.

The specification language supports:

- `coords`: state variables, including any scratch or stage variables. Declare
  coordinates and `params` before `pre`, `post`, or `program`. Parameters are
  constant with respect to state differentiation.
- Assertions: `dq`, `dp`, `d(expression)`, addition/subtraction, scalar
  multiplication, and `^` for wedge product. For multiple canonical pairs use,
  for example, `dq1 ^ dp1 + dq2 ^ dp2`.
- Scalar expressions: `+ - * /`, `**` for powers, parentheses, integer/decimal
  literals, and functions registered in `parse.py:FUNCTIONS`. These include
  `sin`, `cos`, `tan`, `exp`, `log`, `sqrt`, hyperbolic functions, inverse trig
  functions, `atan2`/`arctan2`, and `sinc`/`sincu`. Decimal literals become exact
  rationals; scientific notation is not supported. Write multiplication
  explicitly, and use `**`, never `^`, for scalar powers.
- Commands: `skip`, scalar or simultaneous tuple assignment, implicit `solve`,
  sequencing, `while^N { ... }`, and `for^T { ... }`.
- Layout: `#` comments and newlines or semicolons between statements. Put a
  newline or semicolon before the closing brace after an assignment. Parentheses
  allow multiline expressions and tuples.

Loop bounds in `.spec` files are nonnegative integer literals. Both loop forms
are implemented by repeated pullback; `for` does not automatically bind or
increment an index. For an assertion family indexed by a symbolic parameter,
use `Verifier.verify_for` through the Python API (see `test_conformal_for`).

An implicit midpoint example uses stage variables as follows:

```text
coords q p Q P
params h
pre  dq ^ dp
post dq ^ dp
program {
    (Q, P) := (q, p)
    (Q, P) := solve (a, b) -> (b/2, -sin(a)/2)
    q := q + h*P
    p := p - h*sin(Q)
}
```

Here `solve (a, b) -> (G1, G2)` denotes the fixed-point equations
`(Q_new, P_new) = (Q_old, P_old) + h*(G1, G2)`, with `a, b` representing the new
values inside `G`. The parser requires a parameter named `h`. Choose distinct
unknown names that do not collide with coordinates or parameters. The verifier
uses implicit differentiation; it does not numerically solve these equations.

## Using and extending the Python API

Save Python drivers in the project root (or add that directory to `PYTHONPATH`)
and execute them with `sage -python driver.py`:

```python
from sage.all import sin
from symplectic_verifier import Context, Assign, Seq, Verifier

ctx = Context(["q", "p"], ["h"])
q, p = ctx.coords
h = ctx.param("h")
omega = ctx.omega0([("q", "p")])
program = Seq(Assign(["p"], [p - h*sin(q)]),
              Assign(["q"], [q + h*p]))
ok, pullback, residual = Verifier(ctx).verify(omega, program, omega)
assert ok, residual
```

`parse.parse_spec(text)` returns `(ctx, pre, program, post)`;
`parse.verify_spec(text, verbose=False)` returns `(ok, residual)`.
`Verifier.verify(pre, program, post, verbose=False)` additionally returns the
computed pullback. Use `Verifier.pull(program, form)` when you need the pullback
without the final equality check.

To extend the project:

1. **Add an integrator or Hamiltonian:** create a `.spec` under `examples/`, using
   existing textbook cases as templates. New combinations of assignments,
   solves, and bounded loops need no changes to the core. Check a known positive
   example and a deliberately perturbed negative example.
2. **Add a benchmark:** edit `benchmark.py:registry()` and use `_tb`, `_cf`, or
   `_spec` to reference the file. Non-scaling specifications are not automatically
   discovered. The runner currently treats any child exit 0 as `PASS`; the
   `Bench.expect` label does not implement expected-failure assertions. Put
   negative regressions in a test driver that asserts rejection and exits 0.
3. **Extend the dimension sweep:** edit `gen_scaling.py` (`_LADDER`, `DIMS`,
   `BUILDERS`, and `coords_line` as appropriate). Scaling filenames follow
   `<variant>_<method>_chain<NN>.spec`, where variant is `coupled` or `regular`.
   Update the runner's method labels/order when adding a method. Reported
   physical dimension is `D = 2*N`; implicit methods also have stage coordinates.
4. **Add a scalar function:** register a Sage symbolic callable in
   `parse.py:FUNCTIONS`. Confirm that differentiation, substitution, and
   simplification behave correctly for its domain. Add a parser regression
   using a known identity and a case that should fail.
5. **Add a command or form operation:** extend the command classes and
   `Verifier.pull` in `symplectic_verifier.py`, then the tokenizer/parser if it
   needs textual syntax. Define its mathematical pullback first. Add regression
   checks to the core suite and call them from `main()`; add parser demos with
   explicit expected outcomes. Preserve simultaneous assignment semantics and
   the reverse traversal used for sequential pullback.

Use exact Sage arithmetic in Python drivers. Unlike `.spec` literals, a plain
Python expression such as `1/3` becomes a float before Sage sees it; use `QQ(1)/3`
or `SR(1)/3` when constructing rational constants.

## Reproducing the evaluation

The current registry has 108 cases: 20 textbook, 80 scaling, 2 conformal, and
6 bespoke cases. Some other specifications in `examples/` are standalone
examples and are not registered.

```bash
python3 benchmark.py --list
python3 benchmark.py --only conformal --timeout 120 --verbose
python3 benchmark.py --only textbook gl-2 harmonic --timeout 120
python3 benchmark.py --timeout 120 > benchmark-results.txt 2>&1
```

Filters are case-insensitive substrings of the family and display name; multiple
filters are combined with AND. For example, use `gl-2` to match the displayed
Gauss–Legendre name, rather than the filename prefix `gl2`.
Set `SAGE=/absolute/path/to/sage` if the executable is not on `PATH`.

Each case runs in a fresh process; timings include startup, parsing, verification,
and printed output. The default budget is 120 seconds per case. A timeout kills
the process group, including algebra subprocesses. Scaling runs skip larger
dimensions for a method/variant after its first timeout. The suite can take many
minutes; no full-suite runtime or memory bound is established here.

Statuses are `PASS`, `FAIL`, `TIMEOUT`, `ERROR`, and `SKIP`. The runner exits 1
for `FAIL` or `ERROR`, but **timeouts and skips alone still yield exit 0**. Read
the result table, not just the shell exit status. The implicit midpoint
double-pendulum case is included as a difficult case expected to time out on
the original evaluation setup; it is not converted into a passing verdict.

To regenerate inputs, run these from the project root:

```bash
python3 gen_scaling.py
sage -python gen_respa_md.py
sage -python gen_respa_ff.py
```

`gen_scaling.py` deletes existing `.spec` files in
`examples/benchmarks/scaling/` before writing its sweep. The RESPA generators
overwrite their respective example files; `gen_respa_ff.py` performs verification
before writing and may be expensive. Regeneration is unnecessary for the
included inputs.

For scaling data and plotting:

```bash
python3 benchmark.py --only scaling --timeout 120 --plot results/scaling.pdf
python3 plot_scaling.py results/scaling.csv results/scaling.pdf 120
```

The runner writes a sibling CSV and attempts plotting with a Python interpreter
that can import Matplotlib. Set `PLOT_PYTHON=/path/to/python` to select one.
If plotting is unavailable, the CSV is retained. Use a separate compatible
Python/Matplotlib environment if the system installation has a NumPy ABI error.

**Current plot limitation:** the runner emits variant names `coupled` and
`regular`, while `plot_scaling.py:PANELS` currently expects `Coupled Oscillators`
and `Independent Oscillators`, matching the bundled reference CSV. With newly
generated CSV data, it falls back to one panel and can omit a variant. To obtain
both panels from fresh data, change `PANELS` to
`[("coupled", "Coupled Oscillators"), ("regular", "Independent Oscillators")]`.
The plotter also intentionally omits Störmer–Verlet and plots recorded times
without separating timeout/failure statuses; consult the CSV and result table.

## Scope of the results

`VERIFIED` means the implemented symbolic simplifier reduced every coefficient
of the pullback residual to zero. This relies on SageMath's symbolic algebra;
it is not a proof certificate checked by a separate proof assistant. A residual
that does not simplify to zero is a failure to verify and may need mathematical
inspection; symbolic simplification is not complete for arbitrary expressions.

The implementation does not discharge smoothness, branch, nonzero-denominator,
or implicit-solution existence/uniqueness obligations. Implicit differentiation
assumes `I - h*J_xi` is invertible, with state-independent step size. Interpret
results on domains where the expressions and those assumptions hold. In
particular, the parser implements `sincu(t)` as `sin(t)/t` without implementing
its removable value at zero.

The exact-bend and Tao double-pendulum scripts print `CERTIFIED SYMPLECTIC`, but
they perform **numerical sampling, not symbolic proofs**. They use 200-bit
arithmetic and a `1e-40` residual threshold, with fixed seeds: exact bend uses
seed 7 and 8 samples per residual coefficient; Tao uses seed 11 and 6 phase-space
samples. A passing sample check provides numerical evidence only. Keep these
results separate from symbolic `VERIFIED` results when reporting the artifact.

## Validation of these instructions

On October 9, 2026, the following checks passed in the available Linux environment
using SageMath 9.5 and Python 3.10.12, with `DOT_SAGE` pointing to a writable
temporary directory:

- The core self-test suite and all parser demos.
- All six textbook harmonic benchmarks with a 30-second per-case budget.
- The exact-bend numerical sampling check with a 30-second budget.
- Both specification examples and the Python API example in this README.
- The packaging command's file coverage for all 108 registered benchmarks.

The full benchmark suite and a fresh Conda installation were not run as part of
this documentation check. Plot rendering was not validated: the available system
Matplotlib failed to import because of a NumPy ABI mismatch. The variant-name
mismatch described above was identified by inspecting the runner and plotter.
