# Symplectic verifier artifact

This is artifact for the POPL paper A SymPL Program Logic for Geometric Numerical Integrators. 

## File Description

| Files | Purpose | Purpose |
| --- | --- | --- |
| `symplectic_verifier.py` | Forms, command representation, pullbacks, verification, built-in tests | Always |
| `parse.py` | Specification parser, CLI, and parser demos | Using `.spec` input |
| `examples/**/*.spec` and `examples/*.spec` | Example programs and benchmark inputs | Demonstrating or evaluating the verifier |
| `benchmark.py` | Timed subprocess runner with per-case timeouts | Reproducing the evaluation |
| `gen_scaling.py` | Regenerates the scaling specifications | Changing or regenerating the dimension sweep |
| `gen_respa_md.py`, `gen_respa_ff.py` | Generate the molecular-dynamics RESPA examples | Regenerating those examples |
| `plot_scaling.py` | Renders scaling CSV data with Matplotlib | Producing plots |
| Existing `.csv` and plot `.pdf` files | Saved evaluation outputs | Optional reference results; not inputs to verification |

## Environment and installation

Use a complete **SageMath installation**, including its symbolic algebra
dependencies. Matplotlib is needed only for plotting.

### Without Conda (Ubuntu 22.04 / Debian 12)

Install SageMath and plotting dependencies using the system package manager:

```bash
sudo apt-get update
sudo apt-get install -y sagemath python3-matplotlib python3-venv
```

SageMath packages are available for [Ubuntu 22.04](https://packages.ubuntu.com/jammy/sagemath)
and [Debian 12](https://packages.debian.org/bookworm/sagemath). On Ubuntu, if
`sagemath` is unavailable, run `sudo add-apt-repository universe`, then repeat
the commands above. Other OS releases may have different package availability.

Check the installation from the project root:

```bash
sage --version
sage -python -c 'from sage.all import SR; print(SR(1))'
python3 -c 'import matplotlib; print(matplotlib.__version__)'
sage -python parse.py examples/benchmarks/textbook/symplectic_euler_pendulum.spec
```

The last command should print `RESULT : VERIFIED`. No project build is needed;
run verifier commands with `sage -python`.

### With Conda

With Conda available, one installation route is:

```bash
conda create -n symplectic-artifact -c conda-forge sage matplotlib
conda activate symplectic-artifact
sage --version
sage -python -c 'from sage.all import SR; print(SR(1))'
```

If these do not work, see the [official SageMath installation instructions](https://doc.sagemath.org/html/en/installation/conda). 

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

## Reproducing the evaluation
All the programs from the paper can be found in /examples. 
To see all reproduce Table 1 of the paper run the following three commands. 

```bash
python3 benchmark.py --only conformal --timeout 120 
python3 benchmark.py --only textbook  --timeout 120
python3 benchmark.py --only bespoke  --timeout 120

```
In the bespoke benchmark, Exact Bend and Multipole Kick correspond to SciBMad Integrators. 

To regenerate inputs, run these from the project root. We do not recommend you do this, since it has been done already, and it might delete all examples:

```bash
python3 gen_scaling.py
sage -python gen_respa_md.py
sage -python gen_respa_ff.py
```

`gen_scaling.py` deletes existing `.spec` files in
`examples/benchmarks/scaling/` before writing its sweep. The RESPA generators
overwrite their respective example files. 

Create a results folder. Then, for scaling data and plotting:

```bash
python3 benchmark.py --only scaling --timeout 120 --plot results/scaling.pdf
```
This will reproduce Figure 10 from the paper.

## Writing a new Program
The following is how you would create a new integrator and its specification.  Params is non coordinate variables. 
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


Save this as `my_integrator.spec` and run
`sage -python parse.py my_integrator.spec`:


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

The project can be extended in the following ways:

1. **Add an integrator or Hamiltonian:** create a `.spec` under `examples/`, using
   existing textbook cases as templates. New combinations of assignments,
   solves, and bounded loops need no changes to the core. Check a known positive
   example and a deliberately perturbed negative example.
2. **Add a scalar function:** register a Sage symbolic callable in
   `parse.py:FUNCTIONS`. Confirm that differentiation, substitution, and
   simplification behave correctly for its domain.
3. **Add a command or form operation:** extend the command classes and
   `Verifier.pull` in `symplectic_verifier.py`, then the tokenizer/parser if it
   needs textual syntax. Define its mathematical pullback first. Add parser demos with
   explicit expected outcomes. 

Use exact Sage arithmetic in Python drivers. Unlike `.spec` literals, a plain
Python expression such as `1/3` becomes a float.  Use `QQ(1)/3`
or `SR(1)/3` when constructing rational constants.
