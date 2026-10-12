"""
symplectic_verifier.py
======================

An industrial-strength verifier for symplectic integrators, implementing the
program logic of `main.tex` ("A Program Logic for Symplectic Integration").

The exact Hoare judgment

        |= {A} C {B}      iff      [[C]]^* [[B]] = [[A]]

is decided by computing the pullback of the postcondition `B` backwards through
the program `C` (substitution = pullback, Lemma "Substitution equals pullback")
and checking, term by term, that the result equals the precondition `A`.

Assertions are differential forms on the state space R^D with coordinates
x^1,...,x^D.  A k-form is stored as a dictionary

        { (i_1,...,i_k) : coefficient }          (strictly increasing indices)

where each coefficient is a SageMath symbolic expression that may depend on the
coordinates and on external parameters (step size h, loop index t, ...).

Supported commands (Section 2 of the paper):

    Skip()                              skip
    Assign(targets, exprs)              simultaneous explicit  x^r := e
    Solve(targets, xi, G, h)           simultaneous implicit  x^r := solve(x^r, G)
    Seq(c1, c2, ...)                   sequencing  C1; C2; ...
    While(N, body)                     while^N  (iterate body N times)
    For(T, body)                       for^T    (with explicit loop index)

Run with the SageMath interpreter:

    sage -python symplectic_verifier.py          # runs the self-test suite
"""

from itertools import count

from sage.all import SR, var, diff, matrix, identity_matrix, vector


_solution_symbol_ids = count()


def _fresh_solution_symbol():
    """Allocate a solution placeholder with a previously unused Sage name.

    Freshness is per pullback occurrence, not per target, Solve object, or
    Verifier. Also avoid names already registered by users: Sage and its CAS
    interfaces identify named symbols globally. Keep the registration while
    returned forms may still contain the symbol.
    """
    while True:
        name = "_sv_phi_{}".format(next(_solution_symbol_ids))
        if name not in SR.symbols:
            return SR.var(name)


def _simplify(coeff):
    """Simplify a scalar coefficient toward a canonical (ideally zero) form.

    Deliberately avoids ``Expression.simplify_full``: its ``simplify_rectform``
    pass rewrites every inverse-trig atom (arcsin / arctan / atan2 / ...) into
    logarithmic ``rectform``.  That pass both fails outright on ``atan2`` and
    blows the Maxima/ECL memory limit on the nested inverse-trig produced by
    pulling back through maps such as the exact sector bend.  This rectform-free
    chain still crushes the trig/radical identities the pullbacks generate
    (bare inverse-trig angles cancel as symbols; sines/cosines of them reduce
    via the trig identities), while keeping the expressions small."""
    try:
        return coeff.simplify_trig().simplify_rational().canonicalize_radical()
    except Exception:
        return coeff.simplify_full()


# ---------------------------------------------------------------------------
#  Multi-index / wedge bookkeeping
# ---------------------------------------------------------------------------

def _merge_indices(I, J):
    """Wedge two increasing multi-indices.

    Returns (sign, merged_tuple); sign is 0 when the indices share a coordinate
    (so the wedge vanishes), otherwise +/-1 with the merged tuple sorted
    increasingly.
    """
    if set(I) & set(J):
        return 0, None
    arr = list(I) + list(J)
    sign = 1
    # insertion sort, tracking the parity of the permutation
    for a in range(1, len(arr)):
        key = arr[a]
        b = a - 1
        while b >= 0 and arr[b] > key:
            arr[b + 1] = arr[b]
            b -= 1
            sign = -sign
        arr[b + 1] = key
    return sign, tuple(arr)


# ---------------------------------------------------------------------------
#  Context: coordinates + parameters
# ---------------------------------------------------------------------------

class Context:
    """The state space R^D with named coordinates and named parameters."""

    def __init__(self, coord_names, param_names=()):
        self.coord_names = list(coord_names)
        self.D = len(coord_names)
        self.coords = [var(n) for n in coord_names]
        self.index = {n: i for i, n in enumerate(coord_names)}
        self.params = {n: var(n) for n in param_names}

    def idx(self, name):
        """Index of a coordinate, accepting either its name or an int."""
        if isinstance(name, str):
            return self.index[name]
        return int(name)

    def coord(self, name):
        """The coordinate symbol (a scalar 0-assertion atom)."""
        return self.coords[self.idx(name)]

    def param(self, name):
        return self.params[name]

    # ---- assertion constructors ----------------------------------------

    def scalar(self, expr):
        """A 0-form (scalar atom) f."""
        return Form(self, {(): SR(expr)})

    def dx(self, name):
        """The coordinate 1-form dx^i = d(x^i)."""
        return Form(self, {(self.idx(name),): SR(1)})

    def zero(self):
        return Form(self, {})

    def omega0(self, pairs):
        """Canonical symplectic form sum_i dq_i ^ dp_i for given (q_i, p_i)."""
        w = self.zero()
        for q, p in pairs:
            w = w + self.dx(q).wedge(self.dx(p))
        return w


# ---------------------------------------------------------------------------
#  Differential forms (assertions)
# ---------------------------------------------------------------------------

class Form:
    """A differential form / assertion over a Context.

    Internally a dict from increasing multi-indices to symbolic coefficients.
    """

    def __init__(self, ctx, terms=None):
        self.ctx = ctx
        self.terms = {}
        if terms:
            for key, coeff in terms.items():
                self._add(key, coeff)

    # ---- internal helpers ----------------------------------------------

    def _add(self, key, coeff):
        key = tuple(key)
        coeff = SR(coeff)
        new = self.terms.get(key, SR(0)) + coeff
        if new.is_trivial_zero():
            self.terms.pop(key, None)
        else:
            self.terms[key] = new

    def copy(self):
        f = Form(self.ctx)
        f.terms = dict(self.terms)
        return f

    # ---- algebra -------------------------------------------------------

    def __add__(self, other):
        res = self.copy()
        for key, coeff in other.terms.items():
            res._add(key, coeff)
        return res

    def __sub__(self, other):
        res = self.copy()
        for key, coeff in other.terms.items():
            res._add(key, -coeff)
        return res

    def scale(self, c):
        c = SR(c)
        return Form(self.ctx, {key: c * coeff for key, coeff in self.terms.items()})

    def wedge(self, other):
        res = Form(self.ctx)
        for k1, c1 in self.terms.items():
            for k2, c2 in other.terms.items():
                sign, merged = _merge_indices(k1, k2)
                if sign == 0:
                    continue
                res._add(merged, sign * c1 * c2)
        return res

    def d(self):
        """Exterior derivative  d(sum f_I dx^I) = sum (df_I/dx^k) dx^k ^ dx^I."""
        res = Form(self.ctx)
        for key, coeff in self.terms.items():
            for k in range(self.ctx.D):
                dc = diff(coeff, self.ctx.coords[k])
                if dc.is_trivial_zero():
                    continue
                sign, merged = _merge_indices((k,), key)
                if sign == 0:
                    continue
                res._add(merged, sign * dc)
        return res

    # ---- pullback ------------------------------------------------------

    def pullback(self, scalar_subs, dimages):
        """Pullback under a smooth map phi: R^D -> R^D.

        `scalar_subs` is the substitution {x^j -> phi_j} applied to scalar
        coefficients; `dimages[j]` is the 1-form image of dx^j, namely d(phi_j).
        This is exactly  phi^*(f dx^I) = (f o phi) (phi^*dx^{i_1}) ^ ... .
        """
        res = Form(self.ctx)
        for key, coeff in self.terms.items():
            new_coeff = coeff.subs(scalar_subs) if scalar_subs else coeff
            w = Form(self.ctx, {(): new_coeff})
            for j in key:
                w = w.wedge(dimages[j])
            res = res + w
        return res

    # ---- simplification / comparison -----------------------------------

    def simplified(self):
        res = Form(self.ctx)
        for key, coeff in self.terms.items():
            res._add(key, _simplify(coeff))
        return res

    def residual(self, other):
        """Per-term residual (self - other), fully simplified; drops zeros."""
        diff_form = self - other
        out = {}
        for key, coeff in diff_form.terms.items():
            s = _simplify(coeff)
            if not s.is_zero():
                out[key] = s
        return out

    def equals(self, other):
        return len(self.residual(other)) == 0

    # ---- display -------------------------------------------------------

    def __repr__(self):
        if not self.terms:
            return "0"
        parts = []
        for key in sorted(self.terms, key=lambda t: (len(t), t)):
            coeff = _simplify(self.terms[key])
            if not key:
                parts.append("({})".format(coeff))
            else:
                basis = "^".join("d{}".format(self.ctx.coord_names[i]) for i in key)
                parts.append("({})*{}".format(coeff, basis))
        return " + ".join(parts)


# ---------------------------------------------------------------------------
#  Commands
# ---------------------------------------------------------------------------

class Command:
    pass


class Skip(Command):
    def __repr__(self):
        return "skip"


class Assign(Command):
    """Simultaneous explicit assignment  x^targets := exprs."""

    def __init__(self, targets, exprs):
        self.targets = list(targets)          # coordinate names or indices
        self.exprs = [SR(e) for e in exprs]
        assert len(self.targets) == len(self.exprs)

    def __repr__(self):
        return "({}) := ({})".format(
            ", ".join(map(str, self.targets)),
            ", ".join(map(str, self.exprs)))


class Solve(Command):
    """Simultaneous implicit assignment  x^targets := solve(x^targets, G).

    The fixed point xi = x^targets + h*G(xi, x, pi) is solved abstractly:
    we introduce symbols phi_i for the solution components and obtain the
    Jacobian D_x(phi) by implicit differentiation (Theorem "IFT").

    Parameters
    ----------
    targets : coordinate names/indices being solved for (length m)
    xi      : list of m sage symbols denoting the unknowns of G
    G       : list of m sage expressions G_i(xi, x, pi)
    h       : the step-size symbol h(pi)
    phi     : must be None. Caller-supplied solution symbols are rejected;
              the verifier allocates fresh placeholders for each occurrence.
    """

    def __init__(self, targets, xi, G, h, phi=None):
        if phi is not None:
            raise ValueError("phi must be None; implicit solution symbols are "
                             "allocated freshly by the verifier")
        self.targets = list(targets)
        self.xi = list(xi)
        self.G = [SR(g) for g in G]
        self.h = SR(h)
        m = len(self.targets)
        assert len(self.xi) == m and len(self.G) == m

    def __repr__(self):
        return "({}) := solve(...)".format(", ".join(map(str, self.targets)))


class Seq(Command):
    def __init__(self, *cmds):
        flat = []
        for c in cmds:
            if isinstance(c, Seq):
                flat.extend(c.cmds)
            else:
                flat.append(c)
        self.cmds = flat

    def __repr__(self):
        return "; ".join(map(repr, self.cmds))


class While(Command):
    """while^N C : iterate the body exactly N times."""

    def __init__(self, N, body):
        self.N = int(N)
        self.body = body

    def __repr__(self):
        return "while^{} ({})".format(self.N, self.body)


class For(Command):
    """for^T C : iterate exactly T times, carrying a loop index symbol."""

    def __init__(self, T, body):
        self.T = int(T)
        self.body = body

    def __repr__(self):
        return "for^{} ({})".format(self.T, self.body)


# ---------------------------------------------------------------------------
#  The verifier
# ---------------------------------------------------------------------------

class Verifier:
    def __init__(self, ctx):
        self.ctx = ctx

    # ---- pullback through a single (non-loop) command ------------------

    def _identity_dimages(self):
        return [self.ctx.dx(j) for j in range(self.ctx.D)]

    def pull_assign(self, cmd, form):
        ctx = self.ctx
        tgt = [ctx.idx(t) for t in cmd.targets]
        scalar_subs = {ctx.coords[r]: e for r, e in zip(tgt, cmd.exprs)}
        dimages = self._identity_dimages()
        for r, e in zip(tgt, cmd.exprs):
            # phi^* dx^r = d(e) = sum_k (de/dx^k) dx^k
            dimages[r] = ctx.scalar(e).d()
        return form.pullback(scalar_subs, dimages)

    def pull_solve(self, cmd, form):
        ctx = self.ctx
        m = len(cmd.targets)
        tgt = [ctx.idx(t) for t in cmd.targets]

        # If G does not depend on the unknown, the fixed-point equation is an
        # explicit simultaneous assignment. Preserve its exact values instead
        # of introducing unconstrained placeholders (in particular G=0 is skip).
        unknowns = set(cmd.xi)
        if all(not unknowns.intersection(g.variables()) for g in cmd.G):
            return self.pull_assign(Assign(cmd.targets, [
                ctx.coords[r] + cmd.h*g for r, g in zip(tgt, cmd.G)
            ]), form)

        # Each dynamic occurrence denotes a different solution value, even
        # when the same Solve object is reused or unrolled in a loop. These
        # independent placeholders give a conservative equality check: an
        # identity must hold for all their values, without assuming relations
        # between solutions at different states. They are not state functions
        # for subsequent differentiation of an unsimplified returned Form.
        phi = [_fresh_solution_symbol() for _ in tgt]

        # Jacobians of G, evaluated at the fixed point xi = phi.
        at_phi = {cmd.xi[l]: phi[l] for l in range(m)}

        # J_xi[i][l] = dG_i/dxi_l ,  J_x[i][k] = dG_i/dx^k   (then xi -> phi)
        J_xi = matrix(SR, m, m,
                      lambda i, l: diff(cmd.G[i], cmd.xi[l]).subs(at_phi))
        J_x = matrix(SR, m, ctx.D,
                     lambda i, k: diff(cmd.G[i], ctx.coords[k]).subs(at_phi))

        # selector E_r (m x D): row i has a 1 in column tgt[i]
        E = matrix(SR, m, ctx.D, lambda i, k: 1 if k == tgt[i] else 0)

        # implicit differentiation:  (I - h J_xi) D_x phi = E + h J_x
        A = identity_matrix(SR, m) - cmd.h * J_xi
        Dphi = A.inverse() * (E + cmd.h * J_x)    # m x D Jacobian of phi

        # scalar substitution x^r -> phi_r
        scalar_subs = {ctx.coords[tgt[i]]: phi[i] for i in range(m)}

        # differential images: dx^r -> d(phi_r) = sum_k Dphi[i,k] dx^k
        dimages = self._identity_dimages()
        for i in range(m):
            row = ctx.zero()
            for k in range(ctx.D):
                c = Dphi[i, k]
                if not c.is_trivial_zero():
                    row = row + ctx.dx(k).scale(c)
            dimages[tgt[i]] = row

        return form.pullback(scalar_subs, dimages)

    def pull(self, cmd, form):
        """Pullback `form` backwards through `cmd` :  [[cmd]]^* form."""
        if isinstance(cmd, Skip):
            return form
        if isinstance(cmd, Assign):
            return self.pull_assign(cmd, form)
        if isinstance(cmd, Solve):
            return self.pull_solve(cmd, form)
        if isinstance(cmd, Seq):
            # [[C1;C2]]^* = [[C1]]^* o [[C2]]^*  -> process right to left
            f = form
            for c in reversed(cmd.cmds):
                f = self.pull(c, f)
            return f
        if isinstance(cmd, While):
            f = form
            for _ in range(cmd.N):
                f = self.pull(cmd.body, f)
            return f
        if isinstance(cmd, For):
            f = form
            for _ in range(cmd.T):
                f = self.pull(cmd.body, f)
            return f
        raise TypeError("unknown command: {!r}".format(cmd))

    # ---- top-level judgments -------------------------------------------

    def verify(self, pre, cmd, post, verbose=True):
        """Check the exact Hoare triple  |= {pre} cmd {post}.

        Returns (ok, computed_pullback, residual_dict).
        """
        computed = self.pull(cmd, post).simplified()
        residual = computed.residual(pre)
        ok = len(residual) == 0
        if verbose:
            print("  program   : {}".format(cmd))
            print("  post  B   : {}".format(post.simplified()))
            print("  [[C]]^* B : {}".format(computed))
            print("  pre   A   : {}".format(pre.simplified()))
            if ok:
                print("  RESULT    : VERIFIED  ([[C]]^* B == A)")
            else:
                print("  RESULT    : FAILED")
                for key, coeff in residual.items():
                    basis = "^".join("d{}".format(self.ctx.coord_names[i]) for i in key) or "1"
                    print("     residual {} : {}".format(basis, coeff))
        return ok, computed, residual

    def verify_for(self, A_of_t, body, T, t, verbose=True):
        """Check the parametric for-loop premise once, in the loop index `t`.

        Given an assertion family A_of_t(t) and body C, check
            |= {A_t} C {A_{t+1}}   for the symbolic index t,
        which (Theorem "Exact for-loop") yields  |= {A_0} for^T C {A_T}.
        Returns (ok, residual_dict).
        """
        post = A_of_t(t + 1)
        pre = A_of_t(t)
        computed = self.pull(body, post).simplified()
        residual = computed.residual(pre)
        ok = len(residual) == 0
        if verbose:
            print("  for-loop premise  |= {A_t} body {A_{t+1}}  (symbolic t)")
            print("  A_t       : {}".format(pre))
            print("  A_(t+1)   : {}".format(post))
            print("  pullback  : {}".format(computed))
            if ok:
                print("  RESULT    : VERIFIED for all t  =>  |= {{A_0}} for^{} C {{A_{}}}".format(T, T))
            else:
                print("  RESULT    : FAILED")
                for key, coeff in residual.items():
                    print("     residual {} : {}".format(key, coeff))
        return ok, residual


# ===========================================================================
#  Self-test suite
# ===========================================================================

def _hr(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def test_form_algebra():
    _hr("TEST: form algebra (antisymmetry, d.d = 0)")
    ctx = Context(["x", "y", "z"])
    dx, dy, dz = ctx.dx("x"), ctx.dx("y"), ctx.dx("z")

    assert (dx.wedge(dy) + dy.wedge(dx)).equals(ctx.zero()), "wedge not alternating"
    assert dx.wedge(dx).equals(ctx.zero()), "dx^dx != 0"
    print("  alternating wedge: OK")

    x, y, z = ctx.coords
    f = ctx.scalar(x ** 2 * y + sin_z(z))
    assert f.d().d().equals(ctx.zero()), "d.d != 0 on 0-form"
    alpha = ctx.dx("x").scale(x * y) + ctx.dx("z").scale(z ** 3)
    assert alpha.d().d().equals(ctx.zero()), "d.d != 0 on 1-form"
    print("  d o d = 0       : OK")


def sin_z(z):
    from sage.all import sin
    return sin(z)


def test_symplectic_euler():
    _hr("TEST: symplectic Euler on the pendulum  (should VERIFY)")
    ctx = Context(["q", "p"], ["h"])
    q, p, h = ctx.coord("q"), ctx.coord("p"), ctx.param("h")
    from sage.all import sin

    # p := p - h sin q ;  q := q + h p
    prog = Seq(Assign(["p"], [p - h * sin(q)]),
               Assign(["q"], [q + h * p]))
    w0 = ctx.omega0([("q", "p")])
    ok, _, _ = Verifier(ctx).verify(w0, prog, w0)
    assert ok, "symplectic Euler should be symplectic"


def test_explicit_euler():
    _hr("TEST: explicit (forward) Euler on the pendulum  (should FAIL)")
    ctx = Context(["q", "p"], ["h"])
    q, p, h = ctx.coord("q"), ctx.coord("p"), ctx.param("h")
    from sage.all import sin

    # simultaneous  (q, p) := (q + h p, p - h sin q)
    prog = Assign(["q", "p"], [q + h * p, p - h * sin(q)])
    w0 = ctx.omega0([("q", "p")])
    ok, _, _ = Verifier(ctx).verify(w0, prog, w0)
    assert not ok, "explicit Euler is NOT symplectic; verifier must reject"


def test_implicit_midpoint():
    _hr("TEST: implicit midpoint on the pendulum (solve)  (should VERIFY)")
    # State (q, p).  Stage values (Q, P) solved implicitly, then update.
    # Midpoint:  Q = q + (h/2) P,  P = p - (h/2) sin Q,
    #            q := q + h P,      p := p - h sin Q.
    ctx = Context(["q", "p", "Q", "P"], ["h"])
    q, p, Q, P = ctx.coords
    h = ctx.param("h")
    from sage.all import sin

    # init stage registers to the current point
    c_init = Assign(["Q", "P"], [q, p])

    # solve  (Q, P) := solve(...,  G)
    # fixed point  xi = x^r + h G(xi),   with x^r = (Q, P) = (q, p) after init.
    # We want Q = q + (h/2) P  ->  G_Q = P/2 ;  P = p - (h/2) sin Q -> G_P = -sin(Q)/2
    xiQ, xiP = var("xiQ"), var("xiP")
    c_solve = Solve(["Q", "P"], [xiQ, xiP], [xiP / 2, -sin(xiQ) / 2], h)

    # update
    c_up = Seq(Assign(["q"], [q + h * P]),
               Assign(["p"], [p - h * sin(Q)]))

    prog = Seq(c_init, c_solve, c_up)
    w0 = ctx.omega0([("q", "p")])
    ok, _, _ = Verifier(ctx).verify(w0, prog, w0)
    assert ok, "implicit midpoint should be symplectic"


def test_gl2_pendulum():
    _hr("TEST: Gauss-Legendre 2-stage (GL-2) on the pendulum  (should VERIFY)")
    # The worked example of main.tex.
    from sage.all import sqrt, sin, QQ, Rational

    ctx = Context(["q", "p", "Q1", "P1", "Q2", "P2"], ["h"])
    q, p, Q1, P1, Q2, P2 = ctx.coords
    h = ctx.param("h")

    s3 = sqrt(SR(3))
    a11 = Rational((1, 4))
    a12 = Rational((1, 4)) - s3 / 6
    a21 = Rational((1, 4)) + s3 / 6
    a22 = Rational((1, 4))
    b1 = Rational((1, 2))
    b2 = Rational((1, 2))

    # init  Q1,P1,Q2,P2 := q,p,q,p
    c_init = Assign(["Q1", "P1", "Q2", "P2"], [q, p, q, p])

    # solve  with unknowns xi = (xi0,xi1,xi2,xi3) ~ (Q1,P1,Q2,P2)
    xi0, xi1, xi2, xi3 = var("xi0 xi1 xi2 xi3")
    G = [a11 * xi1 + a12 * xi3,
         -(a11 * sin(xi0) + a12 * sin(xi2)),
         a21 * xi1 + a22 * xi3,
         -(a21 * sin(xi0) + a22 * sin(xi2))]
    c_solve = Solve(["Q1", "P1", "Q2", "P2"], [xi0, xi1, xi2, xi3], G, h)

    # update
    c_up = Seq(Assign(["q"], [q + h * (b1 * P1 + b2 * P2)]),
               Assign(["p"], [p - h * (b1 * sin(Q1) + b2 * sin(Q2))]))

    prog = Seq(c_init, c_solve, c_up)
    w0 = ctx.omega0([("q", "p")])
    ok, _, _ = Verifier(ctx).verify(w0, prog, w0)
    assert ok, "GL-2 should be symplectic"


def test_boris_volume():
    _hr("TEST: Boris pusher (full 3D) preserves the volume form (6-form)")
    # Charged particle in an EM field; the 3D Boris integrator is
    # volume-preserving but not canonically symplectic.  We check both.
    # Every intermediate quantity (a, v^-, t, rho, s, v', v^+) is computed by
    # the program in auxiliary registers (the augmented-state device).
    ctx = Context(
        ["x1", "x2", "x3", "v1", "v2", "v3", "a",
         "vm1", "vm2", "vm3", "t1", "t2", "t3", "rho",
         "s1", "s2", "s3", "vp1", "vp2", "vp3", "vu1", "vu2", "vu3"],
        ["dt", "q", "m", "E1", "E2", "E3", "B1", "B2", "B3"])

    g = ctx.coord
    x1, x2, x3 = g("x1"), g("x2"), g("x3")
    v1, v2, v3 = g("v1"), g("v2"), g("v3")
    a = g("a")
    vm1, vm2, vm3 = g("vm1"), g("vm2"), g("vm3")
    t1, t2, t3 = g("t1"), g("t2"), g("t3")
    rho = g("rho")
    s1, s2, s3 = g("s1"), g("s2"), g("s3")
    vp1, vp2, vp3 = g("vp1"), g("vp2"), g("vp3")
    vu1, vu2, vu3 = g("vu1"), g("vu2"), g("vu3")
    dt, q, m = ctx.param("dt"), ctx.param("q"), ctx.param("m")
    E1, E2, E3 = ctx.param("E1"), ctx.param("E2"), ctx.param("E3")
    B1, B2, B3 = ctx.param("B1"), ctx.param("B2"), ctx.param("B3")

    prog = Seq(
        Assign(["a"], [dt * q / (2 * m)]),
        Assign(["vm1", "vm2", "vm3"], [v1 + a * E1, v2 + a * E2, v3 + a * E3]),
        Assign(["t1", "t2", "t3"], [a * B1, a * B2, a * B3]),
        Assign(["rho"], [1 + t1 ** 2 + t2 ** 2 + t3 ** 2]),
        Assign(["s1", "s2", "s3"], [2 * t1 / rho, 2 * t2 / rho, 2 * t3 / rho]),
        Assign(["vp1", "vp2", "vp3"],
               [vm1 + vm2 * t3 - vm3 * t2,
                vm2 + vm3 * t1 - vm1 * t3,
                vm3 + vm1 * t2 - vm2 * t1]),
        Assign(["vu1", "vu2", "vu3"],
               [vm1 + vp2 * s3 - vp3 * s2,
                vm2 + vp3 * s1 - vp1 * s3,
                vm3 + vp1 * s2 - vp2 * s1]),
        Assign(["v1", "v2", "v3"], [vu1 + a * E1, vu2 + a * E2, vu3 + a * E3]),
        Assign(["x1", "x2", "x3"], [x1 + dt * v1, x2 + dt * v2, x3 + dt * v3]),
    )

    # volume form  dx1^dx2^dx3^dv1^dv2^dv3  is preserved
    vol = (ctx.dx("x1").wedge(ctx.dx("x2")).wedge(ctx.dx("x3"))
           .wedge(ctx.dx("v1")).wedge(ctx.dx("v2")).wedge(ctx.dx("v3")))
    ok_vol, _, _ = Verifier(ctx).verify(vol, prog, vol)
    assert ok_vol, "Boris should preserve phase-space volume"

    # canonical 2-form  sum_i dxi ^ dvi  is NOT preserved
    w2 = (ctx.dx("x1").wedge(ctx.dx("v1"))
          + ctx.dx("x2").wedge(ctx.dx("v2"))
          + ctx.dx("x3").wedge(ctx.dx("v3")))
    ok_w2, _, _ = Verifier(ctx).verify(w2, prog, w2, verbose=False)
    assert not ok_w2, "Boris is not canonically symplectic; verifier must reject"
    print("  canonical 2-form sum dxi^dvi : correctly REJECTED (not symplectic)")


def test_conformal_for():
    _hr("TEST: conformal-symplectic for-loop  (indexed assertions)")
    # A dissipative one-step map with kappa = e^{-lambda h}:
    #   q := q + h p ;  p := e^{-lambda h} (p - h sin q)
    # Pullback of omega0 = dq^dp gives e^{-lambda h} omega0  (conformal).
    # With A_t = e^{t lambda h} omega0 we get  |= {A_t} C {A_{t+1}}.
    from sage.all import sin, exp

    ctx = Context(["q", "p"], ["h", "lam", "t"])
    q, p = ctx.coord("q"), ctx.coord("p")
    h, lam, t = ctx.param("h"), ctx.param("lam"), ctx.param("t")

    body = Seq(Assign(["q"], [q + h * p]),
               Assign(["p"], [exp(-lam * h) * (p - h * sin(q))]))

    def A_of_t(tt):
        return ctx.omega0([("q", "p")]).scale(exp(tt * lam * h))

    ok, _ = Verifier(ctx).verify_for(A_of_t, body, T=10, t=t)
    assert ok, "conformal for-loop premise should hold for all t"


def main():
    test_form_algebra()
    test_symplectic_euler()
    test_explicit_euler()
    test_implicit_midpoint()
    test_gl2_pendulum()
    test_boris_volume()
    test_conformal_for()
    _hr("ALL TESTS PASSED")


if __name__ == "__main__":
    main()
