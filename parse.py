"""
parse.py
========

A frontend parser for the symplectic verifier.  It turns a concrete textual
specification into the `Context`, precondition/postcondition `Form`s and the
`Command` program consumed by `symplectic_verifier.Verifier`.

Specification syntax
--------------------

    coords q p Q1 P1 Q2 P2        # state coordinates (define the dimension)
    params h                      # external parameters (step size, etc.)

    pre  dq ^ dp                  # precondition  (a differential form)
    post dq ^ dp                  # postcondition (a differential form)

    program {
        (Q1, P1, Q2, P2) := (q, p, q, p)            # simultaneous explicit
        (Q1, P1, Q2, P2) := solve (a, b, c, d) -> ( # simultaneous implicit
            (1/4)*b + (1/4 - sqrt(3)/6)*d,
            -((1/4)*sin(a) + (1/4 - sqrt(3)/6)*sin(c)),
            (1/4 + sqrt(3)/6)*b + (1/4)*d,
            -((1/4 + sqrt(3)/6)*sin(a) + (1/4)*sin(c))
        )
        q := q + h*((1/2)*P1 + (1/2)*P2)
        p := p - h*((1/2)*sin(Q1) + (1/2)*sin(Q2))
    }

Expression / form notation
--------------------------
  * scalars: + - * /, ** for powers, parentheses, and the smooth functions
    sin cos tan exp log sqrt sinh cosh tanh arctan arcsin(=asin)
    arccos(=acos) sincu(=sinc, unnormalized sin t / t) and the two-argument
    arctan2(=atan2);  functions may take a comma-separated argument list;
  * `dq`, `dp`, ... are the coordinate 1-forms d(q), d(p);
  * `^` is the wedge product (NOTE: powers use `**`, never `^`);
  * `d( expr )` is the exterior derivative of a form;
  * `+`/`-` combine forms; a scalar times a form scales it.

Commands
--------
  * explicit:   x := e            or   (x1,...,xk) := (e1,...,ek)
  * implicit:   (x1,...,xm) := solve (u1,...,um) -> (G1,...,Gm)
                (the unknowns u_i name the components of the fixed point;
                 the step size defaults to the parameter named `h`)
  * loops:      while^N { ... }        for^T { ... }
  * statements are separated by newlines or `;`.

Run:
    sage -python parse.py spec_file        # parse + verify a spec
    sage -python parse.py                  # run the built-in demos
"""

from sage.all import SR, var, Integer, QQ
from sage.all import (sin, cos, tan, exp, log, sqrt, sinh, cosh, tanh,
                      arctan, arcsin, arccos, arctan2)

from symplectic_verifier import (
    Context, Form, Verifier,
    Skip, Assign, Solve, Seq, While, For,
)


def _sincu(t):
    """Unnormalized sinc, sin(t)/t (smooth, with sincu(0) = 1).  Matches
    BeamTracking.jl's `sincu`.  Kept symbolic; never evaluated at t = 0."""
    return sin(t) / t


# Smooth functions usable in expressions.  Each entry is a plain callable, so
# multi-argument functions (e.g. arctan2) are dispatched by parsed arity.
FUNCTIONS = {
    "sin": sin, "cos": cos, "tan": tan, "exp": exp, "log": log,
    "sqrt": sqrt, "sinh": sinh, "cosh": cosh, "tanh": tanh,
    "arctan": arctan, "arcsin": arcsin, "asin": arcsin,
    "arccos": arccos, "acos": arccos,
    "arctan2": arctan2, "atan2": arctan2,            # two-argument
    "sincu": _sincu, "sinc": _sincu,
}


# ---------------------------------------------------------------------------
#  Tokenizer
# ---------------------------------------------------------------------------

class Token:
    __slots__ = ("kind", "val", "pos")

    def __init__(self, kind, val, pos):
        self.kind = kind     # NAME, NUM, OP, NEWLINE, EOF
        self.val = val
        self.pos = pos

    def __repr__(self):
        return "{}({!r})".format(self.kind, self.val)


# multi-character operators, matched longest-first
_MULTI_OPS = [":=", "->", "**"]
_SINGLE_OPS = set("+-*/^(){},[]")


def tokenize(text):
    """Lex `text`.  Newlines inside ( ) or [ ] are suppressed so that
    multi-line solve kernels and tuples work; `;` acts like a newline."""
    toks = []
    i, n = 0, len(text)
    depth = 0                     # nesting depth of () and []
    line = 1

    def emit(kind, val, pos):
        toks.append(Token(kind, val, pos))

    while i < n:
        c = text[i]

        # comments
        if c == "#":
            while i < n and text[i] != "\n":
                i += 1
            continue

        # newlines / statement separators
        if c == "\n":
            if depth == 0 and toks and toks[-1].kind != "NEWLINE":
                emit("NEWLINE", "\\n", i)
            line += 1
            i += 1
            continue
        if c == ";":
            if depth == 0 and toks and toks[-1].kind != "NEWLINE":
                emit("NEWLINE", ";", i)
            i += 1
            continue

        if c in " \t\r":
            i += 1
            continue

        # identifiers
        if c.isalpha() or c == "_":
            j = i + 1
            while j < n and (text[j].isalnum() or text[j] == "_"):
                j += 1
            emit("NAME", text[i:j], i)
            i = j
            continue

        # numbers (integer or decimal)
        if c.isdigit() or (c == "." and i + 1 < n and text[i + 1].isdigit()):
            j = i + 1
            while j < n and (text[j].isdigit() or text[j] == "."):
                j += 1
            emit("NUM", text[i:j], i)
            i = j
            continue

        # multi-char operators
        matched = False
        for op in _MULTI_OPS:
            if text.startswith(op, i):
                emit("OP", op, i)
                i += len(op)
                matched = True
                break
        if matched:
            continue

        # single-char operators
        if c in _SINGLE_OPS:
            if c in "([":
                depth += 1
            elif c in ")]":
                depth = max(0, depth - 1)
            emit("OP", c, i)
            i += 1
            continue

        raise SyntaxError("unexpected character {!r} at offset {}".format(c, i))

    emit("EOF", None, n)
    return toks


# ---------------------------------------------------------------------------
#  Parser
# ---------------------------------------------------------------------------

class Parser:
    def __init__(self, text):
        self.toks = tokenize(text)
        self.pos = 0
        self.ctx = None
        self.coord_names = []
        self.param_names = []
        self.namespace = {}        # name -> sage symbol (coords + params + locals)

    # ---- token stream helpers ------------------------------------------

    def peek(self):
        return self.toks[self.pos]

    def next(self):
        t = self.toks[self.pos]
        self.pos += 1
        return t

    def at_op(self, *vals):
        t = self.peek()
        return t.kind == "OP" and t.val in vals

    def at_name(self, *vals):
        t = self.peek()
        return t.kind == "NAME" and (not vals or t.val in vals)

    def eat_op(self, val):
        t = self.next()
        if not (t.kind == "OP" and t.val == val):
            raise SyntaxError("expected {!r}, got {!r} at offset {}".format(val, t.val, t.pos))
        return t

    def skip_newlines(self):
        while self.peek().kind == "NEWLINE":
            self.next()

    # ---- top level ------------------------------------------------------

    def parse_spec(self):
        """Parse a full specification; returns (ctx, pre, program, post)."""
        pre = post = program = None
        while True:
            self.skip_newlines()
            t = self.peek()
            if t.kind == "EOF":
                break
            if t.kind != "NAME":
                raise SyntaxError("expected a directive, got {!r}".format(t.val))
            kw = t.val
            if kw == "coords":
                self.next()
                self.coord_names = self._name_list_until_newline()
            elif kw in ("params", "param"):
                self.next()
                self.param_names = self._name_list_until_newline()
            elif kw == "pre":
                self.next()
                self._ensure_ctx()
                pre = self.parse_form()
                self._expect_newline()
            elif kw == "post":
                self.next()
                self._ensure_ctx()
                post = self.parse_form()
                self._expect_newline()
            elif kw == "program":
                self.next()
                self._ensure_ctx()
                program = self.parse_block()
            else:
                raise SyntaxError("unknown directive {!r}".format(kw))

        if program is None:
            raise SyntaxError("specification has no `program { ... }` block")
        if pre is None or post is None:
            raise SyntaxError("specification needs both `pre` and `post`")
        return self.ctx, pre, program, post

    def _name_list_until_newline(self):
        names = []
        while self.peek().kind == "NAME":
            names.append(self.next().val)
        self._expect_newline()
        return names

    def _expect_newline(self):
        t = self.peek()
        if t.kind in ("NEWLINE", "EOF"):
            if t.kind == "NEWLINE":
                self.next()
        else:
            raise SyntaxError("expected end of line, got {!r} at offset {}".format(t.val, t.pos))

    def _ensure_ctx(self):
        if self.ctx is None:
            if not self.coord_names:
                raise SyntaxError("`coords` must be declared before pre/post/program")
            self.ctx = Context(self.coord_names, self.param_names)
            for nm in self.coord_names:
                self.namespace[nm] = self.ctx.coord(nm)
            for nm in self.param_names:
                self.namespace[nm] = self.ctx.param(nm)

    # ---- program / statements ------------------------------------------

    def parse_block(self):
        self.eat_op("{")
        cmds = []
        while True:
            self.skip_newlines()
            if self.at_op("}"):
                self.next()
                break
            if self.peek().kind == "EOF":
                raise SyntaxError("unclosed `{`")
            cmds.append(self.parse_statement())
        if len(cmds) == 1:
            return cmds[0]
        return Seq(*cmds)

    def parse_statement(self):
        if self.at_name("while", "for"):
            kw = self.next().val
            self.eat_op("^")
            n = self.next()
            if n.kind != "NUM":
                raise SyntaxError("expected loop bound after {}^".format(kw))
            bound = int(n.val)
            body = self.parse_block()
            return While(bound, body) if kw == "while" else For(bound, body)
        if self.at_name("skip"):
            self.next()
            self._expect_newline()
            return Skip()
        return self.parse_assignment()

    def parse_assignment(self):
        targets = self.parse_lhs()

        self.eat_op(":=")

        if self.at_name("solve"):
            cmd = self.parse_solve(targets)
            self._expect_newline()
            return cmd

        # explicit assignment
        if len(targets) == 1:
            expr = self._as_scalar(self.parse_form())
            self._expect_newline()
            return Assign(targets, [expr])
        else:
            exprs = self.parse_paren_scalar_list()
            if len(exprs) != len(targets):
                raise SyntaxError("assignment arity mismatch: {} targets vs {} expressions"
                                  .format(len(targets), len(exprs)))
            self._expect_newline()
            return Assign(targets, exprs)

    def parse_lhs(self):
        if self.at_op("("):
            self.next()
            names = [self._expect_coord_name()]
            while self.at_op(","):
                self.next()
                names.append(self._expect_coord_name())
            self.eat_op(")")
            return names
        return [self._expect_coord_name()]

    def _expect_coord_name(self):
        t = self.next()
        if t.kind != "NAME" or t.val not in self.ctx.index:
            raise SyntaxError("expected a coordinate name, got {!r}".format(t.val))
        return t.val

    def parse_solve(self, targets):
        self.next()                                   # consume 'solve'
        # unknown symbols
        self.eat_op("(")
        unknown_names = [self._next_name("unknown")]
        while self.at_op(","):
            self.next()
            unknown_names.append(self._next_name("unknown"))
        self.eat_op(")")
        if len(unknown_names) != len(targets):
            raise SyntaxError("solve: {} unknowns for {} targets"
                              .format(len(unknown_names), len(targets)))
        self.eat_op("->")

        # bind unknowns into the namespace while parsing G
        xi = [var(nm) for nm in unknown_names]
        saved = {nm: self.namespace.get(nm) for nm in unknown_names}
        for nm, sym in zip(unknown_names, xi):
            self.namespace[nm] = sym
        try:
            G = self.parse_paren_scalar_list()
        finally:
            for nm in unknown_names:
                if saved[nm] is None:
                    self.namespace.pop(nm, None)
                else:
                    self.namespace[nm] = saved[nm]

        if len(G) != len(targets):
            raise SyntaxError("solve: {} kernel components for {} targets"
                              .format(len(G), len(targets)))
        if "h" not in self.ctx.params:
            raise SyntaxError("solve requires a parameter named `h` (the step size)")
        return Solve(targets, xi, G, self.ctx.param("h"))

    def _next_name(self, what):
        t = self.next()
        if t.kind != "NAME":
            raise SyntaxError("expected {} name, got {!r}".format(what, t.val))
        return t.val

    def parse_paren_scalar_list(self):
        self.eat_op("(")
        exprs = [self._as_scalar(self.parse_form())]
        while self.at_op(","):
            self.next()
            exprs.append(self._as_scalar(self.parse_form()))
        self.eat_op(")")
        return exprs

    # ---- expression / form grammar -------------------------------------
    #   form  := wedge (('+'|'-') wedge)*
    #   wedge := mul ('^' mul)*
    #   mul   := unary (('*'|'/') unary)*
    #   unary := ('+'|'-') unary | power
    #   power := primary ('**' unary)?
    #   primary := NUM | dName | d(form) | func(args) | name | '(' form ')'

    def parse_form(self):
        left = self.parse_wedge()
        while self.at_op("+", "-"):
            op = self.next().val
            right = self.parse_wedge()
            left = left + right if op == "+" else left - right
        return left

    def parse_wedge(self):
        left = self.parse_mul()
        while self.at_op("^"):
            self.next()
            right = self.parse_mul()
            left = left.wedge(right)
        return left

    def parse_mul(self):
        left = self.parse_unary()
        while self.at_op("*", "/"):
            op = self.next().val
            right = self.parse_unary()
            left = self._mul(left, right) if op == "*" else self._div(left, right)
        return left

    def parse_unary(self):
        if self.at_op("+"):
            self.next()
            return self.parse_unary()
        if self.at_op("-"):
            self.next()
            return self.parse_unary().scale(SR(-1))
        return self.parse_power()

    def parse_power(self):
        base = self.parse_primary()
        if self.at_op("**"):
            self.next()
            exp_form = self.parse_unary()
            b = self._as_scalar(base)
            e = self._as_scalar(exp_form)
            return self.ctx.scalar(b ** e)
        return base

    def parse_primary(self):
        t = self.peek()

        if t.kind == "NUM":
            self.next()
            if "." in t.val:
                return self.ctx.scalar(SR(QQ(t.val)))
            return self.ctx.scalar(SR(Integer(t.val)))

        if self.at_op("("):
            self.next()
            f = self.parse_form()
            self.eat_op(")")
            return f

        if t.kind == "NAME":
            name = t.val

            # exterior derivative  d( ... )
            if name == "d" and self.toks[self.pos + 1].kind == "OP" \
                    and self.toks[self.pos + 1].val == "(":
                self.next()                 # 'd'
                self.next()                 # '('
                inner = self.parse_form()
                self.eat_op(")")
                return inner.d()

            # smooth function application (one or more comma-separated args)
            if name in FUNCTIONS and self.toks[self.pos + 1].kind == "OP" \
                    and self.toks[self.pos + 1].val == "(":
                self.next()                 # function name
                self.eat_op("(")
                args = [self._as_scalar(self.parse_form())]
                while self.at_op(","):
                    self.next()
                    args.append(self._as_scalar(self.parse_form()))
                self.eat_op(")")
                return self.ctx.scalar(FUNCTIONS[name](*args))

            # coordinate differential  dq, dQ1, ...
            if name.startswith("d") and len(name) > 1 and name[1:] in self.ctx.index:
                self.next()
                return self.ctx.dx(name[1:])

            # known symbol (coordinate / parameter / bound unknown)
            if name in self.namespace:
                self.next()
                return self.ctx.scalar(self.namespace[name])

            raise SyntaxError(
                "unknown name {!r} at offset {} (not a coord/param/function; "
                "did you forget to declare it?)".format(name, t.pos))

        raise SyntaxError("unexpected token {!r} at offset {}".format(t.val, t.pos))

    # ---- scalar/form helpers -------------------------------------------

    def _as_scalar(self, form):
        keys = [k for k in form.terms if k != ()]
        if keys:
            raise SyntaxError("expected a scalar expression but got a "
                              "degree-{} form".format(len(keys[0])))
        return form.terms.get((), SR(0))

    def _is_scalar(self, form):
        return all(k == () for k in form.terms)

    def _mul(self, a, b):
        if self._is_scalar(a):
            return b.scale(self._as_scalar(a))
        if self._is_scalar(b):
            return a.scale(self._as_scalar(b))
        raise SyntaxError("cannot `*` two forms of positive degree; use `^` (wedge)")

    def _div(self, a, b):
        if not self._is_scalar(b):
            raise SyntaxError("cannot divide by a form of positive degree")
        return a.scale(SR(1) / self._as_scalar(b))


# ---------------------------------------------------------------------------
#  Convenience API
# ---------------------------------------------------------------------------

def parse_spec(text):
    """Parse a specification string -> (ctx, pre, program, post)."""
    return Parser(text).parse_spec()


def verify_spec(text, verbose=True):
    """Parse and verify a specification string.  Returns (ok, residual)."""
    ctx, pre, program, post = parse_spec(text)
    ok, _computed, residual = Verifier(ctx).verify(pre, program, post, verbose=verbose)
    return ok, residual


# ---------------------------------------------------------------------------
#  Demos / CLI
# ---------------------------------------------------------------------------

DEMOS = {
    "symplectic Euler (pendulum)": r"""
        coords q p
        params h
        pre  dq ^ dp
        post dq ^ dp
        program {
            p := p - h*sin(q)
            q := q + h*p
        }
    """,

    "explicit Euler (pendulum, NOT symplectic)": r"""
        coords q p
        params h
        pre  dq ^ dp
        post dq ^ dp
        program {
            (q, p) := (q + h*p, p - h*sin(q))
        }
    """,

    "implicit midpoint (pendulum, solve)": r"""
        coords q p Q P
        params h
        pre  dq ^ dp
        post dq ^ dp
        program {
            (Q, P) := (q, p)
            (Q, P) := solve (a, b) -> ( b/2, -sin(a)/2 )
            q := q + h*P
            p := p - h*sin(Q)
        }
    """,

    "Gauss-Legendre 2-stage (pendulum)": r"""
        coords q p Q1 P1 Q2 P2
        params h
        pre  dq ^ dp
        post dq ^ dp
        program {
            (Q1, P1, Q2, P2) := (q, p, q, p)
            (Q1, P1, Q2, P2) := solve (a, b, c, e) -> (
                (1/4)*b + (1/4 - sqrt(3)/6)*e,
                -((1/4)*sin(a) + (1/4 - sqrt(3)/6)*sin(c)),
                (1/4 + sqrt(3)/6)*b + (1/4)*e,
                -((1/4 + sqrt(3)/6)*sin(a) + (1/4)*sin(c))
            )
            q := q + h*((1/2)*P1 + (1/2)*P2)
            p := p - h*((1/2)*sin(Q1) + (1/2)*sin(Q2))
        }
    """,

    "Boris pusher 3D (volume form dx1^dx2^dx3^dv1^dv2^dv3 preserved)": r"""
        coords x1 x2 x3 v1 v2 v3 a vm1 vm2 vm3 t1 t2 t3 rho s1 s2 s3 vp1 vp2 vp3 vu1 vu2 vu3
        params dt q m E1 E2 E3 B1 B2 B3
        pre  dx1 ^ dx2 ^ dx3 ^ dv1 ^ dv2 ^ dv3
        post dx1 ^ dx2 ^ dx3 ^ dv1 ^ dv2 ^ dv3
        program {
            a := dt*q/(2*m)
            (vm1, vm2, vm3) := (v1 + a*E1, v2 + a*E2, v3 + a*E3)
            (t1, t2, t3) := (a*B1, a*B2, a*B3)
            rho := 1 + t1**2 + t2**2 + t3**2
            (s1, s2, s3) := (2*t1/rho, 2*t2/rho, 2*t3/rho)
            (vp1, vp2, vp3) := (vm1 + vm2*t3 - vm3*t2,
                                vm2 + vm3*t1 - vm1*t3,
                                vm3 + vm1*t2 - vm2*t1)
            (vu1, vu2, vu3) := (vm1 + vp2*s3 - vp3*s2,
                                vm2 + vp3*s1 - vp1*s3,
                                vm3 + vp1*s2 - vp2*s1)
            (v1, v2, v3) := (vu1 + a*E1, vu2 + a*E2, vu3 + a*E3)
            (x1, x2, x3) := (x1 + dt*v1, x2 + dt*v2, x3 + dt*v3)
        }
    """,

    "Boris pusher 3D (canonical 2-form sum dxi^dvi, NOT symplectic)": r"""
        coords x1 x2 x3 v1 v2 v3 a vm1 vm2 vm3 t1 t2 t3 rho s1 s2 s3 vp1 vp2 vp3 vu1 vu2 vu3
        params dt q m E1 E2 E3 B1 B2 B3
        pre  dx1 ^ dv1 + dx2 ^ dv2 + dx3 ^ dv3
        post dx1 ^ dv1 + dx2 ^ dv2 + dx3 ^ dv3
        program {
            a := dt*q/(2*m)
            (vm1, vm2, vm3) := (v1 + a*E1, v2 + a*E2, v3 + a*E3)
            (t1, t2, t3) := (a*B1, a*B2, a*B3)
            rho := 1 + t1**2 + t2**2 + t3**2
            (s1, s2, s3) := (2*t1/rho, 2*t2/rho, 2*t3/rho)
            (vp1, vp2, vp3) := (vm1 + vm2*t3 - vm3*t2,
                                vm2 + vm3*t1 - vm1*t3,
                                vm3 + vm1*t2 - vm2*t1)
            (vu1, vu2, vu3) := (vm1 + vp2*s3 - vp3*s2,
                                vm2 + vp3*s1 - vp1*s3,
                                vm3 + vp1*s2 - vp2*s1)
            (v1, v2, v3) := (vu1 + a*E1, vu2 + a*E2, vu3 + a*E3)
            (x1, x2, x3) := (x1 + dt*v1, x2 + dt*v2, x3 + dt*v3)
        }
    """,

    "Stormer-Verlet via while^3 (should stay symplectic)": r"""
        coords q p
        params h
        pre  dq ^ dp
        post dq ^ dp
        program {
            while^3 {
                p := p - (h/2)*sin(q)
                q := q + h*p
                p := p - (h/2)*sin(q)
            }
        }
    """,
}


def run_demos():
    for title, spec in DEMOS.items():
        print("\n" + "=" * 72)
        print("DEMO:", title)
        print("=" * 72)
        ok, _ = verify_spec(spec)
        expect_fail = "NOT symplectic" in title
        assert ok != expect_fail, "demo {!r} gave unexpected result".format(title)
    print("\n" + "=" * 72)
    print("ALL DEMOS BEHAVED AS EXPECTED")
    print("=" * 72)


def main(argv):
    if len(argv) > 1:
        with open(argv[1]) as fh:
            text = fh.read()
        ok, _ = verify_spec(text)
        return 0 if ok else 1
    run_demos()
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv))
