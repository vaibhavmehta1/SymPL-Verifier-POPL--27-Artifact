"""Implicit-solve regressions: sage -python -m unittest discover -s tests -p test_implicit_solve.py -v"""

import unittest

try:
    from sage.all import SR, var, sin
except ImportError:
    SR = None
else:
    from parse import verify_spec
    from symplectic_verifier import Context, Assign, Solve, Seq, While, For, Verifier


@unittest.skipIf(SR is None, "run with sage -python")
class ImplicitSolveTests(unittest.TestCase):
    def test_successive_solves_reject_nonsymplectic_map_with_exact_residual(self):
        ok, residual = verify_spec('''
coords q p r
params h
pre dq ^ dp
post dq ^ dp
program {
    q := solve (a) -> (0)
    r := q
    q := solve (a) -> (1)
    r := q - r
    p := (1+r)*p
}
''', verbose=False)
        self.assertFalse(ok)
        self.assertEqual(residual, {(0, 1): var('h')})

    def test_user_parameter_named_phi_q_is_not_a_solution(self):
        ok, residual = verify_spec('''
coords q
params h phi_q
pre phi_q
post q
program {
    q := solve (a) -> (0)
}
''', verbose=False)
        self.assertFalse(ok)
        self.assertEqual(residual, {(): var('q') - var('phi_q')})

    def test_unresolved_nonlinear_solve_cannot_alias_user_symbols(self):
        ctx = Context(['q', 'phi_q'], ['h', '_sv_phi_0'])
        q, collision = ctx.coords
        a = var('a')
        cmd = Solve(['q'], [a], [a**2], ctx.param('h'))
        for candidate in (collision, ctx.param('_sv_phi_0')):
            ok, _, _ = Verifier(ctx).verify(ctx.scalar(candidate), cmd,
                                            ctx.scalar(q), verbose=False)
            self.assertFalse(ok)

    def test_reused_nonlinear_solve_command_has_distinct_values(self):
        ctx = Context(['q', 'r'], ['h'])
        q, r = ctx.coords
        a = var('a')
        cmd = Solve(['q'], [a], [a**2], ctx.param('h'))
        program = Seq(cmd, Assign(['r'], [q]), cmd, Assign(['r'], [q-r]))
        ok, _, _ = Verifier(ctx).verify(ctx.zero(), program, ctx.scalar(r),
                                        verbose=False)
        self.assertFalse(ok)

    def test_each_loop_iteration_gets_distinct_solution_values(self):
        ctx = Context(['q', 'r'], ['h'])
        q, r = ctx.coords
        a = var('a')
        body = Seq(Assign(['r'], [q]),
                   Solve(['q'], [a], [a**2], ctx.param('h')))
        for loop in (While, For):
            with self.subTest(loop=loop.__name__):
                ok, _, _ = Verifier(ctx).verify(ctx.zero(), loop(2, body),
                                                ctx.scalar(q-r), verbose=False)
                self.assertFalse(ok)

    def test_identity_and_state_only_kernels_are_exact_simultaneous_maps(self):
        ctx = Context(['q', 'p'], ['h'])
        q, p = ctx.coords
        h = ctx.param('h')
        a, b = var('a b')
        verifier = Verifier(ctx)
        post = ctx.scalar(q*p) + ctx.omega0([('q', 'p')])
        identity = Solve(['q', 'p'], [a, b], [0, 0], h)
        self.assertTrue(verifier.verify(post, identity, post, verbose=False)[0])
        implicit = Solve(['q', 'p'], [a, b], [p, -sin(q)], h)
        explicit = Assign(['q', 'p'], [q+h*p, p-h*sin(q)])
        self.assertTrue(verifier.pull(implicit, post).equals(verifier.pull(explicit, post)))

    def test_nonlinear_midpoint_remains_symplectic_when_repeated(self):
        ctx = Context(['q', 'p'], ['h'])
        q, p = ctx.coords
        a, b = var('a b')
        cmd = Solve(['q', 'p'], [a, b],
                    [(p+b)/2, -sin((q+a)/2)], ctx.param('h'))
        omega = ctx.omega0([('q', 'p')])
        verifier = Verifier(ctx)
        for program in (cmd, Seq(cmd, cmd), While(2, cmd)):
            self.assertTrue(verifier.verify(omega, program, omega, verbose=False)[0])

    def test_caller_supplied_phi_cannot_bypass_freshness(self):
        ctx = Context(['q'], ['h', 'alias'])
        a = var('a')
        with self.assertRaisesRegex(ValueError, 'phi'):
            Solve(['q'], [a], [a**2], ctx.param('h'), phi=[ctx.param('alias')])


if __name__ == '__main__':
    unittest.main()
