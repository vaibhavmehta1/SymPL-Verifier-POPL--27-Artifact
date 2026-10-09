# Implicit midpoint rule -- H1: harmonic oscillator.
#
# Separable H = p^2/2 + q^2/2, V'(q) = q.  The midpoint values (Q,P) solve the
# implicit stage equations
#       Q = q + (h/2) P,      P = p - (h/2) Q,
# then the state advances with the midpoint slopes
#       q := q + h P,         p := p - h Q.
# The stage fixed point xi = x^r + h G(xi) has kernel G = (P/2, -Q/2); the
# verifier handles the implicit solve via the implicit-function theorem.
# Symplectic, 2nd-order.

coords q p Q P
params h

pre  dq ^ dp
post dq ^ dp

program {
    (Q, P) := (q, p)
    (Q, P) := solve (a, b) -> ( b/2, -a/2 )
    q := q + h*P
    p := p - h*Q
}
