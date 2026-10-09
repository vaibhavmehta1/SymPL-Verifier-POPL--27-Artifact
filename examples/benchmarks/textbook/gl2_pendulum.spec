# Gauss-Legendre 2-stage (GL-2, 4th-order) -- H2: pendulum.
#
# Separable H = p^2/2 - cos q, V'(q) = sin q.  Same Butcher tableau as the
# harmonic case; stage forces use sin.  Symplectic, 4th-order.  (Same map as
# examples/gl2.spec, included so the ladder is self-contained.)

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
