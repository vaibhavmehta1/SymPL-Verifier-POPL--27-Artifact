# Gauss-Legendre 2-stage (GL-2, 4th-order) -- H1: harmonic oscillator.
#
# Separable H = p^2/2 + q^2/2, V'(q) = q.  Butcher tableau
#       a11 = 1/4,             a12 = 1/4 - sqrt(3)/6,
#       a21 = 1/4 + sqrt(3)/6, a22 = 1/4,            b1 = b2 = 1/2.
# Stage values (Q1,P1,Q2,P2) solve  Qi = q + h sum_j a_ij P_j,
#                                    Pi = p - h sum_j a_ij V'(Qj),
# then  q := q + h(b1 P1 + b2 P2),  p := p - h(b1 V'(Q1) + b2 V'(Q2)).
# Symplectic, 4th-order.

coords q p Q1 P1 Q2 P2
params h

pre  dq ^ dp
post dq ^ dp

program {
    (Q1, P1, Q2, P2) := (q, p, q, p)
    (Q1, P1, Q2, P2) := solve (a, b, c, e) -> (
        (1/4)*b + (1/4 - sqrt(3)/6)*e,
        -((1/4)*a + (1/4 - sqrt(3)/6)*c),
        (1/4 + sqrt(3)/6)*b + (1/4)*e,
        -((1/4 + sqrt(3)/6)*a + (1/4)*c)
    )
    q := q + h*((1/2)*P1 + (1/2)*P2)
    p := p - h*((1/2)*Q1 + (1/2)*Q2)
}
