# Gauss-Legendre 2-stage (GL-2) applied to the pendulum  q' = p, p' = -sin q.
# Verifies that one step preserves the canonical symplectic form omega0 = dq ^ dp.

coords q p Q1 P1 Q2 P2
params h

pre  dq ^ dp
post dq ^ dp

program {
    # initialise the stage registers to the current point
    (Q1, P1, Q2, P2) := (q, p, q, p)

    # solve the GL-2 stage equations  xi = x^r + h*G(xi)
    (Q1, P1, Q2, P2) := solve (a, b, c, e) -> (
        (1/4)*b + (1/4 - sqrt(3)/6)*e,
        -((1/4)*sin(a) + (1/4 - sqrt(3)/6)*sin(c)),
        (1/4 + sqrt(3)/6)*b + (1/4)*e,
        -((1/4 + sqrt(3)/6)*sin(a) + (1/4)*sin(c))
    )

    # final update
    q := q + h*((1/2)*P1 + (1/2)*P2)
    p := p - h*((1/2)*sin(Q1) + (1/2)*sin(Q2))
}
