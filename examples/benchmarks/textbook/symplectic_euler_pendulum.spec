# Symplectic (semi-implicit) Euler -- H2: pendulum.
#
# Separable H(q,p) = p^2/2 - cos q, so V'(q) = sin q.  Momentum kick with the
# old position, position drift with the new momentum:
#       p := p - h sin q
#       q := q + h p
# Symplectic, 1st-order.  Middle rung: a transcendental (non-polynomial) force.

coords q p
params h

pre  dq ^ dp
post dq ^ dp

program {
    p := p - h*sin(q)
    q := q + h*p
}
