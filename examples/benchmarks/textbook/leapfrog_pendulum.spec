# Leapfrog / velocity-Verlet (kick-drift-kick) -- H2: pendulum.
#
# Separable H = p^2/2 - cos q, V'(q) = sin q.  Half kick, drift, half kick:
#       p := p - (h/2) sin q
#       q := q + h p
#       p := p - (h/2) sin q
# Symplectic, 2nd-order.  (Same map as examples/leapfrog.spec; included here so
# the benchmark ladder is self-contained.)

coords q p
params h

pre  dq ^ dp
post dq ^ dp

program {
    p := p - (h/2)*sin(q)
    q := q + h*p
    p := p - (h/2)*sin(q)
}
