# Stoermer-Verlet (drift-kick-drift) -- H2: pendulum.
#
# Separable H = p^2/2 - cos q, V'(q) = sin q:
#       q := q + (h/2) p
#       p := p - h sin q
#       q := q + (h/2) p
# Symplectic, 2nd-order.  (Same map as examples/stormer_verlet.spec.)

coords q p
params h

pre  dq ^ dp
post dq ^ dp

program {
    q := q + (h/2)*p
    p := p - h*sin(q)
    q := q + (h/2)*p
}
