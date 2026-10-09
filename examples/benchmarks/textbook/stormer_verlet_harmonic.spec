# Stoermer-Verlet (drift-kick-drift) -- H1: harmonic oscillator.
#
# Adjoint arrangement of leapfrog.  Separable H = p^2/2 + q^2/2, V'(q) = q:
#       q := q + (h/2) p
#       p := p - h q
#       q := q + (h/2) p
# Symplectic, 2nd-order.

coords q p
params h

pre  dq ^ dp
post dq ^ dp

program {
    q := q + (h/2)*p
    p := p - h*q
    q := q + (h/2)*p
}
