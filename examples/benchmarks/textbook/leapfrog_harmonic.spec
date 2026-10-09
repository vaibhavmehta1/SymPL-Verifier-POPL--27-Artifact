# Leapfrog / velocity-Verlet (kick-drift-kick) -- H1: harmonic oscillator.
#
# Separable H = p^2/2 + q^2/2, V'(q) = q.  Half momentum kick, full position
# drift, half momentum kick:
#       p := p - (h/2) q
#       q := q + h p
#       p := p - (h/2) q
# Symmetric composition of shears: symplectic and 2nd-order.

coords q p
params h

pre  dq ^ dp
post dq ^ dp

program {
    p := p - (h/2)*q
    q := q + h*p
    p := p - (h/2)*q
}
