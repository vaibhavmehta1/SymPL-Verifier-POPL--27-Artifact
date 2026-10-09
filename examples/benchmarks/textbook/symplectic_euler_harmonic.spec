# Symplectic (semi-implicit) Euler -- H1: harmonic oscillator.
#
# Separable H(q,p) = p^2/2 + q^2/2, so V'(q) = q.  The symplectic Euler step
# kicks the momentum with the OLD position, then drifts the position with the
# NEW momentum:
#       p := p - h V'(q)
#       q := q + h p
# Each substep is a shear (determinant 1), so the map is symplectic and
# 1st-order accurate.  Simplest rung of the benchmark ladder (linear flow).

coords q p
params h

pre  dq ^ dp
post dq ^ dp

program {
    p := p - h*q
    q := q + h*p
}
