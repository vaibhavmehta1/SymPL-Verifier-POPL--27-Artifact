# Stoermer-Verlet integrator, position form (drift-kick-drift).
#
# The adjoint arrangement of leapfrog: a half position drift, a full momentum
# kick, and a half position drift:
#       q := q + (h/2) p
#       p := p - h V'(q)
#       q := q + (h/2) p
# Separable pendulum H = p^2/2 - cos q, so V'(q) = sin q.  Symplectic and
# 2nd-order.  Same method family as examples/leapfrog.spec, differing only by
# the DKD vs KDK ordering of the symmetric splitting.
#
# Refs: C. Stoermer (1907); L. Verlet, Phys. Rev. 159, 98 (1967);
#       Hairer, Lubich & Wanner, "Geometric Numerical Integration" (2nd ed.,
#         Springer 2006), Secs. I.3 and VI.3 (the Stoermer-Verlet method).

coords q p
params h

pre  dq ^ dp
post dq ^ dp

program {
    q := q + (h/2)*p
    p := p - h*sin(q)
    q := q + (h/2)*p
}
