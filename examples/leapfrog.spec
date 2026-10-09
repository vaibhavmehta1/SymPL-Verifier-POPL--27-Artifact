# Leapfrog integrator (velocity-Verlet / kick-drift-kick form).
#
# Separable Hamiltonian H(q,p) = p^2/2 + V(q); here the pendulum V(q) = -cos q,
# so V'(q) = sin q.  One step is a half momentum kick, a full position drift,
# and a half momentum kick:
#       p := p - (h/2) V'(q)
#       q := q + h p
#       p := p - (h/2) V'(q)
# Each substep is a shear (area-preserving), so the symmetric composition is
# symplectic and 2nd-order accurate.  ("Leapfrog" and "Stoermer-Verlet" are the
# same 2nd-order method; this is its KDK arrangement, examples/stormer_verlet.spec
# is the adjoint DKD arrangement.)
#
# Refs: L. Verlet, Phys. Rev. 159, 98 (1967);
#       Hairer, Lubich & Wanner, "Geometric Numerical Integration" (2nd ed.,
#         Springer 2006), Sec. I.3.1 (leapfrog = Stoermer-Verlet);
#       Leimkuhler & Reich, "Simulating Hamiltonian Dynamics" (CUP 2004), Ch. 2.

coords q p
params h

pre  dq ^ dp
post dq ^ dp

program {
    p := p - (h/2)*sin(q)
    q := q + h*p
    p := p - (h/2)*sin(q)
}
