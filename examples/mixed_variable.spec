# Mixed-variable symplectic (Wisdom-Holman) integrator.
#
# Split  H = H_0(q,p) + eps * V(q),  integrate the integrable part H_0 EXACTLY
# in its own variables and apply the perturbation as a Cartesian kick -- the
# "mixed variable" idea (Keplerian variables for the drift, Cartesian for the
# kick).  Here the integrable part is the harmonic oscillator
#       H_0 = (om/2)(q^2 + p^2),
# whose exact flow is a rotation by angle om*h (standing in for the exact Kepler
# drift of the celestial-mechanics original), and the perturbation is
# V(q) = -cos q, so V'(q) = sin q.  Second-order kick-drift-kick step:
#       p := p - (h/2) eps V'(q)
#       (q,p) := rotation(om*h) (q,p)        # exact integrable drift
#       p := p - (h/2) eps V'(q)
# The rotation has determinant 1 (symplectic) and the kicks are shears, so the
# composition is symplectic.
#
# Refs: J. Wisdom & M. Holman, "Symplectic maps for the n-body problem",
#         Astron. J. 102, 1528-1538 (1991);
#       P. Saha & S. Tremaine, "Symplectic integrators for solar system
#         dynamics", Astron. J. 104, 1633 (1992).

coords q p
params h om eps

pre  dq ^ dp
post dq ^ dp

program {
    p := p - (h/2)*eps*sin(q)
    (q, p) := (q*cos(om*h) + p*sin(om*h), -q*sin(om*h) + p*cos(om*h))
    p := p - (h/2)*eps*sin(q)
}
