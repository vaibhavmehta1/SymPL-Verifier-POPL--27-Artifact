# Wisdom-Holman (mixed-variable) symplectic integrator for the KEPLER problem.
#
# Split  H = H_Kepler(L) + eps * R(l),  advance the Kepler part by its EXACT
# flow (the drift) and kick the perturbation, composed symmetrically (KDK).
#
# The Kepler Hamiltonian in Delaunay action-angle variables (L, l) is
#       H_Kepler = -mu^2 / (2 L^2),        L = sqrt(mu a),  l = mean anomaly,
# which depends only on the action L.  Its EXACT time-h flow is therefore the
# uniform advance of the mean anomaly at the mean motion
#       n = dH_Kepler/dL = mu^2 / L^3 = sqrt(mu / a^3)        (Kepler's 3rd law):
#       l := l + n h,   L unchanged.
# This is the exact Kepler drift -- the action-angle linearization sidesteps
# Kepler's equation E - e sin E = M entirely (no transcendental solve).
#
# The perturbation here is a toy disturbing function R(l) = cos l (a function of
# the mean anomaly), whose exact half-step flow is the action shear
#       L := L + (h/2) eps sin l,   l unchanged.
# Step (kick-drift-kick):
#       L := L + (h/2) eps sin l           # half perturbation kick
#       l := l + (mu^2 / L^3) h            # exact Kepler drift
#       L := L + (h/2) eps sin l           # half perturbation kick
# Drift and kicks are each shears (Jacobian determinant 1), so the step is
# symplectic; the verifier confirms it preserves dL ^ dl exactly.
#
# Honest caveat on "mixed variable": the celebrated WH map keeps the *kick* in
# Cartesian coordinates (where a position-only interaction gives an elementary
# momentum shear) and the *drift* exact-Keplerian -- the conversion between the
# two is Kepler's equation and is NOT elementary.  Working entirely in Delaunay
# variables (as here) keeps every substep elementary while using the genuine
# exact Kepler drift; the perturbation is then a function of the angle rather
# than a Cartesian position.
#
# Refs: J. Wisdom & M. Holman, "Symplectic maps for the n-body problem",
#         Astron. J. 102, 1528-1538 (1991);
#       T. Y. Petrosky, Phys. Lett. A 117, 328 (1986) and B. V. Chirikov &
#         V. V. Vecheslavov (1986) -- the action-angle "Kepler map";
#       A. Morbidelli, "Modern Celestial Mechanics" (2002), Delaunay variables.

coords L l
params h mu eps

pre  dL ^ dl
post dL ^ dl

program {
    L := L + (h/2)*eps*sin(l)
    l := l + (mu**2 / L**3)*h
    L := L + (h/2)*eps*sin(l)
}
