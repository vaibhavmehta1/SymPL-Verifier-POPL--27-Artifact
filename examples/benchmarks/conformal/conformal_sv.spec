# Conformal symplectic Stoermer-Verlet (dissipative leapfrog).
#
# Conformal Hamiltonian system (Franca, Jordan & Vidal, NeurIPS 2020, eq. 3):
#       x' = dH/dp,    p' = -dH/dx - gam p,        gam >= 0 the damping.
# Its flow does not preserve the symplectic form but CONTRACTS it exponentially,
#       omega_t = e^{-gam t} omega_0                                     (eq. 8),
# so one step of any conformal symplectic integrator obeys omega_{k+1} =
# e^{-gam h} omega_k.  We verify exactly that one step pulls the canonical form
# back to e^{-gam h} times itself, i.e.  pre = e^{-gam h} (dx ^ dp).
#
# This is the 2nd-order "dissipative leapfrog" of eq. (11), for the separable
# pendulum H = p^2/2 - cos x  (so dH/dx = sin x, dH/dp = p):
#       xh := x + (h/2) e^{-gam h/2} p          (11a)
#       pt := e^{-gam h/2} p - h sin(xh)        (11b)
#       x  := xh + (h/2) pt                     (11c)
#       p  := e^{-gam h/2} pt                   (11d)
# Setting gam = 0 recovers the ordinary symplectic Stoermer-Verlet.  The two
# factors e^{-gam h/2} applied to the momentum compose to the e^{-gam h}
# contraction of the 2-form.  (xh, pt are augmented-state registers holding the
# leapfrog half-step quantities.)
#
#     Conformal Symplectic and
#      Relativistic Optimization", NeurIPS 2020, eqs. (3), (8), (11).

coords x p xh pt
params h gam

pre  exp(-gam*h) * (dx ^ dp)
post dx ^ dp

program {
    xh := x + (h/2)*exp(-gam*h/2)*p
    pt := exp(-gam*h/2)*p - h*sin(xh)
    x  := xh + (h/2)*pt
    p  := exp(-gam*h/2)*pt
}
