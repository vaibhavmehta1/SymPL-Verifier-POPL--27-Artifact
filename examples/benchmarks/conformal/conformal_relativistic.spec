# Conformal symplectic RELATIVISTIC integrator (Relativistic Gradient Descent).
#
# Franca, Jordan & Vidal (NeurIPS 2020) discretize a DISSIPATIVE RELATIVISTIC
# system whose kinetic energy uses the relativistic dispersion relation:
#       H(x,p) = c sqrt(p^2 + m^2 c^2) + f(x),                          (eq. 21)
# with conformal (damped) dynamics  x' = dH/dp,  p' = -dH/dx - gam p.  The
# momentum is normalized by the sqrt factor, so dH/dp = c p / sqrt(p^2+m^2c^2)
# stays bounded even for large p -- a "relativistic" speed limit that stabilizes
# optimization with large gradients.
#
# First-order conformal symplectic integrator (dissipative symplectic Euler,
# eq. 23) for the pendulum potential f(x) = -cos x  (so f'(x) = sin x):
#       p := e^{-gam h} p - h sin x
#       x := x + h c p / sqrt(p^2 + m^2 c^2)
# As c -> infinity this reduces to classical-momentum / heavy-ball descent.
# Being conformal symplectic, one step contracts the canonical 2-form by exactly
# e^{-gam h}:  omega_{k+1} = e^{-gam h} omega_k  (Theorem 3), which we verify as
#       pre = e^{-gam h} (dx ^ dp),   post = dx ^ dp.
#
# Ref: G. Franca, M. I. Jordan & R. Vidal, "Conformal Symplectic and
#      Relativistic Optimization", NeurIPS 2020, eqs. (21), (23), Thm. 3.

coords x p
params h gam c m

pre  exp(-gam*h) * (dx ^ dp)
post dx ^ dp

program {
    p := exp(-gam*h)*p - h*sin(x)
    x := x + h*c*p / sqrt(p**2 + m**2*c**2)
}
