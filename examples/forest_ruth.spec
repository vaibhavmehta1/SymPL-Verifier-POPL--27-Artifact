# Forest-Ruth 4th-order symplectic integrator.
#
# A symmetric triple composition of the 2nd-order Stoermer-Verlet (DKD) leapfrog
# S(s) = [ q += (s/2)p ; p -= s V'(q) ; q += (s/2)p ]  with Yoshida's weights
#       w1 = 1 / (2 - 2^(1/3)),   w0 = -2^(1/3) / (2 - 2^(1/3)),
# satisfying  w0 + 2 w1 = 1.   One step is  S(w1 h) o S(w0 h) o S(w1 h).
# Separable pendulum H = p^2/2 - cos q, so V'(q) = sin q.
#
# Symplecticity holds structurally -- a composition of symplectic leapfrog maps
# is symplectic for ANY weights; these particular w0, w1 are what cancel the
# 3rd-order error terms and raise the order to 4 (the verifier checks the
# symplectic property, not the order).
#
# Refs: E. Forest & R. D. Ruth, "Fourth-order symplectic integration",
#         Physica D 43, 105-117 (1990);
#       H. Yoshida, "Construction of higher order symplectic integrators",
#         Phys. Lett. A 150, 262-268 (1990);
#       R. D. Ruth, IEEE Trans. Nucl. Sci. 30, 2669 (1983).

coords q p w0 w1
params h

pre  dq ^ dp
post dq ^ dp

program {
    w1 := 1/(2 - 2**(1/3))
    w0 := -2**(1/3)/(2 - 2**(1/3))

    # S(w1 h)
    q := q + (w1*h/2)*p
    p := p - (w1*h)*sin(q)
    q := q + (w1*h/2)*p

    # S(w0 h)
    q := q + (w0*h/2)*p
    p := p - (w0*h)*sin(q)
    q := q + (w0*h/2)*p

    # S(w1 h)
    q := q + (w1*h/2)*p
    p := p - (w1*h)*sin(q)
    q := q + (w1*h/2)*p
}
