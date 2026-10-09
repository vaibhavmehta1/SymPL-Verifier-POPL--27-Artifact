# Tao's explicit symplectic integrator -- H1: harmonic oscillator.
#
# Tao (Phys. Rev. E 94, 043303, 2016) makes ANY H(q,p) explicitly integrable by
# splitting on a DOUBLED phase space (q,p,x,y) with symplectic form dq^dp+dx^dy:
#       Hbar = H_A + H_B + om*H_C,
#       H_A = H(q,y),  H_B = H(x,p),  H_C = (1/2)((q-x)^2 + (p-y)^2).
# Each piece has an explicit flow; the copies (x,y) shadow (q,p) and the binding
# strength om couples them.  Second-order Strang step:
#       phi_{H_A}^{h/2} . phi_{H_B}^{h/2} . phi_{om H_C}^{h}
#                       . phi_{H_B}^{h/2} . phi_{H_A}^{h/2}.
#
# Here H = p^2/2 + q^2/2, so dH/dq = q, dH/dp = p.  During phi_{H_A} the pair
# (q,y) is frozen, so dH/dq(q,y)=q and dH/dp(q,y)=y are CONSTANT and the flow is
# exact; likewise (x,p) is frozen during phi_{H_B}.  Symplectic on dq^dp+dx^dy.

coords q p x y
params h om

pre  dq ^ dp + dx ^ dy
post dq ^ dp + dx ^ dy

program {
    # phi_{H_A}^{h/2} :  p -= (h/2) dH/dq(q,y) ;  x += (h/2) dH/dp(q,y)
    (p, x) := (p - (h/2)*q,  x + (h/2)*y)

    # phi_{H_B}^{h/2} :  q += (h/2) dH/dp(x,p) ;  y -= (h/2) dH/dq(x,p)
    (q, y) := (q + (h/2)*p,  y - (h/2)*x)

    # phi_{om H_C}^{h} : exact flow of H_C (rotate (q-x,p-y) by angle 2*om*h)
    (q, p, x, y) := (
        ((q + x) + (q - x)*cos(2*om*h) + (p - y)*sin(2*om*h))/2,
        ((p + y) - (q - x)*sin(2*om*h) + (p - y)*cos(2*om*h))/2,
        ((q + x) - (q - x)*cos(2*om*h) - (p - y)*sin(2*om*h))/2,
        ((p + y) + (q - x)*sin(2*om*h) - (p - y)*cos(2*om*h))/2
    )

    # phi_{H_B}^{h/2}
    (q, y) := (q + (h/2)*p,  y - (h/2)*x)

    # phi_{H_A}^{h/2}
    (p, x) := (p - (h/2)*q,  x + (h/2)*y)
}
