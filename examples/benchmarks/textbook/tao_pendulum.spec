# Tao's explicit symplectic integrator -- H2: pendulum.
#
# Tao (Phys. Rev. E 94, 043303, 2016) doubles phase space to (q,p,x,y), form
# dq^dp+dx^dy, with Hbar = H(q,y) + H(x,p) + om*(1/2)((q-x)^2+(p-y)^2) and takes
# the 2nd-order Strang step
#       phi_{H_A}^{h/2} . phi_{H_B}^{h/2} . phi_{om H_C}^{h}
#                       . phi_{H_B}^{h/2} . phi_{H_A}^{h/2}.
#
# Here H = p^2/2 - cos q, so dH/dq = sin q, dH/dp = p.  (q,y) frozen in phi_{H_A}
# and (x,p) frozen in phi_{H_B}, so each split flow is exact.  Symplectic on
# dq^dp+dx^dy; the copies (x,y) shadow (q,p).

coords q p x y
params h om

pre  dq ^ dp + dx ^ dy
post dq ^ dp + dx ^ dy

program {
    # phi_{H_A}^{h/2} :  p -= (h/2) sin(q) ;  x += (h/2) y
    (p, x) := (p - (h/2)*sin(q),  x + (h/2)*y)

    # phi_{H_B}^{h/2} :  q += (h/2) p ;  y -= (h/2) sin(x)
    (q, y) := (q + (h/2)*p,  y - (h/2)*sin(x))

    # phi_{om H_C}^{h} : exact flow of H_C (rotate (q-x,p-y) by angle 2*om*h)
    (q, p, x, y) := (
        ((q + x) + (q - x)*cos(2*om*h) + (p - y)*sin(2*om*h))/2,
        ((p + y) - (q - x)*sin(2*om*h) + (p - y)*cos(2*om*h))/2,
        ((q + x) - (q - x)*cos(2*om*h) - (p - y)*sin(2*om*h))/2,
        ((p + y) + (q - x)*sin(2*om*h) - (p - y)*cos(2*om*h))/2
    )

    # phi_{H_B}^{h/2}
    (q, y) := (q + (h/2)*p,  y - (h/2)*sin(x))

    # phi_{H_A}^{h/2}
    (p, x) := (p - (h/2)*sin(q),  x + (h/2)*y)
}
