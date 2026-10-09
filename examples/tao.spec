# Tao's explicit symplectic integrator for NONSEPARABLE Hamiltonians.
#
# A nonseparable H(q,p) is made integrable-by-splitting on a DOUBLED phase space
# (q,p,x,y) with the augmented Hamiltonian
#       Hbar = H_A + H_B + om * H_C,
#       H_A = H(q,y),  H_B = H(x,p),  H_C = (1/2)((q-x)^2 + (p-y)^2),
# whose three pieces each have an explicit flow.  The copies start equal,
# (x,y) = (q,p), and the coupling strength om binds them.  Second-order step
# (Strang composition):
#       phi_{H_A}^{h/2} . phi_{H_B}^{h/2} . phi_{om H_C}^{h}
#                       . phi_{H_B}^{h/2} . phi_{H_A}^{h/2}.
#
# Concrete nonseparable example from the paper:  H(a,b) = (1/2)(a^2+1)(b^2+1),
# so  dH/da = a(b^2+1),  dH/db = (a^2+1)b.   Flows (each exact, hence symplectic
# on dq^dp + dx^dy):
#   phi_{H_A}^t :  p -= t*q*(y^2+1) ;        x += t*(q^2+1)*y      (q,y fixed)
#   phi_{H_B}^t :  q += t*(x^2+1)*p ;        y -= t*x*(p^2+1)      (x,p fixed)
#   phi_{om H_C}^t : linear rotation of (q-x, p-y) by angle 2*om*t,
#                    leaving (q+x, p+y) fixed.
#
# Ref: Molei Tao, "Explicit symplectic approximation of nonseparable
#      Hamiltonians: algorithm and long time performance",
#      Phys. Rev. E 94, 043303 (2016).

coords q p x y
params h om

pre  dq ^ dp + dx ^ dy
post dq ^ dp + dx ^ dy

program {
    # phi_{H_A}^{h/2}
    (p, x) := (p - (h/2)*q*(y**2 + 1),  x + (h/2)*(q**2 + 1)*y)

    # phi_{H_B}^{h/2}
    (q, y) := (q + (h/2)*(x**2 + 1)*p,  y - (h/2)*x*(p**2 + 1))

    # phi_{om H_C}^{h} : exact flow of H_C (angle 2*om*h)
    (q, p, x, y) := (
        ((q + x) + (q - x)*cos(2*om*h) + (p - y)*sin(2*om*h))/2,
        ((p + y) - (q - x)*sin(2*om*h) + (p - y)*cos(2*om*h))/2,
        ((q + x) - (q - x)*cos(2*om*h) - (p - y)*sin(2*om*h))/2,
        ((p + y) + (q - x)*sin(2*om*h) - (p - y)*cos(2*om*h))/2
    )

    # phi_{H_B}^{h/2}
    (q, y) := (q + (h/2)*(x**2 + 1)*p,  y - (h/2)*x*(p**2 + 1))

    # phi_{H_A}^{h/2}
    (p, x) := (p - (h/2)*q*(y**2 + 1),  x + (h/2)*(q**2 + 1)*y)
}
