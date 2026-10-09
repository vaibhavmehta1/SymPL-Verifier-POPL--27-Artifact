# Boris algorithm (full 3D) for a charged particle in an electromagnetic field
#
#     x' = v,    v' = (q/m) (E + v x B).
#
# The Boris pusher is the standard particle-in-cell integrator.  It is NOT
# canonically symplectic, but it exactly preserves the phase-space volume form
#
#     dx1 ^ dx2 ^ dx3 ^ dv1 ^ dv2 ^ dv3.
#
# The half-kicks are velocity shears; the magnetic substep is the Rodrigues /
# Cayley form of an exact rotation (orthogonal, determinant 1), with
# s = 2 t / (1 + |t|^2) and t = (q dt / 2m) B.  The program computes every
# intermediate quantity (a, v^-, t, rho, s, v', v^+) in auxiliary registers.

coords x1 x2 x3  v1 v2 v3  a  vm1 vm2 vm3  t1 t2 t3  rho  s1 s2 s3  vp1 vp2 vp3  vu1 vu2 vu3
params dt q m  E1 E2 E3  B1 B2 B3

pre  dx1 ^ dx2 ^ dx3 ^ dv1 ^ dv2 ^ dv3
post dx1 ^ dx2 ^ dx3 ^ dv1 ^ dv2 ^ dv3

program {
    a := dt*q/(2*m)

    # first half electric kick:  v^- = v + a E
    (vm1, vm2, vm3) := (v1 + a*E1, v2 + a*E2, v3 + a*E3)

    # rotation vector t = a B  and  rho = 1 + |t|^2
    (t1, t2, t3) := (a*B1, a*B2, a*B3)
    rho := 1 + t1**2 + t2**2 + t3**2
    (s1, s2, s3) := (2*t1/rho, 2*t2/rho, 2*t3/rho)

    # v' = v^- + v^- x t
    (vp1, vp2, vp3) := (vm1 + vm2*t3 - vm3*t2,
                        vm2 + vm3*t1 - vm1*t3,
                        vm3 + vm1*t2 - vm2*t1)

    # v^+ = v^- + v' x s   (completes the exact rotation)
    (vu1, vu2, vu3) := (vm1 + vp2*s3 - vp3*s2,
                        vm2 + vp3*s1 - vp1*s3,
                        vm3 + vp1*s2 - vp2*s1)

    # second half electric kick:  v = v^+ + a E
    (v1, v2, v3) := (vu1 + a*E1, vu2 + a*E2, vu3 + a*E3)

    # position drift with the updated velocity
    (x1, x2, x3) := (x1 + dt*v1, x2 + dt*v2, x3 + dt*v3)
}
