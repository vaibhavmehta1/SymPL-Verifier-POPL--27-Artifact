# Tao's explicit symplectic integrator -- H4: double pendulum (NONSEPARABLE).
#
# The planar double pendulum (equal masses m, lengths l, gravity g) has the
# nonseparable Hamiltonian (Delta = th1 - th2, D = 1 + sin^2 Delta)
#   H = [p1^2 + 2 p2^2 - 2 p1 p2 cos Delta] / (2 m l^2 D)
#       - m g l (2 cos th1 + cos th2),
# whose kinetic term couples both angles AND both momenta through the trig
# denominator D, so the separable kick-drift methods do not apply.  Implicit
# midpoint does, but its stage solve is CAS-intractable (see
# implicit_midpoint_double_pendulum.spec, a deliberate timeout).
#
# Tao (Phys. Rev. E 94, 043303, 2016) instead makes H EXPLICITLY integrable by
# doubling: shadow the angles (th1,th2) with position copies (X1,X2) and the
# momenta (p1,p2) with momentum copies (Y1,Y2), on the form
#   dth1^dp1 + dth2^dp2 + dX1^dY1 + dX2^dY2.
# Augmented Hbar = H(th,Y) + H(X,p) + om*H_C,
#   H_C = (1/2)((th1-X1)^2+(p1-Y1)^2 + (th2-X2)^2+(p2-Y2)^2).
# In phi_{H_A} the pair (th, Y) is FROZEN, so H's messy derivatives evaluate to
# CONSTANTS and the (p,X) update is an exact linear flow; likewise (X,p) is
# frozen in phi_{H_B}.  Thus every split flow is explicit -- no stage solve --
# even though H is nonseparable.  The needed derivatives:
#   dH/dp1  = (p1 - p2 cosD)/(m l^2 D),   dH/dp2  = (2 p2 - p1 cosD)/(m l^2 D)
#   dH/dth1 =  K + 2 m g l sin th1,       dH/dth2 = -K + m g l sin th2
#   K = dH_kin/dDelta
#     = [2 p1 p2 sinD*D - (p1^2+2p2^2-2p1p2 cosD)(2 sinD cosD)] / (2 m l^2 D^2).
# Second-order Strang step; symplectic on the doubled form.

coords th1 th2 p1 p2  X1 X2 Y1 Y2
params h om m l g

pre  dth1 ^ dp1 + dth2 ^ dp2 + dX1 ^ dY1 + dX2 ^ dY2
post dth1 ^ dp1 + dth2 ^ dp2 + dX1 ^ dY1 + dX2 ^ dY2

program {
    # phi_{H_A}^{h/2} = flow of H(th1,th2,Y1,Y2): (th,Y) frozen, so the momentum
    # copies Y stand in for (p1,p2) in every derivative.  Kick (p1,p2) by
    # dH/dth(th,Y); drift (X1,X2) by dH/dp(th,Y).
    (p1, p2, X1, X2) := (
        p1 - (h/2)*( ( 2*Y1*Y2*sin(th1-th2)*(1 + sin(th1-th2)**2)
                       - (Y1**2 + 2*Y2**2 - 2*Y1*Y2*cos(th1-th2))*(2*sin(th1-th2)*cos(th1-th2)) )
                     / (2*m*l**2*(1 + sin(th1-th2)**2)**2)
                     + 2*m*g*l*sin(th1) ),
        p2 - (h/2)*( -( 2*Y1*Y2*sin(th1-th2)*(1 + sin(th1-th2)**2)
                        - (Y1**2 + 2*Y2**2 - 2*Y1*Y2*cos(th1-th2))*(2*sin(th1-th2)*cos(th1-th2)) )
                      / (2*m*l**2*(1 + sin(th1-th2)**2)**2)
                     + m*g*l*sin(th2) ),
        X1 + (h/2)*(Y1 - Y2*cos(th1-th2)) / (m*l**2*(1 + sin(th1-th2)**2)),
        X2 + (h/2)*(2*Y2 - Y1*cos(th1-th2)) / (m*l**2*(1 + sin(th1-th2)**2))
    )

    # phi_{H_B}^{h/2} = flow of H(X1,X2,p1,p2): (X,p) frozen, so the position
    # copies X stand in for (th1,th2).  Drift (th1,th2) by dH/dp(X,p); kick
    # (Y1,Y2) by dH/dth(X,p).
    (th1, th2, Y1, Y2) := (
        th1 + (h/2)*(p1 - p2*cos(X1-X2)) / (m*l**2*(1 + sin(X1-X2)**2)),
        th2 + (h/2)*(2*p2 - p1*cos(X1-X2)) / (m*l**2*(1 + sin(X1-X2)**2)),
        Y1 - (h/2)*( ( 2*p1*p2*sin(X1-X2)*(1 + sin(X1-X2)**2)
                       - (p1**2 + 2*p2**2 - 2*p1*p2*cos(X1-X2))*(2*sin(X1-X2)*cos(X1-X2)) )
                     / (2*m*l**2*(1 + sin(X1-X2)**2)**2)
                     + 2*m*g*l*sin(X1) ),
        Y2 - (h/2)*( -( 2*p1*p2*sin(X1-X2)*(1 + sin(X1-X2)**2)
                        - (p1**2 + 2*p2**2 - 2*p1*p2*cos(X1-X2))*(2*sin(X1-X2)*cos(X1-X2)) )
                      / (2*m*l**2*(1 + sin(X1-X2)**2)**2)
                     + m*g*l*sin(X2) )
    )

    # phi_{om H_C}^{h} : rotate each shadow pair -- (th1-X1,p1-Y1) and
    # (th2-X2,p2-Y2) -- by angle 2*om*h, leaving the sums fixed.
    (th1, p1, X1, Y1) := (
        ((th1 + X1) + (th1 - X1)*cos(2*om*h) + (p1 - Y1)*sin(2*om*h))/2,
        ((p1 + Y1) - (th1 - X1)*sin(2*om*h) + (p1 - Y1)*cos(2*om*h))/2,
        ((th1 + X1) - (th1 - X1)*cos(2*om*h) - (p1 - Y1)*sin(2*om*h))/2,
        ((p1 + Y1) + (th1 - X1)*sin(2*om*h) - (p1 - Y1)*cos(2*om*h))/2
    )
    (th2, p2, X2, Y2) := (
        ((th2 + X2) + (th2 - X2)*cos(2*om*h) + (p2 - Y2)*sin(2*om*h))/2,
        ((p2 + Y2) - (th2 - X2)*sin(2*om*h) + (p2 - Y2)*cos(2*om*h))/2,
        ((th2 + X2) - (th2 - X2)*cos(2*om*h) - (p2 - Y2)*sin(2*om*h))/2,
        ((p2 + Y2) + (th2 - X2)*sin(2*om*h) - (p2 - Y2)*cos(2*om*h))/2
    )

    # phi_{H_B}^{h/2}
    (th1, th2, Y1, Y2) := (
        th1 + (h/2)*(p1 - p2*cos(X1-X2)) / (m*l**2*(1 + sin(X1-X2)**2)),
        th2 + (h/2)*(2*p2 - p1*cos(X1-X2)) / (m*l**2*(1 + sin(X1-X2)**2)),
        Y1 - (h/2)*( ( 2*p1*p2*sin(X1-X2)*(1 + sin(X1-X2)**2)
                       - (p1**2 + 2*p2**2 - 2*p1*p2*cos(X1-X2))*(2*sin(X1-X2)*cos(X1-X2)) )
                     / (2*m*l**2*(1 + sin(X1-X2)**2)**2)
                     + 2*m*g*l*sin(X1) ),
        Y2 - (h/2)*( -( 2*p1*p2*sin(X1-X2)*(1 + sin(X1-X2)**2)
                        - (p1**2 + 2*p2**2 - 2*p1*p2*cos(X1-X2))*(2*sin(X1-X2)*cos(X1-X2)) )
                      / (2*m*l**2*(1 + sin(X1-X2)**2)**2)
                     + m*g*l*sin(X2) )
    )

    # phi_{H_A}^{h/2}
    (p1, p2, X1, X2) := (
        p1 - (h/2)*( ( 2*Y1*Y2*sin(th1-th2)*(1 + sin(th1-th2)**2)
                       - (Y1**2 + 2*Y2**2 - 2*Y1*Y2*cos(th1-th2))*(2*sin(th1-th2)*cos(th1-th2)) )
                     / (2*m*l**2*(1 + sin(th1-th2)**2)**2)
                     + 2*m*g*l*sin(th1) ),
        p2 - (h/2)*( -( 2*Y1*Y2*sin(th1-th2)*(1 + sin(th1-th2)**2)
                        - (Y1**2 + 2*Y2**2 - 2*Y1*Y2*cos(th1-th2))*(2*sin(th1-th2)*cos(th1-th2)) )
                      / (2*m*l**2*(1 + sin(th1-th2)**2)**2)
                     + m*g*l*sin(th2) ),
        X1 + (h/2)*(Y1 - Y2*cos(th1-th2)) / (m*l**2*(1 + sin(th1-th2)**2)),
        X2 + (h/2)*(2*Y2 - Y1*cos(th1-th2)) / (m*l**2*(1 + sin(th1-th2)**2))
    )
}
