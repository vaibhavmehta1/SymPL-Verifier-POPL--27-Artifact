# Implicit midpoint rule -- H4: double pendulum (2 DOF, NONSEPARABLE, chaotic).
#
# The planar double pendulum (equal masses m, equal lengths l, gravity g) in
# canonical angle/momentum coordinates (th1, th2, p1, p2) has the nonseparable
# Hamiltonian (Delta = th1 - th2)
#
#   H = [p1^2 + 2 p2^2 - 2 p1 p2 cos Delta] / (2 m l^2 (1 + sin^2 Delta))
#       - m g l (2 cos th1 + cos th2).
#
# Its kinetic term couples BOTH angles and BOTH momenta through a trig
# denominator, so the separable kick-drift integrators do not apply -- only a
# fully implicit method does.  The canonical derivatives (all from this one H,
# so the map is genuinely Hamiltonian and the verifier should VERIFY rather than
# reject):
#   dH/dp1  = (p1 - p2 cosD) / (m l^2 D),     D = 1 + sin^2 Delta
#   dH/dp2  = (2 p2 - p1 cosD) / (m l^2 D)
#   dH/dth1 =  K + 2 m g l sin th1,           dH/dth2 = -K + m g l sin th2
#   K = dH_kin/dDelta
#     = [2 p1 p2 sinD' D - (p1^2+2p2^2-2p1p2 cosD)(2 sinD' cosD')] / (2 m l^2 D^2)
#   (sinD' = sin Delta, cosD' = cos Delta).
#
# Implicit midpoint: the midpoint state (T1,T2,R1,R2) solves the stage equations
# x_mid = x + (h/2) f(x_mid) with f = (dH/dp1, dH/dp2, -dH/dth1, -dH/dth2), then
# x := x + h f(x_mid).  The slope values are stored in augmented registers
# (hp1,hp2,ht1,ht2 and the shared kinetic-angle derivative kterm).  Symplectic,
# 2nd-order; preserves dth1 ^ dp1 + dth2 ^ dp2 for the nonseparable flow.
#
# Refs: Hairer, Lubich & Wanner, "Geometric Numerical Integration" (2006),
#       Sec. VI.3 (midpoint is symplectic for any Hamiltonian); the double
#       pendulum is a standard chaotic 2-DOF test system.

coords th1 th2 p1 p2  T1 T2 R1 R2  hp1 hp2 ht1 ht2 kterm
params h m l g

pre  dth1 ^ dp1 + dth2 ^ dp2
post dth1 ^ dp1 + dth2 ^ dp2

program {
    # initialise the midpoint registers to the current state
    (T1, T2, R1, R2) := (th1, th2, p1, p2)

    # solve the implicit midpoint stage  x_mid = x + (h/2) f(x_mid).
    # unknowns (a1,a2,b1,b2) ~ (T1,T2,R1,R2); kernel G = f/2 with f from H above.
    (T1, T2, R1, R2) := solve (a1, a2, b1, b2) -> (
        (b1 - b2*cos(a1 - a2)) / (2*m*l**2*(1 + sin(a1 - a2)**2)),
        (2*b2 - b1*cos(a1 - a2)) / (2*m*l**2*(1 + sin(a1 - a2)**2)),
        -( ( 2*b1*b2*sin(a1-a2)*(1 + sin(a1-a2)**2)
             - (b1**2 + 2*b2**2 - 2*b1*b2*cos(a1-a2))*(2*sin(a1-a2)*cos(a1-a2)) )
           / (2*m*l**2*(1 + sin(a1-a2)**2)**2) ) / 2
        - m*g*l*sin(a1),
          ( ( 2*b1*b2*sin(a1-a2)*(1 + sin(a1-a2)**2)
             - (b1**2 + 2*b2**2 - 2*b1*b2*cos(a1-a2))*(2*sin(a1-a2)*cos(a1-a2)) )
           / (2*m*l**2*(1 + sin(a1-a2)**2)**2) ) / 2
        - (1/2)*m*g*l*sin(a2)
    )

    # midpoint slopes f(x_mid), evaluated at the solved (T1,T2,R1,R2)
    hp1   := (R1 - R2*cos(T1 - T2)) / (m*l**2*(1 + sin(T1 - T2)**2))
    hp2   := (2*R2 - R1*cos(T1 - T2)) / (m*l**2*(1 + sin(T1 - T2)**2))
    kterm := ( ( 2*R1*R2*sin(T1-T2)*(1 + sin(T1-T2)**2)
                 - (R1**2 + 2*R2**2 - 2*R1*R2*cos(T1-T2))*(2*sin(T1-T2)*cos(T1-T2)) )
               / (2*m*l**2*(1 + sin(T1-T2)**2)**2) )
    ht1   := kterm + 2*m*g*l*sin(T1)
    ht2   := -kterm + m*g*l*sin(T2)

    # full midpoint update  x := x + h f(x_mid)
    th1 := th1 + h*hp1
    th2 := th2 + h*hp2
    p1  := p1 - h*ht1
    p2  := p2 - h*ht2
}
