# Gauss-Legendre 2-stage (GL-2, 4th-order) -- H3: Henon-Heiles (2 DOF).
#
# H = (px^2+py^2)/2 + (x^2+y^2)/2 + x^2 y - y^3/3,
#       V_x = x + 2 x y,      V_y = y + x^2 - y^2.
# The HEAVIEST benchmark: an 8-dimensional implicit stage solve.  The two GL-2
# stages each carry a full (x,y,px,py) copy, so the unknowns are
#       (Q1x,Q1y,P1x,P1y, Q2x,Q2y,P2x,P2y),
# matched here by  (a1x,a1y,b1x,b1y, a2x,a2y,b2x,b2y).  Stage equations
#       Qi = q + h sum_j a_ij Pj,   Pi = p - h sum_j a_ij V(Qj),
# with the Butcher tableau  a11=1/4, a12=1/4-sqrt(3)/6, a21=1/4+sqrt(3)/6,
# a22=1/4, b1=b2=1/2.  Symplectic, 4th-order; preserves dx ^ dpx + dy ^ dpy.
# The verifier must invert an 8x8 symbolic Jacobian (I - h J_xi) -- this is the
# case most likely to exercise the 2-minute kill.

coords x y px py  Q1x Q1y P1x P1y  Q2x Q2y P2x P2y
params h

pre  dx ^ dpx + dy ^ dpy
post dx ^ dpx + dy ^ dpy

program {
    (Q1x, Q1y, P1x, P1y, Q2x, Q2y, P2x, P2y) := (x, y, px, py, x, y, px, py)

    (Q1x, Q1y, P1x, P1y, Q2x, Q2y, P2x, P2y) := solve (a1x, a1y, b1x, b1y, a2x, a2y, b2x, b2y) -> (
            # Q1 = q + h (a11 P1 + a12 P2)
            (1/4)*b1x + (1/4 - sqrt(3)/6)*b2x,
            (1/4)*b1y + (1/4 - sqrt(3)/6)*b2y,
            # P1 = p - h (a11 V(Q1) + a12 V(Q2))
            -((1/4)*(a1x + 2*a1x*a1y)        + (1/4 - sqrt(3)/6)*(a2x + 2*a2x*a2y)),
            -((1/4)*(a1y + a1x**2 - a1y**2)  + (1/4 - sqrt(3)/6)*(a2y + a2x**2 - a2y**2)),
            # Q2 = q + h (a21 P1 + a22 P2)
            (1/4 + sqrt(3)/6)*b1x + (1/4)*b2x,
            (1/4 + sqrt(3)/6)*b1y + (1/4)*b2y,
            # P2 = p - h (a21 V(Q1) + a22 V(Q2))
            -((1/4 + sqrt(3)/6)*(a1x + 2*a1x*a1y)       + (1/4)*(a2x + 2*a2x*a2y)),
            -((1/4 + sqrt(3)/6)*(a1y + a1x**2 - a1y**2) + (1/4)*(a2y + a2x**2 - a2y**2))
        )

    x  := x  + h*((1/2)*P1x + (1/2)*P2x)
    y  := y  + h*((1/2)*P1y + (1/2)*P2y)
    px := px - h*((1/2)*(Q1x + 2*Q1x*Q1y)       + (1/2)*(Q2x + 2*Q2x*Q2y))
    py := py - h*((1/2)*(Q1y + Q1x**2 - Q1y**2) + (1/2)*(Q2y + Q2x**2 - Q2y**2))
}
