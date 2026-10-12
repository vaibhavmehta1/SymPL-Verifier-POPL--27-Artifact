# Gauss-Legendre 2-stage (GL-2, 4th-order) -- H4: double pendulum.
#
# Same equal-mass, equal-length Hamiltonian as implicit_midpoint_double_pendulum.spec:
# H = (p1^2 + 2*p2^2 - 2*p1*p2*cos(th1-th2))
#       / (2*m*l^2*(1 + sin(th1-th2)^2)) - m*g*l*(2*cos(th1) + cos(th2)).
# The canonical vector field is f = (H_p1, H_p2, -H_th1, -H_th2).
#
# Y_i = x + h*sum_j A_ij*f(Y_j), x_new = x + h*(f(Y_1)+f(Y_2))/2,
# A = [[1/4, 1/4-sqrt(3)/6], [1/4+sqrt(3)/6, 1/4]], b = [1/2, 1/2].
# Each stage holds (th1, th2, p1, p2). All eight stage components are solved
# together, with registers initialized to x because solve uses target anchors.
# This nonseparable system requires an 8x8 symbolic implicit Jacobian;
# the standard verifier may exhaust the benchmark time budget.

coords th1 th2 p1 p2  Q11 Q12 P11 P12 Q21 Q22 P21 P22
params h m l g

pre  dth1 ^ dp1 + dth2 ^ dp2
post dth1 ^ dp1 + dth2 ^ dp2

program {
    (Q11, Q12, P11, P12, Q21, Q22, P21, P22) := (th1, th2, p1, p2, th1, th2, p1, p2)

    (Q11, Q12, P11, P12, Q21, Q22, P21, P22) := solve (a11, a12, b11, b12, a21, a22, b21, b22) -> (
        (1/4)*((b11 - b12*cos(a11 - a12)) / (m*l**2*(1 + sin(a11 - a12)**2)))
        + (1/4 - sqrt(3)/6)*((b21 - b22*cos(a21 - a22)) / (m*l**2*(1 + sin(a21 - a22)**2))),

        (1/4)*((2*b12 - b11*cos(a11 - a12)) / (m*l**2*(1 + sin(a11 - a12)**2)))
        + (1/4 - sqrt(3)/6)*((2*b22 - b21*cos(a21 - a22)) / (m*l**2*(1 + sin(a21 - a22)**2))),

        (1/4)*(-((2*b11*b12*sin(a11 - a12)*(1 + sin(a11 - a12)**2)
                - (b11**2 + 2*b12**2 - 2*b11*b12*cos(a11 - a12))*2*sin(a11 - a12)*cos(a11 - a12))
               / (2*m*l**2*(1 + sin(a11 - a12)**2)**2)) - 2*m*g*l*sin(a11))
        + (1/4 - sqrt(3)/6)*(-((2*b21*b22*sin(a21 - a22)*(1 + sin(a21 - a22)**2)
                - (b21**2 + 2*b22**2 - 2*b21*b22*cos(a21 - a22))*2*sin(a21 - a22)*cos(a21 - a22))
               / (2*m*l**2*(1 + sin(a21 - a22)**2)**2)) - 2*m*g*l*sin(a21)),

        (1/4)*(((2*b11*b12*sin(a11 - a12)*(1 + sin(a11 - a12)**2)
                - (b11**2 + 2*b12**2 - 2*b11*b12*cos(a11 - a12))*2*sin(a11 - a12)*cos(a11 - a12))
               / (2*m*l**2*(1 + sin(a11 - a12)**2)**2)) - m*g*l*sin(a12))
        + (1/4 - sqrt(3)/6)*(((2*b21*b22*sin(a21 - a22)*(1 + sin(a21 - a22)**2)
                - (b21**2 + 2*b22**2 - 2*b21*b22*cos(a21 - a22))*2*sin(a21 - a22)*cos(a21 - a22))
               / (2*m*l**2*(1 + sin(a21 - a22)**2)**2)) - m*g*l*sin(a22)),

        (1/4 + sqrt(3)/6)*((b11 - b12*cos(a11 - a12)) / (m*l**2*(1 + sin(a11 - a12)**2)))
        + (1/4)*((b21 - b22*cos(a21 - a22)) / (m*l**2*(1 + sin(a21 - a22)**2))),

        (1/4 + sqrt(3)/6)*((2*b12 - b11*cos(a11 - a12)) / (m*l**2*(1 + sin(a11 - a12)**2)))
        + (1/4)*((2*b22 - b21*cos(a21 - a22)) / (m*l**2*(1 + sin(a21 - a22)**2))),

        (1/4 + sqrt(3)/6)*(-((2*b11*b12*sin(a11 - a12)*(1 + sin(a11 - a12)**2)
                - (b11**2 + 2*b12**2 - 2*b11*b12*cos(a11 - a12))*2*sin(a11 - a12)*cos(a11 - a12))
               / (2*m*l**2*(1 + sin(a11 - a12)**2)**2)) - 2*m*g*l*sin(a11))
        + (1/4)*(-((2*b21*b22*sin(a21 - a22)*(1 + sin(a21 - a22)**2)
                - (b21**2 + 2*b22**2 - 2*b21*b22*cos(a21 - a22))*2*sin(a21 - a22)*cos(a21 - a22))
               / (2*m*l**2*(1 + sin(a21 - a22)**2)**2)) - 2*m*g*l*sin(a21)),

        (1/4 + sqrt(3)/6)*(((2*b11*b12*sin(a11 - a12)*(1 + sin(a11 - a12)**2)
                - (b11**2 + 2*b12**2 - 2*b11*b12*cos(a11 - a12))*2*sin(a11 - a12)*cos(a11 - a12))
               / (2*m*l**2*(1 + sin(a11 - a12)**2)**2)) - m*g*l*sin(a12))
        + (1/4)*(((2*b21*b22*sin(a21 - a22)*(1 + sin(a21 - a22)**2)
                - (b21**2 + 2*b22**2 - 2*b21*b22*cos(a21 - a22))*2*sin(a21 - a22)*cos(a21 - a22))
               / (2*m*l**2*(1 + sin(a21 - a22)**2)**2)) - m*g*l*sin(a22))
    )

    # Simultaneous physical update using only the solved stage values.
    (th1, th2, p1, p2) := (
        th1 + (h/2)*(((P11 - P12*cos(Q11 - Q12)) / (m*l**2*(1 + sin(Q11 - Q12)**2)))
                       + ((P21 - P22*cos(Q21 - Q22)) / (m*l**2*(1 + sin(Q21 - Q22)**2)))),

        th2 + (h/2)*(((2*P12 - P11*cos(Q11 - Q12)) / (m*l**2*(1 + sin(Q11 - Q12)**2)))
                       + ((2*P22 - P21*cos(Q21 - Q22)) / (m*l**2*(1 + sin(Q21 - Q22)**2)))),

        p1 + (h/2)*((-((2*P11*P12*sin(Q11 - Q12)*(1 + sin(Q11 - Q12)**2)
                - (P11**2 + 2*P12**2 - 2*P11*P12*cos(Q11 - Q12))*2*sin(Q11 - Q12)*cos(Q11 - Q12))
               / (2*m*l**2*(1 + sin(Q11 - Q12)**2)**2)) - 2*m*g*l*sin(Q11))
                       + (-((2*P21*P22*sin(Q21 - Q22)*(1 + sin(Q21 - Q22)**2)
                - (P21**2 + 2*P22**2 - 2*P21*P22*cos(Q21 - Q22))*2*sin(Q21 - Q22)*cos(Q21 - Q22))
               / (2*m*l**2*(1 + sin(Q21 - Q22)**2)**2)) - 2*m*g*l*sin(Q21))),

        p2 + (h/2)*((((2*P11*P12*sin(Q11 - Q12)*(1 + sin(Q11 - Q12)**2)
                - (P11**2 + 2*P12**2 - 2*P11*P12*cos(Q11 - Q12))*2*sin(Q11 - Q12)*cos(Q11 - Q12))
               / (2*m*l**2*(1 + sin(Q11 - Q12)**2)**2)) - m*g*l*sin(Q12))
                       + (((2*P21*P22*sin(Q21 - Q22)*(1 + sin(Q21 - Q22)**2)
                - (P21**2 + 2*P22**2 - 2*P21*P22*cos(Q21 - Q22))*2*sin(Q21 - Q22)*cos(Q21 - Q22))
               / (2*m*l**2*(1 + sin(Q21 - Q22)**2)**2)) - m*g*l*sin(Q22)))
    )
}

