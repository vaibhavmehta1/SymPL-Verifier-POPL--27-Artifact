# Implicit midpoint rule -- H2: pendulum.
#
# Separable H = p^2/2 - cos q, V'(q) = sin q.  Midpoint values (Q,P) solve
#       Q = q + (h/2) P,      P = p - (h/2) sin Q,
# then
#       q := q + h P,         p := p - h sin Q.
# Stage kernel G = (P/2, -sin(Q)/2).  Symplectic, 2nd-order.

coords q p Q P
params h

pre  dq ^ dp
post dq ^ dp

program {
    (Q, P) := (q, p)
    (Q, P) := solve (a, b) -> ( b/2, -sin(a)/2 )
    q := q + h*P
    p := p - h*sin(Q)
}
