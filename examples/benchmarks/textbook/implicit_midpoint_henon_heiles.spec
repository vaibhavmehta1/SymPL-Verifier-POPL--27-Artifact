# Implicit midpoint rule -- H3: Henon-Heiles (2 DOF, coupled implicit solve).
#
# H = (px^2+py^2)/2 + (x^2+y^2)/2 + x^2 y - y^3/3,
#       V_x = x + 2 x y,      V_y = y + x^2 - y^2.
# The four midpoint values (X,Y,PX,PY) solve the coupled stage equations
#       X  = x  + (h/2) PX,        Y  = y  + (h/2) PY,
#       PX = px - (h/2) V_x(X,Y),  PY = py - (h/2) V_y(X,Y),
# then the state advances with the midpoint slopes
#       x  := x  + h PX,           y  := y  + h PY,
#       px := px - h V_x(X,Y),     py := py - h V_y(X,Y).
# A genuine 4-dimensional implicit solve.  Symplectic, 2nd-order; preserves
# dx ^ dpx + dy ^ dpy.

coords x y px py X Y PX PY
params h

pre  dx ^ dpx + dy ^ dpy
post dx ^ dpx + dy ^ dpy

program {
    (X, Y, PX, PY) := (x, y, px, py)
    (X, Y, PX, PY) := solve (ax, ay, bx, by) -> (
        bx/2,
        by/2,
        -(ax + 2*ax*ay)/2,
        -(ay + ax**2 - ay**2)/2
    )
    x  := x  + h*PX
    y  := y  + h*PY
    px := px - h*(X + 2*X*Y)
    py := py - h*(Y + X**2 - Y**2)
}
