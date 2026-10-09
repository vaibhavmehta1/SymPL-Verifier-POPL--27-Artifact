# Tao's explicit symplectic integrator -- H3: Henon-Heiles (2 DOF).
#
# Tao (Phys. Rev. E 94, 043303, 2016) doubles each DOF.  The original positions
# (x,y) with momenta (px,py) gain shadow copies: position copies (X,Y) and
# momentum copies (PX,PY), on the form dx^dpx+dy^dpy+dX^dPX+dY^dPY.  Augmented
# Hbar = H(x,y,PX,PY) + H(X,Y,px,py) + om*H_C with
#       H_C = (1/2)((x-X)^2+(px-PX)^2 + (y-Y)^2+(py-PY)^2).
# Second-order Strang step (H_A then H_B then H_C then H_B then H_A).
#
# H = (px^2+py^2)/2 + (x^2+y^2)/2 + x^2 y - y^3/3, so
#       V_x = dH/dx = x + 2 x y,   V_y = dH/dy = y + x^2 - y^2,
#       dH/dpx = px,               dH/dpy = py.
# phi_{H_A} freezes (x,y,PX,PY); phi_{H_B} freezes (X,Y,px,py); each split flow
# is exact.  om*H_C rotates each shadow pair.  Symplectic on the doubled form.

coords x y px py  X Y PX PY
params h om

pre  dx ^ dpx + dy ^ dpy + dX ^ dPX + dY ^ dPY
post dx ^ dpx + dy ^ dpy + dX ^ dPX + dY ^ dPY

program {
    # phi_{H_A}^{h/2} = flow of H(x,y,PX,PY): kick (px,py) by V(x,y), drift (X,Y) by (PX,PY)
    (px, py, X, Y) := (
        px - (h/2)*(x + 2*x*y),
        py - (h/2)*(y + x**2 - y**2),
        X + (h/2)*PX,
        Y + (h/2)*PY
    )

    # phi_{H_B}^{h/2} = flow of H(X,Y,px,py): drift (x,y) by (px,py), kick (PX,PY) by V(X,Y)
    (x, y, PX, PY) := (
        x + (h/2)*px,
        y + (h/2)*py,
        PX - (h/2)*(X + 2*X*Y),
        PY - (h/2)*(Y + X**2 - Y**2)
    )

    # phi_{om H_C}^{h} : rotate each shadow pair -- (x-X,px-PX) and (y-Y,py-PY) -- by 2*om*h
    (x, px, X, PX) := (
        ((x + X) + (x - X)*cos(2*om*h) + (px - PX)*sin(2*om*h))/2,
        ((px + PX) - (x - X)*sin(2*om*h) + (px - PX)*cos(2*om*h))/2,
        ((x + X) - (x - X)*cos(2*om*h) - (px - PX)*sin(2*om*h))/2,
        ((px + PX) + (x - X)*sin(2*om*h) - (px - PX)*cos(2*om*h))/2
    )
    (y, py, Y, PY) := (
        ((y + Y) + (y - Y)*cos(2*om*h) + (py - PY)*sin(2*om*h))/2,
        ((py + PY) - (y - Y)*sin(2*om*h) + (py - PY)*cos(2*om*h))/2,
        ((y + Y) - (y - Y)*cos(2*om*h) - (py - PY)*sin(2*om*h))/2,
        ((py + PY) + (y - Y)*sin(2*om*h) - (py - PY)*cos(2*om*h))/2
    )

    # phi_{H_B}^{h/2}
    (x, y, PX, PY) := (
        x + (h/2)*px,
        y + (h/2)*py,
        PX - (h/2)*(X + 2*X*Y),
        PY - (h/2)*(Y + X**2 - Y**2)
    )

    # phi_{H_A}^{h/2}
    (px, py, X, Y) := (
        px - (h/2)*(x + 2*x*y),
        py - (h/2)*(y + x**2 - y**2),
        X + (h/2)*PX,
        Y + (h/2)*PY
    )
}
