# Leapfrog / velocity-Verlet (kick-drift-kick) -- H3: Henon-Heiles (2 DOF).
#
# H = (px^2+py^2)/2 + (x^2+y^2)/2 + x^2 y - y^3/3,
#       V_x = x + 2 x y,      V_y = y + x^2 - y^2.
# Half kick of both momenta, full drift of both positions, half kick again.
# The second kick reads the forces at the DRIFTED positions.  Symplectic,
# 2nd-order; preserves dx ^ dpx + dy ^ dpy.

coords x y px py
params h

pre  dx ^ dpx + dy ^ dpy
post dx ^ dpx + dy ^ dpy

program {
    (px, py) := (px - (h/2)*(x + 2*x*y), py - (h/2)*(y + x**2 - y**2))
    (x, y)   := (x + h*px, y + h*py)
    (px, py) := (px - (h/2)*(x + 2*x*y), py - (h/2)*(y + x**2 - y**2))
}
