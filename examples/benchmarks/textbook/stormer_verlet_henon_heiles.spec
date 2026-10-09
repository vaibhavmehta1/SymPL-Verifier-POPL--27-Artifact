# Stoermer-Verlet (drift-kick-drift) -- H3: Henon-Heiles (2 DOF).
#
# H = (px^2+py^2)/2 + (x^2+y^2)/2 + x^2 y - y^3/3,
#       V_x = x + 2 x y,      V_y = y + x^2 - y^2.
# Half drift of both positions, full kick of both momenta (forces at the
# half-drifted positions), half drift again.  Symplectic, 2nd-order; preserves
# dx ^ dpx + dy ^ dpy.

coords x y px py
params h

pre  dx ^ dpx + dy ^ dpy
post dx ^ dpx + dy ^ dpy

program {
    (x, y)   := (x + (h/2)*px, y + (h/2)*py)
    (px, py) := (px - h*(x + 2*x*y), py - h*(y + x**2 - y**2))
    (x, y)   := (x + (h/2)*px, y + (h/2)*py)
}
