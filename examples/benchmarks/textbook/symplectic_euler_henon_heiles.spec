# Symplectic (semi-implicit) Euler -- H3: Henon-Heiles (2 DOF, coupled).
#
# H(x,y,px,py) = (px^2 + py^2)/2 + (x^2 + y^2)/2 + x^2 y - y^3/3,
# the classic non-integrable 2-DOF system, with forces
#       V_x = dV/dx = x + 2 x y,      V_y = dV/dy = y + x^2 - y^2.
# Symplectic Euler in both degrees of freedom: kick both momenta with the old
# positions, then drift both positions with the new momenta.  Preserves the
# canonical 2-form dx ^ dpx + dy ^ dpy.  Top rung: coupled, multi-DOF.
#
# Ref: M. Henon & C. Heiles, Astron. J. 69, 73 (1964).

coords x y px py
params h

pre  dx ^ dpx + dy ^ dpy
post dx ^ dpx + dy ^ dpy

program {
    (px, py) := (px - h*(x + 2*x*y), py - h*(y + x**2 - y**2))
    (x, y)   := (x + h*px, y + h*py)
}
