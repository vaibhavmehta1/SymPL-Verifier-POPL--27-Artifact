# Thin multipole kick from BeamTracking.jl  (src/kernels/multipole.jl).
#
# multipole_kick! applies the transverse magnetic field of a multipole magnet
# as a momentum kick:
#
#       px -= By(x, y),     py += Bx(x, y).
#
# The field is the Horner / complex-power expansion built by normalized_field:
#
#       By + i*Bx = sum_n (knl_n + i*ksl_n) * (x + i*y)^(n-1) / (n-1)!
#
# where n is the multipole order (n = 1 dipole, 2 quad, 3 sext, ...).  The
# bend-kick-bend driver calls it with `excluding = 1`, i.e. the dipole term
# (n = 1) is handled by the exact bend, so the kick carries n >= 2.  Here we
# take a concrete finite multipole content: normal/skew quadrupole (kn1,ks1),
# sextupole (kn2,ks2) and octupole (kn3,ks3).
#
# Simplifications vs. the Julia kernel (per the request):
#   * no radiation, no spin, no entrance/exit rotation  -> those branches drop;
#   * all particles assumed alive -
#every `vifelse(alive, new, old)` collapses
#     to `new`, and the STATE_LOST momentum checks disappear.
#
# Since x, y are untouched and By, Bx come from a scalar potential
# (dBy/dy + dBx/dx = 0), the kick preserves the full canonical form
#   omega = dx ^ dpx + dy ^ dpy + dz ^ dpz   (z, pz ride along unchanged).

coords x px y py z pz
params kn1 ks1 kn2 ks2 kn3 ks3

pre  dx ^ dpx + dy ^ dpy + dz ^ dpz
post dx ^ dpx + dy ^ dpy + dz ^ dpz

program {
    # By(x,y) and Bx(x,y) from the (x + i y)^(n-1)/(n-1)! expansion, n = 2,3,4
    px := px - ( kn1*x - ks1*y
               + (kn2*(x**2 - y**2))/2 - ks2*x*y
               + (kn3*(x**3 - 3*x*y**2))/6 - (ks3*(3*x**2*y - y**3))/6 )

    py := py + ( ks1*x + kn1*y
               + (ks2*(x**2 - y**2))/2 + kn2*x*y
               + (ks3*(x**3 - 3*x*y**2))/6 + (kn3*(3*x**2*y - y**3))/6 )
}
