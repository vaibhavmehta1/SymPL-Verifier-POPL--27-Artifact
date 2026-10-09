# Exact sector-bend map from BeamTracking.jl  (src/kernels/bend_kick.jl,
# function exact_bend!).  Tracks a particle through a curved, finite-length
# combined-function dipole by the exact solution of the bend Hamiltonian, in
# canonical accelerator coordinates (x, px, y, py, z, pz).
#
# This is the middle step of the bend-kick-bend driver bkb_multipole!:
#       multipole_kick! ;  exact_bend! ;  multipole_kick!
# (the kick is in examples/multipole_kick.spec).  Being an exact Hamiltonian
# flow, exact_bend! must preserve  omega = dx^dpx + dy^dpy + dz^dpz.
#
# Simplifications applied :
#   * all particles alive  -> every vifelse(alive, new, old) becomes `new`,
#     and the STATE_LOST guards (good_momenta, good_arg, good_cond, good_gp)
#     are dropped, so the masked sqrt/asin arguments are taken directly.
#   * sgn = sign(L): L is the element length (a fixed parameter), so sgn is the
#     compile-time constant +1 (assume L > 0).
#   * theta = g*L at the call site, substituted inline.
#   * The numerical-stability branch
#         xi = vifelse(!good_gp | pos_cplus, xi1, xi2)
#     selects between two algebraically equal formulas (xi1 == xi2 since
#     nsq^2 - cpl^2 = gp*alpha).  We verify xi1,
#     which stays regular as gp -> 0.
#
# Needs asin, atan2 (2-arg) and sincu(t)=sin(t)/t, added to parse.py.

coords x px y py z pz  relp pt arg phi1 gp hh cpl spl th sincth coscth alpha cond nsq xi Lcv thp Lp
params g L Kn0 tilde_m beta_0

pre  dx ^ dpx + dy ^ dpy + dz ^ dpz
post dx ^ dpx + dy ^ dpy + dz ^ dpz

program {
    # bend-angle constants (theta = g*L)
    th     := g*L
    sincth := sincu(th)
    coscth := sincu(th/2)*sincu(th/2)/2

    # entrance kinematics
    relp := 1 + pz
    pt   := sqrt(relp**2 - py**2)
    arg  := px/pt
    phi1 := th + asin(arg)
    gp   := Kn0/pt
    hh   := 1 + g*x
    cpl  := cos(phi1)
    spl  := sin(phi1)

    # intersection of the trajectory arc with the exit face
    alpha := 2*hh*spl*L*sincth - gp*(hh*L*sincth)**2
    cond  := cpl**2 + gp*alpha
    nsq   := sqrt(cond)
    xi    := alpha/(nsq + cpl)

    # transformed bend angle and arc length.
    # The kernel uses atan2(xi, -Lcv); here -Lcv = L*sincu(th) + x*sin(th) > 0
    # for L > 0, so atan2(xi, -Lcv) == arctan(xi/(-Lcv)) on the principal branch
    #.  We use arctan so the CAS can rectify the angle.
    Lcv := -(L*sincth + x*sin(th))
    thp := 2*(phi1 - arctan(xi/(-Lcv)))
    Lp  := sqrt(Lcv**2 + xi**2)/sincu(thp/2)

    # exit coordinates (py, pz are invariants of the bend)
    (x, px, y, z) := (
        x*cos(th) - L*L*g*coscth + xi,
        pt*sin(phi1 - thp),
        y + py*Lp/pt,
        z - relp*Lp/pt + L*relp/sqrt(tilde_m**2 + relp**2)/beta_0
    )
}
