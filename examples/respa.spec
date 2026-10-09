# RESPA / r-RESPA reversible multiple-time-stepping integrator.
#
# For a force split into a cheap-but-stiff FAST part and an expensive SLOW part,
#       H = p^2/2 + V_fast(q) + V_slow(q),
# RESPA wraps many small fast leapfrog steps inside one large slow kick:
#       slow half-kick (h/2)
#       repeat N times:  fast leapfrog with inner step dt = h/N
#       slow half-kick (h/2)
# so the expensive slow force is evaluated once per big step h while the stiff
# fast force is resolved on the fine scale dt.  Each kick and drift is a shear
# (Jacobian determinant 1), and the symmetric nesting is reversible, so the
# whole map is symplectic -- which the verifier confirms preserves dq ^ dp.
#
# Concrete system:  V_fast = (k/2) q^2  (stiff harmonic, F_fast = -k q),
#                   V_slow = -cos q     (soft pendulum,  F_slow = -sin q).
# Here N = 4 inner steps, dt = h/4  (so dt/2 = h/8); the structure is identical
# for any fixed N.
#
# Ref: M. Tuckerman, B. J. Berne & G. J. Martyna, "Reversible multiple time
#      scale molecular dynamics", J. Chem. Phys. 97, 1990-2001 (1992).

coords q p
params h k

pre  dq ^ dp
post dq ^ dp

program {
    p := p - (h/2)*sin(q)          # slow half-kick  (V_slow' = sin q)

    while^4 {                      # N = 4 inner fast leapfrog steps, dt = h/4
        p := p - (h/8)*k*q         # fast half-kick  ((dt/2) V_fast' = (h/8) k q)
        q := q + (h/4)*p           # drift           (dt = h/4)
        p := p - (h/8)*k*q         # fast half-kick
    }

    p := p - (h/2)*sin(q)          # slow half-kick
}
