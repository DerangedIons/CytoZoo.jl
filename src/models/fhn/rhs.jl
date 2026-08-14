"""
    _fhn_rhs_impl!(du, u, model, x, t, overrides) -> Nothing

FitzHugh–Nagumo right-hand side, shared by both functor dispatches.

    dv/dt = v(1 - v)(v - a) - s - Iₛₜᵢₘ
    ds/dt = e(b·v - c·s - d)

Every parameter is spatially overridable through `_resolve_parameter`, which dispatches
the `overrides === nothing` case away so the non-spatial call emits exactly the five
struct-field loads and no branch.
"""
function _fhn_rhs_impl!(
    du::AbstractVector{T},
    u::AbstractVector{T},
    model,
    x,
    t,
    overrides,
) where {T}
    v = u[1]
    s = u[2]

    a = T(_resolve_parameter(model.a, overrides, Val(:a), x, t))
    b = T(_resolve_parameter(model.b, overrides, Val(:b), x, t))
    c = T(_resolve_parameter(model.c, overrides, Val(:c), x, t))
    d = T(_resolve_parameter(model.d, overrides, Val(:d), x, t))
    e = T(_resolve_parameter(model.e, overrides, Val(:e), x, t))

    # Sign convention: `Istim` is subtracted, so a NEGATIVE amplitude depolarizes —
    # the negative-inward convention `ToRORd` and the rest of the zoo use.
    Istim = T(model.stim(x, t))

    du[1] = v * (T(1) - v) * (v - a) - s - Istim
    du[2] = e * (b * v - c * s - d)
    return nothing
end
