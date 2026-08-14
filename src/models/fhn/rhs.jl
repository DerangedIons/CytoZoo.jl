"""
    _fhn_rhs_impl!(du, u, model, x, t, overrides) -> Nothing

FitzHugh–Nagumo right-hand side, shared by both functor dispatches.

    dv/dt = v(1 - v)(v - a) - s - Iₛₜᵢₘ
    ds/dt = e(b·v - c·s - d)

Every parameter is spatially overridable. The `F !== Nothing` guards are resolved at
compile time, so the non-spatial call (`overrides === nothing`) emits exactly the
five struct-field loads and no branch.
"""
function _fhn_rhs_impl!(
    du::AbstractVector{T},
    u::AbstractVector{T},
    model,
    x,
    t,
    overrides::F,
) where {T, F}
    v = u[1]
    s = u[2]

    a = if F !== Nothing && hasproperty(overrides, :a)
        T(_resolve_spatial(overrides.a, x, t))
    else
        T(model.a)
    end
    b = if F !== Nothing && hasproperty(overrides, :b)
        T(_resolve_spatial(overrides.b, x, t))
    else
        T(model.b)
    end
    c = if F !== Nothing && hasproperty(overrides, :c)
        T(_resolve_spatial(overrides.c, x, t))
    else
        T(model.c)
    end
    d = if F !== Nothing && hasproperty(overrides, :d)
        T(_resolve_spatial(overrides.d, x, t))
    else
        T(model.d)
    end
    e = if F !== Nothing && hasproperty(overrides, :e)
        T(_resolve_spatial(overrides.e, x, t))
    else
        T(model.e)
    end

    # Sign convention: `Istim` is subtracted, so a NEGATIVE amplitude depolarizes —
    # the negative-inward convention `ToRORd` and the rest of the zoo use.
    Istim = T(model.stim(x, t))

    du[1] = v * (T(1) - v) * (v - a) - s - Istim
    du[2] = e * (b * v - c * s - d)
    return nothing
end
