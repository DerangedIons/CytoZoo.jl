"""
    Constant(value)

Spatial function that returns a fixed value, ignoring position and time.
Isbits-safe alternative to `(x, t) -> value`.
"""
struct Constant{T} <: SpatialFunction
    value::T
end
@inline (f::Constant)(x, t) = f.value

"""
    SpatialStep(dim, threshold, below, above)

Step function in one spatial dimension. Returns `below` when `x[dim] <= threshold`,
`above` otherwise.
"""
struct SpatialStep{T} <: SpatialFunction
    dim::Int
    threshold::T
    below::T
    above::T
end
@inline (f::SpatialStep)(x, t) = x[f.dim] > f.threshold ? f.above : f.below

"""
    SpatialGradient(dim, x_start, x_end, val_start, val_end)

Linear interpolation in one spatial dimension. Clamps to `[val_start, val_end]`
outside the `[x_start, x_end]` range.
"""
struct SpatialGradient{T} <: SpatialFunction
    dim::Int
    x_start::T
    x_end::T
    val_start::T
    val_end::T
end
@inline function (f::SpatialGradient)(x, t)
    frac = clamp((x[f.dim] - f.x_start) / (f.x_end - f.x_start), zero(f.val_start), one(f.val_start))
    return f.val_start + frac * (f.val_end - f.val_start)
end
