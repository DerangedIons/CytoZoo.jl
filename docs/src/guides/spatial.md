# Spatial Heterogeneity

Tissue is not uniform. Ion channel conductances vary across the ventricular wall, an
ischaemic region has different kinetics from healthy tissue, and a drug block may only reach
part of the domain. CytoZoo expresses all of this by letting individual parameters become
functions of position and time, without duplicating the model or its parameter vector per
cell.

## The Mechanism: `p`

The DiffEq right-hand side signature `f(du, u, p, t)` carries a parameter slot `p` that
CytoZoo leaves free for exactly this purpose. Pass a [`SpatialContext`](@ref) and the model
takes its spatial path:

```@example spatial
using CytoZoo

model = ToRORd()
u  = default_initial_state(model)
du = similar(u)

p = SpatialContext([1.2, 0.5, 1.8], (
    IKr_Multiplier = (x, t) -> x[1] > 1.5 ? 0.5 : 1.0,
    isHypoxic      = (x, t) -> x[3] > 2.0 ? 1.0 : 0.0,
))

model(du, u, p, 0.0)
du[1]
```

A `SpatialContext` bundles two things:

- `x` — this cell's position, whatever the caller's coordinate convention is.
- `overrides` — a `NamedTuple` whose keys are **parameter names** and whose values say how
  that parameter varies.

Any parameter not named in `overrides` keeps its value from the model's parameter vector.
Names must match the model's parameter names:

```@example spatial
parameter_index(model, :IKr_Multiplier)
```

## Three Kinds of Spatial Function

An override value can be any of three things, resolved internally by the same helper:

**A scalar** — a uniform override, no position dependence. Useful for setting a value once
without mutating the shared parameter vector:

```@example spatial
p_uniform = SpatialContext([0.0, 0.0, 0.0], (IKr_Multiplier = 0.5,))
model(du, u, p_uniform, 0.0)
du[1]
```

**A callable** `(x, t) -> value` — maximum flexibility, CPU-only if it closes over
heap-allocated values:

```@example spatial
gradient_fn = (x, t) -> 0.5 + 0.5 * x[1]
p_closure   = SpatialContext([1.0, 0.0, 0.0], (IKr_Multiplier = gradient_fn,))
model(du, u, p_closure, 0.0)
du[1]
```

**An isbits functor** `<: SpatialFunction` — the GPU-safe option. CytoZoo ships three:

```@example spatial
Constant(0.5)                              # fixed value, ignores x and t
SpatialStep(1, 1.5, 1.0, 0.5)              # dim, threshold, below, above
SpatialGradient(1, 0.0, 2.0, 1.0, 0.5)     # dim, x_start, x_end, val_start, val_end
```

[`SpatialStep`](@ref) returns `below` when `x[dim] <= threshold` and `above` otherwise.
[`SpatialGradient`](@ref) interpolates linearly between `val_start` at `x_start` and
`val_end` at `x_end`, clamping outside that range:

```@example spatial
g = SpatialGradient(1, 0.0, 2.0, 1.0, 0.5)
[g([xi, 0.0, 0.0], 0.0) for xi in (-1.0, 0.0, 1.0, 2.0, 3.0)]
```

## GPU Compatibility

A `SpatialContext` is isbits when its position and every override are isbits. That is the
property a GPU kernel needs:

```@example spatial
using StaticArrays

p_gpu = SpatialContext(
    SVector(1.2, 0.5, 1.8),
    (IKr_Multiplier = SpatialStep(1, 1.5, 1.0, 0.5),),
)

isbitstype(typeof(p_gpu))
```

Compare with the closure version, which is not:

```@example spatial
isbitstype(typeof(p_closure))
```

Use a plain `Vector` for `x` and closures for overrides on CPU; use an `SVector` and
[`SpatialFunction`](@ref) functors when the model has to run on a device. You can also write
your own isbits callable — it need not subtype `SpatialFunction`, which exists for
discoverability rather than dispatch:

```@example spatial
struct RadialBlock{T}
    centre::SVector{3, T}
    radius::T
    inside::T
    outside::T
end

@inline function (f::RadialBlock)(x, t)
    d2 = sum(abs2, SVector{3}(x) .- f.centre)
    return d2 <= f.radius^2 ? f.inside : f.outside
end

blk = RadialBlock(SVector(1.0, 1.0, 1.0), 0.5, 0.2, 1.0)
isbitstype(typeof(blk))
```

## Zero Cost When Unused

The internal right-hand side dispatches on the *type* of the overrides. When that type is
`Nothing`, every spatial branch is eliminated at compile time — a non-spatial call is exactly
as fast as a model that never had spatial support:

```@example spatial
model(du, u, nothing, 0.0)     # no spatial branches compiled at all
du[1]
```

This is why models dispatch on `p` rather than testing it at runtime. There is no branch to
predict and no `nothing` check in the hot loop.

## Tissue Simulation

For actual tissue-level simulation you generally do not construct the `SpatialContext`
yourself — the tissue solver supplies each cell's position from the mesh. With Thunderbolt.jl
loaded, the adapter builds it internally:

```julia
using CytoZoo, Thunderbolt

overrides = (IKr_Multiplier = (x, t) -> x[1] > 2.0 ? 0.5 : 1.0,)
ionic = thunderbolt_model(ToRORd(); overrides)
```

See [Integrations](integrations.md).

## See Also

- [Stimulus](stimulus.md) — the stimulus receives the same `x`, but is configured on the
  model rather than through `overrides`.
- [The Cell Model Interface](interface.md) — why dispatch on `p` rather than a runtime check.
