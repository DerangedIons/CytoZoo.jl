# Integrations

CytoZoo's base package has **zero runtime dependencies**. Everything that needs an external
package arrives through a Julia package extension, loaded automatically when you load the
relevant package.

| Extension | Loads with | Provides |
|---|---|---|
| `SciMLBaseExt` | SciMLBase / OrdinaryDiffEq | `ODEProblem` constructor, [`monitor_history`](@ref) |
| `ForwardDiffExt` | ForwardDiff | primal extraction for `connect` inputs |
| `ThunderboltExt` | Thunderbolt | [`thunderbolt_model`](@ref) tissue adapter |

You never load an extension directly — loading the trigger package is enough.

## SciMLBase: Solving

Loading any SciMLBase-based solver stack adds an `ODEProblem` constructor for any
`AbstractCellModel`:

```@example integ
using CytoZoo, OrdinaryDiffEq

model = ToRORd()
prob  = ODEProblem(model, (0.0, 500.0))
prob.tspan
```

The constructor defaults `u0` to `default_initial_state(model)` and `p` to `nothing`, and
forwards any other keyword to `ODEProblem`:

```julia
ODEProblem(model, tspan; u0 = my_state, p = my_spatial_context, saveat = 1.0)
```

This is also the entry point for solving a [`CoupledModel`](@ref) — a coupling is a model, so
it constructs a problem the same way.

The extension also provides [`monitor_history`](@ref), the post-solve helper for derived
observables. See [Derived Observables](monitors.md).

## ForwardDiff: Implicit Solvers

Implicit solvers pull in ForwardDiff, which loads this extension automatically. It exists for
exactly one purpose: a [`connect`](@ref) edge stages a source state into a receiver's
`Float64` parameter slot, and a dual number cannot be stored there. The extension extracts
the primal so the connect input is frozen within the Newton step.

You never call anything from this extension. But it changes what your results *mean*, and the
consequences are serious enough to have their own page — **derivatives taken through a
`connect` edge are silently wrong**. See [Coupling Limitations](coupling/limitations.md).

`share` edges are unaffected, because they flow through the state vector.

## Thunderbolt: Tissue-Level Simulation

[Thunderbolt.jl](https://github.com/termi-official/Thunderbolt.jl) solves tissue-scale
electrophysiology, and its `MonodomainModel` requires an ionic model subtyping
`Thunderbolt.AbstractIonicModel`. Rather than making CytoZoo depend on Thunderbolt, the
extension defines an adapter, `CytoZooIonicModel`, and exposes it through
[`thunderbolt_model`](@ref):

```julia
using CytoZoo, Thunderbolt

ionic = thunderbolt_model(ToRORd())
# use `ionic` with Thunderbolt.MonodomainModel(...)
```

For spatially heterogeneous tissue, pass `overrides`. You do **not** construct the
[`SpatialContext`](@ref) yourself — the adapter builds one per cell from the mesh position
Thunderbolt supplies:

```julia
overrides = (IKr_Multiplier = (x, t) -> x[1] > 2.0 ? 0.5 : 1.0,)
ionic = thunderbolt_model(ToRORd(); overrides)
```

With `overrides = nothing` (the default) the adapter passes `p = nothing` straight through,
so the non-spatial fast path is used and no spatial machinery is compiled in.

The adapter forwards the state count, transmembrane-potential index, and default initial
state to Thunderbolt's equivalents, so a CytoZoo model behaves as a native ionic model.

Calling [`thunderbolt_model`](@ref) without Thunderbolt loaded raises a `MethodError` — the
function exists as a stub in the base package so it can be documented and referenced, but it
has no methods until the extension loads.

## Writing an Adapter for Another Package

The Thunderbolt extension is the template for the **ext fallback** pattern: use it when an
upstream package cannot depend on CytoZoo. Define a wrapper type in `ext/<Pkg>Ext.jl` that
subtypes whatever the upstream framework requires, and forward the interface methods.

When you *do* control the model package, prefer native adherence instead — depend on CytoZoo
and subtype `AbstractCardiacCellModel` directly. See
[Implementing a Model](implementing_a_model.md).

## Documenting Extension Functions

A note for contributors, since it looks odd on first encounter. The docstrings for
[`thunderbolt_model`](@ref) and [`monitor_history`](@ref) live on **stubs in
`src/CytoZoo.jl`**, not in the extension modules. A docstring attached inside an extension
registers in that module's documentation metadata rather than CytoZoo's, so Documenter's
`Modules = [CytoZoo]` would never find it. Keep new extension entry points following the same
pattern.

## See Also

- [Coupling Limitations](coupling/limitations.md) — what ForwardDiffExt implies for gradients.
- [Spatial Heterogeneity](spatial.md) — the overrides Thunderbolt passes through.
