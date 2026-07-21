# Implementing a Model

This guide covers writing a new cell model against the CytoZoo interface — whether it lives
inside CytoZoo or in your own package.

## A Minimal Complete Model

The smallest thing that satisfies the core contract:

```@example impl
using CytoZoo

struct TwoState{T} <: AbstractCardiacCellModel
    parameters::Vector{T}      # [g_leak, E_leak, tau_w]
end
TwoState() = TwoState([0.1, -80.0, 5.0])

CytoZoo.num_states(::TwoState)                    = 2
CytoZoo.num_parameters(::TwoState)                = 3
CytoZoo.state_names(::TwoState)                   = (:v, :w)
CytoZoo.parameter_names(::TwoState)               = (:g_leak, :E_leak, :tau_w)
CytoZoo.default_initial_state(::TwoState{T}) where {T} = T[-80.0, 0.0]
CytoZoo.transmembrane_potential_index(::TwoState) = 1

CytoZoo.state_index(::TwoState, n::Symbol)     = findfirst(==(n), (:v, :w))
CytoZoo.parameter_index(::TwoState, n::Symbol) = findfirst(==(n), (:g_leak, :E_leak, :tau_w))

function (m::TwoState)(du, u, p, t)
    g, E, τ = m.parameters
    du[1] = -g * (u[1] - E) + u[2]
    du[2] = (0.5 - u[2]) / τ
    return nothing
end
nothing # hide
```

```@example impl
m = TwoState()
du = similar(default_initial_state(m))
m(du, default_initial_state(m), nothing, 0.0)
du
```

Note `state_index` and `parameter_index` return `nothing` for an unknown name here, because
`findfirst` does. Preserve that if you switch to a `Dict` for a large model — use
`get(INDEX, name, nothing)` rather than `INDEX[name]`, which throws a `KeyError`.

It matters most for `parameter_index`: [`couple`](@ref) validates a `connect` target's slot by
testing `parameter_index(receiver, slot) === nothing`, so a throwing implementation turns an
actionable "no parameter slot `:x`" message into a bare `KeyError`. `ToRORd` currently has
this problem.

## The Standard File Layout

Models inside CytoZoo live in `src/models/<name>/` with a fixed structure:

| File | Contents |
|---|---|
| `<Name>.jl` | the struct, constructors, functor, and interface methods |
| `parameters.jl` | `<NAME>_PARAMETER_NAMES`, `<NAME>_PARAM_INDEX`, `_<name>_init_parameters!` |
| `states.jl` | the same pattern for states |
| `rhs.jl` | `_<name>_rhs_impl!` |
| `rush_larsen.jl` | `_<name>_rush_larsen_impl!` (if supported) |
| `monitors.jl` | derived quantities (if supported) |

Splitting the right-hand side into an internal `_<name>_rhs_impl!` keeps the functor methods
thin and lets both the non-spatial and spatial dispatches share one implementation.

## Parameters Live in a Flat Vector

Store parameters as a flat `AbstractVector` field, not as individual struct fields. Two
reasons:

1. **GPU compatibility** — a flat numeric vector with a concrete element type stays
   isbits-friendly.
2. **Coupling** — a [`connect`](@ref) edge writes into a parameter slot, which requires an
   indexable, mutable vector. [`writable_parameters`](@ref) defaults to `model.parameters`,
   so naming the field `parameters` gets connect-participation for free.

Named access is provided by `parameter_index`, so callers never hard-code positions.

## Element-Type Genericity

This is the convention most likely to trip you up, and it matters: a `Float32` model must run
end to end in `Float32`, with no `Float64` intermediates leaking in. That is what keeps the
same source usable on CPU and GPU.

Write the internal right-hand side so the working type is derived from the signature:

```julia
function _mymodel_rhs_impl!(
        du::AbstractVector{T}, u::AbstractVector{T}, parameters::AbstractVector, ...
    ) where {T}
```

Then follow three rules.

**Wrap every numeric literal in `T(...)`.** A bare `Float64` literal promotes the whole
expression back to `Float64`:

```julia
GNa = T(1.7) * T(11.7802)     # correct
GNa = 1.7 * 11.7802           # wrong — Float64 leaks in
```

**Guard integer division.** `Int / Int` yields `Float64` no matter what `T` is, because
`1 / 2 === 0.5`:

```julia
half = T(1) / T(2)            # correct
half = T(0.5)                 # also correct
half = 1 / 2                  # wrong
```

**Known caveat, not a bug.** For a `Float32` base and a non-integer `Float32` exponent, Base
computes `x^y` through a `Float64` `log2`/`exp2` scratch and narrows the result back to
`Float32`. This is unavoidable and does not violate the convention.

!!! tip "Watch the `du`/`u` element types"
    Typing both as `AbstractVector{T}` with the *same* `T`, as `ToRORd` does, forces them to
    match. That rules out Rosenbrock solvers, which compute a time gradient by calling the
    model with a dual-typed `du` against a plain `u`. If you want Rosenbrock methods to work,
    give `du` and `u` independent type parameters.

### Naming Collisions

Physical constants collide with Julia builtins and with the type parameter. The convention
used in `ToRORd`'s right-hand side:

| Quantity | Name |
|---|---|
| Faraday's constant | `F_param` |
| temperature | `T_val` |
| cell type, after spatial resolution | `celltype_val` |

## Adding Spatial Support

Give the internal right-hand side an `overrides::F` type parameter and resolve each
overridable parameter through `_resolve_spatial`:

```julia
function _mymodel_rhs_impl!(du::AbstractVector{T}, u::AbstractVector{T},
                            parameters, x, t, overrides::F) where {T, F}
```

Then add both functor dispatches:

```julia
(m::MyModel)(du, u, ::Nothing, t)         = _mymodel_rhs_impl!(du, u, m.parameters, nothing, t, nothing)
(m::MyModel)(du, u, p::SpatialContext, t) = _mymodel_rhs_impl!(du, u, m.parameters, p.x, t, p.overrides)
```

Because the implementation dispatches on `F`, every spatial branch compiles away when
`F === Nothing`. See [Spatial Heterogeneity](spatial.md).

## Optional Capabilities

**Rush-Larsen** — implement [`rush_larsen_step!`](@ref) with the model last, and flip
[`has_rush_larsen`](@ref):

```julia
CytoZoo.has_rush_larsen(::MyModel) = true
CytoZoo.rush_larsen_step!(u_new, u, p, t, dt, m::MyModel) = ...
```

**Monitors** — override all three of [`num_monitors`](@ref), [`monitor_names`](@ref), and
[`monitor_values!`](@ref). See [Derived Observables](monitors.md).

## Checklist for a New Model in CytoZoo

1. Create `src/models/<name>/` with the standard file structure.
2. Define the struct `<: AbstractCardiacCellModel` with a `parameters` field and metadata.
3. Implement `_<name>_rhs_impl!` with `overrides::F` dispatch, using `_resolve_spatial`.
4. Add the interface methods and both functor dispatches.
5. Add `rush_larsen_step!` if applicable.
6. Include it in `src/CytoZoo.jl` and export the type.

## Shipping a Model in Your Own Package

There are two integration patterns, and which one applies depends on whether you control the
package.

### Native Adherence — you own the package

Depend on CytoZoo, subtype `CytoZoo.AbstractCardiacCellModel`, and implement the interface
inside your package. This is the preferred pattern.

[TWorld.jl](https://github.com/DerangedIons/TWorld.jl) is the reference: it defines
`TWorldCellModel{P} <: CytoZoo.AbstractCardiacCellModel` in its own `src/cytozoo_interface.jl`
and exports it. Users write `using TWorld` and get the CytoZoo interface with no extra step,
and `using CytoZoo, TWorld, SomeOtherModel` lets them hot-swap behind one uniform interface.

Stimulus types are owned by CytoZoo, so re-export them (`Stimulus`, `FunctionStimulus`,
`AbstractStimulus`) rather than defining your own.

### Extension Fallback — you do not own the package

When an upstream package cannot take a CytoZoo dependency, CytoZoo can carry a thin adapter in
`ext/<Pkg>Ext.jl` that wraps the upstream type and implements the interface against it.
`ext/ThunderboltExt.jl` is the canonical example. See [Integrations](integrations.md).

## Testing a New Model

Existing practice in the repository, worth mirroring:

- **Correctness** — compare derivatives against a trusted reference implementation at a tight
  relative tolerance (`ToRORd` uses `rtol = 1e-10`).
- **Allocations** — assert the functor allocates nothing in the hot path.
- **Element types** — check a `Float32` model produces `Float32` derivatives throughout.

!!! note "Validate the initial state too, not just the right-hand side"
    A reference-matched right-hand side does not imply a reference-matched initial condition.
    `ToRORd` currently ships a `default_initial_state` that differs substantially from the
    pre-paced steady state its correctness tests use, so simulations started from the default
    need several beats to settle. Test both.

## See Also

- [The Cell Model Interface](interface.md) — the contract in full.
- [API Reference](../reference/api.md) — every signature.
