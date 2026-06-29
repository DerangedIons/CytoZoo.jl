"""
Supertype for all cell models in CytoZoo.

Every concrete model must implement:
- Functor `(model)(du, u, p, t) -> Nothing` — DiffEq-compatible ODE RHS
- [`num_states`](@ref)
- [`num_parameters`](@ref)
- [`transmembrane_potential_index`](@ref)
- [`default_initial_state`](@ref)

# Element-type genericity

A model's RHS must compute in the element type of its state vector, so a `Float32`
model runs end-to-end in `Float32` — no `Float64` intermediates leaking in from
numeric literals — which keeps the same code `Float64`- and GPU-compatible. The
convention (reference implementation: `_torord_rhs_impl!` in
`src/models/torord/rhs.jl`):

- **Derive the working type `T` from the signature.** Write the internal RHS as
  `f!(du::AbstractVector{T}, u::AbstractVector{T}, parameters::AbstractVector, …) where {T}`.
- **Wrap every numeric literal in `T(...)`:** `GNa = T(1.7) * T(11.7802)`, never
  `1.7 * 11.7802`. A bare `Float64` literal promotes the whole expression back to
  `Float64`.
- **Guard integer division:** `Int / Int` returns `Float64` regardless of `T`
  (`1 / 2 === 0.5`). Write `T(1) / T(2)` or `T(0.5)`.
- **Known caveat (not a bug):** for a `Float32` base `x` and non-integer `Float32`
  exponent `y`, Base computes `x^y` through a `Float64` `log2`/`exp2` scratch and
  narrows back to `Float32`. This is unavoidable and does not violate the convention.
"""
abstract type AbstractCellModel end

"""
Supertype for cardiac electrophysiology cell models.
"""
abstract type AbstractCardiacCellModel <: AbstractCellModel end

# ---------------------------------------------------------------------------
# Required interface
# ---------------------------------------------------------------------------

"""
    num_states(model::AbstractCellModel) -> Int

Total number of ODE state variables.
"""
function num_states end

"""
    num_parameters(model::AbstractCellModel) -> Int

Total number of model parameters.
"""
function num_parameters end

"""
    transmembrane_potential_index(model::AbstractCellModel) -> Int

Index of the transmembrane potential (Vm) in the state vector.
"""
function transmembrane_potential_index end

"""
    default_initial_state(model::AbstractCellModel) -> Vector

Default initial condition vector (length `num_states(model)`).
"""
function default_initial_state end

# ---------------------------------------------------------------------------
# Optional interface — name-based access
# ---------------------------------------------------------------------------

"""
    state_index(model::AbstractCellModel, name::Symbol) -> Int

Return the index of state variable `name` in the state vector.
"""
function state_index end

"""
    parameter_index(model::AbstractCellModel, name::Symbol) -> Int

Return the index of parameter `name` in the parameter vector.
"""
function parameter_index end

"""
    state_names(model::AbstractCellModel) -> NTuple{N, Symbol}

Tuple of all state variable names in order.
"""
function state_names end

"""
    parameter_names(model::AbstractCellModel) -> NTuple{N, Symbol}

Tuple of all parameter names in order.
"""
function parameter_names end

# ---------------------------------------------------------------------------
# Optional interface — Rush-Larsen
# ---------------------------------------------------------------------------

"""
    has_rush_larsen(model::AbstractCellModel) -> Bool

Whether the model provides a Rush-Larsen exponential integrator step.
"""
has_rush_larsen(::AbstractCellModel) = false

"""
    rush_larsen_step!(u_new, u, p, t, dt, model::AbstractCellModel) -> Nothing

Single Rush-Larsen exponential integration step. `p` is `nothing` for non-spatial
or `SpatialContext` for per-cell spatial variation. Only available when
`has_rush_larsen(model)` returns `true`.
"""
function rush_larsen_step! end

# ---------------------------------------------------------------------------
# Optional interface — monitors (DERIVED observables)
# ---------------------------------------------------------------------------
#
# Derived/monitored quantities are algebraic functions of the state (e.g. conservation-law
# values like `ATPm = C_A - ADPm`) surfaced as observables. A model opts in by overriding all
# three hooks; the post-solve `monitor_history(sol, model)` helper (SciMLBase extension) walks
# the saved solution and collects them. Defaults make a non-implementing model report zero
# monitors, so the helper returns an empty result rather than erroring.

"""
    num_monitors(model::AbstractCellModel) -> Int

Number of derived/monitored quantities the model can compute. Defaults to `0`; override
alongside [`monitor_names`](@ref) and [`monitor_values!`](@ref) to opt in.
"""
num_monitors(::AbstractCellModel) = 0

"""
    monitor_names(model::AbstractCellModel) -> NTuple{N, Symbol}

Names of the derived/monitored quantities, in the order [`monitor_values!`](@ref) writes
them (length `num_monitors(model)`). Defaults to `()`; mirrors [`state_names`](@ref).
"""
monitor_names(::AbstractCellModel) = ()

"""
    monitor_values!(mon, u, t, model::AbstractCellModel) -> Nothing

Compute derived quantities from state `u` at time `t` and store the `num_monitors(model)`
values in `mon`. Reads parameters from `model` (struct fields), not from a `p` argument.
"""
function monitor_values! end

# ---------------------------------------------------------------------------
# Spatial context — per-cell parameter variation via p
# ---------------------------------------------------------------------------

"""
    SpatialContext(x, overrides)

Per-cell spatial context passed as the `p` argument in `model(du, u, p, t)`.
Bundles the cell's position `x` with a `NamedTuple` of spatial parameter
overrides (scalars, callables, or isbits functors).

GPU-compatible (isbits) when both `X` and every element of `SF` are isbits.
"""
struct SpatialContext{X, SF}
    x::X
    overrides::SF
end

@inline _resolve_spatial(v::Number, x, t) = v
@inline _resolve_spatial(f, x, t) = f(x, t)

"""
    SpatialFunction

Abstract supertype for isbits spatial functor types (e.g., `SpatialStep`,
`SpatialGradient`). Used for discoverability, not dispatch.
"""
abstract type SpatialFunction end

"""
    AbstractStimulus

Supertype for stimulus current models. A subtype must be callable as
`(s::AbstractStimulus)(x, t)` and return the stimulus current directly, where
`x` is a position vector (matching `SpatialFunction`) and `t` is time.

A stimulus used on the non-spatial path (`model(du, u, nothing, t)`) must be
position-independent — it must not dereference `x`, so that `s(nothing, t)`
works. Subtypes that index `x` are only valid under a `SpatialContext`.
"""
abstract type AbstractStimulus end
