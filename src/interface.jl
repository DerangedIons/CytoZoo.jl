"""
Supertype for all cell models in CytoZoo.

Every concrete model must implement:
- Functor `(model)(du, u, p, t) -> Nothing` — DiffEq-compatible ODE RHS
- [`num_states`](@ref)
- [`num_parameters`](@ref)
- [`transmembrane_potential_index`](@ref)
- [`default_initial_state`](@ref)
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
# Optional interface — monitors
# ---------------------------------------------------------------------------

"""
    num_monitors(model::AbstractCellModel) -> Int

Number of derived/monitored quantities the model can compute.
"""
num_monitors(::AbstractCellModel) = 0

"""
    monitor_values!(mon, u, t, model::AbstractCellModel) -> Nothing

Compute derived quantities from state `u` at time `t` and store in `mon`.
"""
function monitor_values! end

# ---------------------------------------------------------------------------
# Optional interface — symbolic system (MTK-backed models)
# ---------------------------------------------------------------------------

"""
    symbolic_system(model::AbstractCellModel)

Return the underlying ModelingToolkit system for symbolic inspection and composition.
Returns `nothing` for models without a symbolic representation.
"""
symbolic_system(::AbstractCellModel) = nothing

"""
    has_symbolic_system(model::AbstractCellModel) -> Bool

Whether the model provides a symbolic MTK system via [`symbolic_system`](@ref).
"""
has_symbolic_system(::AbstractCellModel) = false

# ---------------------------------------------------------------------------
# Spatial context — per-cell parameter variation via p
# ---------------------------------------------------------------------------

"""
    SpatialContext(x, spatial_funcs)

Per-cell spatial context passed as the `p` argument in `model(du, u, p, t)`.
Bundles the cell's position `x` with a `NamedTuple` of spatial parameter
overrides (scalars, callables, or isbits functors).

GPU-compatible (isbits) when both `X` and every element of `SF` are isbits.
"""
struct SpatialContext{X, SF}
    x::X
    spatial_funcs::SF
end

@inline _resolve_spatial(v::Number, x, t) = v
@inline _resolve_spatial(f, x, t) = f(x, t)

"""
    SpatialFunction

Abstract supertype for isbits spatial functor types (e.g., `SpatialStep`,
`PeriodicPulse`). Used for discoverability, not dispatch.
"""
abstract type SpatialFunction end

