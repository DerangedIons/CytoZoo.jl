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
    rush_larsen_step!(u_new, u, t, dt, model::AbstractCellModel) -> Nothing

Single Rush-Larsen exponential integration step. Only available when
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
# Spatial wrapper
# ---------------------------------------------------------------------------

"""
    Spatial(model, spatial_funcs)

Wrapper that adds position-dependent parameter modulation to any cell model.
`spatial_funcs` is a `NamedTuple` of `(x, t) -> value` functions that override
or modulate specific parameters based on spatial coordinates.

The base model is position-independent. This wrapper adds tissue-level spatial
heterogeneity as an opt-in layer.

# Example
```julia
spatial = Spatial(ToRORd(), (
    IKr_Multiplier = (x, t) -> x[1] > 1.5 ? 0.5 : 1.0,
    isHypoxic = (x, t) -> x[1] > 2.0 ? 1.0 : 0.0,
))
```
"""
struct Spatial{M <: AbstractCellModel, F} <: AbstractCellModel
    model::M
    spatial_funcs::F
end

num_states(s::Spatial) = num_states(s.model)
num_parameters(s::Spatial) = num_parameters(s.model)
transmembrane_potential_index(s::Spatial) = transmembrane_potential_index(s.model)
default_initial_state(s::Spatial) = default_initial_state(s.model)
has_rush_larsen(s::Spatial) = has_rush_larsen(s.model)
state_index(s::Spatial, name::Symbol) = state_index(s.model, name)
parameter_index(s::Spatial, name::Symbol) = parameter_index(s.model, name)
state_names(s::Spatial) = state_names(s.model)
parameter_names(s::Spatial) = parameter_names(s.model)
num_monitors(s::Spatial) = num_monitors(s.model)
