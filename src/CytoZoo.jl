module CytoZoo

include("interface.jl")
include("spatial.jl")
include("stimulus.jl")
include("coupling.jl")

export AbstractCellModel, AbstractCardiacCellModel
export SpatialContext, SpatialFunction
export Constant, SpatialStep, SpatialGradient
export AbstractStimulus, Stimulus, FunctionStimulus
export num_states, num_parameters, transmembrane_potential_index, default_initial_state
export has_rush_larsen, rush_larsen_step!
export state_index, parameter_index, state_names, parameter_names, writable_parameters
export num_monitors, monitor_names, monitor_values!
export couple, CoupledModel, Subsystem, share, connect, overwrite

# Stubs for extensions. Docstrings live here rather than on the extension methods so
# `@autodocs Modules = [CytoZoo]` can see them — a docstring attached inside an extension
# module registers in that module's doc META, not CytoZoo's.

"""
    thunderbolt_model(model::AbstractCellModel; overrides = nothing)

Wrap `model` as a `Thunderbolt.AbstractIonicModel` for use with `Thunderbolt.MonodomainModel`
and friends. `overrides` is an optional `NamedTuple` of spatial parameter overrides; the
adapter builds the [`SpatialContext`](@ref) from the mesh position internally.

Requires Thunderbolt to be loaded (implemented in `ext/ThunderboltExt.jl`); calling it
otherwise raises a `MethodError`.
"""
function thunderbolt_model end
export thunderbolt_model

"""
    monitor_history(sol, model::AbstractCellModel) -> (; t, names, values)

Recompute `model`'s derived/monitored quantities (see [`monitor_values!`](@ref)) across a
solved trajectory `sol`. Returns a `NamedTuple` of the time points `t`, the monitor `names`,
and a `values` matrix (rows = monitors in `names` order, columns = time points).

Post-solve by design: the saved `sol.u` is plain (no `Dual`s), so this sidesteps any
implicit-solver primal-extraction concerns. A model with `num_monitors(model) == 0` yields an
empty (`0 × length(sol.t)`) `values` matrix rather than erroring.

Requires a SciMLBase solver stack to be loaded (implemented in `ext/SciMLBaseExt.jl`).
"""
function monitor_history end
export monitor_history

# Models
include("models/torord/ToRORd.jl")
export ToRORd

include("models/fhn/FHN.jl")
# `ParametrizedFHNModel` is an alias of `FHNModel` (Thunderbolt's spelling), not a
# separate type — both are defined in `models/fhn/FHN.jl`.
export FHNModel, ParametrizedFHNModel

end
