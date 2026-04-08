module CytoZoo

include("interface.jl")

export AbstractCellModel, AbstractCardiacCellModel, Spatial
export num_states, num_parameters, transmembrane_potential_index, default_initial_state
export has_rush_larsen, rush_larsen_step!
export state_index, parameter_index, state_names, parameter_names
export num_monitors, monitor_values!
export symbolic_system, has_symbolic_system

# Stubs for extensions
function thunderbolt_model end
export thunderbolt_model

function BeelerReuter end
export BeelerReuter

struct TWorldCellModel{P} <: AbstractCardiacCellModel
    params::P
end
export TWorldCellModel

# Models
include("models/torord/ToRORd.jl")
export ToRORd

end
