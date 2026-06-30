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
export state_index, parameter_index, state_names, parameter_names
export num_monitors, monitor_names, monitor_values!
export couple, CoupledModel, Subsystem, share, connect

# Stubs for extensions
function thunderbolt_model end
export thunderbolt_model
function monitor_history end
export monitor_history

# Models
include("models/torord/ToRORd.jl")
export ToRORd

end
