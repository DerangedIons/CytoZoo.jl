# TODO: Port monitor_values from ArmyHeart trauma_monitors.jl (3090 lines)
# Monitors are derived quantities (individual currents, fluxes, etc.) computed
# from states and parameters. Not required for core ODE integration.

const TORORD_NUM_MONITORS = 0  # Will be 492 once ported
