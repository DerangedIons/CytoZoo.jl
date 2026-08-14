```@meta
CurrentModule = CytoZoo
```

# API Reference

```@index
```

## Types

```@docs
AbstractCellModel
AbstractCardiacCellModel
```

## Core Interface

Implemented by every model.

```@docs
default_initial_state
state_names
transmembrane_potential_index
```

## Atomic-Model Interface

Implemented by any model backed by a flat parameter vector. A composite
([`CoupledModel`](@ref)) leaves `num_parameters` and `parameter_names` undefined.

```@docs
num_states
num_parameters
parameter_names
```

## Named Access

```@docs
state_index
parameter_index
writable_parameters
```

## Rush-Larsen

```@docs
has_rush_larsen
rush_larsen_step!
```

## Derived Observables

```@docs
num_monitors
monitor_names
monitor_values!
monitor_history
```

## Spatial Context

```@docs
SpatialContext
SpatialFunction
Constant
SpatialStep
SpatialGradient
```

## Stimulus

```@docs
AbstractStimulus
Stimulus
FunctionStimulus
```

## Coupling

```@docs
couple
Subsystem
CoupledModel
share
connect
overwrite
```

## Models

```@docs
ToRORd
FHNModel
ParametrizedFHNModel
```

## Extensions

```@docs
thunderbolt_model
```

## Internal

Not public API. Names and behaviour may change without notice.

```@autodocs
Modules = [CytoZoo]
Public  = false
Order   = [:type, :function, :constant]
```
