<p align="center">
  <img src="docs/src/assets/logo.png" alt="CytoZoo.jl" width="200"/>
</p>

<h1 align="center">CytoZoo.jl</h1>

<p align="center">
  <a href="https://kylebeggs.github.io/CytoZoo.jl/stable/"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Stable"></a>
  <a href="https://kylebeggs.github.io/CytoZoo.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Dev"></a>
  <a href="https://github.com/kylebeggs/CytoZoo.jl/actions/workflows/CI.yml?query=branch%3Amain"><img src="https://github.com/kylebeggs/CytoZoo.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="Build Status"></a>
  <a href="https://codecov.io/gh/kylebeggs/CytoZoo.jl"><img src="https://codecov.io/gh/kylebeggs/CytoZoo.jl/branch/main/graph/badge.svg" alt="Coverage"></a>
</p>

A Julia package providing a registry of cardiac cell models with a common functor-based interface. Models work standalone with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) and integrate with [Thunderbolt.jl](https://github.com/termi-official/Thunderbolt.jl) for tissue-level simulation via a package extension.

**Key features:**
- Zero runtime dependencies — all model code is pure Julia arithmetic
- Dual-backend — hand-coded functor models and MTK symbolic models behind the same interface
- Functor interface — models are callable structs compatible with DifferentialEquations.jl
- Rush-Larsen exponential integrator support
- Spatial heterogeneity via the `Spatial` wrapper (zero overhead when unused)
- GPU-friendly design (flat parameter vectors, generic element types)

## Installation

```julia
using Pkg
Pkg.add("CytoZoo")
```

## Quick Start

```julia
using CytoZoo

model = ToRORd()                    # endocardial cell (default)
u     = default_initial_state(model)
du    = similar(u)

model(du, u, nothing, 0.0)         # evaluate the RHS
```

### With OrdinaryDiffEq.jl

```julia
using OrdinaryDiffEq

prob = ODEProblem(model, (0.0, 1000.0))   # uses default initial state
sol  = solve(prob, Tsit5())
```

### Rush-Larsen Step

```julia
u_new = similar(u)
rush_larsen_step!(u_new, u, 0.0, 0.01, model)
```

## Available Models

| Model | States | Parameters | Cell types | Rush-Larsen | Backend |
|-------|--------|------------|------------|-------------|---------|
| `ToRORd` | 65 | 181 | endocardial (0), epicardial (1), midmyocardial (2) | yes | functor |
| `BeelerReuter` | 8 | 7 | — | no | MTK |

```julia
ToRORd()                         # Float64, endocardial
ToRORd(; celltype=1)             # Float64, epicardial
ToRORd(Float32)                  # Float32, endocardial
ToRORd(Float32; celltype=2)      # Float32, midmyocardial
```

`BeelerReuter()` requires loading `MTKCardiacCellModels` and `ModelingToolkit` (see [MTK Integration](#mtk-integration) below).

## Interface

Every model `<: AbstractCardiacCellModel` implements:

```julia
(model)(du, u, p, t)                        # ODE right-hand side (functor)
num_states(model)                            # number of state variables
num_parameters(model)                        # number of parameters
transmembrane_potential_index(model)          # index of Vm in the state vector
default_initial_state(model)                 # initial condition vector
```

Optional:

```julia
has_rush_larsen(model)                       # whether Rush-Larsen is available
rush_larsen_step!(u_new, u, t, dt, model)    # exponential integrator step
state_index(model, :v)                       # state index by name
parameter_index(model, :GNa)                 # parameter index by name
state_names(model)                           # tuple of all state names
parameter_names(model)                       # tuple of all parameter names
symbolic_system(model)                       # MTK system (MTK-backed models only)
has_symbolic_system(model)                   # whether a symbolic system is available
```

## Named Access

```julia
model = ToRORd()

# look up indices by name
vi = state_index(model, :v)
gi = parameter_index(model, :GNa)

# read or modify parameters
model.parameters[gi]       # default GNa value
model.parameters[gi] = 11.0
```

## Spatial Heterogeneity

The `Spatial` wrapper adds position-dependent parameter modulation to any model. Provide a `NamedTuple` of `(x, t) -> value` functions keyed by parameter name:

```julia
spatial = Spatial(ToRORd(), (
    IKr_Multiplier = (x, t) -> x[1] > 1.5 ? 0.5 : 1.0,
    isHypoxic      = (x, t) -> x[3] > 2.0 ? 1.0 : 0.0,
))

x = [1.2, 0.5, 1.8]   # spatial coordinates
spatial(du, u, x, t)
```

When no spatial functions are provided, all spatial branches compile away at zero cost.

## MTK Integration

With [MTKCardiacCellModels](https://github.com/jClugstor/MTKCardiacCellModels) and [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) loaded, MTK-backed models become available:

```julia
using CytoZoo, MTKCardiacCellModels, ModelingToolkit, OrdinaryDiffEq

model = BeelerReuter()              # 8-state Beeler-Reuter 1977 model
prob  = ODEProblem(model, (0.0, 500.0))
sol   = solve(prob, Tsit5())

# access the symbolic system for inspection/composition
sys = symbolic_system(model)
```

MTK-backed models implement the same interface as traditional functor models — `num_states`, `state_names`, `default_initial_state`, etc. all work identically.

## Thunderbolt Integration

With [Thunderbolt.jl](https://github.com/termi-official/Thunderbolt.jl) loaded, convert any CytoZoo model to a `Thunderbolt.AbstractIonicModel`:

```julia
using CytoZoo, Thunderbolt

model = ToRORd()
ionic = thunderbolt_model(model)    # CytoZooIonicModel wrapper
# use `ionic` with Thunderbolt.MonodomainModel(...)
```

This also works with `Spatial`-wrapped models.
