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
- Functor interface — models are callable structs compatible with DifferentialEquations.jl
- Rush-Larsen exponential integrator support
- Spatial heterogeneity via `SpatialContext` in the DiffEq `p` argument (zero overhead when unused)
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
rush_larsen_step!(u_new, u, nothing, 0.0, 0.01, model)
```

## Available Models

| Model | States | Parameters | Cell types | Rush-Larsen |
|-------|--------|------------|------------|-------------|
| `ToRORd` | 65 | 177 | endocardial (0), epicardial (1), midmyocardial (2) | yes |

```julia
ToRORd()                         # Float64, endocardial
ToRORd(; celltype=1)             # Float64, epicardial
ToRORd(Float32)                  # Float32, endocardial
ToRORd(Float32; celltype=2)      # Float32, midmyocardial
```

Additional models live in their own packages that adhere to the interface natively (they depend on CytoZoo and subtype `AbstractCardiacCellModel`) — e.g. [`TWorld.jl`](https://github.com/DerangedIons/TWorld.jl) exposes `TWorldCellModel`. Load the package and drive it behind the same uniform interface (`using CytoZoo, TWorld` lets you hot-swap models).

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
rush_larsen_step!(u_new, u, p, t, dt, model)  # exponential integrator step
state_index(model, :v)                       # state index by name
parameter_index(model, :GNa)                 # parameter index by name
state_names(model)                           # tuple of all state names
parameter_names(model)                       # tuple of all parameter names
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

`SpatialContext` carries per-cell position and spatial parameter overrides through the DiffEq `p` argument:

```julia
model = ToRORd()
u = default_initial_state(model)
du = similar(u)

p = SpatialContext([1.2, 0.5, 1.8], (
    IKr_Multiplier = (x, t) -> x[1] > 1.5 ? 0.5 : 1.0,
    isHypoxic      = (x, t) -> x[3] > 2.0 ? 1.0 : 0.0,
))

model(du, u, p, 0.0)
```

Spatial functions can be scalars (`0.5`), closures (`(x, t) -> ...`), or GPU-safe isbits functors:

```julia
using StaticArrays

p = SpatialContext(
    SVector(1.2, 0.5, 1.8),
    (IKr_Multiplier = SpatialStep(1, 1.5, 1.0, 0.5),),  # step at x[1] = 1.5
)
isbitstype(typeof(p))  # true — GPU compatible
```

When `p = nothing`, all spatial branches compile away at zero cost.

The stimulus current is configured on the model, not through `overrides`. It
is evaluated as `stim(x, t)`, so a position-dependent stimulus is first-class —
pass a `FunctionStimulus` or a custom `AbstractStimulus` subtype:

```julia
# Periodic pulse that fires only where x[1] > 1.5
local_stim = FunctionStimulus((x, t) -> (mod(t, 1000.0) < 1.0 && x[1] > 1.5) ? -53.0 : 0.0)
model = ToRORd(; stim = local_stim)
```

## Thunderbolt Integration

With [Thunderbolt.jl](https://github.com/termi-official/Thunderbolt.jl) loaded, convert any CytoZoo model to a `Thunderbolt.AbstractIonicModel`:

```julia
using CytoZoo, Thunderbolt

model = ToRORd()
ionic = thunderbolt_model(model)    # CytoZooIonicModel wrapper
# use `ionic` with Thunderbolt.MonodomainModel(...)

# With spatial functions — Thunderbolt provides x from the mesh automatically
overrides = (IKr_Multiplier = (x, t) -> x[1] > 2.0 ? 0.5 : 1.0,)
ionic = thunderbolt_model(model; overrides)
```

## Coupling

Compose two or more models into one combined model and solve it by operator splitting. Coupling is a graph: **`Subsystem`** nodes (a model + its inner solver) joined by directed **edges**. Requires `OrdinaryDiffEqOperatorSplitting` loaded; it is a weak dependency, so the base package stays dependency-free.

```julia
using CytoZoo, OrdinaryDiffEqOperatorSplitting, OrdinaryDiffEq

coupled = couple(
    [Subsystem(ModelA(), Tsit5(); name = :A),
     Subsystem(ModelB(), Tsit5(); name = :B)],
    [ share(:A => :d, :B => :x; owner = :A),   # A.d ≡ B.x — A's equation governs the shared state
      connect(:A => :Vm, :B => :Vm_ext) ],     # B reads A's Vm through its :Vm_ext parameter slot
)

prob  = OperatorSplittingProblem(coupled, (0.0, 1000.0))
integ = init(prob, coupled_algorithm(coupled); dt = 0.05, adaptive = false)   # solvers read off the nodes
solve!(integ)
integ.u[state_index(coupled, :d)]    # value of the shared state at the final time
```

Two edge kinds, freely mixed in the edge list:

- **`share`** — two states are the *same* variable (one global slot); the `owner`'s equation governs it (the other's is discarded), while the non-owner still reads the value. Zero authoring change.
- **`connect`** — a directed dataflow edge: a source state is written into a receiver's parameter slot before the receiver steps, so the receiver reads it as an input. The receiver must expose that writable slot. Carries an operation `op`: `overwrite` (default, copy) or `+` to sum several edges into one slot (reset to zero then summed each step).

`CoupledModel` is itself an `AbstractCardiacCellModel`, so couplings nest. See [`examples/coupling_toy.jl`](examples/coupling_toy.jl) for a runnable demo.
