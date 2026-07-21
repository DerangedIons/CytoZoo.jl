<p align="center">
  <img src="docs/src/assets/logo.png" alt="CytoZoo.jl" width="200"/>
</p>

<h1 align="center">CytoZoo.jl</h1>

<p align="center"><em>Cardiac cell models with a common interface — swappable, spatially heterogeneous, and composable.</em></p>

<p align="center">
  <a href="https://derangedions.github.io/CytoZoo.jl/stable/"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Stable"></a>
  <a href="https://derangedions.github.io/CytoZoo.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Dev"></a>
  <a href="https://github.com/DerangedIons/CytoZoo.jl/actions/workflows/CI.yml?query=branch%3Amain"><img src="https://github.com/DerangedIons/CytoZoo.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="Build Status"></a>
  <a href="https://codecov.io/gh/DerangedIons/CytoZoo.jl"><img src="https://codecov.io/gh/DerangedIons/CytoZoo.jl/branch/main/graph/badge.svg" alt="Coverage"></a>
</p>

---

## Why CytoZoo.jl?

Cardiac cell models are usually shipped as one-off scripts: a hard-coded state vector, magic
indices, and no way to reuse the model in a different context. CytoZoo gives every model the
same functor-based interface, so a model is a callable struct you can hand straight to
DifferentialEquations.jl — and swapping one model for another is a one-line change.

It also **composes** models. Join several into a graph and CytoZoo assembles one combined
right-hand side, integrated by a single ODE solver: no operator-splitting error, and a stiff
coupling can use one implicit method across the whole system.

Other things that might matter to you:

- **Zero runtime dependencies** — the base package is pure Julia arithmetic; solving, tissue
  coupling, and autodiff all arrive as package extensions
- **Spatial heterogeneity** via `SpatialContext` in the DiffEq `p` argument, compiled away
  entirely when unused
- **Rush-Larsen** exponential integration for long pacing runs
- **GPU-friendly design** — flat parameter vectors, generic element types, isbits throughout
- **[Thunderbolt.jl](https://github.com/termi-official/Thunderbolt.jl) integration** for
  tissue-level simulation, via an extension

## Installation

```julia
using Pkg
Pkg.add("CytoZoo")
```

Requires Julia 1.10 or later.

## Quick Start

```julia
using CytoZoo, OrdinaryDiffEq

model = ToRORd()                          # endocardial cell, Float64
prob  = ODEProblem(model, (0.0, 1000.0))  # uses the model's default initial state
sol   = solve(prob, FBDF())

# everything is addressable by name — no magic indices
state_names(model)                        # (:v, :jca, :m, :mL, ...)
sol.u[end][state_index(model, :v)]        # transmembrane potential at t = 1000 ms
model.parameters[parameter_index(model, :GNa)] = 11.0
```

Composing two models, with one shared state and one dataflow edge:

```julia
coupled = couple(
    [Subsystem(ModelA(); name = :A),
     Subsystem(ModelB(); name = :B)],
    [share(:A => :d, :B => :x; owner = :A),    # one variable; A's equation governs it
     connect(:A => :v, :B => :v_ext)],         # B reads A's v through a parameter slot
)

sol = solve(ODEProblem(coupled, (0.0, 1000.0)), Tsit5())
```

## Available Models

| Model | States | Parameters | Cell types | Rush-Larsen |
|-------|--------|------------|------------|-------------|
| `ToRORd` | 65 | 177 | endocardial, epicardial, midmyocardial | yes |

Further models live in their own packages and adhere to the interface natively — for example
[TWorld.jl](https://github.com/DerangedIons/TWorld.jl). Load both and hot-swap between them
behind one uniform interface.

> [!TIP]
> The examples above use defaults throughout. Cell type, element type, stimulus waveform,
> spatial overrides, and solver choice are all configurable — see the
> **[Getting Started guide](https://derangedions.github.io/CytoZoo.jl/stable/getting_started)**.

## Documentation

- **[Getting Started](https://derangedions.github.io/CytoZoo.jl/stable/getting_started)** — build, solve, and plot an action potential
- **[The Cell Model Interface](https://derangedions.github.io/CytoZoo.jl/stable/guides/interface)** — the tiered contract every model satisfies
- **[Coupling](https://derangedions.github.io/CytoZoo.jl/stable/guides/coupling/)** — composing models, and the limitations to know about first
- **[Spatial Heterogeneity](https://derangedions.github.io/CytoZoo.jl/stable/guides/spatial)** — per-cell parameter variation for tissue work
- **[Quick Reference](https://derangedions.github.io/CytoZoo.jl/stable/guides/quickref)** — the cheat sheet
- **[API Reference](https://derangedions.github.io/CytoZoo.jl/stable/reference/api)** — full function documentation

## Contributing

Contributions are welcome:

- Report bugs or request features on [GitHub Issues](https://github.com/DerangedIons/CytoZoo.jl/issues)
- Submit pull requests

Part of the [DerangedIons](https://github.com/DerangedIons) organization.
