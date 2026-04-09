# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is CytoZoo?

A Julia package providing a registry of cardiac cell models with a common functor-based interface. Models can be traditional hand-coded functors or MTK symbolic systems — both implement the same interface. Works standalone with DifferentialEquations.jl and integrates with Thunderbolt.jl and MTKCardiacCellModels via package extensions. Zero runtime dependencies in the base package.

## Commands

```bash
# Run tests
julia --project=. -e "using Pkg; Pkg.test()"

# Format (Blue style)
julia --project=. -e "using JuliaFormatter; format(\"src\")"

# Load and quick-test
julia --project=. -e "using CytoZoo; m = ToRORd(); du = similar(default_initial_state(m)); m(du, default_initial_state(m), nothing, 0.0); println(du[1])"
```

## Architecture

### Interface (src/interface.jl)

Type hierarchy: `AbstractCellModel` → `AbstractCardiacCellModel` → concrete models (e.g., `ToRORd`).

**Functor-first design** — no separate `cell_rhs!`. The model IS a callable struct that dispatches on `p`:
```julia
(model::ToRORd)(du, u, ::Nothing, t)          # non-spatial (parameters on struct)
(model::ToRORd)(du, u, p::SpatialContext, t)   # spatial (per-cell variation via p)
```

Required interface: functor, `num_states`, `num_parameters`, `transmembrane_potential_index`, `default_initial_state`.

Optional: `has_rush_larsen`/`rush_larsen_step!` (signature: `rush_larsen_step!(u_new, u, p, t, dt, model)`), `state_index`/`parameter_index` (Symbol-keyed Dict lookup), `monitor_values!`, `symbolic_system`/`has_symbolic_system` (MTK-backed models only).

`SpatialContext{X, SF}` carries per-cell position (`x`) and spatial parameter overrides (`spatial_funcs` NamedTuple) through the DiffEq `p` argument. Spatial functions can be scalars, callables, or isbits functors (`<: SpatialFunction`). The internal RHS uses `_resolve_spatial` to handle all three and dispatches on `spatial_funcs::F` — when `F === Nothing`, all spatial branches compile away (zero overhead). GPU-compatible when all fields are isbits.

### Model layout (src/models/torord/)

Each model lives in its own directory with a standard file structure:
- `ToRORd.jl` — struct, constructors, functor, interface methods
- `parameters.jl` — `TORORD_PARAMETER_NAMES` tuple, `TORORD_PARAM_INDEX` Dict, `_torord_init_parameters!`
- `states.jl` — same pattern for states
- `rhs.jl` — `_torord_rhs_impl!(du, u, parameters, celltype, x, t, spatial_funcs::F) where {T, F}`
- `rush_larsen.jl` — `_torord_rush_larsen_impl!` (same signature + `dt`)
- `monitors.jl` — derived quantities (stubbed, TODO: port from ArmyHeart)

Parameters are stored as flat vectors on the struct for GPU compatibility. Named access via `parameter_index(model, :GNa)`. The type parameter `T` is the element type of the vectors; `F` is the spatial_funcs type (Nothing or NamedTuple).

Naming collision avoidance in the RHS: Faraday's constant → `F_param`, temperature → `T_val`, celltype → `celltype_val` (after spatial function resolution).

### Spatial context (src/spatial.jl)

GPU-safe isbits spatial functor types for use with `SpatialContext`: `Constant`, `SpatialStep`, `SpatialGradient`, `PeriodicPulse`. All `<: SpatialFunction`. Users can define custom isbits callables `f(x, t) -> T` for GPU compatibility, or use closures for CPU-only simulations.

### Extensions

Five package extensions, all via weak dependencies:

**SciMLBaseExt** (`ext/SciMLBaseExt.jl`) — loaded when OrdinaryDiffEq/SciMLBase is available. Adds `ODEProblem(model, tspan; u0=..., p=...)` convenience constructor for any `AbstractCellModel`.

**ThunderboltExt** (`ext/ThunderboltExt.jl`) — `MonodomainModel` requires `ION <: Thunderbolt.AbstractIonicModel`. The extension defines `CytoZooIonicModel{M, SF} <: Thunderbolt.AbstractIonicModel` as an adapter with an optional `spatial_funcs` field. Users call `thunderbolt_model(model; spatial_funcs=nothing)` (stub in base, implemented in ext). The extension constructs `SpatialContext(x, spatial_funcs)` from the mesh position internally.

**TWorldExt** (`ext/TWorldExt.jl`) — loaded when TWorld is available. Implements the full CytoZoo interface for `TWorldCellModel{P}` (92 states), including Rush-Larsen support with task-local workspace. Accepts `SpatialContext` in the functor but spatial_funcs threading to TWorld internals is pending.

**ThunderboltTWorldExt** (`ext/ThunderboltTWorldExt.jl`) — loaded when both Thunderbolt and TWorld are available. Adds `cell_rhs!` overloads for `CytoZooIonicModel{<:TWorldCellModel}`, constructing `SpatialContext` from mesh position.

**MTKCardiacCellModelsExt** (`ext/MTKCardiacCellModelsExt.jl`) — loaded when MTKCardiacCellModels + ModelingToolkit + SciMLBase are available. Defines `MTKCardiacModel{S,Prob} <: AbstractCardiacCellModel` which wraps an MTK-compiled `ODEProblem` and implements the full CytoZoo interface. Contains the `BeelerReuter()` model (8 states, 7 parameters) built from MTKCardiacCellModels components. The compiled model is cached per session to avoid re-compilation. MTK-backed models expose their symbolic system via `symbolic_system(model)`.

### Adding a new model

**Traditional (hand-coded functor):**

1. Create `src/models/<name>/` with the standard file structure
2. Define struct `<: AbstractCardiacCellModel` with `parameters::T` and metadata fields
3. Implement the internal `_<name>_rhs_impl!` with `spatial_funcs::F where {T, F}` dispatch using `_resolve_spatial` for spatial parameter resolution
4. Add interface methods (functor with `p::Nothing` and `p::SpatialContext` dispatches, num_states, etc.)
5. Add `rush_larsen_step!` with `p` argument if applicable
6. Include in `src/CytoZoo.jl` and export
7. Add Thunderbolt dispatch in `ext/ThunderboltExt.jl`

**MTK-backed (symbolic):**

1. Add model constructor stub in `src/CytoZoo.jl` (`function ModelName end; export ModelName`)
2. Implement model in `ext/MTKCardiacCellModelsExt.jl` using MTKCardiacCellModels components
3. Call `_build_mtk_model()` helper with the simplified system, ODEProblem, and Vm symbol
4. Cache the compiled model in a `Ref{Any}(nothing)` to avoid re-compilation

### Testing

Correctness tests compare CytoZoo output against ArmyHeart reference values (embedded in `test/test_torord_correctness.jl`) at `rtol=1e-10`. Performance tests verify zero allocations on the functor. SciMLBase extension tests verify `ODEProblem(model, tspan)` + `solve`. MTK extension tests are conditional — skipped when MTKCardiacCellModels is unavailable (not registered in General). Source models for cross-validation live at `~/dev/ArmyHeart/` and `~/.julia/dev/TWorld/`.
