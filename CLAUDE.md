# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is CytoZoo?

A Julia package providing a registry of cardiac cell models with a common functor-based interface. Models work standalone with DifferentialEquations.jl and integrate with Thunderbolt.jl via a package extension. Zero runtime dependencies in the base package — all model code is pure Julia arithmetic.

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

**Functor-first design** — no separate `cell_rhs!`. The model IS a callable struct:
```julia
(model::ToRORd)(du, u, p, t) -> Nothing  # p is unused; parameters live on the struct
```

Required interface: functor, `num_states`, `num_parameters`, `transmembrane_potential_index`, `default_initial_state`.

Optional: `has_rush_larsen`/`rush_larsen_step!`, `state_index`/`parameter_index` (Symbol-keyed Dict lookup), `monitor_values!`.

`Spatial{M,F}` wrapper adds position-dependent parameter modulation via `NamedTuple` of `(x, t) -> value` functions. The internal RHS dispatches on `spatial_funcs::F` — when `F === Nothing`, all spatial branches compile away (zero overhead).

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

### Thunderbolt extension (ext/ThunderboltExt.jl)

`MonodomainModel` requires `ION <: Thunderbolt.AbstractIonicModel`. Since CytoZoo can't depend on Thunderbolt, the extension defines `CytoZooIonicModel{M} <: Thunderbolt.AbstractIonicModel` as an adapter. Users call `thunderbolt_model(model)` (stub in base, implemented in ext).

### Adding a new model

1. Create `src/models/<name>/` with the standard file structure
2. Define struct `<: AbstractCardiacCellModel` with `parameters::T` and metadata fields
3. Implement the internal `_<name>_rhs_impl!` with `spatial_funcs::F where {T, F}` dispatch
4. Add interface methods (functor, num_states, etc.)
5. Add `Spatial` wrapper support
6. Include in `src/CytoZoo.jl` and export
7. Add Thunderbolt dispatch in `ext/ThunderboltExt.jl`

### Testing

Correctness tests compare CytoZoo output against ArmyHeart reference values (embedded in `test/test_torord_correctness.jl`) at `rtol=1e-10`. Performance tests verify zero allocations on the functor. Source models for cross-validation live at `~/dev/ArmyHeart/` and `~/.julia/dev/TWorld/`.
