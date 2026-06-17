# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is CytoZoo?

A Julia package providing a registry of cardiac cell models with a common functor-based interface. Models are hand-coded callable structs implementing the interface. Works standalone with DifferentialEquations.jl and integrates with Thunderbolt.jl and TWorld.jl via package extensions. Zero runtime dependencies in the base package.

## Commands

```bash
# Run tests
julia --project=. -e "using Pkg; Pkg.test()"

# Format (Runic) — install once: julia -e 'using Pkg; Pkg.Apps.add("Runic")'
runic --inplace src test          # or: julia -m Runic --inplace src test

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

Optional: `has_rush_larsen`/`rush_larsen_step!` (signature: `rush_larsen_step!(u_new, u, p, t, dt, model)`), `state_index`/`parameter_index` (Symbol-keyed Dict lookup). The `num_monitors`/`monitor_values!` hooks exist but are internal (unexported) until a model ships a real monitor implementation.

`SpatialContext{X, SF}` carries per-cell position (`x`) and spatial parameter overrides (`overrides` NamedTuple) through the DiffEq `p` argument. Spatial functions can be scalars, callables, or isbits functors (`<: SpatialFunction`). The internal RHS uses `_resolve_spatial` to handle all three and dispatches on `overrides::F` — when `F === Nothing`, all spatial branches compile away (zero overhead). GPU-compatible when all fields are isbits.

### Model layout (src/models/torord/)

Each model lives in its own directory with a standard file structure:
- `ToRORd.jl` — struct, constructors, functor, interface methods
- `parameters.jl` — `TORORD_PARAMETER_NAMES` tuple, `TORORD_PARAM_INDEX` Dict, `_torord_init_parameters!`
- `states.jl` — same pattern for states
- `rhs.jl` — `_torord_rhs_impl!(du, u, parameters, celltype, x, t, overrides::F) where {T, F}`
- `rush_larsen.jl` — `_torord_rush_larsen_impl!` (same signature + `dt`)
- `monitors.jl` — derived quantities (stubbed, TODO: port from ArmyHeart)

Parameters are stored as flat vectors on the struct for GPU compatibility. Named access via `parameter_index(model, :GNa)`. The type parameter `T` is the element type of the vectors; `F` is the overrides type (Nothing or NamedTuple).

Naming collision avoidance in the RHS: Faraday's constant → `F_param`, temperature → `T_val`, celltype → `celltype_val` (after spatial function resolution).

### Spatial context (src/spatial.jl)

GPU-safe isbits spatial functor types for use with `SpatialContext`: `Constant`, `SpatialStep`, `SpatialGradient`. All `<: SpatialFunction`. Users can define custom isbits callables `f(x, t) -> T` for GPU compatibility, or use closures for CPU-only simulations.

### Coupling (src/coupling.jl + ext/CouplingExt.jl)

Compose two or more `AbstractCellModel`s into one combined model solved by operator splitting (`OrdinaryDiffEqOperatorSplitting`, OS). Coupling is expressed as a **graph**: `Subsystem` nodes (a model + its inner solver + a `name`) joined by directed **edges**. `CoupledModel <: AbstractCardiacCellModel`, so couplings nest. The base `src/coupling.jl` is pure Julia (no OS dep): `couple(nodes, edges)` builds the node NamedTuple (keyed by node name), partitions `edges` by type, validates names, and precomputes the global-state layout (per-component `solution_indices`, canonical names, ICs, operator order); it implements the layout/query interface but is **not** a direct functor (it throws — splitting ≠ a single RHS). OS-dependent solving lives in `ext/CouplingExt.jl` (OS is a weakdep, the ext triggers with SciMLBase): `OperatorSplittingProblem(::CoupledModel, tspan)` and `coupled_algorithm(cm; scheme)` (reads each node's solver and orders them to match the internal operator order — no `inner` arg; errors if a node has no solver).

`Subsystem(model, alg = nothing; name = gensym(:subsystem))` — `alg` may be `nothing` for pure layout/query work; multi-node graphs whose edges reference nodes must pass explicit `name=` (the `gensym` default is only ergonomic for a single node). First node = primary component (bare state names + `vm_index`).

Two edge kinds (freely mixed in the edge list):
- **`share`** (`share(:A => :d, :B => :x; owner = :A)`) — `A.d` and `B.x` are one variable in a single global slot. **Hard-discard**: only the owner's equation drives the slot; the non-owner's operator zeroes the shared state's derivative (freezing it through the substep), so it reads the value but never writes it. Zero authoring change.
- **`connect`** (`connect(:A => :Vm, :B => :Vm_ext)`) — a directed dataflow edge: before `B` steps, an OS synchronizer (`forward_sync_external!`) writes `A`'s `Vm` into `B`'s parameter slot `:Vm_ext`, which `B`'s functor reads. The receiver must expose a writable `parameters` slot (the only authoring change coupling imposes). Carries an operation `op`: `overwrite` (default, copy) or `+` (sum all `+` edges into one slot — the slot is reset to zero then summed each step, a cross-sectional sum, not a running total); other ops are rejected. `_synchronizer` partitions edges by op into homogeneous `overwrites`/`adds` lists, so the per-substep write in `forward_sync_external!` has no dynamic dispatch or allocation.

Layout naming: the first component's states keep bare names; others are prefixed (`:B_y`); shared slots take the owner's name (or an explicit `name=`). Design + decisions in `coupling-redesign.md`; runnable demo in `examples/coupling_toy.jl`.

### Native adherence vs. ext fallback

Two integration patterns for model packages:

1. **Native adherence (packages we own)** — the model package depends on CytoZoo and declares its types as `<: CytoZoo.AbstractCardiacCellModel`, implementing the interface methods inside the model package. Reference example: `DerangedIons/TWorld.jl` defines `TWorldCellModel{P} <: CytoZoo.AbstractCardiacCellModel` in `src/cytozoo_interface.jl` and exports it. User writes `using TWorld` and gets the CytoZoo interface for free; `using CytoZoo, TWorld, OtherModel` lets them hot-swap behind a uniform interface.

2. **Ext fallback (third-party packages)** — when the upstream package can't take a CytoZoo dependency, CytoZoo writes a thin adapter in `ext/<Pkg>Ext.jl` that wraps the upstream type and implements the interface. The current `ThunderboltExt.jl` is the canonical example.

### Extensions

Three package extensions:

**SciMLBaseExt** (`ext/SciMLBaseExt.jl`) — loaded when OrdinaryDiffEq/SciMLBase is available. Adds `ODEProblem(model, tspan; u0=..., p=...)` convenience constructor for any `AbstractCellModel`.

**CouplingExt** (`ext/CouplingExt.jl`) — loaded when `OrdinaryDiffEqOperatorSplitting` (+ SciMLBase) is available. Adds operator-splitting solving for `CoupledModel` (see Coupling above).

**ThunderboltExt** (`ext/ThunderboltExt.jl`) — Thunderbolt's `MonodomainModel` requires `ION <: Thunderbolt.AbstractIonicModel`. The extension defines `CytoZooIonicModel{M, SF} <: Thunderbolt.AbstractIonicModel` as an adapter with an optional `overrides` field. Users call `thunderbolt_model(model; overrides=nothing)` (stub in base, implemented in ext). The extension constructs `SpatialContext(x, overrides)` from the mesh position internally.

### Stimulus

`AbstractStimulus` (interface.jl) is the supertype for stimulus current models; the contract is a callable `(s)(x, t) -> current` returning the full `Istim`. `x` is a position vector (matching `SpatialFunction`); a stimulus used on the non-spatial path must ignore `x` so `s(nothing, t)` works. Spatial dependence is first-class — a stimulus may index `x`. Built-ins: `Stimulus{T}` (closure-free isbits periodic pulse — amplitude/period/duration/start — for GPU and Rush-Larsen) and `FunctionStimulus{F}` (wraps an arbitrary `(x, t)` function for biphasic/S1–S2/ramps; isbits iff `F` is). Models call `stim(x, t)` directly. All are owned by CytoZoo and re-exported by model packages that adhere natively (e.g., TWorld).

### Adding a new model

1. Create `src/models/<name>/` with the standard file structure
2. Define struct `<: AbstractCardiacCellModel` with `parameters::T` and metadata fields
3. Implement the internal `_<name>_rhs_impl!` with `overrides::F where {T, F}` dispatch using `_resolve_spatial` for spatial parameter resolution
4. Add interface methods (functor with `p::Nothing` and `p::SpatialContext` dispatches, num_states, etc.)
5. Add `rush_larsen_step!` with `p` argument if applicable
6. Include in `src/CytoZoo.jl` and export

### Testing

Correctness tests compare CytoZoo output against ArmyHeart reference values (embedded in `test/test_torord_correctness.jl`) at `rtol=1e-10`. Performance tests verify zero allocations on the functor. SciMLBase extension tests verify `ODEProblem(model, tspan)` + `solve`. TWorld tests are conditional — skipped when TWorld is unavailable; they exercise the native-adherence path (`using TWorld` exposes `TWorldCellModel`). Source models for cross-validation live at `~/dev/ArmyHeart/` and `~/dev/TWorld/`.
