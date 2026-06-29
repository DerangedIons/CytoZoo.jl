# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is CytoZoo?

A Julia package providing a registry of cardiac cell models with a common functor-based interface. Models are hand-coded callable structs implementing the interface. Works standalone with DifferentialEquations.jl and integrates with Thunderbolt.jl via a package extension and with TWorld.jl via native interface adherence. Zero runtime dependencies in the base package.

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

### Coupling (src/coupling.jl)

Compose two or more `AbstractCellModel`s into one combined model with a shared global state. Coupling is expressed as a **graph**: `Subsystem` nodes (a model + a `name`) joined by directed **edges**. `CoupledModel <: AbstractCardiacCellModel` and is a **real functor** `(cm)(dU, U, p, t)` assembled from the submodels, so couplings nest *and* a coupling is solved monolithically with a single ODE solver. The base `src/coupling.jl` is pure Julia (no solver dep): `couple(nodes, edges)` builds the node NamedTuple (keyed by node name), partitions `edges` by type, validates names, precomputes the global-state layout (per-component `solution_indices`, canonical names, ICs, operator order), and builds the per-component **execution plan** (a concretely-typed tuple of `CompEntry` stored on the `CoupledModel`) that the functor walks: each submodel writes into its `view` of the shared `dU`/`U`, connect inputs are read live from `U` (`_connect!`), and non-owned shared-slot derivatives are zeroed (`frozen`) so the owner's write — last, by `operator_order` — wins. **Solve path:** `ODEProblem(cm, tspan)` + `solve(prob, alg)` (via `ext/SciMLBaseExt.jl`) — one solver over the whole system, no splitting error, stiff-capable. `num_parameters` throws (no single parameter vector). The plan-building helpers `_frozen_indices`/`_connect_plan` take raw `(components, shares/connects, ck)` so the plan can be assembled at `couple()` time.

`Subsystem(model; name = gensym(:subsystem))` — multi-node graphs whose edges reference nodes must pass explicit `name=` (the `gensym` default is only ergonomic for a single node). First node = primary component (bare state names + `vm_index`).

Two edge kinds (freely mixed in the edge list):
- **`share`** (`share(:A => :d, :B => :x; owner = :A)`) — `A.d` and `B.x` are one variable in a single global slot. **Hard-discard**: only the owner's equation drives the slot; the non-owner's contribution to that slot's derivative is zeroed (its `frozen` local indices), so it reads the value but never writes it. The functor orders the owner last (`operator_order`) so its write into the shared slot wins. Zero authoring change.
- **`connect`** (`connect(:A => :Vm, :B => :Vm_ext)`) — a directed dataflow edge: before `B`'s equations run, `A`'s `Vm` is written into `B`'s parameter slot `:Vm_ext` (via `_connect!`, read **live from `U`** each eval), which `B`'s functor reads. The receiver must expose a writable `parameters` slot (the only authoring change coupling imposes). Carries an operation `op`: `overwrite` (default, copy) or `+` (sum all `+` edges into one slot — reset to zero then summed each eval, a cross-sectional sum, not a running total); other ops are rejected. `_connect_plan` partitions edges by op into homogeneous `overwrites`/`adds` lists, so the per-eval write has no dynamic dispatch or allocation. Under an implicit solver the connect read passes through `_connect_value` — identity in base, but `ext/ForwardDiffExt.jl` extracts the primal of a `Dual` so the input is frozen within the Newton step (a `Dual` is never stored in the receiver's `Float64` slot; `share` is unaffected — it flows through `U`).

Layout naming: the first component's states keep bare names; others are prefixed (`:B_y`); shared slots take the owner's name (or an explicit `name=`). Design + decisions in `handoffs/2026-06-25-1515-coupling-monolithic-rhs.md`; runnable demo in `examples/coupling_toy.jl`.

### Derived observables (monitors)

DERIVED-mode quantities — algebraic functions of the state (e.g. conservation laws like `ATPm = C_A − ADPm`) — are surfaced as observables via the optional monitor hooks in `src/interface.jl`: `num_monitors` (default `0`), `monitor_names` (default `()`, mirrors `state_names`), and `monitor_values!(mon, u, t, model)` (no `p` arg — reads params from the struct). A model opts in by overriding all three. Surfaced **post-solve** by `monitor_history(sol, model)` (`ext/SciMLBaseExt.jl`) → `(; t, names, values::Matrix)` (rows = monitors, cols = `sol.t`); post-solve because the saved `sol.u` is plain (no `Dual`s), so it sidesteps `_connect_value` entirely. `CoupledModel` aggregates: `num_monitors` sums over components, `monitor_names` concatenates with the same non-primary `:<comp>_<name>` prefixing as states, and `monitor_values!` slices each component's own state (`layout.solution_indices[ck]`) — both walk `keys(cm.components)` (declaration order, not operator order) so names and values stay aligned. A 0-monitor model yields a `0×N` matrix without error. **Not ported:** ToRORd's ~492 ArmyHeart monitors (`TORORD_NUM_MONITORS = 0`). **Module switches** (gating whole subsystems) are a separate, unimplemented concern being reframed as graph surgery — see `handoffs/2026-06-29-1259-cytozoo-switches-graph-surgery.md`.

### Native adherence vs. ext fallback

Two integration patterns for model packages:

1. **Native adherence (packages we own)** — the model package depends on CytoZoo and declares its types as `<: CytoZoo.AbstractCardiacCellModel`, implementing the interface methods inside the model package. Reference example: `DerangedIons/TWorld.jl` defines `TWorldCellModel{P} <: CytoZoo.AbstractCardiacCellModel` in `src/cytozoo_interface.jl` and exports it. User writes `using TWorld` and gets the CytoZoo interface for free; `using CytoZoo, TWorld, OtherModel` lets them hot-swap behind a uniform interface.

2. **Ext fallback (third-party packages)** — when the upstream package can't take a CytoZoo dependency, CytoZoo writes a thin adapter in `ext/<Pkg>Ext.jl` that wraps the upstream type and implements the interface. The current `ThunderboltExt.jl` is the canonical example.

### Extensions

Package extensions:

**SciMLBaseExt** (`ext/SciMLBaseExt.jl`) — loaded when OrdinaryDiffEq/SciMLBase is available. Adds the `ODEProblem(model, tspan; u0=..., p=...)` convenience constructor for any `AbstractCellModel` (the default coupled-solve entry point for a `CoupledModel`) and the post-solve `monitor_history(sol, model)` helper for DERIVED observables (see Derived observables above).

**ForwardDiffExt** (`ext/ForwardDiffExt.jl`) — loaded when `ForwardDiff` is available (implicit solvers pull it in). Overrides `_connect_value(::Dual)` to extract the primal so a `connect` input is frozen within an implicit solver's Newton step instead of being stored as a `Dual` in the receiver's `Float64` parameter slot (see Coupling above).

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
