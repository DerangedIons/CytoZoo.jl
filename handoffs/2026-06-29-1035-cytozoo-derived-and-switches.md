---
slug: cytozoo-derived-and-switches
created: 2026-06-29-1035
status: open
---

# Handoff: add DERIVED-mode observables + module-switch support to CytoZoo

## Goal / why this matters
The ECCMitoRedox G+K coupling (Gauthier ⊗ Kembro cardiac models) is built on CytoZoo's
"socket" abstraction. A readiness review found CytoZoo already supports two of the four socket
modes as first-class primitives — `WIRE` (`connect`) and `STATE`-merge (`share`) — and `PARAM`
works as plain model constants. The two modes with **no framework support** are:

1. **`DERIVED`** — algebraic conservation-law quantities (`ATPm = C_A − ADPm`,
   `NAD = C_PN − NADH`, `NADP = 0.1 − NADPH`) surfaced as **observables** in the solution.
2. **Module switches** — boolean flags gating whole subsystems on/off
   (`redox_on`, `imac_dpsi`, `dynamic_pH`, `CII_dynamic`, …), each with an OFF-invariant that
   doubles as a regression test.

Implement both so the coupling can (a) report derived quantities without hand-rolling them per
script, and (b) stage subsystems on/off with tested invariants per the architecture build order.

## Background & current state
- CytoZoo is a minimalist, solver-agnostic coupling framework. `couple(nodes, edges)` assembles
  one monolithic `f(dU,U,p,t)` over a shared global state vector; submodels implement
  `AbstractCellModel`. Two edge kinds exist: `connect` (directed dataflow → a source state is
  written into a receiver's parameter slot each RHS eval = **WIRE**) and `share` (two states
  become one owner-governed global slot = **STATE-merge**).
- **DERIVED has dormant scaffolding already present**: `num_monitors(model)` and
  `monitor_values!(mon, u, t, model)` are declared but **unexported, unimplemented, and never
  called**. `ToRORd` wires `num_monitors(::ToRORd) = TORORD_NUM_MONITORS` with the constant set
  to `0` and a TODO to port ~492 monitors from ArmyHeart. That heavy port is **not** what
  ECCMitoRedox needs — keep DERIVED lightweight (a handful of conservation quantities).
- **Module switches do not exist in any form.** No model gates terms on a boolean today; all
  fluxes are always computed. Crucially, `CardiacGauthierCell` and `MitoKembro` functors
  **ignore the `p` argument** entirely (`(m)(du,u,p,t) = _rhs!(du,u,m.parameters,...)`), so any
  `p`-based switch mechanism would force a functor rewrite. The dynamics already compute the
  conservation algebra internally, so DERIVED is about *exposing* values, not changing integration.

## Key files / locations
CytoZoo (`pro/dev/CytoZoo`):
- `src/interface.jl:125-139` — the dormant `num_monitors` / `monitor_values!` hook (no name
  accessor yet; no `p` arg in the monitor signature — fine, models hold params on the struct).
- `src/interface.jl:154-181` — `SpatialContext` + `_resolve_spatial`; the *time/space-varying*
  `p` payload. (Deliberately **not** where switches go — see Decisions.)
- `src/coupling.jl` — `CoupledModel`, `CouplingLayout` (`:122-130`), `_compute_layout` (`:241-300`,
  non-primary states prefixed `:<comp>_<state>` ~`:263`), `solution_indices` (global idx per
  component, local order), monolithic RHS `_run!` (`:402-410`, same `p` to every component),
  `_connect!` (`:415-428`), `_connect_value` primal-extraction hook (`:372`).
- `ext/SciMLBaseExt.jl` — `ODEProblem(model, tspan; u0, p=nothing, kwargs...)`; **`kwargs`
  pass through**, so `SavingCallback`/`observed`/`callback` already work with no new plumbing.
- `ext/ForwardDiffExt.jl` — `_connect_value(::Dual) = _connect_value(value(x))`; the primal-
  extraction pattern needed only if monitors are computed *mid-solve* under an implicit solver.
- `src/CytoZoo.jl:8-15` — export list. `src/models/torord/{ToRORd.jl,rhs.jl,monitors.jl}` —
  reference model + `num_monitors` stub. `test/test_coupling.jl:5-30,35-82` — mock-model patterns.
- `Project.toml` `[weakdeps]`/`[extensions]` — `SciMLBaseExt`, `ForwardDiffExt`, `ThunderboltExt`.

Model packages (consumers, separate repos under `~/dev`):
- `MitoKembro/src/rhs.jl` (`_mitokembro_rhs!`, 25 states, ~87 params; redox/CII/pH machinery lives
  here — most switches gate *these* terms, e.g. shunt source `du[15]/du[16]`, `V_IMAC` in `du[3]`,
  `Hm` ODE `du[5]`).
- `CardiacGauthierCell.jl/src/rhs.jl` (`_gauthier_rhs!`, 76 states).

Architecture spec (`ECCMitoRedox/architecture/`): `coupling-master-doc.md` §2 (socket modes), §3
(module-switch table + OFF-invariants), §6 (DERIVED examples); `feedback.md`/`feedforward.md`.

## Decisions & conclusions
- **DERIVED = activate the existing monitor hook, surface post-solve.** Add a `monitor_names`
  accessor (the hook has none), export the three functions, add `CoupledModel` glue, and a
  post-solve `monitor_history(sol, model)` helper. **Use post-solve, not a callback**: `sol.u` is
  plain `Float64` (Duals never reach the saved solution), so it sidesteps `_connect_value`
  entirely, and recompute cost for conservation algebra is trivial. A `SavingCallback` variant is
  only worth it for expensive monitors (e.g. ToRORd's 834 intermediates) — **mention, don't build.**
- **Module switches live as static config on the model struct — NOT through `p`, NOT by extending
  `SpatialContext`.** Reasons: (1) switches are set-once-per-run config, constant in time/space —
  `SpatialContext` is the wrong abstraction and overloading it conflates two concerns; (2) the G/K
  functors ignore `p` anyway, but already read `model.parameters`/struct fields, so struct config
  flows through the channel the RHS uses; (3) precedent — `ToRORd` already carries static config as
  struct fields (`celltype`, `stim`). Therefore **most switch work is in the model packages**, not
  CytoZoo. CytoZoo's role is thin: an introspection method + documenting the wire-gating idiom.
- **Wire-gating switches need no new infra.** A switch that toggles a coupling wire (e.g.
  `cyto_ions_dynamic`) is just edge inclusion: `edges = cyto_ions_dynamic ? (connect(...),) : ()`.
  Codify the idiom; don't build machinery for it.

## What's left / next steps
**DERIVED (CytoZoo):**
1. `src/interface.jl`: add `monitor_names(model)::NTuple{N,Symbol}` (mirror `state_names`). Keep
   `monitor_values!(mon, u, t, model)` signature.
2. `src/CytoZoo.jl`: export `num_monitors`, `monitor_names`, `monitor_values!`.
3. `src/coupling.jl`: implement for `CoupledModel` — `num_monitors(cm)` = Σ over components;
   `monitor_names(cm)` = component monitor names with the same non-primary prefixing as states;
   `monitor_values!(mon, U, t, cm)` loops components, slices `U[layout.solution_indices[ck]]`,
   calls each component's `monitor_values!` into its segment of `mon`, concatenates.
4. `ext/SciMLBaseExt.jl`: add `monitor_history(sol, model) -> (t, names::NTuple, values::Matrix)`
   (rows = monitors, cols = `sol.t`). Iterate `sol.u`, call `monitor_values!` per step.
5. `test/test_monitors.jl` (new): mock model with a known law `mon = C − u[i]`; assert
   `monitor_history` recovers it across a solve; assert coupled-model names are prefixed and
   per-component values match.

**Module switches (CytoZoo — thin):**
6. `src/interface.jl`: add optional `subsystem_switches(model)::NamedTuple` (default `(;)`); export it.
7. Document (docstring / `coupling.jl` header or README) the wire-gating idiom (edge inclusion).

**Module switches (model packages — the real work; track as a follow-up issue per repo):**
8. `MitoKembroModel` (and `GauthierCell` where a switch applies): add a typed `switches`
   `NamedTuple{Bool}` field, default = architecture minimal config (§3: `CII_dynamic=ON`,
   `cyto_ions_dynamic=ON`, `dynamic_pH=OFF`, `dynamic_Pi=OFF`, `acid_base=OFF`; redox per build
   order). Functor forwards `model.switches` into `_rhs!`; gate the §3 terms.
9. OFF-invariant regression tests (§3 table): each switch OFF reproduces its baseline
   (`shunt=0 ⇒ reproduces Gauthier`, `cyto_ions_dynamic=OFF ⇒ reproduces Kembro`, etc.), reusing
   the §13.1 G1-isolation harness pattern.

Suggested order: 1→5 (DERIVED end-to-end in CytoZoo, fully testable with mocks) → 6,7 (switch
introspection + idiom) → 8,9 (model-package gating, one subsystem at a time per the build order).

## Gotchas / constraints
- **Monitor signature has no `p`.** `monitor_values!(mon, u, t, model)` reads params from
  `model.parameters` (or struct fields). Fine for G/K/ToRORd. Don't add `p` — it'd diverge from the
  declared hook and the coupled `monitor_values!(mon, U, t, cm)` form.
- **Coupled `u` is the global merged vector.** A component's monitor must be computed from its own
  local slice via `layout.solution_indices[ck]` (global indices in local order), not the whole `U`.
- **Don't compute monitors in a callback unless forced.** Mid-solve under an implicit solver, `u`
  carries `Dual`s; you'd need `_connect_value`-style primal extraction (`coupling.jl:372`,
  `ext/ForwardDiffExt.jl`). Post-solve `monitor_history` avoids this — saved `sol.u` is `Float64`.
- **Element-type genericity (interface.jl:11-28).** Any new RHS-side code (switch gating in model
  packages) must compute in the state eltype `T` — wrap literals `T(...)`, guard `Int/Int`. Monitor
  buffers should be `zeros(eltype(u), num_monitors(model))`.
- **`CoupledModel` has no parameter vector** — `num_parameters(::CoupledModel)` throws by design
  (`coupling.jl:441`). Don't try to stuff switches into a coupled param vector; they belong on each
  component model.
- **Path C from the scouting (separate `FlagsContext` through `p`) and the `SpatialContext`-
  extension path were both rejected** — see Decisions. If you revisit, the blocker is that the G/K
  functors discard `p`.
- No secrets involved in this work. The model-package steps (8,9) touch *other* repos
  (`~/dev/MitoKembro`, `~/dev/CardiacGauthierCell.jl`) — open tracking issues there rather than
  editing them from a CytoZoo branch.
