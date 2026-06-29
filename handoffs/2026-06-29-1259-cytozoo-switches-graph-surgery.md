---
slug: cytozoo-switches-graph-surgery
created: 2026-06-29-1259
status: open
---

# Handoff: reframe CytoZoo module switches as graph surgery

## Goal / why this matters
Design module-switch support in CytoZoo around a **graph-surgery** model — a switch
*cuts off a whole subgraph* of the coupled system — rather than the original "boolean
flags gating individual RHS terms" framing. This unifies couplings and switches under one
vocabulary and makes the OFF-invariant regression tests fall out naturally. This is a
design task: produce the abstraction + interface, then (separately) implement it.

## Background & current state
- The predecessor handoff (`handoffs/2026-06-29-1035-cytozoo-derived-and-switches.md`)
  scoped two missing socket modes: **DERIVED** observables and **module switches**.
- **DERIVED is being implemented now** (monitor hooks + post-solve `monitor_history`); it
  is independent of this work and not blocked by it.
- **Module switches were parked** (this doc). The original framing — a typed
  `switches::NamedTuple{Bool}` field on each model, the functor forwarding it into `_rhs!`,
  and the RHS doing `if switches.redox_on; du[i] += flux; end` — was rejected as the design
  target. It is imperative, scattered across the RHS, hand-rolled per model package, and
  reduces CytoZoo's role to a flat-`NamedTuple{Bool}` introspection accessor
  (`subsystem_switches`). That accessor was **deliberately not implemented** so we don't
  commit to a likely-wrong abstraction.

## Decisions & conclusions

### The reframe (Kyle)
Both the couplings we declare (`share`/`connect`) **and** a model's own internal dynamics
are **graphs**:
- A model *is* a graph — its states are nodes, the RHS dependencies linking them are edges.
- Coupling two models = taking two disjoint graphs and joining them with new edges.
- The unified coupled system is therefore one graph.

A **module switch is graph surgery: cutting off / disabling a whole subgraph** (a set of
nodes + their in/out edges) as a first-class operation — not ad-hoc per-term `if` guards.
Wire-gating (cut a *coupling* edge, e.g. `cyto_ions_dynamic`) and term-gating (cut an
*internal* subgraph, e.g. `redox_on`, `CII_dynamic`, `imac_dpsi`, `dynamic_pH`) collapse
into the **same operation at different scopes**.

### Why it's better than the term-flag framing
- **One vocabulary**: couplings *add* edges; switches *remove* subgraphs. Symmetric.
- **Declarative, not imperative**: "disable subsystem X" at the graph/config level instead
  of N scattered, easy-to-desync `if` guards in the RHS.
- **OFF-invariants fall out for free**: cutting a cleanly-added subgraph ⇒ reproduces the
  baseline model. The §3 regression criteria (`shunt=0 ⇒ reproduces Gauthier`,
  `cyto_ions_dynamic=OFF ⇒ reproduces Kembro`) become "remove the added subgraph, recover
  the original graph."

### Invariant that still holds from the original handoff
Switches are **set-once-per-run config, constant in time and space**. Whatever the
representation, it is resolved at `couple()` / construction time — **NOT** carried through
`p` / `SpatialContext` (that is the time/space-varying payload, the wrong abstraction). The
graph surgery should happen at build time so the compiled RHS only contains live subgraphs
(ideally zero runtime cost).

## Key files / locations
- `src/coupling.jl` — the current (coarse) graph machinery to extend or mirror at finer
  granularity: `Subsystem` nodes (whole models), `share`/`connect` edges (between models),
  `couple()`, `_compute_layout` (global-state layout + non-primary `Symbol(ck, :_, name)`
  prefixing), the `CompEntry` execution plan, and the `_run!` monolithic functor. **Study
  the build-time compile-away pattern** (`overrides::Nothing` ⇒ spatial branches vanish,
  zero overhead) as the model for zero-cost surgery.
- `src/interface.jl` — where a "toggleable subgraphs" accessor would live; also the
  `SpatialContext`/`p` payload that switches must avoid.
- `handoffs/2026-06-29-1035-cytozoo-derived-and-switches.md` — the superseded term-flag
  framing (§"Module switches"); its referenced §3 OFF-invariant table is still the
  acceptance criteria for the eventual implementation.
- ECCMitoRedox `architecture/coupling-master-doc.md` §2 (socket modes), §3 (switch table +
  OFF-invariants), §13.1 (G1-isolation harness pattern for the regression tests).
- Model packages (consumers, separate repos under `~/dev`): `MitoKembro/src/rhs.jl`
  (`_mitokembro_rhs!` — redox/CII/pH machinery; the internal subgraphs/fluxes to be cut,
  e.g. shunt source `du[15]/du[16]`, `V_IMAC` in `du[3]`, `Hm` ODE `du[5]`),
  `CardiacGauthierCell.jl/src/rhs.jl` (`_gauthier_rhs!`).

## What's left / next steps
1. **Decide granularity.** CytoZoo's graph is coarse (nodes = whole models; edges = between
   models). Cutting a subgraph *inside* a model needs finer granularity. Choose: does
   CytoZoo gain an internal graph representation, or do switches stay model-package-owned
   with CytoZoo providing only a declare/toggle *protocol*?
2. **Decide ownership.** The model package knows its internal structure; CytoZoo can't
   enumerate a model's subgraphs unless the model declares named, toggleable subsystems.
   Likely: a model declares its toggleable subgraphs; CytoZoo provides the vocabulary.
3. **Decide mechanism.** Options: (a) RHS reads a flag and zeros a subgraph's contributions
   (imperative but organized by subgraph); (b) a structural transform that yields a reduced
   model (fewer states — true surgery, but heavier); (c) an `enabled::Bool` on graph
   elements resolved at build time. Reuse the compile-away pattern for zero runtime cost.
4. **Unify wire-gating + term-gating** under one "cut this graph element" operation if the
   chosen representation allows it.
5. **Design the introspection surface** that replaces the rejected flat-`NamedTuple{Bool}`
   `subsystem_switches` (probably "list toggleable subgraphs / their state").
6. **Implement + OFF-invariant regression tests** (§3 table), reusing the §13.1 G1-isolation
   harness pattern. Model-package gating lands in those repos (track as issues there per the
   global rule — don't edit them from a CytoZoo branch).

## Gotchas / constraints
- **Do not route switches through `p`/`SpatialContext`** — they are static, set-once-per-run
  config; `p` is the time/space-varying payload. Resolve at `couple()`/construction.
- The Gauthier/Kembro functors currently **ignore `p`** entirely but already read
  `model.parameters` / struct fields — so static struct config flows through the channel the
  RHS uses. A `p`-based mechanism would force a functor rewrite.
- **Don't reintroduce the flat `subsystem_switches(model)::NamedTuple` accessor** as the
  answer; it was rejected pending this redesign.
- Any RHS-side gating added in model packages must compute in the state eltype `T` (wrap
  literals `T(...)`, guard `Int/Int`) — element-type genericity (`interface.jl` doc).
- No secrets involved. Model-package work touches *other* repos (`~/dev/MitoKembro`,
  `~/dev/CardiacGauthierCell.jl`) — open tracking issues there, don't edit from a CytoZoo
  branch.
