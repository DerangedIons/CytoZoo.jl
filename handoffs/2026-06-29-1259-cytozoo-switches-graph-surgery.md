---
slug: cytozoo-switches-graph-surgery
created: 2026-06-29-1259
status: in-progress
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

## Resolution (2026-06-29, branch `feat/module-switches`)

Design questions answered, CytoZoo-side machinery implemented + smoke-verified (23/23). Full
rationale in the plan; the answers:

1. **Granularity → protocol, not an internal graph.** CytoZoo provides a declare/toggle protocol;
   it does NOT model each model's internal flux DAG (the hand-coded flat-index RHS is the design).
   The only subgraph CytoZoo sees is a coupling edge — exactly the one wire-gate.
2. **Ownership → model declares, CytoZoo provides vocabulary** (same shape as the monitor hooks).
3. **Mechanism → compile-time, type-parameter compile-away.** `SwitchConfig{names, vals}` encodes
   resolved on/off state in the type domain; `is_on(sw, Val(:name))` constant-folds so disabled
   subgraphs are eliminated (the `overrides::Nothing` pattern). Rejected the reduced-model
   structural transform; noted the `frozen` idiom is insufficient (OFF also kills transport fluxes).
4. **Unify wire + term → unified vocabulary, dispatched by owner.** One `couple(...; switches=(…))`
   namespace; term/state → component (`with_switches`), wire → edge inclusion. No single execution
   path (that would need the rejected internal graph).
5. **Introspection → `Switch` descriptor (`default`/`requires`/`members`/`gates`/`off_recovers`) +
   `toggleable_subgraphs`;** the flat-Bool view demoted to a derived `switch_state`.
6. **Regression → OFF-invariants via the G1-isolation harness;** `off_recovers` tags each baseline.
   Physiological harness stays in the consumer (ECCMitoRedox).

**Implemented (CytoZoo, this branch):** `src/interface.jl` — `Switch`, `SwitchConfig`, `is_on`,
`switch_state`, `toggleable_subgraphs`, `with_switches`, `resolve_switches` (all exported);
`src/coupling.jl` — `enabled`/`switch` on `share`/`connect`, `switches=` distribution +
`_route_component_switch` + disabled-edge dropping in `couple`, `CoupledModel` aggregation of
`toggleable_subgraphs`/`switch_state`. Tests: `test/test_switches.jl`. Docs: CLAUDE.md "Module
switches". PRs back into `feat/coupling-interface`.

**Still open (consumer repos — track as issues, don't edit from a CytoZoo branch):** MitoKembro &
CardiacGauthierCell add the `S` switch type param + per-subgraph `is_on(sw, Val(:…)) ? flux :
zero(T)` guards and implement `toggleable_subgraphs`/`switch_state`/`with_switches` from §3 (incl.
`requires`/`members`/`off_recovers`). OFF-invariant regression tests + the G1-isolation harness
land in ECCMitoRedox (§13.1).

## Scope correction (2026-06-30, branch `feat/module-switches`)

Re-reading the ECCMitoRedox `architecture/` docs against the implementation (current focus:
**feedforward** G→K) surfaced a scope mismatch, and the Resolution machinery above was **pared
back to wire-gating only**. Rationale:

- The §3 switch table is **two things**: one **wire-gate** (`cyto_ions_dynamic`, the only
  feedforward switch — cuts the Cai/Nai *coupling wires*) and ~nine **internal term/state gates**
  (redox/CII/pH, all *inside MitoKembro's RHS* — the **feedback** stage).
- **Feedforward uses none of the term/state vocabulary.** Its one switch is a wire-gate, and the
  coupling graph already does wire-gating via `enabled`/`switch` on `share`/`connect` (one switch
  name can gate both Cai and Nai wires).
- The term/state vocabulary was **unvalidated** — no model implements the gating yet (the
  "Still open" consumer work above never landed), so exporting `Switch`/`SwitchConfig`/`is_on`/
  `resolve_switches`/`with_switches`/`toggleable_subgraphs` would lock in a likely-wrong
  abstraction — the exact failure this handoff set out to avoid.

**Removed (this branch):** the entire component term/state-switch vocabulary from
`src/interface.jl`; `_route_component_switch` + the component-routing half of `_apply_switches` +
the `CoupledModel` `toggleable_subgraphs`/`switch_state` aggregation from `src/coupling.jl`; the
exports from `src/CytoZoo.jl`; the vocabulary testsets from `test/test_switches.jl`.

**Kept:** edge wire-gating (`enabled`/`switch` on `share`/`connect`, `_set_enabled`,
disabled-edge dropping, `couple(...; switches=...)` toggling edges by name, now erroring on an
unknown key). This is what feedforward needs.

**Returns later:** the term/state-gate vocabulary comes back **validated by a real consumer**
when the feedback stage lands — answers 3/5/6 above still hold as the design; only the *timing*
was wrong (built before a consumer existed). Decisions captured in
`~/.claude/plans/we-are-adding-features-recursive-melody.md`.

**Next (follow-up branch):** wire feedforward G→K — two `connect` edges (`Cai`: G `u[12]` → K Cai
param; `Nai`: G `u[10]` → K Nai param), both `switch=:cyto_ions_dynamic`; ATPi/ADPi stay
G-native. Consumer issue: MitoKembro exposes `Cai`/`Nai` as writable parameter slots. G1-isolation
harness (`switches=(cyto_ions_dynamic=false,)` ⇒ standalone Kembro).

## Final decision (2026-06-30) — no switch primitive; wire-gating removed

Re-framing the whole effort at the **API/capability level** (not the Gauthier/Kembro instance)
settled it: the architecture's needs are the four **variable roles** (state / parameter /
derived / input), and the CytoZoo primitives **are** the role changes (`share` = state↔state
merge, `connect` = parameter→input, monitors = derived, `couple` = compose). A "module switch" is **not a coupling-API primitive**:

- Gating a coupling **element** (edge/node) ≡ **whether you include it when composing**. The
  OFF-invariant is "compose without the edge" (`cyto_ions_dynamic ? [edges] : []`) — no machinery.
- Gating a model's **internal** fluxes/states is invisible to the coupling graph → model-package
  concern (an optional declare/accept protocol), added only when a real model implements it.

So the kept wire-gating (`switch=`/`enabled=` on `share`/`connect` + `couple(...; switches=...)`)
was **removed**: redundant with composition, and premature (toy-only, no consumer) — the same
failure that retired the term/state vocabulary. The follow-up feedforward wiring is therefore just
the two `connect` edges (no `switch=`/`switches=`); the OFF-invariant is the edge-less `couple`.

**Removed (this branch, on top of the earlier pare-back):** `src/coupling.jl` reverted to HEAD
(drops `enabled`/`switch` fields, `_apply_switches`, `_set_enabled`, the `switches=` kwarg,
disabled-edge filtering); `test/test_switches.jl` deleted + dropped from `runtests.jl`. **Docs:**
README gains a "Variable roles" capability table; CLAUDE.md reframed ("Variable roles — and why
there's no switch primitive").

**Feedback foot-gun flagged (not built):** some feedback couplings need two models to each
contribute a term to the *same* shared state's derivative (A's core eq + B's extra flux). Current
`share` is hard-discard (one owner governs) and can't sum contributions — don't bake "one
governing equation per shared state" in as a permanent invariant. Resolve later via an additive/
contributory `share` or flux-injection through `connect`'s `+` op, when a consumer needs it.

Decisions captured in `~/.claude/plans/we-are-adding-features-cuddly-papert.md`.
