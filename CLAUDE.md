# CytoZoo.jl

Uniform interface for cardiac cell models, plus a graph-based coupling system that composes
several models into one monolithically-solved ODE.

## Key files
- interface contract and monitor hooks: `src/interface.jl`
- coupling graph (`couple`, `share`, `connect`, `CoupledModel`): `src/coupling.jl`
- state clamps (`ClampedCell`, `clamp_states`, `base_model`): `src/clamp.jl`
- spatial context and GPU-safe spatial functors: `src/spatial.jl`
- stimulus types: `src/stimulus.jl`
- models, one directory each: `src/models/<name>/` (`<Name>.jl`, `parameters.jl`, `states.jl`, `rhs.jl`, `rush_larsen.jl`, `monitors.jl`)
- extensions: `ext/SciMLBaseExt.jl` (ODEProblem + `monitor_history`), `ext/ForwardDiffExt.jl`, `ext/ThunderboltExt.jl`

## Traps
- **ForwardDiff derivatives taken *through* a `connect` edge silently lose the coupling term.** `ext/ForwardDiffExt.jl` extracts the primal of a `Dual` so the input is frozen within a Newton step — correct fixed point, approximate Jacobian. Consequence: sensitivity analysis and gradient-based parameter fitting across a `connect` are **wrong, not merely slow**. `share` is unaffected (it flows through `U`). Fix, if it becomes load-bearing, is a `PreallocationTools.DiffCache` on the receiver's parameter vector — not built.
- **A `CoupledModel` is not thread-safe to solve concurrently.** `couple` deepcopies each connect receiver, but a single `cm` still stages into shared scratch — deepcopy the `cm` per trajectory in an `EnsembleProblem` `prob_func`.
- **A `connect` source resolves against states first, then monitors.** A DERIVED quantity (conservation-law complement, etc.) is wired by declaring it as a monitor on the model that owns the law — no per-edge transform exists or is wanted. Monitor sources are resolved in a single pre-pass before any staging, which is why a monitor source may **not** also receive a `connect` edge (rejected at `couple` time; it would lag by one eval). Wiring one monitor computes that model's whole monitor vector.
- **`share` defaults to hard-discard (one owner governs the slot, non-owners' derivatives zeroed); `op = +` makes a non-owner *contributory*, adding its derivative instead.** Contributions sum, the sum is cross-sectional (every accumulating slot is reset once per eval, never a running total), a class may mix both kinds, and an accumulating class imposes **no operator ordering** — so mutual contribution is legal where hard-discard shares reject it as a cycle. The owner still supplies the slot's initial value; a contributor's is ignored. Mechanism is save/restore around each member's write, not a re-ordering — see `_run!` / `_accumulated_indices` in `src/coupling.jl`. `inject(...)` (a monitor-sourced flux edge) was deliberately NOT built: it would inherit the monitor-receive ban and lose the term under ForwardDiff.
- **A component RHS must WRITE every `du` slot, never `+=` into one it has not assigned.** `du` arrives dirty. Under a contributory share the coupling adds each write onto a running total, so this bug becomes a silent double-count instead of an obvious `NaN`. CytoZoo deliberately does not pre-zero component blocks to hide it.
- **A `ClampedCell` on a coupling component only drops that component's own write.** The clamp zeroes `du` inside the component's view, and `_add_back!` then adds the other members' saved contributions, so a state under a contributory `share` (`op = +`) keeps moving. Wrap the `CoupledModel` to hold it. Under hard-discard the question is moot (non-owners are frozen anyway). Rush-Larsen is the other place a clamp cannot work through `du`: `rush_larsen_step!` writes `u_new` directly, so the wrapper restores the held entries after the base step.
- **Switching is composition — there is no switch primitive.** Turn a capability on/off by composing with or without a `Subsystem` node or an edge. No construction-time `switch=` kwarg exists (a prototype was removed as redundant).
- `num_parameters` / `parameter_names` are deliberately **left undefined** on `CoupledModel`, so calling one is a truthful `MethodError` and `hasmethod` reports `false`. Don't "fix" that by adding a throwing method.
- `state_index` / `parameter_index` must return `nothing` for an unknown name — `couple`'s validation depends on that convention.
- Multi-node coupling graphs need explicit `Subsystem(model; name=...)`; the `gensym` default only works for a single node.
- In RHS bodies, watch the name collisions: Faraday's constant is `F_param`, temperature `T_val`, celltype `celltype_val`.
- ToRORd's ~492 ArmyHeart monitors are **not ported** (`TORORD_NUM_MONITORS = 0`); `src/models/torord/monitors.jl` is a stub.

## Working on coupling
Drive coupling-API changes with `examples/coupling_mwe.jl` (+ `examples/coupling_mwe.md`) — a small
toy pair exercising the full coupling taxonomy — **not** the 101-state ECCMitoRedox model.
Rationale for the architecture (why MTK was dropped, why operator splitting was measured and
abandoned for a monolithic RHS) is in `docs/adr/0001-coupling-architecture.md`; the canonical
variable-roles table is in `docs/src/guides/coupling/index.md`; a code-level tour is in
`handoffs/2026-07-16-0745-coupling-infrastructure-tour.md`.

## Docs
Documenter + DocumenterVitepress. Set `VITEPRESS_DEV=1` to skip the Node stage and run Documenter
only — fast iteration that still catches `@example` failures, broken `@ref`s, and `checkdocs=:exports`.
Examples are live `@example` blocks with one namespace per page (Documenter shares state per page,
not across pages). Every export must appear in `docs/src/reference/api.md` or the build fails.
Prose belongs in the docs; the README is a taster only. Keep `docs/adr/` and
`docs/src/reference/design.md` in sync.

## Testing
Correctness is checked against ArmyHeart reference values embedded in
`test/test_torord_correctness.jl` at `rtol=1e-10`. Performance tests assert zero allocations on the
functor. TWorld tests are conditional and skipped when TWorld is unavailable. Cross-validation
sources live at `~/dev/ArmyHeart/` and `~/dev/TWorld/`.
