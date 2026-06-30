---
slug: coupling-monolithic-rhs
created: 2026-06-25-1515
status: open
---

# Handoff: switch coupling from operator-splitting to a monolithic single-RHS

## Goal / why this matters

Make `CoupledModel` a **real functor** — assemble one combined `f(dU, U, p, t)` from the
submodels and solve it with a single ODE solver — instead of integrating it by operator
splitting (OS). A prototype comparison (below) showed monolithic is dramatically more accurate,
zero-alloc, and robust on stiff systems where the OS path diverges. This is the productionization
step: move the prototype's assembly logic into `src` as the default solve path, keeping OS as an
opt-in backend.

## Background & current state

Today a `CoupledModel` is **not** a functor: `(::CoupledModel)(du,u,p,t)` throws
(`src/coupling.jl:342`), and it's solved by stepping each component with its own inner solver,
exchanging data at substep boundaries. Two OS backends exist: `ext/CouplingExt.jl` (drives the
unregistered `OrdinaryDiffEqOperatorSplitting`) and `ext/SplitSolveExt.jl` (a hand-rolled
adaptive Lie–Trotter loop, written *because* OS adaptive sub-solves are "broken for stiff
subsystems" — see its header comment).

A **non-invasive prototype** was built and run: `examples/coupling_mono_vs_os.jl`. It assembles a
single RHS from an existing `cm` using only public-ish internals (`cm.layout.operator_order`,
`cm.layout.solution_indices`, `CytoZoo._frozen_indices`, `CytoZoo._connect_plan`), wraps it in an
`AssembledRHS` callable, and solves via `ODEProblem(f, u0, tspan, nothing)`. `src/` and all tests
are **untouched** — nothing is committed to the design yet.

### Prototype results (the motivation)

1. **Accuracy** — `connect` toy, error vs analytic `acc(2)=1−exp(−2)`, same `Tsit5` + fixed `dt`
   (isolates splitting error):
   - dt=0.1 → OS 4.25e-2 vs monolithic 1.88e-10
   - dt=0.05 → OS 2.14e-2 vs monolithic 4.73e-12
   - dt=0.025 → OS 1.08e-2 vs monolithic 1.31e-13
   - OS error halves with `dt` (textbook O(dt) Lie–Trotter splitting error); monolithic has none,
     so it sits at `Tsit5`'s 5th-order truncation — ~10 orders of magnitude better.
2. **Allocations** — one assembled-RHS call = **0 bytes**, matching the OS per-component
   operators. Zero-alloc on the contiguous-block case (range-views into the shared `dU`/`U`).
3. **Stiffness** — stiff `share` pair (K=500): monolithic `Rodas5P` is **stable**
   (`v(5)=0.6067`, max|v|=10), OS-LTG with explicit `Tsit5` inner → **NaN / −Inf** at dt=0.05 and
   0.01. Monolithic lets you put one A-stable implicit solver on the whole coupled system; OS
   can't.

## Key files / locations

- `src/coupling.jl:342` — the throwing functor to **replace** with the assembled RHS.
  `num_parameters` (`:345`) stays "no single vector".
- `src/coupling.jl` reusable helpers: `_compute_layout` (`:233`), `solution_indices`,
  `_frozen_indices` (`:410`), `_connect_plan` (`:429`), `operator_order` (`_operator_order`,
  `:305`).
- `examples/coupling_mono_vs_os.jl` — the prototype. The `CompEntry` / `AssembledRHS` / `_run!` /
  `_connect!` code is the assembly logic to port into `src` verbatim (it's written to move).
- `ext/SciMLBaseExt.jl` — already defines `ODEProblem(model, tspan; p=nothing)`; once the functor
  exists, `ODEProblem(cm, tspan)` + `solve` works with no new solve code.
- `ext/CouplingExt.jl`, `ext/SplitSolveExt.jl` — keep as the opt-in OS backends. They wrap each
  *component* in `ComponentOperator`, never the `CoupledModel` functor, so adding the functor
  doesn't touch them.
- `test/test_coupling.jl` — contains an assertion that `cm(...)` **throws** `ArgumentError`. This
  test must be **updated** (the functor now works).
- Design doc / approved plan: `/Users/kylebeggs/.claude/plans/it-occured-to-me-radiant-yao.md`
  (full head-to-head comparison and rationale).

## Decisions & conclusions

- Monolithic is the recommended **default** solve path; OS stays as an opt-in backend for the one
  thing it genuinely buys — per-component solver heterogeneity (different method/tolerance per
  subsystem). Most cell-internal coupling is stiff-everywhere, so that's a minor loss.
- The assembly is correct via **overlapping views into the single `dU`/`U` + owner-last
  `operator_order` + frozen-zeroing**. No gather/scatter buffers, no per-component state — the
  earlier "preallocated scratch" worry does not apply. The non-owner writes then zeroes its own
  shared-slot contribution; the owner writes last into the same slot.
- `connect` inputs are read **live from `U`** during the single eval (instantaneous), not staged
  into parameter slots between substeps (lagged) as OS does.
- Tissue scale (Thunderbolt monodomain) is a *later, non-driving* factor — but note: a monolithic
  `CoupledModel` is a callable, so it could serve as a Thunderbolt `AbstractIonicModel`, which the
  throwing OS version structurally cannot. One-way door worth keeping in mind.

## What's left / next steps

1. Port `CompEntry`/`AssembledRHS`/`_run!`/`_connect!` from `examples/coupling_mono_vs_os.jl` into
   `src/coupling.jl`, replacing the throw at `:342` with the assembled functor. Build the
   per-component execution plan **once at `couple()` time** (store it on `CoupledModel`, or
   build it concretely-typed) so the functor specializes and stays zero-alloc — don't rebuild the
   plan per call. Resolve the `couple()`-time chicken-and-egg by refactoring `_frozen_indices` /
   `_connect_plan` to take `(components, shares/connects, ck)` with thin `cm`-based wrappers.
2. Update `test/test_coupling.jl` — the "functor throws" assertion becomes "functor evaluates the
   assembled RHS". Add monolithic-vs-OS cross-validation tests (accuracy at coarse `dt`,
   zero-alloc, stiff stability) mirroring the prototype.
3. **Resolve the `connect`-under-implicit AD issue** (the one real caveat): `connect` writes a
   source state into the receiver's `Float64` `.parameters` slot inside the RHS, which breaks
   ForwardDiff under an implicit solver (a `Dual` can't be stored in a `Float64` slot). `share`
   is unaffected (everything flows through `U`). Options: eltype-generic param slots, *or* pass
   the coupling input through `p` instead of mutating `.parameters`, *or* freeze connect inputs
   within a Newton step.
4. Decide whether `ODEProblem(cm, tspan)` + `solve` is the documented entry point and update
   CLAUDE.md / README coupling sections (currently say "solved by operator splitting … not a
   direct functor").

## Gotchas / constraints

- **Don't break the OS path.** Both ext backends must keep working; they read `cm` fields and
  wrap components, so adding the functor is safe — but run `test/test_coupling_ext.jl` and
  `test/test_splitsolve_ext.jl` after.
- The prototype passes `p=nothing` through to each submodel functor (toy models ignore `p`; real
  models dispatch on `::Nothing`/`::SpatialContext`). Preserve that — `ODEProblem(cm,tspan)`
  defaults `p=nothing` in `SciMLBaseExt`.
- `share` requires **owner-last** ordering for correctness in the monolithic functor (last write
  into the shared slot wins). Use `cm.layout.operator_order`, not declaration order.
- For zero-alloc, derive each component's block as a `UnitRange` when its `solution_indices` are
  contiguous (the common case); fall back to an index-vector view otherwise.
- `examples/` env: `SciMLBase` is only a transitive dep — `using SciMLBase` fails. Use the
  `ODEProblem`/`solve` re-exported by `OrdinaryDiffEq`.
- No secrets involved in this work.
