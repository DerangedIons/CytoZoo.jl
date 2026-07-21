# ADR 0001 — Coupling architecture: monolithic single-RHS, not operator splitting

**Status:** accepted · **Date:** 2026-07-20 · **Supersedes:** `coupling-redesign.md` (deleted)

Architecture Decision Records live outside `docs/src/`, so the ADR files themselves are not
part of the rendered documentation site. This ADR is the **archival record** — full context,
supersession history, and the reasoning as it stood at decision time. A user-facing summary of
the same decisions is published on the site at `docs/src/reference/design.md`; keep the two in
sync when a decision changes.

## Context

CytoZoo was repositioned from "a registry of cardiac cell models" to "an interface package for
cell models **and for coupling them together**." Two decisions drove the coupling design, and
both had their rationale recorded only in `coupling-redesign.md`, which this branch deletes.
This ADR preserves it.

The driving example: model A has states `(b, c, d)`, model B has `(x, y)`, and physically
`d ≡ x`. The user keeps one side's equation and discards the other, and the result must solve
with the standard `f(du, u, p, t)` form. Read-only cross-references (B reads A's `Vm`) are also
in scope.

Priorities, in order: **(1)** high performance — within ~1% of uncoupled component performance,
no allocations in the hot path; **(2)** a simple *authoring* API — someone writing a single
component writes today's `f(du, u, p, t)` and learns nothing new for the share case;
**(3)** a simple *coupling* API.

## Decision 1 — Drop ModelingToolkit.jl

The MTK-backed `BeelerReuter` model and the `MTKCardiacCellModelsExt` extension were removed,
along with the `symbolic_system` / `has_symbolic_system` interface hooks, which lost their only
consumer. Removal scope: `ext/MTKCardiacCellModelsExt/`, the `MTKCardiacCellModels` and
`ModelingToolkit` weakdeps/extensions entries, the `BeelerReuter` stub and export, and
`test/test_mtk_ext.jl`.

**Consequence:** CytoZoo models are hand-coded callable structs, and the package keeps zero
runtime dependencies. Nothing in the current tree records that MTK was ever a design axis, which
is precisely why it is written down here.

## Decision 2 — Operator splitting was built, evaluated, and abandoned

Coupling was **originally built on operator splitting** (OS), via
`OrdinaryDiffEqOperatorSplitting.jl`. The reasoning at the time: Lie–Trotter and Strang splitting
have well-established convergence theory, each operator can use its own inner solver, and shared
state lives in a single global vector so shared variables *are* the same memory (no callback
dance keeping two copies consistent). `PartitionedIntegrators.jl` was considered as an
alternative and rejected as a worse fit for intra-cell sub-model coupling, where shared variables
are literally the same physical quantity.

Two OS backends existed: `ext/CouplingExt.jl` (driving the unregistered
`OrdinaryDiffEqOperatorSplitting`) and `ext/SplitSolveExt.jl`, a hand-rolled adaptive Lie–Trotter
loop written *because* OS adaptive sub-solves were broken for stiff subsystems.

**This was abandoned** after a non-invasive prototype (`examples/coupling_mono_vs_os.jl`)
compared it against assembling one combined RHS. The measured results:

| Criterion | Operator splitting | Monolithic single-RHS |
|---|---|---|
| Error vs analytic, `dt = 0.1` | 4.25e-2 | 1.88e-10 |
| Error vs analytic, `dt = 0.05` | 2.14e-2 | 4.73e-12 |
| Error vs analytic, `dt = 0.025` | 1.08e-2 | 1.31e-13 |
| Stiff `share` pair (K = 500) | NaN / −Inf with explicit inner solver | stable under `Rodas5P` |
| Allocations per RHS call | 0 B | 0 B |

OS error halves with `dt` — textbook O(dt) Lie–Trotter splitting error. The monolithic RHS has
*no* splitting error at all, so it sits at the integrator's own truncation order: roughly ten
orders of magnitude better on the same solver and step size, at identical allocation cost.

**Decision:** `CoupledModel` is a real functor `(cm)(dU, U, p, t)` assembled from its submodels
and solved by a single ODE solver. The OS backends (`ext/CouplingExt.jl`,
`ext/SplitSolveExt.jl`, `test/test_coupling_ext.jl`) were deleted rather than kept as an opt-in
path — no consumer justified maintaining a second, strictly less accurate solve path.

## Consequences

- One solver over the whole system: no splitting error, and a stiff coupling can use one
  implicit method across all components.
- `OrdinaryDiffEqOperatorSplitting` is not a dependency. The base package stays dependency-free;
  the solve path needs only a SciMLBase stack, via `ext/SciMLBaseExt.jl`.
- Because every component's derivative is written into one shared `dU`, the "which equation
  governs a shared slot" question is resolved by *ordering* (owner last) plus zeroing the
  non-owners' contributions, rather than by a splitting schedule.
- Higher-order (Strang) splitting, listed as future work in the original design, is moot.

## Still-open questions carried forward from the original design

- **Additive / contributory share.** The current `share` is hard-discard: one owner governs, and
  two models cannot each contribute a term to the same shared derivative. This is Modelica's
  *through* (flow) variable, where contributions sum, as opposed to the *across* variable that
  `share` implements. Deferred until the ECCMitoRedox port needs it; see socket 5 in
  `examples/coupling_mwe.md`.
- **GPU compatibility.** Not precluded — the design keeps structs isbits-friendly — but not
  shipped or tested.
- **Non-ODE components.** DAE / algebraic-constraint coupling remains out of scope.
- **Nested coupling.** Supported (`CoupledModel <: AbstractCardiacCellModel`), but the layout
  collapsing logic deserves review against a real use case.

## References

- Full pre-deletion text: `git show main:coupling-redesign.md` and
  `git show 25512c6:handoffs/2026-06-25-1515-coupling-monolithic-rhs.md`.
- Code-level tour: `handoffs/2026-07-16-0745-coupling-infrastructure-tour.md`.
- Coupling taxonomy and API backlog: `examples/coupling_mwe.md`.
