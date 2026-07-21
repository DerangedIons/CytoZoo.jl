# Design Notes

Why the architecture looks the way it does. The archival records, with full context and
supersession history, are the Architecture Decision Records in `docs/adr/` in the repository;
this page carries the summary that users and contributors actually need.

## Models Are Hand-Coded Structs, Not Symbolic

CytoZoo previously carried a ModelingToolkit.jl-backed path, with `symbolic_system` /
`has_symbolic_system` interface hooks and an MTK-based model. It was removed.

**Consequence:** models are hand-coded callable structs, and the base package keeps **zero
runtime dependencies**. Everything requiring an external package — solving, tissue coupling,
automatic differentiation — arrives through a package extension.

The trade is deliberate. A symbolic layer buys automatic simplification and code generation,
but costs a heavyweight dependency, compile latency, and control over the exact arithmetic in
the hot loop. For a package whose value is a *registry of validated models with predictable
performance*, hand-written right-hand sides matching a reference implementation to `1e-10`
are worth more than symbolic manipulation.

## Coupling Is Monolithic, Not Operator Splitting

This is the load-bearing decision, and it was made by measurement rather than taste.

### What was tried

Coupling was **originally built on operator splitting**, via
`OrdinaryDiffEqOperatorSplitting.jl`. The reasoning at the time was sound: Lie–Trotter and
Strang splitting have well-established convergence theory, each operator can use its own
inner solver, and shared state lives in one global vector so shared variables genuinely are
the same memory rather than two copies kept in sync by callbacks.

Two splitting backends existed. The second was a hand-rolled adaptive Lie–Trotter loop,
written because the library's adaptive sub-solves were broken for stiff subsystems — itself a
signal about the maintenance burden.

### What the measurements showed

A non-invasive prototype compared splitting against assembling one combined right-hand side:

| Criterion | Operator splitting | Monolithic single-RHS |
|---|---|---|
| Error vs. analytic, `dt = 0.1` | 4.25e-2 | 1.88e-10 |
| Error vs. analytic, `dt = 0.05` | 2.14e-2 | 4.73e-12 |
| Error vs. analytic, `dt = 0.025` | 1.08e-2 | 1.31e-13 |
| Stiff `share` pair (K = 500) | NaN / −Inf with an explicit inner solver | stable under `Rodas5P` |
| Allocations per RHS call | 0 B | 0 B |

The splitting error halves as `dt` halves — textbook first-order Lie–Trotter behaviour. The
monolithic right-hand side has *no splitting error at all*, so its accuracy is bounded only
by the integrator's own truncation order: roughly ten orders of magnitude better on the same
solver at the same step size, at identical allocation cost.

### The decision

[`CoupledModel`](@ref) is a real functor `(cm)(dU, U, p, t)` assembled from its submodels and
integrated by a single ODE solver. The splitting backends were **deleted** rather than kept as
an opt-in path — no consumer justified maintaining a second, strictly less accurate solve
path.

### Consequences

- One solver spans the whole system: no splitting error, and a stiff coupling can use one
  implicit method across every component.
- `OrdinaryDiffEqOperatorSplitting` is not a dependency. The solve path needs only a
  SciMLBase stack, through `ext/SciMLBaseExt.jl`.
- Because every component writes into one shared `dU`, "which equation governs a shared slot"
  is resolved by **ordering plus zeroing** — owner last, non-owner contributions zeroed —
  rather than by a splitting schedule. That is where [`share`](@ref)'s hard-discard semantics
  come from.
- Higher-order (Strang) splitting, once listed as future work, is moot.

## Priorities That Shaped the Coupling API

In order:

1. **Performance** — within roughly 1% of uncoupled component performance, with no
   allocations in the hot path.
2. **A simple *authoring* API** — someone writing a single component writes today's
   `f(du, u, p, t)` and learns nothing new in order to participate in a `share`.
3. **A simple *coupling* API.**

Priority 2 explains why `share` requires no authoring change at all, and why `connect`'s
receiver contract ([`parameter_index`](@ref) plus [`writable_parameters`](@ref)) is the *only*
thing coupling asks a model author to add.

## Switching Is Composition, Not a Primitive

A `switch=` / `switches=` API was prototyped and **removed as redundant**. Composing with or
without a subsystem or an edge already expresses every on/off case, and it makes the
"OFF invariant" — that disabling a module recovers the validated baseline — automatic rather
than something to be tested.

Composition therefore resolves at `couple()` time and never through `p`, which carries
time- and space-varying payload rather than structural configuration. See
[Patterns Cookbook](../guides/coupling/patterns.md).

## Open Questions

Carried forward deliberately, not overlooked:

- **Additive / contributory share.** The current `share` is hard-discard, so two models cannot
  each contribute a term to one shared derivative. This is Modelica's *through* (flow)
  variable, as opposed to the *across* variable `share` implements. Deferred until a consumer
  needs it — and explicitly **not** to be treated as a permanent invariant.
- **GPU compatibility.** Not precluded; the design keeps structs isbits-friendly. Not shipped
  or tested.
- **Non-ODE components.** DAE and algebraic-constraint coupling remain out of scope.
- **Nested coupling.** Supported, since `CoupledModel <: AbstractCardiacCellModel`, but the
  layout-collapsing logic deserves review against a real use case.

See [Coupling Limitations](../guides/coupling/limitations.md) for what these mean in practice.

## See Also

- [Coupling Internals](internals.md) — how the monolithic right-hand side is assembled.
- [Coupling Overview](../guides/coupling/index.md) — the user-facing model.
