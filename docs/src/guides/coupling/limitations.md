# Coupling Limitations

Everything on this page is a real constraint of the current implementation. Two of them can
produce **silently wrong results**, so read those before building anything serious on
`connect`.

## Derivatives Through a `connect` Edge Are Wrong

This is the sharpest limitation in the package.

**What happens.** Under an implicit solver, ForwardDiff threads dual numbers through the
global state vector `U`. A `connect` edge has to stage a source state into the receiver's
parameter vector, which holds plain `Float64`s — storing a dual there would error. The
`ForwardDiff` extension therefore extracts the *primal* of the dual before writing it, so the
connect input is frozen at its current value for the duration of the Newton step.

**Consequence 1 — approximate Jacobian.** The solver still converges to the correct fixed
point, because the staged value is refreshed on every evaluation. But the Jacobian it uses
omits the coupling term, so it is approximate. For a loosely-coupled system this costs a few
extra Newton iterations. For a *tightly*-coupled stiff `connect`, convergence can degrade
badly.

**Consequence 2 — gradients are incorrect, not merely inaccurate.** Any ForwardDiff
derivative taken *through* a connect edge loses that coupling term entirely. Concretely:

- Sensitivity analysis across a `connect` gives wrong sensitivities.
- Gradient-based parameter fitting across a `connect` optimises the wrong objective gradient.
- The failure is **silent** — no error, no warning, just a number that is not the derivative
  you asked for.

**What to do about it.**

- [`share`](share.md) is completely unaffected, because it flows through `U` like any other
  state. If a coupling can be expressed as a share, and you need derivatives, use a share.
- If you only need to *solve*, `connect` is fine — the trajectory is correct.
- If you need derivatives through a `connect`, the standard fix is to route the receiver's
  parameter vector through a `PreallocationTools.DiffCache` so it can hold duals. **This is
  not built.**

## Thread Safety

A `CoupledModel` stages `connect` inputs into a private copy of each receiver, held on the
coupled model itself. That copy is mutable scratch, so **one `CoupledModel` cannot be solved
concurrently from several threads** — the trajectories would scribble over each other's
staged inputs.

`couple` deepcopies each receiver at construction, so building a coupling never mutates the
model instances you passed in. The hazard is only in sharing one `CoupledModel` across
threads.

For a threaded parameter sweep, give each trajectory its own copy:

```julia
prob_func = (prob, i, repeat) -> remake(prob; f = deepcopy(coupled))
solve(EnsembleProblem(prob; prob_func), Tsit5(), EnsembleThreads(); trajectories = 100)
```

Full in-RHS reentrancy — routing connect inputs through the receiver's per-evaluation `p`
rather than a stored vector — is a tracked follow-up, not a current property.

## Contributions to a Shared State Cannot Sum

Today's [`share`](share.md) is **hard-discard**: exactly one owner's equation governs the
slot, and every other contribution is zeroed.

Some feedback couplings need something different — two models each contributing a *term* to
the same derivative:

```math
\frac{dw}{dt} = \underbrace{f_{\text{core}}(w, \ldots)}_{\text{owner}} \;+\; \underbrace{J_{\text{module}}(\ldots)}_{\text{another subsystem}}
```

In Modelica's vocabulary, the current `share` implements an *across* variable (both sides see
one value), whereas this needs a *through* (flow) variable, where contributions sum. There is
no way to express it today.

Candidate designs, none implemented:

```julia
share(:D => :w, :X => :w_flux; owner = :D, op = +)   # contributory share
inject(:X => :J, :D => :w)                           # a new edge kind
```

!!! note "Do not assume one equation per shared state is permanent"
    "A shared state has exactly one governing equation" is a property of today's
    implementation, not a design invariant. It is expected to change when a consumer needs
    additive coupling, so avoid building anything that depends on it.

## `connect` Sources Must Be States

A `connect` source is resolved with [`state_index`](@ref), so it must name a real state. A
derived quantity — one computed algebraically from state, such as `ADP = C_A - ATP` — cannot
be a connect source, even though it is exactly the kind of thing one model wants to hand
another.

The workaround is to promote the quantity to an actual state, which means integrating a
variable that a conservation law already determines. Supporting derived sources would mean
resolving a monitor value on each evaluation.

## Other Constraints

**No DAE or algebraic-constraint coupling.** Every component must be an ODE. Coupling that
introduces an algebraic constraint between components is out of scope.

**Nested coupling is supported but lightly exercised.** `CoupledModel <:
AbstractCardiacCellModel`, so a coupling can be a component of another coupling, and the
layout logic collapses accordingly. It has not been validated against a demanding real use
case.

**Composition is construction-time only.** Structure is fixed when `couple` returns. Nothing
about the graph can change during a solve, and nothing routes through `p`.

**GPU execution is not shipped.** The design keeps structs isbits-friendly, but no coupled
model has been run or tested on a device.

## Summary

| Limitation | Severity | Affects |
|---|---|---|
| Derivatives through `connect` are wrong | **silent incorrectness** | `connect` only |
| Newton convergence with tight stiff `connect` | performance | `connect` only |
| One `CoupledModel` is not thread-safe | crash / corruption | `connect` receivers |
| Shared-state contributions cannot sum | expressiveness | `share` |
| `connect` sources must be states | expressiveness | `connect` |
| No DAE coupling | scope | all |

## See Also

- [Design Notes](../../reference/design.md) — why the architecture is shaped this way.
- [Coupling Internals](../../reference/internals.md) — the mechanism behind these constraints.
