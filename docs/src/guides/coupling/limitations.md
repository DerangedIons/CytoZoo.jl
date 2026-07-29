# Coupling Limitations

Everything on this page is a real constraint of the current implementation. Two of them can
produce **silently wrong results** — derivatives through a `connect` edge, and a model that
accumulates into its own `du` under a contributory `share` — so read those before building
anything serious on either edge kind.

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
coupled model itself, and resolves any monitor sources into scratch vectors held the same way.
That is all mutable scratch, so **one `CoupledModel` cannot be solved concurrently from several
threads** — the trajectories would scribble over each other's staged inputs.

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

## A Contributor's Initial Value Is Discarded

Not a limitation of expressiveness — [`share(…; op = +)`](share.md) sums contributions to a
shared derivative, so "one equation per shared state" is no longer a constraint. What survives
is smaller and easy to trip over.

A contributory member supplies a *term in the owner's equation*, not a variable of its own. The
shared slot's initial value comes from the **owner**, exactly as under hard-discard, so a
contributor's `default_initial_state` entry for that state is ignored. If your module's local
initial value for the contributed state matters, it is not a contribution — it is a separate
state that should not be shared.

## A Component Must Fully Write Its Own `du`

The interface has always expected a component to *write* every derivative slot it owns rather
than accumulate into one. `dU` arrives from the solver dirty, so a model that does
`du[i] += …` without first assigning `du[i]` was already reading garbage.

An accumulating slot sharpens the consequence. Because the coupling adds each member's write
onto a running total, a model that accumulates into its own `du` sees a *plausible partial sum*
there rather than obvious garbage — so the bug turns into a silent double-count instead of an
obvious `NaN`. Assign first, then modify:

```julia
du[1] = f(u, p, t)      # assign — never `du[1] += ...` as the first write
du[1] -= g(u, p, t)     # fine once the slot is initialised
```

CytoZoo deliberately does **not** pre-zero each component's block to paper over this: it would
repair a broken model on its shared states while leaving its private states broken, and add
writes to the hot path.

## Derived Sources Must Be Declared, and Cost the Whole Monitor Vector

A `connect` source may name a state or a monitor, so a derived quantity such as
`ADP = C_A - ATP` is wirable (see [Connect Edges](connect.md)). Three constraints come with it.

**The source model must declare the observable.** A derived quantity that no model exposes
through [`monitor_names`](@ref) cannot be wired. For a model you own this is a one-line
addition; for a third-party model it means wrapping it in an adapter, the same pattern
`ext/ThunderboltExt.jl` uses.

**Wiring one monitor computes them all.** [`monitor_values!`](@ref) fills a model's entire
monitor vector, once per evaluation per sourcing component. A model with a handful of monitors
costs nothing measurable; a model with hundreds pays for all of them to deliver one. Splitting
that into per-monitor access is a possible future addition, not a current property.

**A monitor source cannot also receive a `connect` edge.** The monitor pre-pass runs before the
component walk stages any parameter, so a monitor that read a staged slot would see the
*previous* evaluation's value. Since "does this monitor read that slot" is not knowable
statically, `couple()` rejects the overlap outright rather than risk a silent one-evaluation
lag. Source the monitor from a component that receives no edges, or split the model.

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
| A contributor's initial value is discarded | correctness trap | `share(…; op = +)` |
| A model that accumulates into its own `du` double-counts | **silent incorrectness** | `share(…; op = +)` |
| A derived source costs its model's whole monitor vector | performance | `connect` |
| A monitor source cannot also receive an edge | expressiveness | `connect` |
| No DAE coupling | scope | all |

## See Also

- [Design Notes](../../reference/design.md) — why the architecture is shaped this way.
- [Coupling Internals](../../reference/internals.md) — the mechanism behind these constraints.
