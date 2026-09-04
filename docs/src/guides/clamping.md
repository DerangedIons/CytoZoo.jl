# Clamping States

A sensitivity or mechanistic experiment usually needs one quantity **held** while the rest of
the cell responds: a voltage clamp, an ion clamp (`[Na]i`, `[Ca]i`), an `[ATP]i` clamp. That
is what a [`ClampedCell`](@ref) does. It wraps any model, calls the base right-hand side, and
zeroes `du` at the held indices.

Holding a state is what separates a *whole-cell* perturbation from an *organelle-local* one.
Raising a parameter that an organelle reads, while leaving the shared ion free, lets the
change leak back out through the sarcolemma and recruit the very pathway the experiment
means to exclude. Clamping the shared state severs that path.

## Freeze Plus Seed

A clamp is two halves that have to agree:

- **the hold** — `du[i] = 0`, so the state cannot move;
- **the seed** — `u[i] = value`, the level it is held at.

A value in `u` without the hold decays back to equilibrium. A hold without the seed pins the
state at baseline. `ClampedCell` carries both, so setting one sets the other:

```@example clamp
using CytoZoo

model = ToRORd()
clamped = ClampedCell(model, (1,), (-20.0,))     # Vm held at -20 mV

u0 = default_initial_state(clamped)              # already seeded
du = similar(u0)
clamped(du, u0, nothing, 0.0)

(vm = u0[1], dvm = du[1])
```

The two-argument form holds a state at the base model's own initial value — the natural
choice for "keep this from drifting" rather than "impose a level":

```@example clamp
default_initial_state(ClampedCell(model, (1,)))[1]
```

## Clamping by Name

[`clamp_states`](@ref) resolves state names and returns the model together with a seeded
state vector:

```@example clamp
c, u = clamp_states(model; v = -20.0, nai = 15.0)
(u[state_index(model, :v)], u[state_index(model, :nai)])
```

Names resolve against [`state_names`](@ref), so a typo names the states that do exist:

```@example clamp
try
    clamp_states(FHNModel(); vm = 0.5)
catch err
    err.msg
end
```

Clamp values are converted to `eltype(u0)`, so a `Float32` model stays `Float32` end to end.

### Multi-Segment Protocols

Pass a state vector to continue from where a previous segment ended. The seed is re-applied
to the state you hand in, so the level can change from segment to segment while the hold
stays:

```julia
c1, u1 = clamp_states(model; nai = 20.0)               # segment 1: elevated
sol1 = solve(ODEProblem(c1, u1, (0.0, 60_000.0), nothing), FBDF())

c2, u2 = clamp_states(model, sol1.u[end]; nai = 7.5)   # segment 2: back to baseline
sol2 = solve(ODEProblem(c2, u2, (0.0, 60_000.0), nothing), FBDF())
```

This is the only time-varying clamp there is. A clamp whose value is a function of `t` was
deliberately not built: segmenting the solve keeps the held level a property of the model
object rather than a hidden schedule inside it, and it matches how CytoZoo handles switching
everywhere else — [compose, don't parameterize](coupling/patterns.md).

## A Clamped State Is Not Conserved

The clamp implicitly injects or removes whatever flux is needed to hold the state, exactly as
a voltage clamp sources current. Mass balance for that species no longer holds — by
construction, not by accident. Fluxes computed *from* the held state remain correct and are
usually the point of the experiment; a budget summed over the held species is not.

## What the Wrapper Forwards

Everything in the interface: states, parameters, named lookup, monitors,
[`writable_parameters`](@ref), and Rush-Larsen. A clamped model is a drop-in for the model it
wraps.

```@example clamp
(num_states(c), transmembrane_potential_index(c), state_index(c, :nai), num_monitors(c))
```

Rush-Larsen needs care and gets it: that scheme writes `u_new` directly instead of going
through `du`, so zeroing derivatives would never reach it. [`rush_larsen_step!`](@ref) on a
`ClampedCell` restores each held state after the base step:

```@example clamp
u_new = similar(u0)
rush_larsen_step!(u_new, u0, nothing, 0.0, 0.01, clamped)
(held = u_new[1] == u0[1], others_moved = u_new[2:end] != u0[2:end])
```

What it cannot forward is a model's *own* API — nothing generic knows those names exist.
Unwrap with [`base_model`](@ref):

```@example clamp
base_model(c).celltype
```

## Clamping and Coupling

A clamp composes with [`couple`](@ref) from either side, and the two positions mean different
things.

**Wrapping a component** holds that component's own states, indexed in its local state
vector. Its seeded initial condition flows into the coupled layout, so there is no global
`u0` to patch by hand.

**Wrapping the [`CoupledModel`](@ref)** holds a state of the whole system, named the way the
coupling names it — a bare name for the primary component, `:<component>_<state>` for the
others.

!!! warning "Under a contributory share, clamp the coupling"
    A component's clamp zeroes only that component's own write. A contributory share
    (`op = +`) adds the other members' derivatives back afterwards, so the shared slot keeps
    moving. Wrap the `CoupledModel` to hold it.

    Under a hard-discard share the distinction does not arise: the owner governs the slot and
    every other member's derivative is discarded, so clamping the owner (or the coupling)
    holds it either way.

```julia
cm = couple(nodes, edges)
held, U0 = clamp_states(cm; mito_cai = 0.1)     # holds the global slot, whoever writes it
```

## See Also

- [The Cell Model Interface](interface.md) — what a wrapper has to forward, and why.
- [Coupling Overview](coupling/index.md) — variable roles, and where a held state sits among
  them.
- [Patterns Cookbook](coupling/patterns.md) — composition as the switching mechanism.
