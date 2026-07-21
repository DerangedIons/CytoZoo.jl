# Coupling Overview

Cell models rarely stay standalone. An electrophysiology model gets bolted to a mitochondrial
energetics model; a signalling module gets switched on to study its effect. CytoZoo composes
models into one combined system that a single ODE solver integrates.

## Coupling Is a Graph

A coupling is **nodes joined by directed edges**:

- A [`Subsystem`](@ref) is a node — a model plus a `name` used to reference it in edges.
- An edge relates one node's variable to another's. There are two kinds,
  [`share`](@ref) and [`connect`](@ref), and you can mix them freely in one edge list.

[`couple`](@ref) takes the nodes and edges and returns a [`CoupledModel`](@ref):

```julia
coupled = couple(
    [Subsystem(ModelA(); name = :A),
     Subsystem(ModelB(); name = :B)],
    [share(:A => :d, :B => :x; owner = :A),
     connect(:A => :d, :B => :d_ext)],
)
```

Pass an explicit `name` to every node whose edges reference it. The `gensym` default exists
only so a single-node graph is ergonomic.

## One Right-Hand Side, One Solver

The important property: a `CoupledModel` is **a real functor** `(cm)(dU, U, p, t)` assembled
from its submodels. It is not a scheduler that alternates between component solvers.

That means a coupled system is solved monolithically:

```julia
prob = ODEProblem(coupled, (0.0, 1000.0))
sol  = solve(prob, Rodas5P())
```

There is **no operator-splitting error** — the accuracy is whatever the integrator itself
delivers — and a stiff coupling can use one implicit method across every component at once.
CytoZoo originally built coupling on operator splitting, measured it against this approach,
and deleted it; the numbers are in [Design Notes](../../reference/design.md).

Because `CoupledModel <: AbstractCardiacCellModel`, a coupling is itself a model, so
couplings nest.

## Variable Roles

This is the concept the rest of the coupling documentation rests on. Every variable in a
combined system has a **role**, and *coupling changes roles*:

| Role | Meaning | Representation |
|---|---|---|
| **state** | an integrated variable, has a `du/dt` | a global state slot |
| **parameter** | a fixed input | a slot in a model's parameter vector |
| **derived** | an algebraic function of state (e.g. a conservation law) | the monitor hooks |
| **input** | driven live by another model on every evaluation | a [`connect`](@ref) edge |

The coupling primitives *are* the role changes:

- [`share`](@ref) merges two **states** into one global slot.
- [`connect`](@ref) turns a **parameter** into an **input**.
- The monitor hooks surface a **derived** quantity — see [Derived Observables](../monitors.md).
- [`couple`](@ref) composes.

Because the primitives are exactly the role changes, no additional primitive is needed to
express a coupling architecture; you are always choosing which role a variable should play in
the combined system.

## Switching Is Composition

There is **no switch primitive**, no construction-time `switch=` keyword, and no
module-switch protocol. The only way to turn an optional capability on or off is to compose
with or without it — either a whole subsystem, or a single edge:

```julia
nodes = [Subsystem(Core(); name = :C)]
redox_on && push!(nodes, Subsystem(Redox(); name = :R))

coupled = couple(nodes, edges)
```

Omitting the subsystem recovers the baseline exactly, which is what makes an "OFF invariant"
free rather than something you have to test for. Composition resolves at `couple()` time, and
never through `p` — the DiffEq parameter slot carries time- and space-varying payload, not
structural configuration.

A model's *own* internal on/off flag is a private implementation detail of that model, not a
CytoZoo concept.

See [Patterns Cookbook](patterns.md) for module-, edge-, and state-versus-parameter switching
demonstrated live.

## A First Coupled Model

Two toy models. `ModelA` owns a decaying variable `d`; `ModelB` has a state `x` that is
physically the same quantity, plus a `y` driven by it:

```@example coupling
using CytoZoo, OrdinaryDiffEq

struct ModelA <: AbstractCardiacCellModel end
CytoZoo.num_states(::ModelA)                    = 2
CytoZoo.state_names(::ModelA)                   = (:c, :d)
CytoZoo.default_initial_state(::ModelA)         = [5.0, 1.0]
CytoZoo.state_index(::ModelA, n::Symbol)        = findfirst(==(n), (:c, :d))
CytoZoo.transmembrane_potential_index(::ModelA) = 2
function (::ModelA)(du, u, p, t)
    du[1] = -0.2 * u[1]     # c
    du[2] = -u[2]           # d — the shared variable, owned by A
    return nothing
end

struct ModelB <: AbstractCardiacCellModel end
CytoZoo.num_states(::ModelB)                    = 2
CytoZoo.state_names(::ModelB)                   = (:x, :y)
CytoZoo.default_initial_state(::ModelB)         = [0.0, 0.0]
CytoZoo.state_index(::ModelB, n::Symbol)        = findfirst(==(n), (:x, :y))
CytoZoo.transmembrane_potential_index(::ModelB) = 1
function (::ModelB)(du, u, p, t)
    du[1] = 3.0             # x — discarded under the share; A owns the slot
    du[2] = u[1] - u[2]     # y — driven by the shared value
    return nothing
end
nothing # hide
```

Declaring `A.d` and `B.x` to be one variable, owned by `A`:

```@example coupling
coupled = couple(
    [Subsystem(ModelA(); name = :A),
     Subsystem(ModelB(); name = :B)],
    [share(:A => :d, :B => :x; owner = :A)],
)

state_names(coupled), num_states(coupled)
```

Four local states collapsed to three global slots, because `d` and `x` became one.

```@example coupling
sol = solve(ODEProblem(coupled, (0.0, 2.0)), Tsit5(); dt = 0.01, adaptive = false)

(shared_d = sol.u[end][state_index(coupled, :d)], exact = exp(-2.0))
```

`A`'s equation `dd/dt = -d` governs the shared slot, so it decays to `exp(-2)`. `B`'s
`dx/dt = 3.0` was discarded entirely — and `B`'s own `y` still saw the correct shared value
while it integrated.

## Layout and Naming

`couple` computes the global layout once, at construction:

- The **first node is the primary component**. Its states keep their bare names and supply
  the coupled model's transmembrane-potential index.
- Every other component's states are prefixed with the node name — `B`'s `y` becomes `:B_y`.
- A shared slot takes the **owner's** name, unless you pass an explicit `name=` to `share`.

```@example coupling
state_names(coupled)
```

All edge names are validated at `couple()` time, so a typo is a clear error at construction
rather than a wrong answer at solve time.

## Where to Go Next

- [Share Edges](share.md) — merging two states into one.
- [Connect Edges](connect.md) — driving one model from another's state.
- [Patterns Cookbook](patterns.md) — the coupling taxonomy, with what is and is not
  expressible.
- [Limitations](limitations.md) — **read this before doing anything involving derivatives or
  threads.**
