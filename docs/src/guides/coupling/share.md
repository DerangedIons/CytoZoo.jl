# Share Edges

A `share` edge declares that a state in one subsystem and a state in another are **the same
physical quantity** — one variable, one global slot, one governing equation.

The motivating case: an electrophysiology model and a metabolic model both carry intracellular
ATP. They are not two quantities that happen to track each other; they are one quantity that
both models refer to.

```@setup share
using CytoZoo, OrdinaryDiffEq

struct ModelA <: AbstractCardiacCellModel end
CytoZoo.num_states(::ModelA)                    = 2
CytoZoo.state_names(::ModelA)                   = (:c, :d)
CytoZoo.default_initial_state(::ModelA)         = [5.0, 1.0]
CytoZoo.state_index(::ModelA, n::Symbol)        = findfirst(==(n), (:c, :d))
CytoZoo.transmembrane_potential_index(::ModelA) = 2
function (::ModelA)(du, u, p, t)
    du[1] = -0.2 * u[1]
    du[2] = -u[2]
    return nothing
end

struct ModelB <: AbstractCardiacCellModel end
CytoZoo.num_states(::ModelB)                    = 2
CytoZoo.state_names(::ModelB)                   = (:x, :y)
CytoZoo.default_initial_state(::ModelB)         = [0.0, 0.0]
CytoZoo.state_index(::ModelB, n::Symbol)        = findfirst(==(n), (:x, :y))
CytoZoo.transmembrane_potential_index(::ModelB) = 1
function (::ModelB)(du, u, p, t)
    du[1] = 3.0
    du[2] = u[1] - u[2]
    return nothing
end

struct ModelC <: AbstractCardiacCellModel end
CytoZoo.num_states(::ModelC)                    = 1
CytoZoo.state_names(::ModelC)                   = (:z,)
CytoZoo.default_initial_state(::ModelC)         = [0.0]
CytoZoo.state_index(::ModelC, n::Symbol)        = findfirst(==(n), (:z,))
CytoZoo.transmembrane_potential_index(::ModelC) = 1
function (::ModelC)(du, u, p, t)
    du[1] = -7.0
    return nothing
end
```

## Spelling

```julia
share(:A => :d, :B => :x; owner = :A)
```

Each argument is a `component => state` pair. `owner` must name one of the two components,
and decides which equation survives. An optional `name =` overrides the slot's global name,
which otherwise comes from the owner's state name.

## Hard-Discard Semantics

Only the owner's equation drives the shared slot. Every non-owner's contribution to that
slot's derivative is **zeroed** — the non-owner reads the value, but never writes it.

`ModelA` says `dd/dt = -d`; `ModelB` says `dx/dt = 3.0`. Under a share owned by `A`, only the
first survives:

```@example share
coupled = couple(
    [Subsystem(ModelA(); name = :A),
     Subsystem(ModelB(); name = :B)],
    [share(:A => :d, :B => :x; owner = :A)],
)

sol = solve(ODEProblem(coupled, (0.0, 2.0)), Tsit5(); dt = 0.01, adaptive = false)

(shared = sol.u[end][state_index(coupled, :d)],
 exact  = exp(-2.0),
 note   = "B's dx/dt = 3.0 was discarded")
```

Crucially, `B` still *reads* the shared value while integrating its own states — `y` is
driven by it and follows the owner's trajectory:

```@example share
(y = sol.u[end][state_index(coupled, :B_y)],)
```

This is what makes the authoring story clean: **neither model needs modification**. `ModelB`
was written as a standalone model with its own `x` equation, and it participates in a share
without a single line changing. It simply finds that `x` now moves the way `A` says it does.

## How the Discard Works

Mechanically, two things happen on every evaluation:

1. Each component writes its derivatives into its own view of the shared global `dU`.
2. Non-owner contributions to shared slots are zeroed, and the components are ordered so the
   **owner writes last** — its value is therefore final.

That ordering is computed once, at `couple()` time, and stored on the model. There is no
per-step bookkeeping. Details in [Coupling Internals](../../reference/internals.md).

## Sharing Across Three or More Components

Shares resolve as **equivalence classes**, not pairs. To put one quantity in three
components, fan the edges out from the owner:

```@example share
tri = couple(
    [Subsystem(ModelA(); name = :A),
     Subsystem(ModelB(); name = :B),
     Subsystem(ModelC(); name = :C)],
    [share(:A => :d, :B => :x; owner = :A),
     share(:A => :d, :C => :z; owner = :A)],
)

state_names(tri), num_states(tri)
```

Five local states collapse to three global slots: `c`, the shared `d`, and `B_y`. Both
`B.x` and `C.z` are now the same variable as `A.d`, and both non-owner equations are
discarded:

```@example share
sol3 = solve(ODEProblem(tri, (0.0, 2.0)), Tsit5(); dt = 0.01, adaptive = false)
(shared = sol3.u[end][state_index(tri, :d)], exact = exp(-2.0))
```

The grouping does not depend on the order the nodes were declared in.

!!! warning "Fan out from the owner — do not chain"
    Writing `share(:A => :d, :B => :x; owner = :A)` followed by
    `share(:B => :x, :C => :z; owner = :B)` is rejected. Both edges land in the same
    equivalence class, but they name two different owners, and a shared slot has exactly one
    governing equation. `couple` raises an error at construction rather than silently picking
    one:

```@example share
try
    couple(
        [Subsystem(ModelA(); name = :A),
         Subsystem(ModelB(); name = :B),
         Subsystem(ModelC(); name = :C)],
        [share(:A => :d, :B => :x; owner = :A),
         share(:B => :x, :C => :z; owner = :B)],
    )
catch e
    println(sprint(showerror, e))
end
```

## Naming the Slot

By default the shared slot takes the owner's state name. When neither name is the natural one
for the combined system, set it explicitly:

```@example share
named = couple(
    [Subsystem(ModelA(); name = :A),
     Subsystem(ModelB(); name = :B)],
    [share(:A => :d, :B => :x; owner = :A, name = :ATP)],
)

state_names(named)
```

## Requirements on the Models

To take part in a share, a model must implement [`state_index`](@ref) — that is how `couple`
resolves a name to a position in the component's local state vector. Nothing else is needed,
and no parameter machinery is involved.

Because a share flows entirely through the global state vector `U`, it is **unaffected** by
the automatic-differentiation caveat that applies to [`connect`](@ref) edges. If you need a
coupling that survives sensitivity analysis, prefer a share. See [Limitations](limitations.md).

## Known Gap: Contributions Cannot Sum

Today's share is hard-discard — exactly one equation governs. Some feedback couplings instead
need two models to *each contribute a term* to the same state's derivative: the owner's core
equation plus another module's extra flux. That is not expressible yet. See
[Limitations](limitations.md).

## See Also

- [Connect Edges](connect.md) — the other edge kind.
- [Coupling Internals](../../reference/internals.md) — how the discard and ordering are built.
