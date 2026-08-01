# Share Edges

A `share` edge declares that a state in one subsystem and a state in another are **the same
physical quantity** — one variable, one global slot, governed by the owner's equation alone or
summed with contributions from the other members.

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
share(:A => :d, :B => :x; owner = :A)            # hard discard (default)
share(:A => :d, :B => :x; owner = :A, op = +)    # contributory
```

Each argument is a `component => state` pair. `owner` must name one of the two components, and
supplies the governing equation, the canonical name and the initial value. `op` says what becomes
of the **non-owner's** derivative — [`overwrite`](@ref) (the default) discards it, `+` adds it to
the owner's. An optional `name =` overrides the slot's global name, which otherwise comes from the
owner's state name.

## Hard-Discard Semantics

This is the default. Only the owner's equation drives the shared slot; every non-owner's
contribution to that slot's derivative is **zeroed** — the non-owner reads the value, but never
writes it. For the additive alternative, see *Contributory Shares* below.

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
per-step bookkeeping. A contributory class works differently — each member saves the slot's
running value across its own write and adds it back, which needs no ordering at all. Details in
[Coupling Internals](../../reference/internals.md).

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

## Contributory Shares: `op = +`

Hard discard is not the only option. Some feedback couplings need two models to *each contribute
a term* to the same state's derivative — the owner's core equation plus another module's extra
flux:

```math
\frac{dw}{dt} = \underbrace{f_{\text{core}}(w, \ldots)}_{\text{owner}} \;+\; \underbrace{J_{\text{module}}(\ldots)}_{\text{contributor}}
```

Pass `op = +` and the non-owner's derivative is **added** instead of dropped. In Modelica's
vocabulary, the default `share` is an *across* variable (both sides see one value) and this is a
*through* (flow) variable, where contributions sum.

`ModelB` says `dx/dt = 3.0`. Owned by `A` with `op = +`, that 3.0 joins `A`'s `dd/dt = -d`
instead of being discarded, so the slot obeys `dd/dt = -d + 3` and relaxes to 3 rather than 0:

```@example share
contrib = couple(
    [Subsystem(ModelA(); name = :A),
     Subsystem(ModelB(); name = :B)],
    [share(:A => :d, :B => :x; owner = :A, op = +)],
)

sol = solve(ODEProblem(contrib, (0.0, 5.0)), Tsit5(); dt = 0.01, adaptive = false)

(shared = sol.u[end][state_index(contrib, :d)],
 exact  = 3 + (1 - 3) * exp(-5.0),
 note   = "B's dx/dt = 3.0 was added, not discarded")
```

Neither model changed. A module contributes a term by carrying it as an ordinary state whose
derivative *is* the flux; whether that term is summed or dropped is decided at the edge.

### What Contributory Shares Give You

**Several contributors sum.** Fan any number of `op = +` edges at one slot and the derivative is
the owner's equation plus every contribution.

**Mixing is fine.** One class may hold the owner, hard-discard members, and contributors at once.
`op` is per edge and describes that edge's non-owner endpoint.

**The sum is cross-sectional, not cumulative.** Every accumulating slot is reset once per
evaluation, so it holds one evaluation's contributions — never a running total across steps.

**No ordering is imposed.** Hard discard needs the owner to write last; an accumulating class
needs no precedence at all, because every member saves and restores the slot around its own
write. Two components may therefore contribute to *each other's* slots — a topology that
hard-discard shares reject as a cyclic operator order.

**The owner still supplies the initial value.** A contributor is a term in someone else's
equation, not a variable of its own, so its own initial value for that state is ignored.

Two things to know before using it: the contributor's initial value being dropped, and the
requirement that a component fully *write* its own `du` rather than accumulate into it. Both are
in [Limitations](limitations.md).

## `gain` — One Quantity, Two Bases

Two models often mean the same physical quantity but express it on different bases: a matrix
volume against a whole-cell volume, mM against µM, a pool normalised differently. A share
asserts the two states *are* one variable, and `gain` is where that conversion lives.

```julia
# The host cell carries NADH on a 10 mM whole-cell pool; the mitochondrial model's kinetics
# are written against its own 1 mM matrix pool.
share(:Cell => :NADH, :Mito => :NADHm; owner = :Cell, op = +, gain = 1/10)
```

**The non-owner reads `gain * slot`. The owner reads the slot untouched.**

**A contributor's derivative is added unscaled.** `gain` converts an input, not a flux. That
asymmetry is deliberate: in practice the two bases differ by a pool normalisation rather than a
unit conversion, and the contributed flux is already in the owner's frame — the contributing
model computes it from rate constants that were fitted there. If you want a true change of
variables, with the contribution scaled by `1/gain` as well, express it in the contributing
model. One keyword cannot mean both, and quietly choosing one would corrupt a calibrated
coupling without failing.

A gain belongs to the edge's non-owner *endpoint*, exactly as `op` does, so every edge naming
that endpoint must declare the same gain — a state reads its slot through one scaling or the
build fails. `gain` must be finite and non-zero; to carry nothing across an edge, omit the edge.

`gain = 1` is the default and costs nothing: no gained slots are recorded and the
save/scale/restore compiles away.

Unlike a [`connect`](connect.md) input, a gained share still flows through the global state
vector, so the coupling term survives in the Jacobian under ForwardDiff.

## See Also

- [Connect Edges](connect.md) — the other edge kind.
- [Coupling Internals](../../reference/internals.md) — how the discard and ordering are built.
