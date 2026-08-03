# Connect Edges

A `connect` edge is a **directed dataflow link**: one component's state is written into
another component's parameter slot before that component's equations run, so the receiver
reads it as a live input.

Where [`share`](@ref) says "these are the same variable", `connect` says "this model *reads*
that model". It is the right tool when a quantity genuinely belongs to one model and merely
influences another — extracellular potential driving a cell, or a cytosolic concentration
that a mitochondrial module responds to but does not own.

```@setup connect
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

struct Hold <: AbstractCardiacCellModel end
CytoZoo.num_states(::Hold)                    = 1
CytoZoo.state_names(::Hold)                   = (:z,)
CytoZoo.default_initial_state(::Hold)         = [2.0]
CytoZoo.state_index(::Hold, n::Symbol)        = findfirst(==(n), (:z,))
CytoZoo.transmembrane_potential_index(::Hold) = 1
(::Hold)(du, u, p, t) = (du[1] = 0.0; nothing)
```

## The Receiver Contract

This is the one authoring change coupling asks for. To *receive* a connect edge, a model must:

1. Implement [`parameter_index`](@ref) so the target slot can be named.
2. Expose [`writable_parameters`](@ref) — the vector the coupling writes into, which the
   model's own functor reads. It defaults to `model.parameters`, so a model with a
   `parameters` field gets this for free.

```@example connect
struct Reader <: AbstractCardiacCellModel
    parameters::Vector{Float64}     # slot 1 = :d_ext, written each evaluation
end
Reader() = Reader([0.0])

CytoZoo.num_states(::Reader)                    = 1
CytoZoo.state_names(::Reader)                   = (:acc,)
CytoZoo.default_initial_state(::Reader)         = [0.0]
CytoZoo.state_index(::Reader, n::Symbol)        = findfirst(==(n), (:acc,))
CytoZoo.transmembrane_potential_index(::Reader) = 1
CytoZoo.parameter_index(::Reader, n::Symbol)    = n === :d_ext ? 1 : nothing

# `acc` integrates whatever the coupling stages into slot 1
(m::Reader)(du, u, p, t) = (du[1] = m.parameters[1]; nothing)
nothing # hide
```

Note `parameter_index` returns `nothing` for an unknown name. That convention is what lets
`couple` report a mistyped slot as an actionable error.

Both requirements are checked at `couple()` time, not at solve time.

## A First Connect

`Reader.acc` integrates `A.d`, which decays as `exp(-t)`:

```@example connect
coupled = couple(
    [Subsystem(ModelA(); name = :A),
     Subsystem(Reader(); name = :R)],
    [connect(:A => :d, :R => :d_ext)],
)

sol = solve(ODEProblem(coupled, (0.0, 2.0)), Tsit5(); dt = 0.01, adaptive = false)

(acc   = sol.u[end][state_index(coupled, :R_acc)],
 exact = 1 - exp(-2.0))
```

The source value is read **live from the global state vector** on every evaluation, not
cached from the start of a step, so the receiver always sees the current value.

Note also that no state was merged: `connect` leaves the state count alone, unlike `share`.

```@example connect
state_names(coupled), num_states(coupled)
```

## Operations: `overwrite` and `+`

Every edge carries an operation deciding how the value lands in the slot.

[`overwrite`](@ref) is the default — a plain copy:

```julia
connect(:A => :d, :R => :d_ext)                 # same as op = overwrite
```

`+` sums every contributing edge into one slot. The slot is reset to zero and re-summed on
each evaluation, so this is a **cross-sectional sum over sources**, not a running total over
time:

```@example connect
summed = couple(
    [Subsystem(ModelA(); name = :A),
     Subsystem(Hold();   name = :H),
     Subsystem(Reader(); name = :R)],
    [connect(:A => :d, :R => :d_ext; op = +),
     connect(:H => :z, :R => :d_ext; op = +)],
)

sol2 = solve(ODEProblem(summed, (0.0, 2.0)), Tsit5(); dt = 0.01, adaptive = false)

# slot = exp(-t) + 2, so acc = (1 - exp(-2)) + 4
(acc = sol2.u[end][state_index(summed, :R_acc)],
 exact = (1 - exp(-2.0)) + 4.0)
```

The rules for combining edges into one slot are strict, and violations are rejected at
`couple()` time rather than producing a silently wrong answer:

| Edges targeting one slot | Result |
|---|---|
| a single `overwrite` | fine |
| several `+` | fine — all summed |
| several `overwrite` | **rejected** — which one wins is ambiguous |
| `overwrite` mixed with `+` | **rejected** — ordering would decide the answer |
| any other operation | **rejected** — only `overwrite` and `+` are supported |

Because edges are partitioned by operation once at construction, the per-evaluation write has
no dynamic dispatch and no allocation.

## Derived Sources

A connect source does not have to be a state. It can also be a **monitor** — a derived quantity
the source model computes algebraically from its own state (see
[Derived Observables](../monitors.md)).

This is what wires a quantity that a conservation law determines. A model tracking a conserved
pool integrates one side and derives the other: cytosolic ADP is not a state anywhere, it is
`C_A - ATP`. Declaring it as a monitor makes it wirable:

```@example connect
const C_A = 8.0

struct Pool <: AbstractCardiacCellModel end
CytoZoo.num_states(::Pool)                    = 1
CytoZoo.state_names(::Pool)                   = (:atp,)
CytoZoo.default_initial_state(::Pool)         = [3.0]
CytoZoo.state_index(::Pool, n::Symbol)        = findfirst(==(n), (:atp,))
CytoZoo.transmembrane_potential_index(::Pool) = 1
(::Pool)(du, u, p, t) = (du[1] = -u[1]; nothing)

# `adp` is DERIVED — the law lives here, in the model that owns the pool
CytoZoo.num_monitors(::Pool)  = 1
CytoZoo.monitor_names(::Pool) = (:adp,)
CytoZoo.monitor_values!(mon, u, t, ::Pool) = (mon[1] = C_A - u[1]; nothing)

derived = couple(
    [Subsystem(Pool();   name = :P),
     Subsystem(Reader(); name = :R)],
    [connect(:P => :adp, :R => :d_ext)],     # a monitor, not a state
)

state_names(derived), num_states(derived)
```

Wiring `adp` added no state — it is recomputed from `atp` on every evaluation, so it always
agrees with the pool. The alternative, promoting it to an integrated variable with
`d(adp)/dt = -d(atp)/dt`, would add a state the law already determines and let round-off drift
break the closure.

Keeping the law inside the model matters for more than tidiness. It reads that model's live
parameters, so changing `C_A` changes the wire; a constant restated at the edge would silently
go stale. And the same declaration serves both roles — `monitor_history` plots `adp` post-solve
whether or not anything is wired to it.

A monitor can depend on several states, on `t`, and on the model's parameters — anything
`monitor_values!` can compute.

Three things to know:

- **Resolution order.** Names resolve against `state_names` first, then `monitor_names`. A
  component declaring the same name in both is rejected at `couple()` time as ambiguous.
- **Cost.** `monitor_values!` computes a model's *whole* monitor vector, once per evaluation per
  sourcing component. Wiring one monitor of a model with hundreds pays for all of them.
- **A monitor source cannot also receive an edge.** Monitors are computed in a single pass
  before any parameter is staged, so a monitor reading a staged slot would see the previous
  evaluation's value. `couple()` rejects the overlap rather than lagging silently.

`share` is unaffected by any of this: it merges *states*, and a monitor has no derivative to
own. Naming one as a share endpoint is an error.

## Receivers Are Copied

`couple` **deepcopies every connect receiver**. Staging an input mutates a private copy, never
the model instance you passed in:

```@example connect
r = Reader()
cm = couple([Subsystem(ModelA(); name = :A), Subsystem(r; name = :R)],
            [connect(:A => :d, :R => :d_ext)])

solve(ODEProblem(cm, (0.0, 2.0)), Tsit5(); dt = 0.01, adaptive = false)

r.parameters[1]     # still 0.0 — your instance was not touched
```

The flip side is that a single `CoupledModel` owns mutable scratch, so it is **not safe to
solve concurrently from several threads**. See [Limitations](limitations.md).

## Scaling a Socket with `gain`

A socket often needs a unit or basis conversion — the two models agree on what the quantity
*is* and disagree on how it is expressed. `gain` keeps that conversion at the edge instead of
restating it inside either model:

```julia
# The host's whole-cell respiration drives the mitochondrial model's matrix-basis ROS shunt;
# 0.18 is the respiratory-complex density ratio between the two frames.
connect(:Cell => :V_O2, :Mito => :VNO_ext; gain = 0.18)
```

The source value is scaled before it enters the slot. With `op = +`, each edge is scaled by its
own gain and then summed. `gain` must be finite and non-zero; `gain = 1` is the default and is
free.

## Two Serious Caveats

Before using `connect` in anger, read [Limitations](limitations.md). In brief:

1. Under an implicit solver, a connect input is frozen to its primal value within each Newton
   step. The fixed point is correct but the Jacobian is approximate, so a tightly-coupled
   stiff connect may converge poorly.
2. **Any derivative taken through a connect edge silently loses the coupling term.**
   Sensitivity analysis and gradient-based parameter fitting across a `connect` produce wrong
   answers, not merely inaccurate ones.

Neither applies to [`share`](@ref), which flows through the state vector.

## See Also

- [Share Edges](share.md) — merging states instead of piping values.
- [Patterns Cookbook](patterns.md) — which real coupling patterns map to `connect`.
- [Limitations](limitations.md) — the caveats above, in full.
