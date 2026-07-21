# Getting Started

This page walks through everything you need to run a cardiac cell model: building one,
evaluating its right-hand side, solving it, reading states and parameters by name, and
switching cell types. Every code block on this page runs during the documentation build,
so the output you see is real.

Install the package first (see [the home page](index.md)), then:

```@example gs
using CytoZoo
```

## Your First Model

A model is a callable struct. Constructing one gives you a fully populated parameter set:

```@example gs
model = ToRORd()
```

`ToRORd` is a modified O'Hara–Rudy ventricular myocyte model with 65 states and 177
parameters. You can ask any model about its own layout:

```@example gs
num_states(model), num_parameters(model)
```

```@example gs
state_names(model)[1:6]
```

The first state is the transmembrane potential. Rather than hard-coding index `1`, ask:

```@example gs
transmembrane_potential_index(model)
```

Every model ships a default initial condition:

```@example gs
u0 = default_initial_state(model)
u0[1:6]
```

## Evaluating the Right-Hand Side

There is no separate `cell_rhs!` function. The model **is** the right-hand side: it is
callable with the standard DifferentialEquations.jl signature `f(du, u, p, t)`.

```@example gs
du = similar(u0)
model(du, u0, nothing, 0.0)
du[1]      # dV/dt in mV/ms
```

The default stimulus is firing at `t = 0`, so this derivative includes the stimulus current
as well as the ionic currents.

The third argument is the DiffEq parameter slot `p`. Passing `nothing` selects the
non-spatial path, where all parameters come from the model struct. Passing a
[`SpatialContext`](@ref) instead lets parameters vary per cell — see
[Spatial Heterogeneity](guides/spatial.md).

## Solving

With a SciMLBase solver stack loaded, you get an `ODEProblem` constructor for free:

```@example gs
using OrdinaryDiffEq

prob = ODEProblem(model, (0.0, 1000.0))
sol  = solve(prob, FBDF(); abstol = 1e-8, reltol = 1e-8)
sol.retcode
```

Cardiac cell models are stiff — the gating variables span timescales from microseconds to
seconds — so an implicit solver is the right default. `FBDF` reaches the same answer as the
explicit `Tsit5` here in roughly 1,600 steps rather than 174,000.

!!! warning "Rosenbrock methods do not work with `ToRORd`"
    `Rodas5P`, `Rodas4`, `Rosenbrock23` and other Rosenbrock methods need a time gradient,
    which they obtain by differentiating the right-hand side with respect to `t`. That calls
    the model with a dual-typed `du` against a plain `u`, and `ToRORd`'s internal right-hand
    side is typed `(du::AbstractVector{T}, u::AbstractVector{T})` — the two must share an
    element type — so the call raises a `MethodError`. Use a BDF or SDIRK method such as
    `FBDF` instead.

`ODEProblem(model, tspan)` uses `default_initial_state(model)` unless you pass `u0`, and
sets `p = nothing` unless you pass a `SpatialContext`.

The default stimulus fires a 1 ms, −53 µA/µF pulse every 1000 ms starting at `t = 0`, so
this single-beat solve captures one full action potential.

## Plotting an Action Potential

Reading the potential trace out of the solution and plotting it:

```@example gs
using CairoMakie

vi = transmembrane_potential_index(model)

fig = Figure(size = (700, 380))
ax  = Axis(fig[1, 1]; xlabel = "time (ms)", ylabel = "Vₘ (mV)",
           title = "ToRORd endocardial action potential")
lines!(ax, sol.t, [u[vi] for u in sol.u]; linewidth = 2)
fig
```

For long multi-beat simulations, the Rush-Larsen exponential integrator is usually a better
fit than an implicit solver — see [Rush-Larsen Integration](guides/rush_larsen.md).

## Named Access

Hard-coding state and parameter indices is fragile. Every model maps names to indices:

```@example gs
state_index(model, :v), state_index(model, :CaMKt)
```

```@example gs
parameter_index(model, :GNa)
```

Parameters live in a plain, mutable, flat vector on the struct, so tuning one is direct
assignment:

```@example gs
gi = parameter_index(model, :GNa)
model.parameters[gi]
```

```@example gs
model.parameters[gi] = 11.0
model.parameters[gi]
```

The flat-vector layout is deliberate: it keeps models isbits-friendly and therefore
GPU-compatible, and it is what lets a `connect` coupling edge write into a parameter slot.
See [The Cell Model Interface](guides/interface.md).

## Changing Cell Type

Ventricular heterogeneity is a constructor keyword — `0` endocardial, `1` epicardial,
`2` midmyocardial:

```@example gs
endo = ToRORd(; celltype = 0)
epi  = ToRORd(; celltype = 1)
mid  = ToRORd(; celltype = 2)

(endo.celltype, epi.celltype, mid.celltype)
```

Because all three satisfy the same interface, one function drives any of them — this is the
hot-swapping property the interface exists to provide:

```@example gs
function peak_potential(m)
    s = solve(ODEProblem(m, (0.0, 400.0)), FBDF(); abstol = 1e-8, reltol = 1e-8)
    vi = transmembrane_potential_index(m)
    return maximum(u[vi] for u in s.u)
end

[peak_potential(m) for m in (endo, epi, mid)]
```

## Element Types

Passing a number type as the first positional argument builds the whole model in that type,
with no `Float64` intermediates leaking into the right-hand side:

```@example gs
m32 = ToRORd(Float32; celltype = 1)
eltype(m32.parameters), eltype(default_initial_state(m32))
```

This matters for GPU work, where `Float32` throughput is much higher. The genericity
convention models must follow is described in
[Implementing a Model](guides/implementing_a_model.md).

## Current Limitations

Worth knowing before you build on this:

- **`ToRORd` is the only model shipped in-tree.** Other models live in separate packages that
  adhere to the interface natively; see the [Model Catalog](reference/models.md).
- **`default_initial_state(ToRORd())` is not a pre-paced steady state.** It sits well away
  from equilibrium — the intrinsic dV/dt there is roughly +240 mV/ms — so the first beat of a
  simulation is not physiologically meaningful. Pace for several beats and take the final
  state as your initial condition before measuring anything. Note this affects the *initial
  condition* only: the right-hand side itself is validated against the ArmyHeart reference
  implementation to a relative tolerance of 1e-10.
- **No shipped model exposes monitors yet.** The derived-observable hooks are implemented and
  tested, but `ToRORd`'s monitor set has not been ported, so `num_monitors(ToRORd()) == 0`.
  See [Derived Observables](guides/monitors.md).
- **GPU execution is designed for but not tested.** The types are isbits-friendly and the
  arithmetic is element-type generic, but no GPU test suite exists.
- **Two coupling patterns are not yet expressible** — additive contributed flux into a shared
  derivative, and `connect` edges sourced from a derived quantity. See
  [Coupling Limitations](guides/coupling/limitations.md).
- **Derivatives taken through a `connect` edge are wrong**, not merely approximate. If you
  plan to do sensitivity analysis or gradient-based fitting across a coupling, read
  [Coupling Limitations](guides/coupling/limitations.md) first.

## Next Steps

- [The Cell Model Interface](guides/interface.md) — the contract every model satisfies, and
  what each tier of it buys you.
- [Coupling Overview](guides/coupling/index.md) — composing several models into one system.
- [Spatial Heterogeneity](guides/spatial.md) — per-cell parameter variation for tissue work.
- [Quick Reference](guides/quickref.md) — the cheat sheet.
