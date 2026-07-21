# Rush-Larsen Integration

Rush-Larsen is the standard time-stepping scheme for cardiac cell models, and CytoZoo models
can implement it as an optional part of the interface. This page explains what it buys you
and how to drive it.

## Why Not Just Use an Implicit Solver?

Cardiac cell models are stiff, which normally argues for an implicit method. But they are
stiff in a specific, exploitable way: most of the stiffness lives in the **gating variables**,
each of which obeys a linear relaxation equation

```math
\frac{dg}{dt} = \frac{g_\infty(V) - g}{\tau_g(V)}
```

Holding `V` fixed over a step, that equation has an exact solution. Rush-Larsen integrates
each gate with that exponential update instead of a generic numerical one:

```math
g_{n+1} = g_\infty(V_n) + \bigl(g_n - g_\infty(V_n)\bigr)\, e^{-\Delta t / \tau_g(V_n)}
```

The result is stable at step sizes an explicit method could never take, while remaining fully
explicit — no Newton iteration, no Jacobian, no linear solve.

That last point is the practical argument. An implicit solver drives the right-hand side with
trial states produced by Newton iteration, and those trial states need not be physical. The
ToRORd right-hand side takes logarithms and fractional powers of physiological quantities —
concentrations, in particular — which are positive in any real state but can go negative in a
Newton iterate. Rush-Larsen never evaluates the model off the physical manifold, so the
problem cannot arise.

## Checking for Support

Rush-Larsen is optional, so query before you use it:

```@example rl
using CytoZoo

model = ToRORd()
has_rush_larsen(model)
```

[`has_rush_larsen`](@ref) defaults to `false`, so a model that has not implemented the scheme
answers correctly rather than erroring.

## Taking a Step

The signature puts the model last, matching the mutating-first convention of the rest of the
interface:

```julia
rush_larsen_step!(u_new, u, p, t, dt, model)
```

`u_new` is written in place; `p` is `nothing` for the non-spatial path or a
[`SpatialContext`](@ref) for per-cell variation, exactly as in the functor.

```@example rl
u     = default_initial_state(model)
u_new = similar(u)

rush_larsen_step!(u_new, u, nothing, 0.0, 0.01, model)
u_new[1]
```

Unlike the functor, which writes *derivatives*, this writes the **next state**. Stepping
forward means swapping the buffers.

## A Full Simulation

Rush-Larsen has no adaptive controller — you choose `dt` and march. For ToRORd, `dt = 0.01`
ms is the conventional choice:

```@example rl
vm_idx = transmembrane_potential_index(model)
dt     = 0.01
n_steps = round(Int, 2000.0 / dt)

u     = default_initial_state(model)
u_new = similar(u)

t_hist = Vector{Float64}(undef, n_steps + 1)
v_hist = Vector{Float64}(undef, n_steps + 1)
t_hist[1] = 0.0
v_hist[1] = u[vm_idx]

for i in 1:n_steps
    t = (i - 1) * dt
    rush_larsen_step!(u_new, u, nothing, t, dt, model)
    u .= u_new
    t_hist[i + 1] = t + dt
    v_hist[i + 1] = u[vm_idx]
end

extrema(v_hist)
```

```@example rl
using CairoMakie

fig = Figure(size = (800, 380))
ax  = Axis(fig[1, 1]; xlabel = "time (ms)", ylabel = "Vₘ (mV)",
           title = "ToRORd — 2 beats, Rush-Larsen at dt = $dt ms")
lines!(ax, t_hist, v_hist; linewidth = 1.5)
fig
```

The two beats are driven by the model's own periodic stimulus at a 1000 ms basic cycle
length; nothing external paces it. See [Stimulus](stimulus.md).

!!! note "The first beat is not physiological"
    `default_initial_state(ToRORd())` is not a pre-paced steady state, so the opening
    transient is an artefact of the initial condition. Discard the first several beats before
    measuring action potential duration or any other biomarker.

## Choosing Between the Two

| | Rush-Larsen | Implicit solver (`FBDF`) |
|---|---|---|
| Step control | fixed `dt` you choose | adaptive, error-controlled |
| Cost per step | one right-hand-side evaluation | Jacobian + linear solve |
| Robustness | never leaves the physical manifold | Newton iterates may go non-physical |
| Interpolation | none — you get what you store | dense output via `sol(t)` |
| Best for | long pacing runs, tissue coupling, GPU | single beats, tight tolerances, stiff couplings |

Rush-Larsen is also what makes an isbits stimulus worthwhile: with no allocation and no
solver state, the whole step is GPU-kernel shaped.

For coupled systems, note that [`couple`](@ref) builds a `CoupledModel` whose functor is a
standard right-hand side — coupled systems are solved through the ODE path, not through
Rush-Larsen. See [Coupling Overview](coupling/index.md).

## See Also

- [Getting Started](../getting_started.md) — the implicit-solver path, and why Rosenbrock
  methods do not work with `ToRORd`.
- [Implementing a Model](implementing_a_model.md) — adding Rush-Larsen to your own model.
