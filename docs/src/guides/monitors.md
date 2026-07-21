# Derived Observables (Monitors)

Not every quantity you want out of a simulation is a state. Individual ionic currents, fluxes,
and conservation-law partners are **algebraic functions of the state** — integrating them
would be wrong as well as wasteful. CytoZoo surfaces these as *monitors*.

The canonical case is a conservation law. If a pool obeys `ATPm + ADPm = C_A`, then only one
of the two is a state; the other is derived. Making both states would add a redundant equation
and let numerical drift break the invariant.

!!! warning "No shipped model implements monitors yet"
    The hooks below are implemented, tested, and ready to use in your own models, but
    `ToRORd` does not expose any — `num_monitors(ToRORd())` returns `0`. Its roughly 492
    monitors from the ArmyHeart reference implementation have not been ported. Every example
    on this page therefore uses a toy model.

## The Three Hooks

A model opts in by overriding all three:

```julia
num_monitors(model)                   # default: 0
monitor_names(model)                  # default: ()
monitor_values!(mon, u, t, model)     # writes num_monitors(model) values into `mon`
```

The defaults mean a model that opts out still answers gracefully rather than erroring:

```@example mon
using CytoZoo

num_monitors(ToRORd()), monitor_names(ToRORd())
```

Note the argument order of [`monitor_values!`](@ref): the model comes **last**, and there is
no `p` argument. Monitors read parameters from the model struct.

## A Model With Monitors

A two-state pool where `a` is integrated and `b = C - a` is derived:

```@example mon
struct Pool <: AbstractCardiacCellModel
    C::Float64
end

CytoZoo.num_states(::Pool)                    = 2
CytoZoo.state_names(::Pool)                   = (:v, :a)
CytoZoo.default_initial_state(m::Pool)        = [-80.0, 0.8 * m.C]
CytoZoo.state_index(::Pool, n::Symbol)        = findfirst(==(n), (:v, :a))
CytoZoo.transmembrane_potential_index(::Pool) = 1

function (m::Pool)(du, u, p, t)
    du[1] = -0.5 * (u[1] + 80.0)
    du[2] = -0.2 * u[2]
    return nothing
end

# --- the monitor triple ---
CytoZoo.num_monitors(::Pool)  = 2
CytoZoo.monitor_names(::Pool) = (:b, :total)
function CytoZoo.monitor_values!(mon, u, t, m::Pool)
    mon[1] = m.C - u[2]        # the conservation partner
    mon[2] = u[2] + mon[1]     # should equal C at all times
    return nothing
end
nothing # hide
```

```@example mon
pool = Pool(5.0)
num_monitors(pool), monitor_names(pool)
```

## Reading Monitors After a Solve

[`monitor_history`](@ref) recomputes every monitor across a solved trajectory:

```@example mon
using OrdinaryDiffEq

sol = solve(ODEProblem(pool, (0.0, 5.0)), Tsit5())
h   = monitor_history(sol, pool)

(names = h.names, size = size(h.values), n_times = length(h.t))
```

It returns a `NamedTuple` of the time points `t`, the monitor `names`, and a `values` matrix
with **rows = monitors** (in `names` order) and **columns = time points**:

```@example mon
h.values[1, 1:5]      # `b` over the first five saved times
```

The conservation law holds throughout, which is exactly the sort of invariant monitors make
cheap to check:

```@example mon
extrema(h.values[2, :])    # `total` — should be 5.0 everywhere
```

### Why Post-Solve?

Monitors are computed *after* the solve, by walking the saved solution, rather than during
it. That is deliberate: the saved `sol.u` contains plain numbers, never dual numbers, so
`monitor_history` works under any solver without touching the implicit-solver machinery or
the primal-extraction concerns described in
[Coupling Limitations](coupling/limitations.md).

The cost is that monitors are only available at the time points the solver saved.

A model with no monitors yields an empty result rather than an error:

```@example mon
empty_h = monitor_history(sol, ToRORd())
size(empty_h.values)
```

## Monitors on a Coupled Model

A [`CoupledModel`](@ref) aggregates its components' monitors automatically:

- [`num_monitors`](@ref) sums over components.
- [`monitor_names`](@ref) concatenates them, prefixing non-primary components exactly as
  states are prefixed — `:<component>_<name>`.
- [`monitor_values!`](@ref) slices each component's own portion of the global state and hands
  it to that component's implementation.

```@example mon
second = Pool(9.0)

coupled = couple(
    [Subsystem(pool;   name = :P),
     Subsystem(second; name = :Q)],
)

monitor_names(coupled), num_monitors(coupled)
```

```@example mon
csol = solve(ODEProblem(coupled, (0.0, 5.0)), Tsit5())
ch   = monitor_history(csol, coupled)

(names = ch.names, totals = (ch.values[2, end], ch.values[4, end]))
```

Each component's monitors are computed from that component's own state slice, so the two
pools report their own conservation totals — 5.0 and 9.0.

Names and values stay aligned because both walk the components in **declaration order**,
which is not necessarily the order the operators run in.

## Switching

Monitors are not a switching mechanism. To turn an optional capability on or off, compose
with or without the subsystem or edge — see [Patterns Cookbook](coupling/patterns.md).

## See Also

- [Coupling Overview](coupling/index.md) — the variable-role table monitors belong to.
- [Implementing a Model](implementing_a_model.md) — adding monitors to your own model.
