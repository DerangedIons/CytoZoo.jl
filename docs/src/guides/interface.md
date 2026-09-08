# The Cell Model Interface

CytoZoo is, above all, an *interface*. Models are hand-coded callable structs that satisfy a
common contract, which is what makes them interchangeable, composable, and usable from any
SciML solver without an adapter layer.

## Functor-First Design

There is no `cell_rhs!` function to implement. **The model is the right-hand side.** A model
is a callable struct with exactly the signature DifferentialEquations.jl expects:

```julia
(model)(du, u, p, t) -> Nothing
```

This has a practical consequence: a model can be handed straight to `ODEProblem` with no
wrapper, and `p` — the DiffEq parameter slot — stays free to carry per-evaluation payload.

Models dispatch on `p` rather than branching on it:

```julia
(model::ToRORd)(du, u, ::Nothing, t)           # non-spatial: parameters from the struct
(model::ToRORd)(du, u, p::SpatialContext, t)   # spatial: per-cell variation through p
```

Because this is dispatch and not a runtime check, the non-spatial path compiles with every
spatial branch removed. Using `p = nothing` costs nothing at all relative to a model that
never had spatial support. See [Spatial Heterogeneity](spatial.md).

## The Contract Is Tiered

A flat "required vs. optional" split under-describes this interface, because what a model
must implement depends on what you intend to do with it. There are four tiers.

### Core — every model

```julia
(model)(du, u, p, t)                    # the ODE right-hand side
default_initial_state(model)            # initial condition vector
state_names(model)                      # tuple of state names, in state-vector order
transmembrane_potential_index(model)    # index of Vₘ in the state vector
```

| Method | Returns | Purpose |
|---|---|---|
| functor | `Nothing` | writes derivatives into `du` |
| [`default_initial_state`](@ref) | `Vector` | length must match the state count |
| [`state_names`](@ref) | `NTuple{N, Symbol}` | names in order; drives coupling layout |
| [`transmembrane_potential_index`](@ref) | `Int` | so callers never hard-code `1` |

### Atomic models — anything with a flat parameter vector

```julia
num_states(model)
num_parameters(model)
parameter_names(model)
```

A **composite** model — one built by [`couple`](@ref) — deliberately leaves
[`num_parameters`](@ref) and [`parameter_names`](@ref) *undefined*. Its parameters live on
its individual components; there is no single flat vector to return. Calling one therefore
raises a `MethodError`, and `hasmethod` correctly reports `false`:

```@example iface
using CytoZoo
hasmethod(num_parameters, Tuple{CoupledModel})
```

This is a deliberate design choice. A method that exists only to throw would lie to
`hasmethod` and to any code doing capability detection; leaving it undefined tells the truth.

### Coupling participants — only to take part in a `couple`

```julia
state_index(model, :v)          # share owners and connect sources
parameter_index(model, :GNa)    # connect receivers: names the target slot
writable_parameters(model)      # connect receivers: defaults to model.parameters
```

By convention [`state_index`](@ref) and [`parameter_index`](@ref) return `nothing` for an
unknown name rather than throwing.

This matters most for [`parameter_index`](@ref): [`couple`](@ref) validates a `connect`
target's slot by testing `parameter_index(receiver, slot) === nothing`, so a receiver whose
implementation *throws* on an unknown name produces a bare `KeyError` instead of the intended
"no parameter slot `:x` on component `:B`" message. State names are validated against
[`state_names`](@ref) instead, so `state_index` is less load-bearing here — but implement both
the same way.

`couple` never calls `num_states` or `num_parameters` on a component — the global layout is
derived from `state_names` plus `default_initial_state`. That is why the atomic tier is not
required for participation.

### Optional

```julia
has_rush_larsen(model)                        # default: false
rush_larsen_step!(u_new, u, p, t, dt, model)

num_monitors(model)                           # default: 0
monitor_names(model)                          # default: ()
monitor_values!(mon, u, t, model)
```

The two optional groups are all-or-nothing. Rush-Larsen means implementing
[`rush_larsen_step!`](@ref) *and* returning `true` from [`has_rush_larsen`](@ref); monitors
mean overriding all three hooks. The defaults are chosen so a model that opts out still
answers the query gracefully — [`monitor_history`](@ref) on a zero-monitor model returns an
empty result rather than erroring.

### Wrappers

A model can also be *wrapped*. A wrapper satisfies the contract by forwarding it to the model
inside, and changes one thing — [`ClampedCell`](@ref) forwards everything and zeroes `du` at
the states it holds. Because the tiers above are all a wrapper needs, one implementation
covers every model in the zoo, including a [`CoupledModel`](@ref). See
[Clamping States](clamping.md).

A wrapper can forward only what the interface names. Reach a model's own accessors through
[`base_model`](@ref), which unwraps to the model underneath and is the identity on anything
that is not a wrapper.

## Hot-Swapping

The payoff of the contract is that a function written against it drives any model. Here is a
minimal second model — a two-state toy — alongside `ToRORd`:

```@example iface
struct Decay <: AbstractCardiacCellModel
    rate::Float64
end

CytoZoo.num_states(::Decay)                     = 2
CytoZoo.state_names(::Decay)                    = (:v, :w)
CytoZoo.default_initial_state(::Decay)          = [-85.0, 1.0]
CytoZoo.transmembrane_potential_index(::Decay)  = 1

function (m::Decay)(du, u, p, t)
    du[1] = -m.rate * (u[1] + 85.0)
    du[2] = -m.rate * u[2]
    return nothing
end
nothing # hide
```

One function, written once against the interface, summarises either:

```@example iface
function describe(m)
    u = default_initial_state(m)
    vi = transmembrane_potential_index(m)
    return (states = length(state_names(m)), resting_Vm = u[vi])
end

describe.((ToRORd(), Decay(0.5)))
```

Nothing in `describe` knows which model it received. Swapping models is a one-line change at
the call site, and packages that adhere to the interface natively — such as
[TWorld.jl](https://github.com/DerangedIons/TWorld.jl), whose `TWorldCellModel` subtypes
`CytoZoo.AbstractCardiacCellModel` directly — drop into the same function:

```julia
using CytoZoo, TWorld

describe(TWorldCellModel(; celltype = 0, sex = 2))
```

## Type Hierarchy

```
AbstractCellModel
└── AbstractCardiacCellModel
    ├── ToRORd
    ├── CoupledModel
    └── (models defined in other packages)
```

Subtype [`AbstractCardiacCellModel`](@ref) for cardiac electrophysiology models.
[`AbstractCellModel`](@ref) is the root, reserved for cell models that are not cardiac
electrophysiology — nothing in the coupling machinery requires the cardiac subtype.

Note that [`CoupledModel`](@ref) is itself an `AbstractCardiacCellModel`. A coupling is a
model, so couplings nest.

## Where to Go Next

- [Implementing a Model](implementing_a_model.md) — writing a new model against this contract.
- [Coupling Overview](coupling/index.md) — what the coupling tier unlocks.
- [API Reference](../reference/api.md) — full signatures.
