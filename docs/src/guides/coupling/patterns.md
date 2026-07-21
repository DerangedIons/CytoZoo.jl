# Patterns Cookbook

This page catalogues the coupling patterns that arise when composing real cell models, and
shows how each maps onto CytoZoo's primitives. Every pattern marked ✅ runs live below; the
two marked ❌ are not expressible yet and are covered in [Limitations](limitations.md).

| Pattern | What changes role | Primitive | Status |
|---|---|---|---|
| Feedforward wire | state → parameter | `connect` | ✅ |
| Adopt-native / drop duplicate state | state + state → one state | `share` | ✅ |
| Module on/off | whole subsystem included or omitted | composition | ✅ |
| Edge on/off | single edge included or omitted | composition | ✅ |
| State ↔ parameter flip | parameter ⇄ state | composition | ✅ |
| Derived observable | state → derived | monitor hooks | ✅ |
| Additive contributed flux | two equations into one derivative | — | ❌ |
| Derived-source wire | derived → parameter | — | ❌ |

## The Running Example

A driver `D` producing two signals, and a responder `R` that reads them:

```@example patterns
using CytoZoo, OrdinaryDiffEq

# Driver: a fast signal `u` and a slow one `v`
struct Driver <: AbstractCardiacCellModel end
CytoZoo.num_states(::Driver)                    = 2
CytoZoo.state_names(::Driver)                   = (:u, :v)
CytoZoo.default_initial_state(::Driver)         = [1.0, 4.0]
CytoZoo.state_index(::Driver, n::Symbol)        = findfirst(==(n), (:u, :v))
CytoZoo.transmembrane_potential_index(::Driver) = 1
function (::Driver)(du, u, p, t)
    du[1] = -u[1]           # u = exp(-t)
    du[2] = -0.1 * u[2]     # v, slow
    return nothing
end

# Responder: reads `p_u` (slot 1) and `h` (slot 2) as inputs
struct Responder <: AbstractCardiacCellModel
    parameters::Vector{Float64}
end
Responder() = Responder([0.0, 1.0])      # p_u = 0, h held at 1.0

CytoZoo.num_states(::Responder)                    = 1
CytoZoo.state_names(::Responder)                   = (:y,)
CytoZoo.default_initial_state(::Responder)         = [0.0]
CytoZoo.state_index(::Responder, n::Symbol)        = findfirst(==(n), (:y,))
CytoZoo.transmembrane_potential_index(::Responder) = 1
function CytoZoo.parameter_index(::Responder, n::Symbol)
    n === :p_u   && return 1
    n === :h_ext && return 2
    return nothing
end
(m::Responder)(du, u, p, t) = (du[1] = m.parameters[1] - u[1] + 0.1 * m.parameters[2]; nothing)
nothing # hide
```

Two optional subsystems, each a candidate for switching on or off:

```@example patterns
# An optional module that nothing reads — a pure leaf
struct Redox <: AbstractCardiacCellModel end
CytoZoo.num_states(::Redox)                    = 1
CytoZoo.state_names(::Redox)                   = (:z,)
CytoZoo.default_initial_state(::Redox)         = [0.5]
CytoZoo.state_index(::Redox, n::Symbol)        = findfirst(==(n), (:z,))
CytoZoo.transmembrane_potential_index(::Redox) = 1
(::Redox)(du, u, p, t) = (du[1] = -0.3 * u[1]; nothing)

# An optional subsystem that promotes `h` from parameter to state
struct HPool <: AbstractCardiacCellModel end
CytoZoo.num_states(::HPool)                    = 1
CytoZoo.state_names(::HPool)                   = (:h,)
CytoZoo.default_initial_state(::HPool)         = [3.0]
CytoZoo.state_index(::HPool, n::Symbol)        = findfirst(==(n), (:h,))
CytoZoo.transmembrane_potential_index(::HPool) = 1
(::HPool)(du, u, p, t) = (du[1] = -0.05 * u[1]; nothing)

solve_to(cm, T = 2.0) =
    solve(ODEProblem(cm, (0.0, T)), Tsit5(); dt = 0.01, adaptive = false)
nothing # hide
```

## ✅ Feedforward Wire

The most common pattern: one model's state drives another's input. `D.u` becomes `R`'s
`p_u`.

```@example patterns
wired = couple(
    [Subsystem(Driver();    name = :D),
     Subsystem(Responder(); name = :R)],
    [connect(:D => :u, :R => :p_u)],
)

sol = solve_to(wired)
(y = sol.u[end][state_index(wired, :R_y)],)
```

This is a one-way influence. `R` sees `u`; `D` is entirely unaffected by `R`.

## ✅ Adopt-Native (Drop the Duplicate State)

When both models carry the same physical quantity, one model's version wins and the other's
is discarded. This is [`share`](share.md):

```julia
share(:D => :a, :R => :e; owner = :D)
```

`R` keeps reading the quantity, but `D`'s equation governs it. Use this when the two states
*are* the same thing; use a wire when one model merely reads the other.

## ✅ Module On/Off Is Composition

There is no switch. To turn the redox module off, compose without it:

```@example patterns
without_redox = couple(
    [Subsystem(Driver(); name = :D), Subsystem(Responder(); name = :R)],
    [connect(:D => :u, :R => :p_u)],
)

with_redox = couple(
    [Subsystem(Driver(); name = :D), Subsystem(Responder(); name = :R),
     Subsystem(Redox();  name = :X)],
    [connect(:D => :u, :R => :p_u)],
)

(off = state_names(without_redox), on = state_names(with_redox))
```

The **OFF invariant** — that switching a module off recovers the baseline exactly — is free,
because "off" literally means the module is not there. `R`'s trajectory is bit-identical:

```@example patterns
y_off = solve_to(without_redox).u[end][state_index(without_redox, :R_y)]
y_on  = solve_to(with_redox).u[end][state_index(with_redox, :R_y)]

(y_off, y_on, identical = y_off == y_on)
```

Building the node list conditionally is the idiom:

```julia
nodes = [Subsystem(Driver(); name = :D), Subsystem(Responder(); name = :R)]
redox_enabled && push!(nodes, Subsystem(Redox(); name = :X))
coupled = couple(nodes, edges)
```

## ✅ Edge On/Off Is Also Composition

The same mechanism at finer granularity — include or omit a single edge. Without the wire,
`R` runs on its own default parameter value and is fully isolated from `D`:

```@example patterns
unwired = couple(
    [Subsystem(Driver(); name = :D), Subsystem(Responder(); name = :R)],
    [],                                   # no edges at all
)

(wired   = solve_to(wired).u[end][state_index(wired, :R_y)],
 unwired = solve_to(unwired).u[end][state_index(unwired, :R_y)])
```

## ✅ State ↔ Parameter Role Flip

A quantity held as a fixed parameter in one configuration can become a live state in another
— again by composing with or without a subsystem, plus its edge.

Held as a parameter (`h = 1.0`, baked into `Responder`'s defaults):

```@example patterns
h_as_param = couple(
    [Subsystem(Driver(); name = :D), Subsystem(Responder(); name = :R)],
    [connect(:D => :u, :R => :p_u)],
)
solve_to(h_as_param).u[end][state_index(h_as_param, :R_y)]
```

Promoted to a state by adding the `HPool` subsystem and wiring it in:

```@example patterns
h_as_state = couple(
    [Subsystem(Driver(); name = :D), Subsystem(Responder(); name = :R),
     Subsystem(HPool();  name = :H)],
    [connect(:D => :u, :R => :p_u),
     connect(:H => :h, :R => :h_ext)],
)

(states = state_names(h_as_state),
 y      = solve_to(h_as_state).u[end][state_index(h_as_state, :R_y)])
```

Same responder, no code change — only the composition differs. `h` now integrates its own
equation and `R` follows it.

## ✅ Derived Observables

Conservation laws and other algebraic functions of state are not states and should not be
integrated. They are surfaced through the monitor hooks and recomputed after the solve — see
[Derived Observables](../monitors.md).

## ❌ Additive Contributed Flux

Some feedback couplings need two models to each add a term to the *same* state's derivative:
the owner's core equation plus another module's extra flux. Today's `share` is hard-discard
— exactly one equation governs — so this cannot be expressed. See
[Limitations](limitations.md).

## ❌ Derived-Source Wire

A `connect` source must resolve to a state. Sourcing one from a derived quantity — wiring
`ADP = C_A - ATP` into another model's parameter slot — is not supported; you would have to
promote the derived quantity to a real state first. See [Limitations](limitations.md).

## See Also

- [Limitations](limitations.md) — the two ❌ rows, and the automatic-differentiation caveat.
- `examples/coupling_mwe.jl` in the repository — the full taxonomy against a toy pair, with
  executable specifications for the patterns that do not exist yet.
