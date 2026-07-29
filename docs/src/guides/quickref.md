# Quick Reference

A cheat sheet. Code here is illustrative and not executed — see the guides for worked,
runnable versions.

## Constructing Models

```julia
ToRORd()                          # Float64, endocardial, default stimulus
ToRORd(; celltype = 1)            # epicardial   (0 = endo, 1 = epi, 2 = mid)
ToRORd(Float32)                   # Float32 throughout
ToRORd(Float32; celltype = 2)     # Float32, midmyocardial
ToRORd(; stim = Stimulus(; amplitude = -40.0, period = 500.0))
```

## Inspecting a Model

```julia
num_states(model)                       # 65 for ToRORd
num_parameters(model)                   # 177 — atomic models only
state_names(model)                      # (:v, :jca, :m, ...)
parameter_names(model)                  # (:GNa, ...)
transmembrane_potential_index(model)    # 1
default_initial_state(model)            # Vector of length num_states

state_index(model, :v)                  # 1
parameter_index(model, :GNa)            # 14
model.parameters[parameter_index(model, :GNa)] = 11.0
```

## Evaluating and Solving

```julia
du = similar(u)
model(du, u, nothing, t)                # non-spatial RHS
model(du, u, spatial_context, t)        # spatial RHS

prob = ODEProblem(model, (0.0, 1000.0))                 # needs OrdinaryDiffEq
prob = ODEProblem(model, tspan; u0 = u, saveat = 1.0)
sol  = solve(prob, FBDF(); abstol = 1e-8, reltol = 1e-8)
```

## Rush-Larsen

```julia
has_rush_larsen(model)                          # true for ToRORd
rush_larsen_step!(u_new, u, nothing, t, dt, model)   # model comes LAST
```

## Stimulus

```julia
Stimulus()                                       # -53, period 1000 ms, 1 ms wide
Stimulus(; amplitude = -80.0, period = 500.0, duration = 2.0, start = 50.0)
Stimulus(Float32)                                # typed for GPU / Float32

FunctionStimulus((x, t) -> ...)                  # arbitrary waveform
struct MyStim <: AbstractStimulus end            # custom, isbits
(s::MyStim)(x, t) = ...                          # must return the FULL current
```

## Spatial

```julia
SpatialContext(x, overrides)                     # overrides is a NamedTuple, or nothing

# three flavours of override value
(IKr_Multiplier = 0.5,)                                  # scalar
(IKr_Multiplier = (x, t) -> x[1] > 1.5 ? 0.5 : 1.0,)     # callable (CPU)
(IKr_Multiplier = SpatialStep(1, 1.5, 1.0, 0.5),)        # isbits (GPU-safe)

Constant(value)
SpatialStep(dim, threshold, below, above)
SpatialGradient(dim, x_start, x_end, val_start, val_end)
```

## Coupling

```julia
Subsystem(model; name = :A)                      # a graph node

share(:A => :d, :B => :x; owner = :A)            # one variable, owner's equation wins
share(:A => :d, :B => :x; owner = :A, op = +)    # ...or the other's derivative is ADDED to it
share(:A => :d, :B => :x; owner = :A, name = :ATP)   # name the global slot

connect(:A => :v, :B => :v_ext)                  # dataflow: state -> parameter slot
connect(:A => :adp, :B => :adp_ext)              # source may be a monitor (derived) too
connect(:A => :j, :B => :flux; op = +)           # sum several sources into one slot

coupled = couple(nodes, edges)
sol = solve(ODEProblem(coupled, tspan), Tsit5())

state_index(coupled, :d)                         # primary component: bare name
state_index(coupled, :B_y)                       # others: prefixed
```

Turning something off means composing without it — there is no switch:

```julia
nodes = [Subsystem(Core(); name = :C)]
enabled && push!(nodes, Subsystem(Module(); name = :M))
```

## Monitors

```julia
num_monitors(model)                              # 0 unless overridden
monitor_names(model)                             # ()
monitor_values!(mon, u, t, model)                # model LAST, no `p`

h = monitor_history(sol, model)                  # post-solve
h.t; h.names; h.values                           # rows = monitors, cols = time points
```

## Extensions

```julia
using OrdinaryDiffEq     # -> ODEProblem constructor, monitor_history
using Thunderbolt        # -> thunderbolt_model(model; overrides)
# ForwardDiff loads automatically with implicit solvers
```

## Common Errors and Fixes

**`MethodError` calling the model inside a Rosenbrock solver**
`Rodas5P`, `Rodas4`, `Rosenbrock23` need a time gradient, which calls the model with a
dual-typed `du` against a plain `u`. `ToRORd`'s internal right-hand side requires `du` and
`u` to share an element type. Use `FBDF` or another BDF/SDIRK method.

**`MethodError` from `num_parameters(coupled_model)`**
Deliberate. A composite has no single flat parameter vector, so `num_parameters` and
`parameter_names` are left undefined. Query the components instead.

**`couple` rejects your share edges as having different owners**
You chained them (`:A`–`:B`, then `:B`–`:C`). Fan every edge out from the owner instead:
`share(:A => :d, :B => :x; owner = :A)` **and** `share(:A => :d, :C => :z; owner = :A)`.

**`couple` rejects several `connect` edges into one slot**
More than one `overwrite`, or `overwrite` mixed with `+`, is ambiguous. Use `op = +` on all
of them if you want a sum.

**A `connect` receiver is rejected at `couple` time**
The receiver must implement `parameter_index` for the slot *and* expose
`writable_parameters` (which defaults to a `parameters` field).

**`KeyError` instead of a readable message from `couple`**
`couple` validates a `connect` slot with `parameter_index(receiver, slot) === nothing`, so a
receiver whose `parameter_index` *throws* on an unknown name leaks a bare `KeyError`. Use
`get(INDEX, name, nothing)`, not `INDEX[name]`. Connecting into a `ToRORd` with a mistyped
slot currently hits this.

**Gradients across a coupling are wrong**
Expected, and silent. Derivatives through a `connect` edge lose the coupling term. Use
`share` where you need differentiability. See
[Coupling Limitations](coupling/limitations.md).

**Corrupted results from a threaded sweep**
One `CoupledModel` is not thread-safe. `deepcopy` it per trajectory.

**The first beat looks wrong**
`default_initial_state(ToRORd())` is not a pre-paced steady state. Pace for several beats and
start from the final state.

## Performance Tips

- **Keep the functor allocation-free.** Write into `du` and never build temporaries; this is
  asserted by the test suite for `ToRORd`.
- **Pass `p = nothing` when you have no spatial variation.** Dispatch removes every spatial
  branch at compile time, so there is nothing to skip at run time.
- **Use isbits everywhere for GPU.** `SVector` positions, `SpatialFunction` functors rather
  than closures, and a typed `Stimulus(Float32)`. Check with `isbitstype(typeof(p))`.
- **Prefer `FBDF` over explicit methods for stiff single-cell runs** — roughly 100× fewer
  steps than `Tsit5` on `ToRORd`.
- **Prefer Rush-Larsen for long pacing runs**, where fixed steps and no linear solve win.
- **Edge operations are resolved at `couple` time**, so a coupled right-hand side has no
  dynamic dispatch — build the coupling once and reuse it.
