# Model Catalog

## Models Shipped With CytoZoo

| Model | States | Parameters | Cell types | Rush-Larsen | Monitors |
|-------|--------|------------|------------|-------------|----------|
| [`ToRORd`](@ref) | 65 | 177 | endocardial, epicardial, midmyocardial | yes | none yet |
| [`FHNModel`](@ref) | 2 | 5 | — | no | none |

### `ToRORd`

A modified O'Hara–Rudy human ventricular myocyte model with mechanical coupling,
trauma/hypoxia effects, and transmural cellular heterogeneity.

```@example models
using CytoZoo

model = ToRORd()
(states = num_states(model), parameters = num_parameters(model))
```

#### Constructors

```julia
ToRORd()                            # Float64, endocardial, default stimulus
ToRORd(; celltype = 1)              # epicardial
ToRORd(Float32)                     # Float32 element type throughout
ToRORd(Float32; celltype = 2)       # Float32, midmyocardial
ToRORd(Vector{Float64})             # specify the parameter vector type (GPU paths)
ToRORd(; stim = Stimulus(; amplitude = -40.0, period = 500.0))
```

| Keyword | Default | Meaning |
|---|---|---|
| `celltype` | `0` | `0` endocardial, `1` epicardial, `2` midmyocardial |
| `stim` | `Stimulus()` | any [`AbstractStimulus`](@ref); −53 amplitude, 1000 ms period, 1 ms wide |

The first positional argument sets the element type. Passing a `Number` type builds a
`Vector` of that type; passing an `AbstractVector` type builds that container, which is the
hook for GPU array types.

#### Fields

```@example models
(celltype = model.celltype, stim = model.stim, nparams = length(model.parameters))
```

#### States

65 states, with the transmembrane potential first:

```@example models
transmembrane_potential_index(model), state_names(model)[1:10]
```

Look up any state by name rather than position:

```@example models
state_index(model, :CaMKt)
```

#### Parameters

177 parameters in a flat, mutable vector:

```@example models
parameter_names(model)[1:8]
```

```@example models
gi = parameter_index(model, :GNa)
model.parameters[gi]
```

Parameters commonly used as spatial overrides include `IKr_Multiplier` and `isHypoxic` — see
[Spatial Heterogeneity](../guides/spatial.md).

#### Integration

Rush-Larsen is available and is the conventional choice for long pacing runs:

```@example models
has_rush_larsen(model)
```

!!! warning "Rosenbrock solvers do not work with `ToRORd`"
    `ToRORd`'s internal right-hand side is typed `(du::AbstractVector{T}, u::AbstractVector{T})`,
    requiring both to share an element type. Rosenbrock methods (`Rodas5P`, `Rodas4`,
    `Rosenbrock23`) compute a time gradient by calling the model with a dual-typed `du`
    against a plain `u`, which raises a `MethodError`. Use `FBDF` or another BDF/SDIRK method,
    or Rush-Larsen.

#### Known Issues

- **Monitors are not ported.** `num_monitors(ToRORd()) == 0`. The reference implementation
  defines roughly 492 derived quantities (individual currents and fluxes); none are available
  yet. See [Derived Observables](../guides/monitors.md).
- **`state_index` and `parameter_index` throw on unknown names.** Both index into a `Dict`
  directly, so they raise a `KeyError` rather than returning `nothing` as the interface
  convention specifies. The visible consequence is that a `connect` edge naming a mistyped
  `ToRORd` parameter slot fails with a bare `KeyError` instead of `couple`'s actionable
  message. Share and connect edges naming a mistyped *state* are unaffected, since those are
  validated against [`state_names`](@ref).
- **`default_initial_state` is not a pre-paced steady state.** It differs substantially from
  the steady state used in the correctness tests, and the intrinsic dV/dt there is roughly
  +240 mV/ms, so the opening transient of a simulation is not physiological. Pace for several
  beats and take the final state before measuring biomarkers. The right-hand side itself is
  validated against the ArmyHeart reference at `rtol = 1e-10`; this affects the initial
  condition only.

### `FHNModel`

The FitzHugh–Nagumo excitable-medium model — a dimensionless caricature of an action
potential with one cubic fast variable and one linear slow recovery variable.

```@example models_fhn
using CytoZoo

model = FHNModel()
(states = num_states(model), parameters = num_parameters(model))
```

It is a *test* model, not a physiological one: `v` runs roughly over `[0, 1]` rather than
millivolts and time is in arbitrary units. Use it to exercise a solver, a coupling graph, or
a tissue framework's spatial machinery without paying for 65-state kinetics.

#### Equations

```math
\begin{aligned}
\frac{dv}{dt} &= v(1 - v)(v - a) - s - I_\text{stim} \\
\frac{ds}{dt} &= e\,(b\,v - c\,s - d)
\end{aligned}
```

#### Constructors

```julia
FHNModel()                            # Float64, published parameters, stimulus off
FHNModel(Float32)                     # Float32 element type throughout
FHNModel(; a = 0.15)                  # override one parameter
FHNModel(; stim = Stimulus(; amplitude = -0.5, duration = 1.0))
```

#### Parameters

```@example models_fhn
parameter_names(model), (model.a, model.b, model.c, model.d, model.e)
```

| Name | Default | Meaning |
|---|---|---|
| `a` | `0.1` | excitation threshold; the middle root of the cubic |
| `b` | `0.5` | recovery gain on `v` |
| `c` | `1.0` | recovery decay on `s` |
| `d` | `0.0` | recovery offset |
| `e` | `0.01` | recovery/excitation timescale ratio |

All five are spatially overridable — see [Spatial Heterogeneity](../guides/spatial.md).

#### States

```@example models_fhn
state_names(model), default_initial_state(model)
```

`(0, 0)` is a genuine steady state, not an un-paced approximation: a solve started there
stays there until something stimulates it.

#### The Stimulus Defaults to Zero

```@example models_fhn
model.stim
```

A tissue framework normally injects the stimulus as a source term in its diffusion half, and
a nonzero default here would silently add a second one. Pass an explicit
[`Stimulus`](@ref) for single-cell use. The sign convention matches the rest of the zoo: `Istim` is *subtracted*
from `dv/dt`, so a **negative** amplitude depolarizes.

```@example models_fhn
firing = FHNModel(; stim = Stimulus(; amplitude = -0.5, duration = 1.0))
du = similar(default_initial_state(firing))
firing(du, default_initial_state(firing), nothing, 0.0)
du
```

#### Parameters Are Immutable Struct Fields

`FHNModel` deliberately breaks the flat-`parameters`-vector convention: its five
parameters are plain fields, which keeps the whole model isbits so it can be captured by value
inside a GPU kernel. The cost is that it cannot receive a [`connect`](@ref) edge —
[`couple`](@ref) rejects that at construction with an actionable message. Build a new model
rather than mutating one.

## Models in Other Packages

Model packages that adhere to the interface natively depend on CytoZoo and subtype
`CytoZoo.AbstractCardiacCellModel` directly. Load the package and drive the model through the
same uniform interface — no adapter, no registration step.

| Package | Type | Notes |
|---|---|---|
| [TWorld.jl](https://github.com/DerangedIons/TWorld.jl) | `TWorldCellModel` | cell type and sex-specific variants |

```julia
using CytoZoo, TWorld

model = TWorldCellModel(; celltype = 0, sex = 2)
state_names(model)          # same interface as any CytoZoo model
```

Because every model satisfies the same contract, `using CytoZoo, TWorld, SomeOtherModel` lets
you hot-swap between them behind one function. See
[The Cell Model Interface](../guides/interface.md).

## Adding a Model

See [Implementing a Model](../guides/implementing_a_model.md) for both the in-tree recipe and
the two integration patterns for external packages.
