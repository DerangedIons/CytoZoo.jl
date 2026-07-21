# Model Catalog

## Models Shipped With CytoZoo

| Model | States | Parameters | Cell types | Rush-Larsen | Monitors |
|-------|--------|------------|------------|-------------|----------|
| [`ToRORd`](@ref) | 65 | 177 | endocardial, epicardial, midmyocardial | yes | none yet |

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
