include("parameters.jl")
include("states.jl")
include("rhs.jl")

"""
    ParametrizedFHNModel{T, S}

FitzHugh–Nagumo excitable-medium model (2 states, 5 parameters).

    dv/dt = v(1 - v)(v - a) - s - Iₛₜᵢₘ
    ds/dt = e(b·v - c·s - d)

A dimensionless caricature of an action potential: a cubic fast variable `v` and a
linear slow recovery variable `s`. It is not calibrated to any physiological units —
`v` runs roughly over `[0, 1]` rather than millivolts, and time is in arbitrary units
set by `e`. Its value is as a *test* model: two states, closed-form-ish behaviour, and
a right-hand side cheap enough that a tissue solve is dominated by the diffusion
operator rather than the kinetics.

# Constructors
```julia
ParametrizedFHNModel()                            # Float64, published parameters, stimulus off
ParametrizedFHNModel(Float32)                     # Float32 throughout
ParametrizedFHNModel(; a = 0.15)                  # override one parameter
ParametrizedFHNModel(; stim = Stimulus(; amplitude = -0.5, duration = 1.0))
```

`FHNModel` is a `Float64` alias for use in signatures and `isa` tests; construct through
`ParametrizedFHNModel`.

# Parameters

| Name | Default | Meaning |
|---|---|---|
| `a` | `0.1` | excitation threshold; the middle root of the cubic |
| `b` | `0.5` | recovery gain on `v` |
| `c` | `1.0` | recovery decay on `s` |
| `d` | `0.0` | recovery offset |
| `e` | `0.01` | recovery/excitation timescale ratio |

All five are spatially overridable through a [`SpatialContext`](@ref).

# Stimulus

Unlike [`ToRORd`](@ref), the default `stim` has **zero amplitude**. A tissue framework
normally supplies the stimulus as a source term in its diffusion half rather than through
the cell model, and a nonzero default would silently add a second one. Pass an explicit
[`Stimulus`](@ref) for single-cell use. The sign convention matches the rest of the zoo:
`Istim` is subtracted from `dv/dt`, so a **negative** amplitude depolarizes.

# Parameters are immutable struct fields

The five parameters are plain fields, not a `parameters` vector, which keeps the whole
model isbits and lets it ride into a GPU kernel by value. The cost is that a
[`connect`](@ref) edge cannot target it — [`couple`](@ref) rejects that at construction
with an actionable message. Build a new model instead of mutating one.
"""
struct ParametrizedFHNModel{T <: Number, S} <: AbstractCardiacCellModel
    a::T
    b::T
    c::T
    d::T
    e::T
    stim::S
end

"""
    FHNModel

`Float64` alias of [`ParametrizedFHNModel`](@ref), matching the name Thunderbolt uses.
It is a type, not a constructor — build models with `ParametrizedFHNModel(...)`.
"""
const FHNModel = ParametrizedFHNModel{Float64}

ParametrizedFHNModel(; kwargs...) = ParametrizedFHNModel(Float64; kwargs...)

function ParametrizedFHNModel(
    ::Type{T};
    a::Real = FHN_DEFAULT_PARAMETERS.a,
    b::Real = FHN_DEFAULT_PARAMETERS.b,
    c::Real = FHN_DEFAULT_PARAMETERS.c,
    d::Real = FHN_DEFAULT_PARAMETERS.d,
    e::Real = FHN_DEFAULT_PARAMETERS.e,
    stim = Stimulus(T; amplitude = 0),
) where {T <: Number}
    return ParametrizedFHNModel(T(a), T(b), T(c), T(d), T(e), stim)
end

# ---------------------------------------------------------------------------
# Required interface
# ---------------------------------------------------------------------------

num_states(::ParametrizedFHNModel) = FHN_NUM_STATES
num_parameters(::ParametrizedFHNModel) = FHN_NUM_PARAMS
transmembrane_potential_index(::ParametrizedFHNModel) = 1

function default_initial_state(::ParametrizedFHNModel{T}) where {T}
    return T[FHN_DEFAULT_INITIAL_STATE...]
end

# ---------------------------------------------------------------------------
# Functor — DiffEq-compatible ODE RHS
# ---------------------------------------------------------------------------

function (model::ParametrizedFHNModel)(du, u, ::Nothing, t)
    _fhn_rhs_impl!(du, u, model, nothing, t, nothing)
    return nothing
end

function (model::ParametrizedFHNModel)(du, u, p::SpatialContext, t)
    _fhn_rhs_impl!(du, u, model, p.x, t, p.overrides)
    return nothing
end

# ---------------------------------------------------------------------------
# Optional interface
# ---------------------------------------------------------------------------

state_names(::ParametrizedFHNModel) = FHN_STATE_NAMES
parameter_names(::ParametrizedFHNModel) = FHN_PARAMETER_NAMES

# `findfirst` returns `nothing` for an unknown name, which is the interface convention
# `couple`'s validation depends on (`ToRORd`'s `Dict` lookups throw instead).
state_index(::ParametrizedFHNModel, name::Symbol) = findfirst(==(name), FHN_STATE_NAMES)
parameter_index(::ParametrizedFHNModel, name::Symbol) =
    findfirst(==(name), FHN_PARAMETER_NAMES)
