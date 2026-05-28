include("parameters.jl")
include("states.jl")
include("rhs.jl")
include("rush_larsen.jl")
include("monitors.jl")

"""
    ToRORd{T, S}

ToRORd trauma-hetero cardiac cell model (65 states, 177 parameters).

A modified O'Hara-Rudy model with mechanical coupling, trauma/hypoxia effects,
and cellular heterogeneity (endo/mid/epi).

# Constructor
```julia
ToRORd()                            # Float64, endocardial, default Stimulus
ToRORd(Float32; celltype=1)         # Float32, epicardial
ToRORd(Vector{Float64})             # specify vector type for GPU
ToRORd(; stim = Stimulus(; amplitude = -40.0, period = 500.0))  # custom stimulus
```

Cell types: `0` = endocardial, `1` = epicardial, `2` = midmyocardial.

The `stim` field is the canonical time-only stimulus (default: -53 mV, period
1000 ms, 1 ms duration). Per-cell variation is layered on via
`SpatialContext(x, (stim = (x, t) -> ...,))` — the spatial override returns
the full `Istim` value (not a multiplier).
"""
struct ToRORd{T <: AbstractVector, S} <: AbstractCardiacCellModel
    parameters::T
    celltype::Int
    stim::S
end

ToRORd(; celltype::Int = 0, stim = Stimulus()) = ToRORd(Float64; celltype, stim)

function ToRORd(::Type{ElT}; celltype::Int = 0, stim = Stimulus(ElT)) where {ElT <: Number}
    p = zeros(ElT, TORORD_NUM_PARAMS)
    _torord_init_parameters!(p)
    return ToRORd(p, celltype, stim)
end

function ToRORd(::Type{VT}; celltype::Int = 0, stim = Stimulus(eltype(VT))) where {VT <: AbstractVector}
    p = zeros(eltype(VT), TORORD_NUM_PARAMS)
    _torord_init_parameters!(p)
    return ToRORd(VT(p), celltype, stim)
end

# ---------------------------------------------------------------------------
# Required interface
# ---------------------------------------------------------------------------

num_states(::ToRORd) = TORORD_NUM_STATES
num_parameters(::ToRORd) = TORORD_NUM_PARAMS
transmembrane_potential_index(::ToRORd) = 1

function default_initial_state(model::ToRORd{T}) where {T}
    u = zeros(eltype(T), TORORD_NUM_STATES)
    _torord_init_state_values!(u)
    return u
end

# ---------------------------------------------------------------------------
# Functor — DiffEq-compatible ODE RHS
# ---------------------------------------------------------------------------

function (model::ToRORd)(du, u, ::Nothing, t)
    _torord_rhs_impl!(du, u, model.parameters, model.celltype, model.stim, nothing, t, nothing)
    return nothing
end

function (model::ToRORd)(du, u, p::SpatialContext, t)
    _torord_rhs_impl!(du, u, model.parameters, model.celltype, model.stim, p.x, t, p.spatial_funcs)
    return nothing
end

# ---------------------------------------------------------------------------
# Optional interface
# ---------------------------------------------------------------------------

has_rush_larsen(::ToRORd) = true

function rush_larsen_step!(u_new, u, ::Nothing, t, dt, model::ToRORd)
    _torord_rush_larsen_impl!(u_new, u, model.parameters, model.celltype, model.stim, nothing, t, dt, nothing)
    return nothing
end

function rush_larsen_step!(u_new, u, p::SpatialContext, t, dt, model::ToRORd)
    _torord_rush_larsen_impl!(u_new, u, model.parameters, model.celltype, model.stim, p.x, t, dt, p.spatial_funcs)
    return nothing
end

num_monitors(::ToRORd) = TORORD_NUM_MONITORS

state_index(::ToRORd, name::Symbol) = TORORD_STATE_INDEX[name]
parameter_index(::ToRORd, name::Symbol) = TORORD_PARAM_INDEX[name]
state_names(::ToRORd) = TORORD_STATE_NAMES
parameter_names(::ToRORd) = TORORD_PARAMETER_NAMES
