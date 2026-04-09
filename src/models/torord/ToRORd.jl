include("parameters.jl")
include("states.jl")
include("rhs.jl")
include("rush_larsen.jl")
include("monitors.jl")

"""
    ToRORd{T}

ToRORd trauma-hetero cardiac cell model (65 states, 181 parameters).

A modified O'Hara-Rudy model with mechanical coupling, trauma/hypoxia effects,
and cellular heterogeneity (endo/mid/epi).

# Constructor
```julia
ToRORd()                        # Float64, endocardial
ToRORd(Float32; celltype=1)     # Float32, epicardial
ToRORd(Vector{Float64})         # specify vector type for GPU
```

Cell types: `0` = endocardial, `1` = epicardial, `2` = midmyocardial.
"""
struct ToRORd{T <: AbstractVector} <: AbstractCardiacCellModel
    parameters::T
    celltype::Int
end

ToRORd(; celltype::Int = 0) = ToRORd(Float64; celltype)

function ToRORd(::Type{ElT}; celltype::Int = 0) where {ElT <: Number}
    p = zeros(ElT, TORORD_NUM_PARAMS)
    _torord_init_parameters!(p)
    return ToRORd(p, celltype)
end

function ToRORd(::Type{VT}; celltype::Int = 0) where {VT <: AbstractVector}
    p = zeros(eltype(VT), TORORD_NUM_PARAMS)
    _torord_init_parameters!(p)
    return ToRORd(VT(p), celltype)
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
    _torord_rhs_impl!(du, u, model.parameters, model.celltype, nothing, t, nothing)
    return nothing
end

function (model::ToRORd)(du, u, p::SpatialContext, t)
    _torord_rhs_impl!(du, u, model.parameters, model.celltype, p.x, t, p.spatial_funcs)
    return nothing
end

# ---------------------------------------------------------------------------
# Optional interface
# ---------------------------------------------------------------------------

has_rush_larsen(::ToRORd) = true

function rush_larsen_step!(u_new, u, ::Nothing, t, dt, model::ToRORd)
    _torord_rush_larsen_impl!(u_new, u, model.parameters, model.celltype, nothing, t, dt, nothing)
    return nothing
end

function rush_larsen_step!(u_new, u, p::SpatialContext, t, dt, model::ToRORd)
    _torord_rush_larsen_impl!(u_new, u, model.parameters, model.celltype, p.x, t, dt, p.spatial_funcs)
    return nothing
end

num_monitors(::ToRORd) = TORORD_NUM_MONITORS

state_index(::ToRORd, name::Symbol) = TORORD_STATE_INDEX[name]
parameter_index(::ToRORd, name::Symbol) = TORORD_PARAM_INDEX[name]
state_names(::ToRORd) = TORORD_STATE_NAMES
parameter_names(::ToRORd) = TORORD_PARAMETER_NAMES

