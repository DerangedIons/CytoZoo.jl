module TWorldExt

using CytoZoo
using TWorld

# ---------------------------------------------------------------------------
# State name extraction from TWorld's IDX_* constants
# ---------------------------------------------------------------------------

const _IDX_PAIRS = let
    pairs = Tuple{Int,Symbol}[]
    for nm in names(TWorld)
        s = string(nm)
        if startswith(s, "IDX_")
            idx = getfield(TWorld, nm)
            if idx isa Integer
                push!(pairs, (idx, Symbol(lowercase(s[5:end]))))
            end
        end
    end
    sort!(pairs; by=first)
    pairs
end

const TWORLD_NUM_STATES = TWorld.N_STATES
const TWORLD_STATE_NAMES = ntuple(i -> _IDX_PAIRS[i][2], TWORLD_NUM_STATES)
const TWORLD_STATE_INDEX = Dict{Symbol,Int}(v => k for (k, v) in _IDX_PAIRS)

# ---------------------------------------------------------------------------
# Parameter name extraction from TWorldParameters fieldnames (filtered)
# ---------------------------------------------------------------------------

const _PARAM_FILTER = Set((:stim_fn, :x_coord, :y_coord, :z_coord))

const TWORLD_PARAM_NAMES = let
    fnames = fieldnames(TWorld.TWorldParameters)
    Tuple(n for n in fnames if n ∉ _PARAM_FILTER)
end

const TWORLD_NUM_PARAMS = length(TWORLD_PARAM_NAMES)
const TWORLD_PARAM_INDEX = Dict{Symbol,Int}(n => i for (i, n) in enumerate(TWORLD_PARAM_NAMES))

# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

function CytoZoo.TWorldCellModel(; kwargs...)
    p = TWorld.TWorldParameters(; kwargs...)
    return CytoZoo.TWorldCellModel(p)
end

# ---------------------------------------------------------------------------
# Required interface
# ---------------------------------------------------------------------------

CytoZoo.num_states(::CytoZoo.TWorldCellModel) = TWORLD_NUM_STATES
CytoZoo.num_parameters(::CytoZoo.TWorldCellModel) = TWORLD_NUM_PARAMS
CytoZoo.transmembrane_potential_index(::CytoZoo.TWorldCellModel) = TWorld.IDX_V

function CytoZoo.default_initial_state(::CytoZoo.TWorldCellModel)
    return TWorld.tworld_initial_conditions()
end

# ---------------------------------------------------------------------------
# Functor — DiffEq-compatible ODE RHS
# ---------------------------------------------------------------------------

function (m::CytoZoo.TWorldCellModel)(du, u, ::Nothing, t)
    TWorld.tworld_ode!(du, u, m.params, t)
    return nothing
end

function (m::CytoZoo.TWorldCellModel)(du, u, p::CytoZoo.SpatialContext, t)
    # TODO: Thread p.spatial_funcs through TWorld.tworld_ode! once TWorld supports them
    TWorld.tworld_ode!(du, u, m.params, t)
    return nothing
end

# ---------------------------------------------------------------------------
# Optional interface — name-based access
# ---------------------------------------------------------------------------

CytoZoo.state_names(::CytoZoo.TWorldCellModel) = TWORLD_STATE_NAMES
CytoZoo.parameter_names(::CytoZoo.TWorldCellModel) = TWORLD_PARAM_NAMES
CytoZoo.state_index(::CytoZoo.TWorldCellModel, name::Symbol) = TWORLD_STATE_INDEX[name]
CytoZoo.parameter_index(::CytoZoo.TWorldCellModel, name::Symbol) = TWORLD_PARAM_INDEX[name]

# ---------------------------------------------------------------------------
# Optional interface — Rush-Larsen
# ---------------------------------------------------------------------------

CytoZoo.has_rush_larsen(::CytoZoo.TWorldCellModel) = true

function _get_tworld_workspace(::Type{T}) where {T}
    key = :tworld_rl_workspace
    ws = get(task_local_storage(), key, nothing)
    if ws === nothing || eltype(ws) !== T
        ws = Vector{T}(undef, TWorld.N_WORKSPACE)
        task_local_storage(key, ws)
    end
    return ws::Vector{T}
end

function CytoZoo.rush_larsen_step!(u_new, u, ::Nothing, t, dt, m::CytoZoo.TWorldCellModel)
    ws = _get_tworld_workspace(eltype(u))
    TWorld.tworld_rl_step!(u_new, u, m.params, t, dt, ws)
    return nothing
end

function CytoZoo.rush_larsen_step!(u_new, u, p::CytoZoo.SpatialContext, t, dt, m::CytoZoo.TWorldCellModel)
    # TODO: Thread p.spatial_funcs through TWorld.tworld_rl_step! once TWorld supports them
    ws = _get_tworld_workspace(eltype(u))
    TWorld.tworld_rl_step!(u_new, u, m.params, t, dt, ws)
    return nothing
end

end
