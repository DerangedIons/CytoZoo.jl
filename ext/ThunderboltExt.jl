module ThunderboltExt

using CytoZoo
using Thunderbolt

struct CytoZooIonicModel{M, SF} <: Thunderbolt.AbstractIonicModel
    model::M
    overrides::SF
end

Thunderbolt.num_states(w::CytoZooIonicModel) = CytoZoo.num_states(w.model)
Thunderbolt.transmembranepotential_index(w::CytoZooIonicModel) =
    CytoZoo.transmembrane_potential_index(w.model)

function Thunderbolt.default_initial_state(w::CytoZooIonicModel, _)
    return CytoZoo.default_initial_state(w.model)
end

function Thunderbolt.cell_rhs!(
    du::AbstractVector{TU},
    u::AbstractVector{TU},
    x::AbstractVector{TX},
    t::TU,
    w::CytoZooIonicModel{<:CytoZoo.AbstractCellModel},
) where {TU, TX}
    p = w.overrides === nothing ? nothing : CytoZoo.SpatialContext(x, w.overrides)
    w.model(du, u, p, t)
    return nothing
end

function Thunderbolt.cell_rhs!(
    du::AbstractVector{TU},
    u::AbstractVector{TU},
    x,
    t::TU,
    w::CytoZooIonicModel{<:CytoZoo.AbstractCellModel},
) where {TU}
    p = w.overrides === nothing ? nothing : CytoZoo.SpatialContext(x, w.overrides)
    w.model(du, u, p, t)
    return nothing
end

CytoZoo.thunderbolt_model(m::CytoZoo.AbstractCellModel; overrides=nothing) =
    CytoZooIonicModel(m, overrides)

end
