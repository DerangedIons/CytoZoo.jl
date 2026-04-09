module ThunderboltExt

using CytoZoo
using Thunderbolt

struct CytoZooIonicModel{M, SF} <: Thunderbolt.AbstractIonicModel
    model::M
    spatial_funcs::SF
end

Thunderbolt.num_states(w::CytoZooIonicModel) = CytoZoo.num_states(w.model)
Thunderbolt.transmembranepotential_index(w::CytoZooIonicModel) =
    CytoZoo.transmembrane_potential_index(w.model)

function Thunderbolt.default_initial_state(w::CytoZooIonicModel, _)
    return CytoZoo.default_initial_state(w.model)
end

# ToRORd — typed x
function Thunderbolt.cell_rhs!(
    du::AbstractVector{TU},
    u::AbstractVector{TU},
    x::AbstractVector{TX},
    t::TU,
    w::CytoZooIonicModel{<:CytoZoo.ToRORd},
) where {TU, TX}
    p = w.spatial_funcs === nothing ? nothing : CytoZoo.SpatialContext(x, w.spatial_funcs)
    w.model(du, u, p, t)
    return nothing
end

# ToRORd — untyped x
function Thunderbolt.cell_rhs!(
    du::AbstractVector{TU},
    u::AbstractVector{TU},
    x,
    t::TU,
    w::CytoZooIonicModel{<:CytoZoo.ToRORd},
) where {TU}
    p = w.spatial_funcs === nothing ? nothing : CytoZoo.SpatialContext(x, w.spatial_funcs)
    w.model(du, u, p, t)
    return nothing
end

CytoZoo.thunderbolt_model(m::CytoZoo.AbstractCellModel; spatial_funcs=nothing) =
    CytoZooIonicModel(m, spatial_funcs)

end
