module ThunderboltExt

using CytoZoo
using Thunderbolt

struct CytoZooIonicModel{M} <: Thunderbolt.AbstractIonicModel
    model::M
end

Thunderbolt.num_states(w::CytoZooIonicModel) = CytoZoo.num_states(w.model)
Thunderbolt.transmembranepotential_index(w::CytoZooIonicModel) =
    CytoZoo.transmembrane_potential_index(w.model)

function Thunderbolt.default_initial_state(w::CytoZooIonicModel, _)
    return CytoZoo.default_initial_state(w.model)
end

# ToRORd without spatial functions
function Thunderbolt.cell_rhs!(
    du::AbstractVector{TU},
    u::AbstractVector{TU},
    x::AbstractVector{TX},
    t::TU,
    w::CytoZooIonicModel{<:CytoZoo.ToRORd},
) where {TU, TX}
    m = w.model
    CytoZoo._torord_rhs_impl!(du, u, m.parameters, m.celltype, x, t, nothing)
    return nothing
end

function Thunderbolt.cell_rhs!(
    du::AbstractVector{TU},
    u::AbstractVector{TU},
    x,
    t::TU,
    w::CytoZooIonicModel{<:CytoZoo.ToRORd},
) where {TU}
    m = w.model
    CytoZoo._torord_rhs_impl!(du, u, m.parameters, m.celltype, x, t, nothing)
    return nothing
end

# Spatial{ToRORd} with spatial functions
function Thunderbolt.cell_rhs!(
    du::AbstractVector{TU},
    u::AbstractVector{TU},
    x::AbstractVector{TX},
    t::TU,
    w::CytoZooIonicModel{<:CytoZoo.Spatial{<:CytoZoo.ToRORd}},
) where {TU, TX}
    m = w.model
    CytoZoo._torord_rhs_impl!(
        du, u, m.model.parameters, m.model.celltype, x, t, m.spatial_funcs,
    )
    return nothing
end

function Thunderbolt.cell_rhs!(
    du::AbstractVector{TU},
    u::AbstractVector{TU},
    x,
    t::TU,
    w::CytoZooIonicModel{<:CytoZoo.Spatial{<:CytoZoo.ToRORd}},
) where {TU}
    m = w.model
    CytoZoo._torord_rhs_impl!(
        du, u, m.model.parameters, m.model.celltype, x, t, m.spatial_funcs,
    )
    return nothing
end

CytoZoo.thunderbolt_model(m::CytoZoo.AbstractCellModel) = CytoZooIonicModel(m)
CytoZoo.thunderbolt_model(m::CytoZoo.Spatial) = CytoZooIonicModel(m)

end
