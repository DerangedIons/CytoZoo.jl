module ThunderboltTWorldExt

using CytoZoo
using Thunderbolt
using TWorld

const CytoZooIonicModel = Base.get_extension(CytoZoo, :ThunderboltExt).CytoZooIonicModel

function Thunderbolt.cell_rhs!(
    du::AbstractVector{TU},
    u::AbstractVector{TU},
    x::AbstractVector{TX},
    t::TU,
    w::CytoZooIonicModel{<:CytoZoo.TWorldCellModel},
) where {TU, TX}
    TWorld.tworld_ode!(du, u, w.model.params, t)
    return nothing
end

function Thunderbolt.cell_rhs!(
    du::AbstractVector{TU},
    u::AbstractVector{TU},
    x,
    t::TU,
    w::CytoZooIonicModel{<:CytoZoo.TWorldCellModel},
) where {TU}
    TWorld.tworld_ode!(du, u, w.model.params, t)
    return nothing
end

end
