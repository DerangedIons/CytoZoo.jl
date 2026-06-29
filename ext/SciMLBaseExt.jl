module SciMLBaseExt

using CytoZoo
using SciMLBase

function SciMLBase.ODEProblem(model::CytoZoo.AbstractCellModel, tspan::Tuple;
                              u0=CytoZoo.default_initial_state(model), p=nothing, kwargs...)
    return SciMLBase.ODEProblem{true}(model, u0, tspan, p; kwargs...)
end

"""
    monitor_history(sol, model::AbstractCellModel) -> (; t, names, values)

Recompute `model`'s derived/monitored quantities (see [`monitor_values!`](@ref)) across a
solved trajectory `sol`. Returns a `NamedTuple` of the time points `t`, the monitor `names`,
and a `values` matrix (rows = monitors in `names` order, columns = time points).

Post-solve by design: the saved `sol.u` is plain (no `Dual`s), so this avoids any
implicit-solver primal-extraction concerns. A model with `num_monitors(model) == 0` yields an
empty (`0 × length(sol.t)`) `values` matrix.
"""
function CytoZoo.monitor_history(sol, model::CytoZoo.AbstractCellModel)
    names = CytoZoo.monitor_names(model)
    nmon = CytoZoo.num_monitors(model)
    T = eltype(eltype(sol.u))
    values = Matrix{T}(undef, nmon, length(sol.t))
    nmon == 0 && return (; t = sol.t, names, values)
    mon = zeros(T, nmon)
    for (j, u) in enumerate(sol.u)
        CytoZoo.monitor_values!(mon, u, sol.t[j], model)
        @views values[:, j] .= mon
    end
    return (; t = sol.t, names, values)
end

end
