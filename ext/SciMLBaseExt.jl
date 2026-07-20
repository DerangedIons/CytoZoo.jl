module SciMLBaseExt

using CytoZoo
using SciMLBase

function SciMLBase.ODEProblem(model::CytoZoo.AbstractCellModel, tspan::Tuple;
                              u0=CytoZoo.default_initial_state(model), p=nothing, kwargs...)
    return SciMLBase.ODEProblem{true}(model, u0, tspan, p; kwargs...)
end

# Docstring lives on the stub in src/CytoZoo.jl so Documenter's `Modules = [CytoZoo]` sees it.
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
