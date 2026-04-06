module SciMLBaseExt

using CytoZoo
using SciMLBase

function SciMLBase.ODEProblem(model::CytoZoo.AbstractCellModel, tspan::Tuple;
                              u0=CytoZoo.default_initial_state(model), p=nothing, kwargs...)
    return SciMLBase.ODEProblem{true}(model, u0, tspan, p; kwargs...)
end

end
