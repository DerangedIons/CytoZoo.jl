# Hot-swap two cardiac cell models behind the CytoZoo interface.
#
# Both `ToRORd()` (native to CytoZoo) and `TWorld.TWorldCellModel()` (native
# adherence from TWorld.jl) are `<: CytoZoo.AbstractCellModel`. The simulation
# code below treats them interchangeably — only the constructor differs.
#
# Run: julia --project=examples examples/hot_swap.jl

using CytoZoo, TWorld, OrdinaryDiffEq

function simulate(model; tspan = (0.0, 100.0))
    prob = ODEProblem(model, tspan)
    return solve(prob, Tsit5(); adaptive = true, dtmax = 0.5)
end

for (label, model) in (
        "ToRORd"          => ToRORd(),
        "TWorldCellModel" => TWorldCellModel(; celltype = 0, sex = 2),
    )
    sol = simulate(model)
    vm_idx = transmembrane_potential_index(model)
    println("$label — final Vm = $(sol.u[end][vm_idx]) mV (over $(length(sol.u)) saved steps)")
end
