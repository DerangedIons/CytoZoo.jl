# Three-beat simulation of the Beeler-Reuter 1977 model (MTK-backed)
#
# The model has a single built-in stimulus at t = 20 ms. Subsequent beats are
# triggered via PresetTimeCallback that injects a brief depolarising current.
# Run: julia --project=examples examples/beeler_reuter_action_potential.jl

using CytoZoo, MTKCardiacCellModels, ModelingToolkit, OrdinaryDiffEq, DiffEqCallbacks, CairoMakie

model = BeelerReuter()
vm_idx = transmembrane_potential_index(model)

BCL = 500.0
stim_times = [BCL + 20.0, 2 * BCL + 20.0]

function stim!(integrator)
    integrator.u[vm_idx] += 40.0
end
cb = PresetTimeCallback(stim_times, stim!)

prob = ODEProblem(model, (0.0, 3 * BCL))
sol = solve(prob, Tsit5(); callback=cb, maxiters=1_000_000)

fig = Figure(size=(800, 400))
ax = Axis(fig[1, 1]; xlabel="Time (ms)", ylabel="Vm (mV)",
    title="Beeler-Reuter — 3 beats (BCL = $(Int(BCL)) ms)")
lines!(ax, sol.t, [u[vm_idx] for u in sol.u])
save("beeler_reuter_action_potential.png", fig)
display(fig)
