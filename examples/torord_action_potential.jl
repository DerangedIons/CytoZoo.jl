# Three-beat simulation of the ToRORd endocardial model using Rush-Larsen
#
# The ToRORd model has a built-in periodic stimulus (BCL = 1000 ms by default)
# and a Rush-Larsen exponential integrator for stable time-stepping.
# Run: julia --project=examples examples/torord_action_potential.jl

using CytoZoo, CairoMakie

model = ToRORd()
vm_idx = transmembrane_potential_index(model)
dt = 0.01  # ms

u = default_initial_state(model)
u_new = similar(u)

n_steps = round(Int, 3000.0 / dt)
t_hist = Vector{Float64}(undef, n_steps + 1)
v_hist = Vector{Float64}(undef, n_steps + 1)
t_hist[1] = 0.0
v_hist[1] = u[vm_idx]

for i in 1:n_steps
    t = (i - 1) * dt
    rush_larsen_step!(u_new, u, t, dt, model)
    u .= u_new
    t_hist[i + 1] = t + dt
    v_hist[i + 1] = u[vm_idx]
end

fig = Figure(size=(800, 400))
ax = Axis(fig[1, 1]; xlabel="Time (ms)", ylabel="Vm (mV)",
    title="ToRORd — 3 beats (BCL = 1000 ms, Rush-Larsen dt = $(dt) ms)")
lines!(ax, t_hist, v_hist)
save("torord_action_potential.png", fig)
display(fig)
