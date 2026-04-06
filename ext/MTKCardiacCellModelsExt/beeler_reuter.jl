function _rounded_pulse(t_c; t_start=20.0, duration=5.0, amplitude=10.0, smoothness=0.1)
    t_end = t_start + duration
    rise = 0.5 * (1 + tanh((t_c - t_start) / smoothness))
    fall = 0.5 * (1 - tanh((t_c - t_end) / smoothness))
    return amplitude * rise * fall
end

@register_symbolic _rounded_pulse(t)

function _build_beeler_reuter()
    _BeelerReuterModel = getfield(MTKCardiacCellModels, :BeelerReuterModel)
    @named br = _BeelerReuterModel(I_stim=_rounded_pulse)
    sys = structural_simplify(br)

    unknowns = ModelingToolkit.unknowns(sys)
    vm_sym = only(filter(u -> endswith(string(u), "φₘ(t)"), unknowns))

    u0_map = Dict(u => 0.0 for u in unknowns)
    for u in unknowns
        s = string(u)
        if endswith(s, "φₘ(t)")
            u0_map[u] = -84.0
        elseif endswith(s, "Caᵢ(t)")
            u0_map[u] = 1e-7
        elseif endswith(s, "m₊y(t)")
            u0_map[u] = 0.011
        elseif endswith(s, "h₊y(t)")
            u0_map[u] = 0.988
        elseif endswith(s, "j₊y(t)")
            u0_map[u] = 0.975
        elseif endswith(s, "d₊y(t)")
            u0_map[u] = 0.003
        elseif endswith(s, "f₊y(t)")
            u0_map[u] = 0.994
        elseif endswith(s, "gate₊y(t)")
            u0_map[u] = 0.0001
        end
    end

    p_map = ModelingToolkit.defaults(sys)
    op = merge(u0_map, p_map)

    prob = ODEProblem(sys, op, (0.0, 1000.0))
    return _build_mtk_model(:BeelerReuter, sys, prob, vm_sym)
end
