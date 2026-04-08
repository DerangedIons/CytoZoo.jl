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

    prob = ODEProblem(sys, [], (0.0, 1000.0))
    return _build_mtk_model(:BeelerReuter, sys, prob, vm_sym)
end
