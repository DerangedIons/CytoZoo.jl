module MTKCardiacCellModelsExt

using CytoZoo
using MTKCardiacCellModels
using MTKCardiacCellModels: AlphaBetaGate, IonChannelConductance, IonCurrent, LeakCurrent,
    SlowInwardCurrent, TimeIndependentK, TimeActivatedOutward, CalciumDynamics,
    StimulationSystem
using ModelingToolkit
using SciMLBase

import ModelingToolkit: t_nounits as t, D_nounits as D

# ---------------------------------------------------------------------------
# MTKCardiacModel — wrapper for MTK-backed cardiac cell models
# ---------------------------------------------------------------------------

struct MTKCardiacModel{S,Prob} <: CytoZoo.AbstractCardiacCellModel
    name::Symbol
    sys::S
    prob::Prob
    state_syms::Vector{Symbol}
    param_syms::Vector{Symbol}
    vm_index::Int
end

# --- CytoZoo interface ---

CytoZoo.num_states(m::MTKCardiacModel) = length(m.state_syms)
CytoZoo.num_parameters(m::MTKCardiacModel) = length(m.param_syms)
CytoZoo.transmembrane_potential_index(m::MTKCardiacModel) = m.vm_index
CytoZoo.default_initial_state(m::MTKCardiacModel) = collect(Float64, m.prob.u0)

CytoZoo.state_names(m::MTKCardiacModel) = Tuple(m.state_syms)
CytoZoo.parameter_names(m::MTKCardiacModel) = Tuple(m.param_syms)

function CytoZoo.state_index(m::MTKCardiacModel, name::Symbol)
    idx = findfirst(==(name), m.state_syms)
    isnothing(idx) && error("State :$name not found in $(m.name)")
    return idx
end

function CytoZoo.parameter_index(m::MTKCardiacModel, name::Symbol)
    idx = findfirst(==(name), m.param_syms)
    isnothing(idx) && error("Parameter :$name not found in $(m.name)")
    return idx
end

function (m::MTKCardiacModel)(du, u, p, t)
    m.prob.f(du, u, m.prob.p, t)
    return nothing
end

CytoZoo.symbolic_system(m::MTKCardiacModel) = m.sys
CytoZoo.has_symbolic_system(::MTKCardiacModel) = true

function SciMLBase.ODEProblem(m::MTKCardiacModel, tspan::Tuple;
                              u0=CytoZoo.default_initial_state(m), kwargs...)
    return remake(m.prob; u0=u0, tspan=tspan, kwargs...)
end

# ---------------------------------------------------------------------------
# Helper: build MTKCardiacModel from a simplified system + ODEProblem
# ---------------------------------------------------------------------------

function _strip_t_suffix(s::Symbol)
    str = string(s)
    endswith(str, "(t)") ? Symbol(str[1:end-3]) : s
end

function _build_mtk_model(name::Symbol, sys, prob, vm_sym)
    state_syms = _strip_t_suffix.(Symbol.(ModelingToolkit.unknowns(sys)))
    param_syms = _strip_t_suffix.(Symbol.(ModelingToolkit.parameters(sys)))
    vm_idx = findfirst(s -> s == _strip_t_suffix(Symbol(vm_sym)), state_syms)
    isnothing(vm_idx) && error("Vm variable $vm_sym not found in system unknowns")
    return MTKCardiacModel(name, sys, prob, state_syms, param_syms, vm_idx)
end

# ---------------------------------------------------------------------------
# BeelerReuter model
# ---------------------------------------------------------------------------

function _rounded_pulse(t_c; t_start=20.0, duration=5.0, amplitude=10.0, smoothness=0.1)
    t_end = t_start + duration
    rise = 0.5 * (1 + tanh((t_c - t_start) / smoothness))
    fall = 0.5 * (1 - tanh((t_c - t_end) / smoothness))
    return amplitude * rise * fall
end

@register_symbolic _rounded_pulse(t)

@component function _BeelerReuterCompartment(; name, Vm, Cm, currents, I_stim, Ca_dyn)
    @parameters Cm = Cm
    Is = [current.Im for current in currents]
    eqs = [D(Vm) ~ -1 / Cm * (sum(Is) - I_stim)]
    System(eqs, t, [Vm, I_stim], [Cm], systems=[currents..., Ca_dyn], name=name)
end

function _build_beeler_reuter()
    @variables V(t) Ca_i(t)
    V = GlobalScope(V)
    Ca_i = GlobalScope(Ca_i)

    # Fast sodium current gates
    @named m = AlphaBetaGate(
        α_expr=(V + 47.0) / (1.0 - exp(-0.1 * (V + 47.0))),
        β_expr=40.0 * exp(-0.056 * (V + 72.0))
    )
    @named h = AlphaBetaGate(
        α_expr=0.126 * exp(-0.25 * (V + 77.0)),
        β_expr=1.7 / (1.0 + exp(-0.082 * (V + 22.5)))
    )
    @named j = AlphaBetaGate(
        α_expr=0.055 * exp(-0.25 * (V + 78)) / (1 + exp(-0.2 * (V + 78))),
        β_expr=0.3 / (1 + exp(-0.1 * (V + 32.0)))
    )
    @named sodium = IonChannelConductance(g_max=4.0, gates=[m, h, j], powers=[3.0, 1.0, 1.0])
    @named sod_current = IonCurrent(Vm=V, gm=sodium, Ex=50.0)

    # Sodium leak
    @named na_leak = LeakCurrent(Vm=V, gL=0.003, EL=50.0)

    # Slow inward current (calcium)
    @named d_gate = AlphaBetaGate(
        α_expr=0.095 * exp(-0.01 * (V - 5.0)) / (1.0 + exp(-0.072 * (V - 5.0))),
        β_expr=0.07 * exp(-0.017 * (V + 44.0)) / (1.0 + exp(0.05 * (V + 44.0)))
    )
    @named f_gate = AlphaBetaGate(
        α_expr=0.012 * exp(-0.008 * (V + 28.0)) / (1.0 + exp(0.15 * (V + 28.0))),
        β_expr=0.0065 * exp(-0.02 * (V + 30.0)) / (1.0 + exp(-0.2 * (V + 30.0)))
    )
    @named slow_inward = IonChannelConductance(g_max=0.09, gates=[d_gate, f_gate])
    @named slow_current = SlowInwardCurrent(Vm=V, Ca_i=Ca_i, gm=slow_inward)

    # Potassium currents
    @named k1_current = TimeIndependentK(Vm=V)

    @named x1 = AlphaBetaGate(
        α_expr=0.0005 * exp(0.083 * (V + 50.0)) / (1.0 + exp(0.057 * (V + 50.0))),
        β_expr=0.0013 * exp(-0.06 * (V + 20.0)) / (1.0 + exp(-0.04 * (V + 20.0)))
    )
    @named x1_current = TimeActivatedOutward(Vm=V, gate=x1, g_x1=0.8)

    # Calcium dynamics
    @named calcium_dyn = CalciumDynamics(Ca_i=Ca_i, I_s=slow_current.Im)

    # Stimulation
    @named stim = StimulationSystem(stim_func=_rounded_pulse)

    # Compose
    @named com = _BeelerReuterCompartment(
        Vm=V, Cm=1.0,
        currents=[sod_current, na_leak, slow_current, k1_current, x1_current],
        I_stim=stim.I_stim,
        Ca_dyn=calcium_dyn
    )
    @named beeler_reuter_system = compose(com, stim)
    sys = structural_simplify(beeler_reuter_system)

    # Initial conditions (Beeler-Reuter 1977 resting state)
    u0_map = Dict(
        Ca_i => 1e-7,
        V => -84.0,
        sod_current.sodium.m.y => 0.011,
        sod_current.sodium.h.y => 0.988,
        sod_current.sodium.j.y => 0.975,
        slow_current.slow_inward.d_gate.y => 0.003,
        slow_current.slow_inward.f_gate.y => 0.994,
        x1_current.x1.y => 0.0001,
    )

    p_map = ModelingToolkit.defaults(sys)
    op = merge(u0_map, p_map)

    prob = ODEProblem(sys, op, (0.0, 1000.0))
    return _build_mtk_model(:BeelerReuter, sys, prob, V)
end

const _BEELER_REUTER_CACHE = Ref{Any}(nothing)

function CytoZoo.BeelerReuter()
    if isnothing(_BEELER_REUTER_CACHE[])
        _BEELER_REUTER_CACHE[] = _build_beeler_reuter()
    end
    return _BEELER_REUTER_CACHE[]
end

end
