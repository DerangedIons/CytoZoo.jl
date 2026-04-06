module MTKCardiacCellModelsExt

using CytoZoo
using MTKCardiacCellModels
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
# Helpers
# ---------------------------------------------------------------------------

function _strip_t_suffix(s::Symbol)
    str = string(s)
    endswith(str, "(t)") ? Symbol(chop(str; tail=3)) : s
end

function _build_mtk_model(name::Symbol, sys, prob, vm_sym)
    state_syms = _strip_t_suffix.(Symbol.(ModelingToolkit.unknowns(sys)))
    param_syms = _strip_t_suffix.(Symbol.(ModelingToolkit.parameters(sys)))
    vm_idx = findfirst(s -> s == _strip_t_suffix(Symbol(vm_sym)), state_syms)
    isnothing(vm_idx) && error("Vm variable $vm_sym not found in system unknowns")
    return MTKCardiacModel(name, sys, prob, state_syms, param_syms, vm_idx)
end

# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

include("MTKCardiacCellModelsExt/beeler_reuter.jl")

CytoZoo.BeelerReuter() = _build_beeler_reuter()

end
