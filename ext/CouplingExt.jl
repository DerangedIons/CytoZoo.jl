module CouplingExt

using CytoZoo
using OrdinaryDiffEqOperatorSplitting
using SciMLBase: ODEFunction
const OS = OrdinaryDiffEqOperatorSplitting

# Each operator runs one component model on its gathered local state slice. The splitting
# solver gathers `u[solution_indices]` into a local vector (in the component's own state
# order), so the component's functor is used unchanged on the non-spatial path.
#
# `frozen` holds the local indices of aliased states this component does NOT own. Their
# derivative is zeroed after the functor runs, so the shared state stays constant through
# this component's substep — it reads the value but never drives it. Only the alias owner's
# operator (no frozen indices for that state) advances the shared slot, giving hard-discard
# semantics: the owner's equation governs, the non-owner's is discarded.
struct ComponentOperator{M, K}
    model::M
    frozen::NTuple{K, Int}
end
function (op::ComponentOperator)(du, u, p, t)
    op.model(du, u, nothing, t)
    @inbounds for i in op.frozen
        du[i] = zero(eltype(du))
    end
    return nothing
end

function CytoZoo.build_split_function(cm::CoupledModel)
    order = cm.layout.operator_order
    fs = Tuple(ODEFunction(ComponentOperator(cm.components[ck], _frozen_indices(cm, ck))) for ck in order)
    indices = Tuple(cm.layout.solution_indices[ck] for ck in order)
    isempty(cm.refs) && return GenericSplitFunction(fs, indices)
    syncs = Tuple(_synchronizer(cm, ck) for ck in order)
    return GenericSplitFunction(fs, indices, syncs)
end

# Local state indices `ck` participates in via an alias but does not own.
function _frozen_indices(cm::CoupledModel, ck::Symbol)
    snames = CytoZoo.state_names(cm.components[ck])
    frozen = Int[]
    for al in cm.aliases
        al.owner === ck && continue
        if al.a === ck
            push!(frozen, findfirst(==(al.a_state), snames))
        elseif al.b === ck
            push!(frozen, findfirst(==(al.b_state), snames))
        end
    end
    return Tuple(frozen)
end

function OS.OperatorSplittingProblem(cm::CoupledModel, tspan; u0 = default_initial_state(cm))
    return OS.OperatorSplittingProblem(CytoZoo.build_split_function(cm), u0, tspan)
end

# Build a splitting algorithm whose inner solvers are ordered to match the model's internal
# operator order (owner-last), so callers never have to know that order. `inner` is one
# solver applied to every component, or a `NamedTuple` of per-component solvers.
function CytoZoo.coupled_algorithm(cm::CoupledModel, inner; scheme = OS.LieTrotterGodunov)
    algs = Tuple(_inner_solver(inner, ck) for ck in cm.layout.operator_order)
    return scheme(algs)
end
_inner_solver(inner, ::Symbol) = inner
_inner_solver(inner::NamedTuple, ck::Symbol) = inner[ck]

# --- read-only cross-references (refs) ---

# Before the receiver steps, copy each source state's current global value into the
# receiver model's parameter slot; the receiver's functor reads that slot as an input.
struct CrossRefSync{P <: AbstractVector}
    params::P                        # receiver's parameter vector (mutated in place)
    writes::Vector{Tuple{Int, Int}}  # (source global state index, destination parameter index)
end

# Override forward_sync_external! itself (not synchronize_solution_with_parameters!, whose
# generic path routes through child.p — our operators carry NullParameters). Types are
# strictly more specific than OS's generic method so there is no dispatch ambiguity.
function OS.forward_sync_external!(parent::OS.OperatorSplittingIntegrator, child::OS.DEIntegrator, sync::CrossRefSync)
    @inbounds for (src, dst) in sync.writes
        sync.params[dst] = parent.u[src]
    end
    return nothing
end
OS.backward_sync_external!(::OS.OperatorSplittingIntegrator, ::OS.DEIntegrator, ::CrossRefSync) = nothing

_parameter_vector(model) = model.parameters

function _synchronizer(cm::CoupledModel, ck::Symbol)
    writes = Tuple{Int, Int}[]
    for r in cm.refs
        r.dst === ck || continue
        src_model = cm.components[r.src]
        src_global = cm.layout.solution_indices[r.src][CytoZoo.state_index(src_model, r.src_state)]
        push!(writes, (src_global, CytoZoo.parameter_index(cm.components[ck], r.dst_slot)))
    end
    isempty(writes) && return OS.NoExternalSynchronization()
    return CrossRefSync(_parameter_vector(cm.components[ck]), writes)
end

end
