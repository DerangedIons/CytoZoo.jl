module CouplingExt

using CytoZoo
using OrdinaryDiffEqOperatorSplitting
using SciMLBase: ODEFunction
const OS = OrdinaryDiffEqOperatorSplitting

# Each operator runs one component model on its gathered local state slice. The splitting
# solver gathers `u[solution_indices]` into a local vector (in the component's own state
# order), so the component's functor is used unchanged on the non-spatial path.
#
# `frozen` holds the local indices of shared states this component does NOT own. Their
# derivative is zeroed after the functor runs, so the shared state stays constant through this
# component's substep — it reads the value but never drives it. Only the share owner's operator
# (no frozen indices for that state) advances the shared slot, giving hard-discard semantics:
# the owner's equation governs, the non-owner's is discarded.
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
    isempty(cm.connects) && return GenericSplitFunction(fs, indices)
    syncs = Tuple(_synchronizer(cm, ck) for ck in order)
    return GenericSplitFunction(fs, indices, syncs)
end

# Local state indices `ck` participates in via a share but does not own.
function _frozen_indices(cm::CoupledModel, ck::Symbol)
    snames = CytoZoo.state_names(cm.components[ck])
    frozen = Int[]
    for sh in cm.shares
        sh.owner === ck && continue
        if sh.a === ck
            push!(frozen, findfirst(==(sh.a_state), snames))
        elseif sh.b === ck
            push!(frozen, findfirst(==(sh.b_state), snames))
        end
    end
    return Tuple(frozen)
end

function OS.OperatorSplittingProblem(cm::CoupledModel, tspan; u0 = default_initial_state(cm))
    return OS.OperatorSplittingProblem(CytoZoo.build_split_function(cm), u0, tspan)
end

# Build a splitting algorithm whose inner solvers (read off each Subsystem node) are ordered
# to match the model's internal operator order (owner-last), so callers never have to know
# that order.
function CytoZoo.coupled_algorithm(cm::CoupledModel; scheme = OS.LieTrotterGodunov)
    algs = Tuple(_node_alg(cm, ck) for ck in cm.layout.operator_order)
    return scheme(algs)
end
function _node_alg(cm::CoupledModel, ck::Symbol)
    alg = cm.algs[ck]
    alg === nothing && throw(ArgumentError(
        "subsystem :$ck has no inner solver; pass one to Subsystem(model, alg), e.g. Tsit5()"))
    return alg
end

# --- connect edges (directed dataflow) ---

# Before the receiver steps, write each source state's current global value into the receiver
# model's parameter slot; the receiver's functor reads that slot as an input. `overwrite` edges
# copy their source into the slot; `+` edges sum their sources (the slot is reset to zero first,
# so it holds the per-step sum, not a running total). Edges are pre-partitioned by op so the
# per-substep path is homogeneous — no dynamic dispatch, no allocation.
struct ConnectSync{P <: AbstractVector}
    params::P                            # receiver's parameter vector (mutated in place)
    overwrites::Vector{Tuple{Int, Int}}  # (source global state idx, dest param idx): slot = source
    adds::Vector{Tuple{Int, Int}}        # (source global state idx, dest param idx): slot = Σ sources
end

# Override forward_sync_external! itself (not synchronize_solution_with_parameters!, whose
# generic path routes through child.p — our operators carry NullParameters). Types are strictly
# more specific than OS's generic method so there is no dispatch ambiguity.
function OS.forward_sync_external!(parent::OS.OperatorSplittingIntegrator, child::OS.DEIntegrator, sync::ConnectSync)
    @inbounds begin
        for (_, dst) in sync.adds
            sync.params[dst] = zero(eltype(sync.params))   # reset before summing the sources
        end
        for (src, dst) in sync.overwrites
            sync.params[dst] = CytoZoo.overwrite(sync.params[dst], parent.u[src])
        end
        for (src, dst) in sync.adds
            sync.params[dst] += parent.u[src]
        end
    end
    return nothing
end
OS.backward_sync_external!(::OS.OperatorSplittingIntegrator, ::OS.DEIntegrator, ::ConnectSync) = nothing

_parameter_vector(model) = model.parameters

function _synchronizer(cm::CoupledModel, ck::Symbol)
    overwrites = Tuple{Int, Int}[]
    adds = Tuple{Int, Int}[]
    for cn in cm.connects
        cn.dst === ck || continue
        src_model = cm.components[cn.src]
        src_global = cm.layout.solution_indices[cn.src][CytoZoo.state_index(src_model, cn.src_state)]
        write = (src_global, CytoZoo.parameter_index(cm.components[ck], cn.dst_slot))
        cn.op === (+) ? push!(adds, write) : push!(overwrites, write)
    end
    (isempty(overwrites) && isempty(adds)) && return OS.NoExternalSynchronization()
    return ConnectSync(_parameter_vector(cm.components[ck]), overwrites, adds)
end

end
