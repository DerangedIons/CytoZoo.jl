# Model coupling.
#
# This file is pure Julia (no solver dependency): it builds a `CoupledModel`, computes the
# global-state layout, and assembles the monolithic single-RHS functor that solves it. Layout,
# naming, interface queries, and the assembled RHS are all testable without a solver; the actual
# solve goes through `ODEProblem(cm, tspan)` + `solve` (ext/SciMLBaseExt.jl).
#
# Coupling is expressed as a graph: `Subsystem` nodes (a model) joined by directed edges. Two
# edge kinds:
#   - `share`   — two states are one variable (one global slot; the `owner`'s equation governs).
#   - `connect` — a directed dataflow edge carrying an operation (copy by default): a source
#     state is written into a receiver's parameter slot before the receiver steps.

"""
    Subsystem(model; name = gensym(:subsystem))

A graph node: a component `model` (any `AbstractCellModel`) plus a `name` used to reference it
in edges. Pass an explicit `name` whenever edges reference the node — the `gensym` default is
only ergonomic for a single-node graph.
"""
struct Subsystem{M}
    name::Symbol
    model::M
end
Subsystem(model; name::Symbol = gensym(:subsystem)) = Subsystem(name, model)

"""
    ShareSpec

One `share` edge: state `a_state` of component `a` and state `b_state` of component `b` are
the same physical quantity, sharing a single global state slot. `owner` (`a` or `b`) wins the
equation (its write into the shared slot is final, by owner-last operator order) and supplies
the canonical name and initial value. Build with [`share`](@ref).
"""
struct ShareSpec
    a::Symbol
    a_state::Symbol
    b::Symbol
    b_state::Symbol
    owner::Symbol
    name::Symbol
end

"""
    share(a => :state_a, b => :state_b; owner, name = <owner's state>)

Declare that `state_a` of component `a` and `state_b` of component `b` are the **same**
quantity (one shared global slot). `a` and `b` are subsystem names (see [`Subsystem`](@ref)).
`owner` must be `a` or `b`; it decides whose equation governs the shared slot (the other's is
hard-discarded) and supplies the canonical name (`name`) and initial value.
"""
function share(
        a_pair::Pair{Symbol, Symbol}, b_pair::Pair{Symbol, Symbol};
        owner::Symbol, name::Symbol = owner === a_pair.first ? a_pair.second : b_pair.second
    )
    a, a_state = a_pair
    b, b_state = b_pair
    owner in (a, b) || throw(ArgumentError("share owner must be :$a or :$b, got :$owner"))
    return ShareSpec(a, a_state, b, b_state, owner, name)
end

"""
    overwrite(old, new) -> new

Default [`connect`](@ref) operation: replace the receiver slot with the source value (copy).
The other supported operation is `+`, which sums all edges targeting a slot.
"""
overwrite(_old, new) = new

"""
    ConnectSpec

One `connect` edge: a directed dataflow link. Before component `dst` steps, the value of
`src`'s state `src_state` (read from the global state vector) is combined into `dst`'s
parameter slot `dst_slot` via the operation `op` (default [`overwrite`](@ref), i.e. copy).
Build with [`connect`](@ref).
"""
struct ConnectSpec{OP}
    src::Symbol
    src_state::Symbol
    dst::Symbol
    dst_slot::Symbol
    op::OP
end

"""
    connect(src => :state, dst => :param_slot; op = overwrite)

Declare a directed dataflow edge so component `dst` reads `src`'s `state` through its
parameter slot `param_slot`. `src` and `dst` are subsystem names. `op` controls how the
source value enters the slot each step:
- [`overwrite`](@ref) (default) — the slot is set to the source value (copy).
- `+` — the slot is summed across all `+` edges targeting it; it is reset to zero before
  summing each step, so it holds the per-step sum of its sources (a cross-sectional sum, not a
  running total over time).

Use one op per slot; mixing `overwrite` and `+` into the same slot is not meaningful. Other
`op`s are not supported.
"""
function connect(src_pair::Pair{Symbol, Symbol}, dst_pair::Pair{Symbol, Symbol}; op = overwrite)
    op in (overwrite, +) || throw(
        ArgumentError(
            "connect: op must be `overwrite` (copy) or `+` (sum into the slot); got $op"
        )
    )
    return ConnectSpec(src_pair.first, src_pair.second, dst_pair.first, dst_pair.second, op)
end

"""
    CouplingLayout

Pre-computed global state layout for a [`CoupledModel`]. Fields:
- `num_states` — number of distinct global slots.
- `solution_indices` — `NamedTuple` mapping each component key to its global index vector,
  in the component's own local state order (the global slots each component writes its `view` of `dU`/`U` into).
- `names` — canonical name per global slot.
- `name_to_index` — reverse lookup.
- `u0` — default initial global state.
- `operator_order` — component keys in operator-application order (owner steps last).
- `vm_index` — global slot of the primary component's transmembrane potential.
"""
struct CouplingLayout{SI <: NamedTuple, U <: AbstractVector}
    num_states::Int
    solution_indices::SI
    names::Vector{Symbol}
    name_to_index::Dict{Symbol, Int}
    u0::U
    operator_order::Vector{Symbol}
    vm_index::Int
end

"""
    CoupledModel <: AbstractCardiacCellModel

Composition of two or more `AbstractCellModel`s into one combined model with a shared global
state vector. `share` edges merge two states into one slot; `connect` edges let one component
read another's state through a parameter slot. Build with [`couple`](@ref).

A `CoupledModel` is a callable `f(dU, U, p, t)` assembled from its submodels, so it is solved
monolithically with a single ODE solver via `ODEProblem(cm, tspan)` + `solve` (requires a
SciMLBase solver stack such as OrdinaryDiffEq). One solver over the full coupled system means no
splitting error and lets a stiff coupling use an implicit method.

It also implements the layout/query interface (`num_states`, [`state_names`](@ref),
[`state_index`](@ref), [`default_initial_state`](@ref),
[`transmembrane_potential_index`](@ref)) so couplings nest.
"""
struct CoupledModel{Cs <: NamedTuple, SS <: Tuple, CS <: Tuple, L <: CouplingLayout, PL <: Tuple} <: AbstractCardiacCellModel
    components::Cs
    shares::SS
    connects::CS
    layout::L
    plan::PL
end

"""
    couple(nodes, edges = ()) -> CoupledModel

Compose `nodes` (a collection of [`Subsystem`](@ref)s) into a `CoupledModel`. `edges` is a
collection mixing [`share`](@ref) and [`connect`](@ref) specs; it is partitioned by type. The
first node is the primary component — its states keep bare names and supply the
transmembrane-potential index. Validates all names and computes the global state layout at
construction time.
"""
function couple(nodes, edges = ())
    isempty(nodes) && throw(ArgumentError("couple requires at least one Subsystem node"))
    subsystems = _subsystems_namedtuple(nodes)
    components = map(s -> s.model, subsystems)
    shares, connects = _split_edges(edges)
    _validate_specs(components, shares, connects)
    layout = _compute_layout(components, shares)
    plan = _build_plan(components, shares, connects, layout)
    return CoupledModel(components, shares, connects, layout, plan)
end

# Build a NamedTuple keyed by each node's name, erroring on duplicates before the build.
function _subsystems_namedtuple(nodes)
    names = Symbol[]
    for n in nodes
        n.name in names && throw(
            ArgumentError(
                "duplicate subsystem name :$(n.name); give each Subsystem a distinct name="
            )
        )
        push!(names, n.name)
    end
    return NamedTuple{Tuple(names)}(Tuple(nodes))
end

# Partition the edge list into share and connect specs by type, erroring on anything else.
function _split_edges(edges)
    shares = ShareSpec[]
    connects = ConnectSpec[]
    for e in edges
        if e isa ShareSpec
            push!(shares, e)
        elseif e isa ConnectSpec
            push!(connects, e)
        else
            throw(
                ArgumentError(
                    "unknown edge type $(typeof(e)); expected a share(...) or connect(...) spec"
                )
            )
        end
    end
    return Tuple(shares), Tuple(connects)
end

# --- validation ---

function _validate_specs(components::NamedTuple, shares::Tuple, connects::Tuple)
    comp_keys = keys(components)
    for sh in shares
        _check_component(comp_keys, sh.a)
        _check_component(comp_keys, sh.b)
        _check_state(components[sh.a], sh.a, sh.a_state)
        _check_state(components[sh.b], sh.b, sh.b_state)
    end
    for cn in connects
        _check_component(comp_keys, cn.src)
        _check_component(comp_keys, cn.dst)
        _check_state(components[cn.src], cn.src, cn.src_state)
        parameter_index(components[cn.dst], cn.dst_slot) === nothing &&
            throw(ArgumentError("connect target :$(cn.dst) has no parameter slot :$(cn.dst_slot)"))
    end
    return nothing
end

_check_component(comp_keys, c::Symbol) =
    c in comp_keys || throw(ArgumentError("unknown component :$c (have $(comp_keys))"))

function _check_state(model, comp::Symbol, s::Symbol)
    s in state_names(model) ||
        throw(ArgumentError("component :$comp has no state :$s"))
    return nothing
end

# --- layout ---

function _compute_layout(components::NamedTuple, shares::Tuple)
    comp_keys = keys(components)
    primary = first(comp_keys)

    share_of = Dict{Tuple{Symbol, Symbol}, ShareSpec}()
    for sh in shares
        share_of[(sh.a, sh.a_state)] = sh
        share_of[(sh.b, sh.b_state)] = sh
    end

    slot_of = Dict{Tuple{Symbol, Symbol}, Int}()
    names = Symbol[]
    u0 = _coupled_initial_eltype(components)[]

    for ck in comp_keys
        model = components[ck]
        snames = state_names(model)
        ic = default_initial_state(model)
        for (li, sname) in enumerate(snames)
            key = (ck, sname)
            sh = get(share_of, key, nothing)
            if sh === nothing
                gname = ck === primary ? sname : Symbol(ck, :_, sname)
                push!(names, gname)
                push!(u0, ic[li])
                slot_of[key] = length(names)
            else
                partner = key == (sh.a, sh.a_state) ? (sh.b, sh.b_state) : (sh.a, sh.a_state)
                if haskey(slot_of, partner)
                    slot_of[key] = slot_of[partner]   # shared slot already created
                else
                    push!(names, sh.name)
                    push!(u0, _share_initial_value(components, sh))
                    slot_of[key] = length(names)
                end
            end
        end
    end

    name_to_index = Dict{Symbol, Int}()
    for (i, nm) in enumerate(names)
        haskey(name_to_index, nm) &&
            throw(ArgumentError("coupled state name collision on :$nm; rename or share to resolve"))
        name_to_index[nm] = i
    end

    solution_indices = NamedTuple{comp_keys}(
        map(comp_keys) do ck
            model = components[ck]
            Int[slot_of[(ck, s)] for s in state_names(model)]
        end
    )

    operator_order = _operator_order(comp_keys, shares)

    primary_model = components[primary]
    vm_index = slot_of[(primary, state_names(primary_model)[transmembrane_potential_index(primary_model)])]

    return CouplingLayout(length(names), solution_indices, names, name_to_index, u0, operator_order, vm_index)
end

_coupled_initial_eltype(components::NamedTuple) =
    promote_type(map(c -> eltype(default_initial_state(c)), values(components))...)

function _share_initial_value(components::NamedTuple, sh::ShareSpec)
    owner_state = sh.owner === sh.a ? sh.a_state : sh.b_state
    model = components[sh.owner]
    return default_initial_state(model)[state_index(model, owner_state)]
end

# Operator-application order: for each share the non-owner steps before the owner, so the
# owner's write into the shared slot is final. Stable topological sort over component keys.
function _operator_order(comp_keys::NTuple{N, Symbol}, shares::Tuple) where {N}
    edges = Tuple{Symbol, Symbol}[]   # (before, after)
    for sh in shares
        nonowner = sh.owner === sh.a ? sh.b : sh.a
        nonowner == sh.owner || push!(edges, (nonowner, sh.owner))
    end
    isempty(edges) && return collect(comp_keys)

    indeg = Dict{Symbol, Int}(c => 0 for c in comp_keys)
    for (_, after) in edges
        indeg[after] += 1
    end
    order = Symbol[]
    remaining = collect(comp_keys)
    while !isempty(remaining)
        ready = filter(c -> indeg[c] == 0, remaining)
        isempty(ready) &&
            throw(ArgumentError("share owners imply a cyclic operator order; restructure with nested couple(...)"))
        next = first(ready)               # stable: first in original order
        push!(order, next)
        deleteat!(remaining, findfirst(==(next), remaining))
        for (before, after) in edges
            before == next && (indeg[after] -= 1)
        end
    end
    return order
end

# --- monolithic single-RHS assembly ---
#
# Assemble one combined f(dU, U, p, t) from the submodels: each component writes into its own
# view of the shared dU/U, connect inputs are read live from U into receiver parameter slots,
# and non-owned shared-slot derivatives are zeroed so the owner's write (last, by operator
# order) governs. The plan is built once at couple() time and stored concretely-typed on the
# CoupledModel so the functor specializes and stays allocation-free.

"""
    CompEntry

One component's pre-resolved execution entry in a [`CoupledModel`](@ref)'s monolithic plan.
`block` is the component's slice of the global state (a `UnitRange` when contiguous, else an
index vector); `frozen` are local indices of shared states it does not own (zeroed after its
functor, so the owner's write wins); `overwrites`/`adds` are connect edges resolved to
`(src_global_index, dst_param_index)`; `params` is the receiver's parameter vector, or `nothing`
when it has no incoming connect edges.
"""
struct CompEntry{M, B, P, F, OW, AD}
    model::M
    block::B
    params::P
    frozen::F
    overwrites::OW
    adds::AD
end

# Identity in base — zero overhead on the Float64 / explicit-solver path. `ext/ForwardDiffExt.jl`
# overrides this for `Dual`s, extracting the primal so a connect input is never stored as a Dual
# in a component's Float64 parameter slot under an implicit solver (it is frozen to its current
# value within the Newton step: correct fixed point, approximate Jacobian).
_connect_value(x) = x

# Build a concretely-typed tuple of entries in operator order (recursion keeps element types
# concrete → the functor specializes and stays allocation-free on the contiguous-block case).
_build_plan(components, shares, connects, layout) =
    _entries(components, shares, connects, layout, layout.operator_order, 1)

_entries(components, shares, connects, layout, order, i) =
    i > length(order) ? () :
    (
        _entry(components, shares, connects, layout, order[i]),
        _entries(components, shares, connects, layout, order, i + 1)...,
    )

function _entry(components, shares, connects, layout, ck)
    si = layout.solution_indices
    idxs = si[ck]
    block = (idxs == first(idxs):last(idxs)) ? (first(idxs):last(idxs)) : copy(idxs)
    frozen = _frozen_indices(components, shares, ck)
    ow_l, ad_l = _connect_plan(components, connects, ck)
    toglobal(e) = (si[e[1]][e[2]], e[3])   # (src_sym, src_local, dst_param) -> (src_global, dst_param)
    overwrites = map(toglobal, ow_l) |> Tuple
    adds = map(toglobal, ad_l) |> Tuple
    params = (isempty(ow_l) && isempty(ad_l)) ? nothing : components[ck].parameters
    return CompEntry(components[ck], block, params, frozen, overwrites, adds)
end

# Evaluate the assembled RHS: walk the plan (fully unrolled), each component staging its connect
# inputs, writing into its view of dU/U, then zeroing its non-owned shared-slot derivatives.
_run!(dU, U, p, t, ::Tuple{}) = nothing
function _run!(dU, U, p, t, plan)
    e = first(plan)
    _connect!(U, e.params, e.overwrites, e.adds)
    e.model(view(dU, e.block), view(U, e.block), p, t)
    @inbounds for i in e.frozen
        dU[e.block[i]] = zero(eltype(dU))
    end
    return _run!(dU, U, p, t, Base.tail(plan))
end

# Stage connect inputs: read source states live from U into the receiver's parameter slots
# (overwrite = copy, + = cross-sectional sum reset to zero then accumulated each eval).
_connect!(U, ::Nothing, ::Tuple{}, ::Tuple{}) = nothing
function _connect!(U, params, overwrites, adds)
    @inbounds begin
        for (_, d) in adds
            params[d] = zero(eltype(params))
        end
        for (s, d) in overwrites
            params[d] = _connect_value(U[s])
        end
        for (s, d) in adds
            params[d] += _connect_value(U[s])
        end
    end
    return nothing
end

# --- AbstractCellModel interface (layout/query methods) ---

num_states(cm::CoupledModel) = cm.layout.num_states
default_initial_state(cm::CoupledModel) = copy(cm.layout.u0)
state_names(cm::CoupledModel) = Tuple(cm.layout.names)
state_index(cm::CoupledModel, name::Symbol) = get(cm.layout.name_to_index, name, nothing)
transmembrane_potential_index(cm::CoupledModel) = cm.layout.vm_index

# Monolithic single-RHS: evaluate every component into the shared dU/U in operator order. This
# is the default solve path — `ODEProblem(cm, tspan)` + `solve` (see ext/SciMLBaseExt.jl).
(cm::CoupledModel)(dU, U, p, t) = (_run!(dU, U, p, t, cm.plan); nothing)
num_parameters(::CoupledModel) = throw(
    ArgumentError("CoupledModel has no single parameter vector; parameters live on each component")
)

# --- coupling-semantics helpers (consumed by the monolithic plan builder) ---

# Local state indices `ck` participates in via a share but does not own (used by the monolithic
# plan to zero the non-owner's contribution to a shared slot, so the owner's write wins).
function _frozen_indices(components::NamedTuple, shares::Tuple, ck::Symbol)
    snames = state_names(components[ck])
    frozen = Int[]
    for sh in shares
        sh.owner === ck && continue
        if sh.a === ck
            push!(frozen, findfirst(==(sh.a_state), snames))
        elseif sh.b === ck
            push!(frozen, findfirst(==(sh.b_state), snames))
        end
    end
    return Tuple(frozen)
end

# Resolve the connect edges targeting `ck` into (overwrites, adds), each a vector of
# (src::Symbol, src_local::Int, dst_param::Int) tuples partitioned by op. The source is named as
# (component, local state index); `_entry` maps it to a global-state index via `solution_indices`
# when building the plan. Partitioning by op keeps the per-eval write homogeneous (no dispatch).
function _connect_plan(components::NamedTuple, connects::Tuple, ck::Symbol)
    overwrites = Tuple{Symbol, Int, Int}[]
    adds = Tuple{Symbol, Int, Int}[]
    for cn in connects
        cn.dst === ck || continue
        entry = (
            cn.src, state_index(components[cn.src], cn.src_state),
            parameter_index(components[ck], cn.dst_slot),
        )
        cn.op === (+) ? push!(adds, entry) : push!(overwrites, entry)
    end
    return overwrites, adds
end
