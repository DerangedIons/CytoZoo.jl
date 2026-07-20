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

Shares are resolved as equivalence **classes**, so a state may appear in more than one edge and
the whole group still collapses to a single slot. To share one quantity across three or more
components, fan the edges out from the owner:

```julia
[share(:A => :d, :B => :x; owner = :A), share(:A => :d, :C => :z; owner = :A)]
```

Every non-owner in the class is frozen, and the result does not depend on the order the nodes
were declared in. *Chaining* edges instead (`:A`–`:B`, then `:B`–`:C`) is an error: both edges
name one class but declare different owners, and a shared slot has exactly one governing
equation.
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

Use one op per slot: a single `overwrite`, or a homogeneous group of `+` edges. Mixing
`overwrite` and `+` (or more than one `overwrite`) into the same slot is rejected at
`couple` time. Other `op`s are not supported.

Under an implicit solver the connect input is frozen to its primal within each Newton step
(see `ext/ForwardDiffExt.jl`): a correct fixed point but an approximate Jacobian that omits
the coupling term. Two consequences:

- A *tightly*-coupled stiff `connect` can see degraded Newton convergence (smaller steps or
  failures).
- Any ForwardDiff derivative taken **through** a connect edge loses that coupling term
  entirely, so sensitivity analysis or gradient-based fitting across a `connect` is
  incorrect — not merely inexact.

`share` is unaffected on both counts, since it flows through the state vector `U`.
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

Pre-computed global state layout for a [`CoupledModel`](@ref). Fields:
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
    components = _copy_connect_receivers(components, connects)
    parent, classes = _share_classes(shares)
    layout = _compute_layout(components, parent, classes)
    plan = _build_plan(components, parent, classes, connects, layout)
    return CoupledModel(components, shares, connects, layout, plan)
end

# Deepcopy every component that receives a `connect` edge so the monolithic RHS scratches a
# private parameter vector (see `_connect!`) rather than mutating the caller's model instance.
# Read-only components are shared by reference. (A single `CoupledModel` is still not safe to
# solve concurrently across threads — deepcopy the `cm` per trajectory in an `EnsembleProblem`
# `prob_func`; full in-RHS reentrancy is tracked as a follow-up.)
function _copy_connect_receivers(components::NamedTuple, connects::Tuple)
    isempty(connects) && return components
    dsts = map(cn -> cn.dst, connects)
    return NamedTuple{keys(components)}(
        map(ck -> ck in dsts ? deepcopy(components[ck]) : components[ck], keys(components))
    )
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
        sh.a === sh.b && throw(
            ArgumentError(
                "share cannot merge two states of the same component :$(sh.a); self-shares are undefined"
            ),
        )
        _check_state(components[sh.a], sh.a, sh.a_state)
        _check_state(components[sh.b], sh.b, sh.b_state)
    end
    for cn in connects
        _check_component(comp_keys, cn.src)
        _check_component(comp_keys, cn.dst)
        _check_state(components[cn.src], cn.src, cn.src_state)
        dst = components[cn.dst]
        hasmethod(parameter_index, Tuple{typeof(dst), Symbol}) || throw(
            ArgumentError(
                "connect target :$(cn.dst) ($(typeof(dst))) must implement `parameter_index` to receive a connect edge"
            ),
        )
        _provides_writable_parameters(dst) || throw(
            ArgumentError(
                "connect target :$(cn.dst) ($(typeof(dst))) has no `parameters` field; override `writable_parameters` or add the field to receive a connect edge"
            ),
        )
        parameter_index(dst, cn.dst_slot) === nothing &&
            throw(ArgumentError("connect target :$(cn.dst) has no parameter slot :$(cn.dst_slot)"))
    end
    _check_connect_op_conflicts(connects)
    return nothing
end

# A connect receiver satisfies `writable_parameters` if it has the default `parameters` field
# or overrides the method beyond the `AbstractCellModel` fallback.
_provides_writable_parameters(model) =
    hasproperty(model, :parameters) ||
    which(writable_parameters, Tuple{typeof(model)}) !== which(writable_parameters, Tuple{AbstractCellModel})

# Reject conflicting ops into one receiver slot: a mix of `overwrite`/`+`, or more than one
# `overwrite` (silent last-wins). Homogeneous `+` fan-in is the supported cross-sectional sum.
function _check_connect_op_conflicts(connects::Tuple)
    counts = Dict{Tuple{Symbol, Symbol}, Tuple{Int, Int}}()
    for cn in connects
        key = (cn.dst, cn.dst_slot)
        no, na = get(counts, key, (0, 0))
        counts[key] = cn.op === (+) ? (no, na + 1) : (no + 1, na)
    end
    for ((dst, slot), (no, na)) in counts
        (no > 1 || (no ≥ 1 && na ≥ 1)) && throw(
            ArgumentError(
                "connect slot :$dst.:$slot has conflicting ops; use a single `overwrite` or only `+` edges"
            ),
        )
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

# --- share equivalence classes ---
#
# A `share` edge asserts two states are one variable, so shares must be resolved as
# equivalence *classes*, not pairwise: chained edges (`A.d ≡ B.x` and `B.x ≡ C.z`) name a
# single variable and must collapse to one global slot. Union-find over `(component, state)`
# keys gives that for chains, fan-out from one state, and groups of any size uniformly.
#
# Each class has exactly one owner — the component whose equation drives the slot; every
# other member is frozen. Two edges in one class naming different owners is a contradiction
# and throws.

const StateKey = Tuple{Symbol, Symbol}

"""
    ShareClass

Resolved metadata for one share equivalence class: the component whose equation governs the
shared slot (`owner`), the owner's local state name (`owner_state`, which supplies the slot's
initial value), and the slot's canonical `name`.
"""
struct ShareClass
    owner::Symbol
    owner_state::Symbol
    name::Symbol
end

# Path-compressed find; a key absent from `parent` is its own root.
function _find_root(parent::Dict{StateKey, StateKey}, k::StateKey)
    p = get(parent, k, k)
    p === k && return k
    r = _find_root(parent, p)
    parent[k] = r
    return r
end

# Build the union-find forest over share endpoints plus one `ShareClass` per class root.
# Every share endpoint is seeded into `parent`, so `haskey(parent, key)` is the membership
# test "does this state participate in a share".
function _share_classes(shares::Tuple)
    parent = Dict{StateKey, StateKey}()
    for sh in shares
        a = (sh.a, sh.a_state)
        b = (sh.b, sh.b_state)
        get!(parent, a, a)
        get!(parent, b, b)
        ra, rb = _find_root(parent, a), _find_root(parent, b)
        ra === rb || (parent[ra] = rb)
    end

    classes = Dict{StateKey, ShareClass}()
    for sh in shares
        r = _find_root(parent, (sh.a, sh.a_state))
        owner_state = sh.owner === sh.a ? sh.a_state : sh.b_state
        prev = get(classes, r, nothing)
        if prev === nothing
            classes[r] = ShareClass(sh.owner, owner_state, sh.name)
        elseif prev.owner !== sh.owner
            throw(
                ArgumentError(
                    "shared states :$(prev.owner).:$(prev.owner_state) and :$(sh.owner).:$owner_state " *
                        "are transitively the same slot but declare different owners; a shared slot has " *
                        "exactly one governing equation. Declare every share in the group against the " *
                        "owning component (e.g. share(:owner => :s, :other => :t; owner = :owner)) " *
                        "rather than chaining them."
                ),
            )
        end
    end
    return parent, classes
end

# --- layout ---

function _compute_layout(components::NamedTuple, parent, classes)
    comp_keys = keys(components)
    primary = first(comp_keys)

    ics = map(default_initial_state, components)   # one IC vector per component, computed once

    slot_of = Dict{StateKey, Int}()
    class_slot = Dict{StateKey, Int}()             # class root -> the one global slot it owns
    names = Symbol[]
    u0 = promote_type(map(eltype, values(ics))...)[]

    for ck in comp_keys
        model = components[ck]
        snames = state_names(model)
        ic = ics[ck]
        length(ic) == length(snames) || throw(
            ArgumentError(
                "component :$ck reports $(length(snames)) state names but its " *
                    "default_initial_state has length $(length(ic)); they must agree"
            ),
        )
        for (li, sname) in enumerate(snames)
            key = (ck, sname)
            if !haskey(parent, key)
                gname = ck === primary ? sname : Symbol(ck, :_, sname)
                push!(names, gname)
                push!(u0, ic[li])
                slot_of[key] = length(names)
            else
                r = _find_root(parent, key)
                slot = get(class_slot, r, nothing)
                if slot === nothing
                    cls = classes[r]
                    push!(names, cls.name)
                    push!(u0, ics[cls.owner][state_index(components[cls.owner], cls.owner_state)])
                    class_slot[r] = length(names)
                    slot_of[key] = length(names)
                else
                    slot_of[key] = slot   # class slot already created
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

    operator_order = _operator_order(comp_keys, components, parent, classes)

    primary_model = components[primary]
    vm_index = slot_of[(primary, state_names(primary_model)[transmembrane_potential_index(primary_model)])]

    return CouplingLayout(length(names), solution_indices, names, name_to_index, u0, operator_order, vm_index)
end

# Operator-application order: every non-owner member of a share class steps before that
# class's owner, so the owner's write into the shared slot is final. Stable topological sort
# over component keys; iterating components/states (not the `parent` dict) keeps it
# deterministic.
function _operator_order(comp_keys::NTuple{N, Symbol}, components::NamedTuple, parent, classes) where {N}
    edges = Tuple{Symbol, Symbol}[]   # (before, after)
    for ck in comp_keys
        for sname in state_names(components[ck])
            key = (ck, sname)
            haskey(parent, key) || continue
            owner = classes[_find_root(parent, key)].owner
            ck === owner || push!(edges, (ck, owner))
        end
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
`(src_global_index, dst_param_index)`; `params` is the receiver's `writable_parameters` vector
(a private deepcopy made by `couple`, so staging never mutates the caller's model), or `nothing`
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
_build_plan(components, parent, classes, connects, layout) =
    _entries(components, parent, classes, connects, layout, layout.operator_order, 1)

_entries(components, parent, classes, connects, layout, order, i) =
    i > length(order) ? () :
    (
        _entry(components, parent, classes, connects, layout, order[i]),
        _entries(components, parent, classes, connects, layout, order, i + 1)...,
    )

function _entry(components, parent, classes, connects, layout, ck)
    si = layout.solution_indices
    idxs = si[ck]
    # UnitRange when contiguous (fast view), else the index vector. The Union here is build-time
    # only — the concrete type is captured by `CompEntry`'s `B`, so the functor stays specialized.
    block = (idxs == first(idxs):last(idxs)) ? (first(idxs):last(idxs)) : copy(idxs)
    frozen = _frozen_indices(components, parent, classes, ck)
    ow_l, ad_l = _connect_plan(components, connects, ck)
    toglobal(e) = (si[e[1]][e[2]], e[3])   # (src_sym, src_local, dst_param) -> (src_global, dst_param)
    overwrites = map(toglobal, ow_l) |> Tuple
    adds = map(toglobal, ad_l) |> Tuple
    params = (isempty(ow_l) && isempty(ad_l)) ? nothing : writable_parameters(components[ck])
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

# --- DERIVED monitors ---
#
# A coupling's monitors are its components' monitors concatenated in component-declaration
# order (`keys(cm.components)`, NOT operator order). Non-primary components' monitor names are
# prefixed `:<comp>_<name>`, matching the state-name prefixing in `_compute_layout`. Each
# component's monitor is computed from its own local state slice (`solution_indices[ck]`, global
# indices in local order), not the whole merged `U`.

num_monitors(cm::CoupledModel) = sum(num_monitors, values(cm.components))

function monitor_names(cm::CoupledModel)
    comp_keys = keys(cm.components)
    primary = first(comp_keys)
    names = Symbol[]
    for ck in comp_keys
        for mn in monitor_names(cm.components[ck])
            push!(names, ck === primary ? mn : Symbol(ck, :_, mn))
        end
    end
    return Tuple(names)
end

function monitor_values!(mon, U, t, cm::CoupledModel)
    si = cm.layout.solution_indices
    offset = 0
    for ck in keys(cm.components)
        model = cm.components[ck]
        nm = num_monitors(model)
        if nm > 0
            monitor_values!(view(mon, (offset + 1):(offset + nm)), view(U, si[ck]), t, model)
        end
        offset += nm
    end
    return nothing
end

# Monolithic single-RHS: evaluate every component into the shared dU/U in operator order. This
# is the default solve path — `ODEProblem(cm, tspan)` + `solve` (see ext/SciMLBaseExt.jl).
(cm::CoupledModel)(dU, U, p, t) = (_run!(dU, U, p, t, cm.plan); nothing)

# `num_parameters` and `parameter_names` are deliberately left UNDEFINED for `CoupledModel`: a
# composite has no single flat parameter vector (parameters live on each component). Defining a
# method that only throws would make `hasmethod(num_parameters, Tuple{CoupledModel})` report a
# capability the type does not have — and this codebase feature-detects with `hasmethod` (see
# `_validate_specs`). Leaving it undefined yields a truthful `MethodError` instead.

# --- coupling-semantics helpers (consumed by the monolithic plan builder) ---

# Local state indices of `ck` that belong to a share class `ck` does not own (the monolithic plan
# zeroes these so the class owner's write into the shared slot wins). Driven by the resolved
# equivalence classes, not by individual edges, so every non-owner member of a multi-way share is
# frozen exactly once.
function _frozen_indices(components::NamedTuple, parent, classes, ck::Symbol)
    snames = state_names(components[ck])
    frozen = Int[]
    for (li, sname) in enumerate(snames)
        key = (ck, sname)
        haskey(parent, key) || continue
        cls = classes[_find_root(parent, key)]
        (cls.owner === ck && cls.owner_state === sname) && continue
        push!(frozen, li)
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
