# Operator-splitting model coupling.
#
# This file is pure Julia (no OrdinaryDiffEqOperatorSplitting dependency): it builds a
# `CoupledModel` and the global-state layout that the OS-based solver in
# `ext/CouplingExt.jl` consumes. Layout, naming, and interface queries are fully testable
# without solving.

"""
    AliasSpec

One alias declaration: state `a_state` of component `a` and state `b_state` of component
`b` are the same physical quantity, sharing a single global state slot. `owner` (`a` or
`b`) wins the equation (enforced by operator order in the splitting solver) and supplies
the canonical name and initial value. Build with [`alias`](@ref).
"""
struct AliasSpec
    a::Symbol
    a_state::Symbol
    b::Symbol
    b_state::Symbol
    owner::Symbol
    name::Symbol
end

"""
    alias(a => :state_a, b => :state_b; owner, name = <owner's state>)

Declare that `state_a` of component `a` and `state_b` of component `b` are the same
quantity (one shared global slot). `a` and `b` are component keys (the `NamedTuple` keys
passed to [`couple`](@ref)). `owner` must be `a` or `b`; it decides whose equation wins
and supplies the canonical name (`name`) and initial value for the shared slot.
"""
function alias(
        a_pair::Pair{Symbol, Symbol}, b_pair::Pair{Symbol, Symbol};
        owner::Symbol, name::Symbol = owner === a_pair.first ? a_pair.second : b_pair.second
    )
    a, a_state = a_pair
    b, b_state = b_pair
    owner in (a, b) || throw(ArgumentError("alias owner must be :$a or :$b, got :$owner"))
    return AliasSpec(a, a_state, b, b_state, owner, name)
end

"""
    CrossRef

One read-only cross-reference: before component `dst` steps, the value of `src`'s state
`src_state` (read from the global state vector) is copied into `dst`'s parameter slot
`dst_slot`. Build with [`crossref`](@ref).
"""
struct CrossRef
    src::Symbol
    src_state::Symbol
    dst::Symbol
    dst_slot::Symbol
end

"""
    crossref(src => :state, dst => :param_slot)

Declare a read-only cross-reference so component `dst` can read `src`'s `state` through
its parameter slot `param_slot`. `src` and `dst` are component keys.
"""
crossref(src_pair::Pair{Symbol, Symbol}, dst_pair::Pair{Symbol, Symbol}) =
    CrossRef(src_pair.first, src_pair.second, dst_pair.first, dst_pair.second)

"""
    CouplingLayout

Pre-computed global state layout for a [`CoupledModel`]. Fields:
- `num_states` — number of distinct global slots.
- `solution_indices` — `NamedTuple` mapping each component key to its global index vector,
  in the component's own local state order (what the splitting solver gathers/scatters).
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

Composition of two or more `AbstractCellModel`s into one combined model with a shared
global state vector. Aliased states occupy a single slot; read-only cross-references let
one component read another's state. Build with [`couple`](@ref) and solve via
`OperatorSplittingProblem(coupled, tspan)` (available once
`OrdinaryDiffEqOperatorSplitting` is loaded).

A `CoupledModel` implements the layout/query interface (`num_states`,
[`state_names`](@ref), [`state_index`](@ref), [`default_initial_state`](@ref),
[`transmembrane_potential_index`](@ref)) so couplings nest, but it is **not** a direct
functor — it is integrated by operator splitting, not a single RHS.
"""
struct CoupledModel{Cs <: NamedTuple, AS <: Tuple, RS <: Tuple, L <: CouplingLayout} <: AbstractCardiacCellModel
    components::Cs
    aliases::AS
    refs::RS
    layout::L
end

"""
    couple(components::NamedTuple; aliases = (), refs = ()) -> CoupledModel

Compose `components` (e.g. `(A = modelA, B = modelB)`) into a `CoupledModel`. `aliases` is
a collection of [`alias`](@ref) specs declaring shared states; `refs` is a collection of
[`crossref`](@ref) specs declaring read-only cross-references. Validates all names and
computes the global state layout at construction time.
"""
function couple(components::NamedTuple; aliases = (), refs = ())
    isempty(components) && throw(ArgumentError("couple requires at least one component"))
    alias_tuple = Tuple(aliases)
    ref_tuple = Tuple(refs)
    _validate_specs(components, alias_tuple, ref_tuple)
    layout = _compute_layout(components, alias_tuple)
    return CoupledModel(components, alias_tuple, ref_tuple, layout)
end

# --- validation ---

function _validate_specs(components::NamedTuple, aliases::Tuple, refs::Tuple)
    comp_keys = keys(components)
    for al in aliases
        _check_component(comp_keys, al.a)
        _check_component(comp_keys, al.b)
        _check_state(components[al.a], al.a, al.a_state)
        _check_state(components[al.b], al.b, al.b_state)
    end
    for r in refs
        _check_component(comp_keys, r.src)
        _check_component(comp_keys, r.dst)
        _check_state(components[r.src], r.src, r.src_state)
        parameter_index(components[r.dst], r.dst_slot) === nothing &&
            throw(ArgumentError("crossref target :$(r.dst) has no parameter slot :$(r.dst_slot)"))
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

function _compute_layout(components::NamedTuple, aliases::Tuple)
    comp_keys = keys(components)
    primary = first(comp_keys)

    alias_of = Dict{Tuple{Symbol, Symbol}, AliasSpec}()
    for al in aliases
        alias_of[(al.a, al.a_state)] = al
        alias_of[(al.b, al.b_state)] = al
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
            al = get(alias_of, key, nothing)
            if al === nothing
                gname = ck === primary ? sname : Symbol(ck, :_, sname)
                push!(names, gname)
                push!(u0, ic[li])
                slot_of[key] = length(names)
            else
                partner = key == (al.a, al.a_state) ? (al.b, al.b_state) : (al.a, al.a_state)
                if haskey(slot_of, partner)
                    slot_of[key] = slot_of[partner]   # shared slot already created
                else
                    push!(names, al.name)
                    push!(u0, _alias_initial_value(components, al))
                    slot_of[key] = length(names)
                end
            end
        end
    end

    name_to_index = Dict{Symbol, Int}()
    for (i, nm) in enumerate(names)
        haskey(name_to_index, nm) &&
            throw(ArgumentError("coupled state name collision on :$nm; rename or alias to resolve"))
        name_to_index[nm] = i
    end

    solution_indices = NamedTuple{comp_keys}(
        map(comp_keys) do ck
            model = components[ck]
            Int[slot_of[(ck, s)] for s in state_names(model)]
        end
    )

    operator_order = _operator_order(comp_keys, aliases)

    primary_model = components[primary]
    vm_index = slot_of[(primary, state_names(primary_model)[transmembrane_potential_index(primary_model)])]

    return CouplingLayout(length(names), solution_indices, names, name_to_index, u0, operator_order, vm_index)
end

_coupled_initial_eltype(components::NamedTuple) =
    promote_type(map(c -> eltype(default_initial_state(c)), values(components))...)

function _alias_initial_value(components::NamedTuple, al::AliasSpec)
    owner_state = al.owner === al.a ? al.a_state : al.b_state
    model = components[al.owner]
    return default_initial_state(model)[state_index(model, owner_state)]
end

# Operator-application order: for each alias the non-owner steps before the owner, so the
# owner's write into the shared slot is final. Stable topological sort over component keys.
function _operator_order(comp_keys::NTuple{N, Symbol}, aliases::Tuple) where {N}
    edges = Tuple{Symbol, Symbol}[]   # (before, after)
    for al in aliases
        nonowner = al.owner === al.a ? al.b : al.a
        nonowner == al.owner || push!(edges, (nonowner, al.owner))
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
            throw(ArgumentError("alias owners imply a cyclic operator order; restructure with nested couple(...)"))
        next = first(ready)               # stable: first in original order
        push!(order, next)
        deleteat!(remaining, findfirst(==(next), remaining))
        for (before, after) in edges
            before == next && (indeg[after] -= 1)
        end
    end
    return order
end

# --- AbstractCellModel interface (layout/query methods) ---

num_states(cm::CoupledModel) = cm.layout.num_states
default_initial_state(cm::CoupledModel) = copy(cm.layout.u0)
state_names(cm::CoupledModel) = Tuple(cm.layout.names)
state_index(cm::CoupledModel, name::Symbol) = get(cm.layout.name_to_index, name, nothing)
transmembrane_potential_index(cm::CoupledModel) = cm.layout.vm_index

# A CoupledModel is solved by operator splitting, not as a single RHS.
(::CoupledModel)(du, u, p, t) = throw(
    ArgumentError("CoupledModel is integrated by operator splitting; build OperatorSplittingProblem(coupled, tspan) and solve (requires OrdinaryDiffEqOperatorSplitting)")
)
num_parameters(::CoupledModel) = throw(
    ArgumentError("CoupledModel has no single parameter vector; parameters live on each component")
)

# Solving entry points — methods are added in ext/CouplingExt.jl, which requires
# OrdinaryDiffEqOperatorSplitting to be loaded.
"""
    build_split_function(cm::CoupledModel)

Build the `GenericSplitFunction` (operators in owner-last order, with the global
`solution_indices`) for `cm`. Requires `OrdinaryDiffEqOperatorSplitting`.
"""
function build_split_function end

"""
    coupled_algorithm(cm::CoupledModel, inner; scheme = LieTrotterGodunov)

Build a splitting algorithm whose inner solvers are ordered to match `cm`'s internal
operator order. `inner` is one solver applied to every component, or a `NamedTuple` of
per-component solvers. Requires `OrdinaryDiffEqOperatorSplitting`.
"""
function coupled_algorithm end
