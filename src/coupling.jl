# Model coupling.
#
# This file is pure Julia (no solver dependency): it builds a `CoupledModel`, computes the
# global-state layout, and assembles the monolithic single-RHS functor that solves it. Layout,
# naming, interface queries, and the assembled RHS are all testable without a solver; the actual
# solve goes through `ODEProblem(cm, tspan)` + `solve` (ext/SciMLBaseExt.jl).
#
# Coupling is expressed as a graph: `Subsystem` nodes (a model) joined by directed edges. Two
# edge kinds:
#   - `share`   — two states are one variable (one global slot; the `owner`'s equation governs,
#     and every other member's derivative is either discarded or added to it).
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
the same physical quantity, sharing a single global state slot. `owner` (`a` or `b`) supplies
the governing equation, the canonical name and the initial value. `op` says what becomes of the
**non-owner's** derivative: [`overwrite`](@ref) (the default) discards it, `+` adds it to the
owner's. Build with [`share`](@ref).
"""
struct ShareSpec{OP}
    a::Symbol
    a_state::Symbol
    b::Symbol
    b_state::Symbol
    owner::Symbol
    name::Symbol
    op::OP
end

"""
    share(a => :state_a, b => :state_b; owner, name = <owner's state>, op = overwrite)

Declare that `state_a` of component `a` and `state_b` of component `b` are the **same**
quantity (one shared global slot). `a` and `b` are subsystem names (see [`Subsystem`](@ref)).
`owner` must be `a` or `b`; it supplies the governing equation, the canonical name (`name`) and
the initial value.

Shares are resolved as equivalence **classes**, so a state may appear in more than one edge and
the whole group still collapses to a single slot. To share one quantity across three or more
components, fan the edges out from the owner:

```julia
[share(:A => :d, :B => :x; owner = :A), share(:A => :d, :C => :z; owner = :A)]
```

The result does not depend on the order the nodes were declared in. *Chaining* edges instead
(`:A`–`:B`, then `:B`–`:C`) is an error: both edges name one class but declare different owners,
and a shared slot has exactly one governing equation. A component may contribute at most one
state to a class — two of its states in one class would silently collapse onto one slot, so
`couple` rejects it.

`op` decides what happens to the **non-owner's** derivative:

- [`overwrite`](@ref) (default) — *hard discard*. The non-owner reads the shared value but never
  writes it; only the owner's equation drives the slot. In Modelica's vocabulary this is an
  *across* variable.
- `+` — *contributory*. The non-owner's derivative is **added** to the owner's, so the slot is
  governed by `dw/dt = f_owner + Σ contributions`. This is a *through* (flow) variable, and it is
  how a module contributes a flux to a core model's equation without either model being rewritten.

The two mix freely within one class: an owner, any number of hard-discard members, and any number
of contributors. Several contributors into one slot sum. The sum is **cross-sectional** — every
accumulating slot is reset once per evaluation, so it holds one evaluation's contributions, never
a running total across evaluations.

An accumulating class imposes **no ordering** on its components (each member saves and restores
the slot around its own write), so two components may contribute to each other's slots — a
topology that hard-discard shares reject as a cyclic operator order.

The owner still supplies the initial value; a contributor's own initial value for its local state
is ignored, exactly as a hard-discard non-owner's is.

Because a share flows entirely through the global state vector `U`, it is unaffected by the
automatic-differentiation caveat that applies to [`connect`](@ref) edges — contributory shares
included.
"""
function share(
        a_pair::Pair{Symbol, Symbol}, b_pair::Pair{Symbol, Symbol};
        owner::Symbol, name::Symbol = owner === a_pair.first ? a_pair.second : b_pair.second,
        op = overwrite
    )
    a, a_state = a_pair
    b, b_state = b_pair
    owner in (a, b) || throw(ArgumentError("share owner must be :$a or :$b, got :$owner"))
    op in (overwrite, +) || throw(
        ArgumentError(
            "share: op must be `overwrite` (the owner's equation governs and the non-owner's is " *
                "discarded) or `+` (the non-owner's derivative is added to the owner's); got $op"
        )
    )
    return ShareSpec(a, a_state, b, b_state, owner, name, op)
end

"""
    overwrite(old, new) -> new

Default [`connect`](@ref) operation: replace the receiver slot with the source value (copy).
The other supported operation is `+`, which sums all edges targeting a slot.

Also [`share`](@ref)'s default `op`, where it names the *hard-discard* semantics: the owner's
equation overwrites the shared slot and the non-owner's derivative is dropped.
"""
overwrite(_old, new) = new

"""
    ConnectSpec

One `connect` edge: a directed dataflow link. Before component `dst` steps, the value of
`src`'s observable `src_state` — a **state** (read from the global state vector) or a
**monitor** (recomputed from state each evaluation) — is combined into `dst`'s parameter slot
`dst_slot` via the operation `op` (default [`overwrite`](@ref), i.e. copy). Build with
[`connect`](@ref).
"""
struct ConnectSpec{OP}
    src::Symbol
    src_state::Symbol
    dst::Symbol
    dst_slot::Symbol
    op::OP
end

"""
    connect(src => :observable, dst => :param_slot; op = overwrite)

Declare a directed dataflow edge so component `dst` reads `src`'s `observable` through its
parameter slot `param_slot`. `src` and `dst` are subsystem names.

The source may be either of `src`'s observables:

- a **state** — read live from the global state vector on every evaluation;
- a **monitor** — a DERIVED quantity (see [`monitor_names`](@ref)), recomputed by `src`'s own
  [`monitor_values!`](@ref) on every evaluation.

A monitor source is how a quantity an algebraic law determines — a conserved pool's complement,
say, `ADP = C_A - ATP` — is wired without promoting it to an integrated state. The law stays
inside the model that owns it, so it reads that model's live parameters and needs no
restatement at the edge. Names resolve against states first, then monitors; a component
declaring the same name in both is rejected at `couple` time.

`op` controls how the source value enters the slot each step:
- [`overwrite`](@ref) (default) — the slot is set to the source value (copy).
- `+` — the slot is summed across all `+` edges targeting it; it is reset to zero before
  summing each step, so it holds the per-step sum of its sources (a cross-sectional sum, not a
  running total over time).

Use one op per slot: a single `overwrite`, or a homogeneous group of `+` edges, mixing state
and monitor sources freely. Mixing `overwrite` and `+` (or more than one `overwrite`) into the
same slot is rejected at `couple` time. Other `op`s are not supported.

Monitor sources cost one `monitor_values!` call per sourcing component per evaluation, which
computes that model's **whole** monitor vector — wiring one monitor of a model with many pays
for all of them. A component that sources a monitor may not itself receive a `connect` edge
(rejected at `couple` time): its monitors are computed before the component walk stages any
parameters, so a monitor reading a staged slot would silently see the previous evaluation's
value.

Under an implicit solver the connect input — state- or monitor-sourced alike — is frozen to its
primal within each Newton step (see `ext/ForwardDiffExt.jl`): a correct fixed point but an
approximate Jacobian that omits the coupling term. Two consequences:

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
state vector. `share` edges merge two states into one slot — the owner's equation governing it
alone, or summed with contributions from members declared `op = +`; `connect` edges let one
component read another's state through a parameter slot. Build with [`couple`](@ref).

A `CoupledModel` is a callable `f(dU, U, p, t)` assembled from its submodels, so it is solved
monolithically with a single ODE solver via `ODEProblem(cm, tspan)` + `solve` (requires a
SciMLBase solver stack such as OrdinaryDiffEq). One solver over the full coupled system means no
splitting error and lets a stiff coupling use an implicit method.

It also implements the layout/query interface (`num_states`, [`state_names`](@ref),
[`state_index`](@ref), [`default_initial_state`](@ref),
[`transmembrane_potential_index`](@ref)) so couplings nest.
"""
struct CoupledModel{Cs <: NamedTuple, SS <: Tuple, CS <: Tuple, L <: CouplingLayout, PL <: Tuple, AS <: Tuple, MP <: Tuple, MS <: AbstractVector} <: AbstractCardiacCellModel
    components::Cs
    shares::SS
    connects::CS
    layout::L
    plan::PL
    additive_slots::AS   # global slots an accumulating class owns; `()` unless a share declares op = +
    monitor_plan::MP
    monitor_scratch::MS
    source_scratch::MS
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
    parent, classes, additive = _share_classes(shares)
    layout = _compute_layout(components, parent, classes)
    monitor_plan, monitor_offsets, num_mon = _build_monitor_plan(components, connects, layout)
    plan = _build_plan(components, parent, classes, additive, connects, layout, monitor_offsets)
    additive_slots = _accumulating_slots(plan)
    T = eltype(layout.u0)
    monitor_scratch = zeros(T, num_mon)
    source_scratch = zeros(T, isempty(monitor_plan) ? 0 : maximum(e -> length(e.block), monitor_plan))
    return CoupledModel(
        components, shares, connects, layout, plan, additive_slots,
        monitor_plan, monitor_scratch, source_scratch,
    )
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
        _resolve_source(components[cn.src], cn.src, cn.src_state)
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
    _check_monitor_source_receivers(components, connects)
    return nothing
end

# Resolve a connect source name against the source component's observables: states first, then
# monitors. Returns `(:state, local_state_index)` or `(:monitor, local_monitor_index)`.
# Membership is tested on the name tuples (not `state_index`, which some models implement by
# throwing on an unknown name — see the `nothing` convention in `src/interface.jl`).
function _resolve_source(model, comp::Symbol, name::Symbol)
    is_state = name in state_names(model)
    mi = findfirst(==(name), monitor_names(model))
    is_state && mi !== nothing && throw(
        ArgumentError(
            "component :$comp declares :$name as both a state and a monitor, so a connect " *
                "source naming it is ambiguous; rename one of them"
        ),
    )
    is_state && return (:state, state_index(model, name))
    mi !== nothing && return (:monitor, mi)
    return throw(
        ArgumentError(
            "component :$comp has no state or monitor :$name (states $(state_names(model)), " *
                "monitors $(monitor_names(model)))"
        ),
    )
end

# A monitor-sourcing component's monitors are computed in one pre-pass, before the component walk
# stages any connect input. If such a component also *received* an edge, a monitor of its that
# read the staged slot would see the previous evaluation's value — a silent one-eval lag. The
# overlap is rejected rather than lag-checked, since "does this monitor read that slot" is not
# statically knowable.
function _check_monitor_source_receivers(components::NamedTuple, connects::Tuple)
    dsts = map(cn -> cn.dst, connects)
    for cn in connects
        _resolve_source(components[cn.src], cn.src, cn.src_state)[1] === :monitor || continue
        cn.src in dsts && throw(
            ArgumentError(
                "component :$(cn.src) sources the monitor :$(cn.src_state) and also receives a " *
                    "connect edge; monitors are computed before any parameter is staged, so a " *
                    "monitor reading a staged slot would lag by one evaluation. Source the " *
                    "monitor from a component that receives no edges, or split the model."
            ),
        )
    end
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
    s in state_names(model) && return nothing
    s in monitor_names(model) && throw(
        ArgumentError(
            "component :$comp declares :$s as a monitor, not a state; `share` merges states and " *
                "a monitor has no derivative to own. Use `connect` to read a derived value."
        ),
    )
    return throw(ArgumentError("component :$comp has no state :$s"))
end

# --- share equivalence classes ---
#
# A `share` edge asserts two states are one variable, so shares must be resolved as
# equivalence *classes*, not pairwise: chained edges (`A.d ≡ B.x` and `B.x ≡ C.z`) name a
# single variable and must collapse to one global slot. Union-find over `(component, state)`
# keys gives that for chains, fan-out from one state, and groups of any size uniformly.
#
# Each class has exactly one owner — the component whose equation drives the slot. Every other
# member is either frozen (`op = overwrite`, its derivative discarded) or contributory
# (`op = +`, its derivative added to the owner's); a class with at least one contributor is
# `additive` and takes the accumulate path. Two edges in one class naming different owners is a
# contradiction and throws, whatever their ops — `op` decides what happens to the non-owners,
# not who owns.

const StateKey = Tuple{Symbol, Symbol}

"""
    ShareClass

Resolved metadata for one share equivalence class: the component whose equation governs the
shared slot (`owner`), the owner's local state name (`owner_state`, which supplies the slot's
initial value), the slot's canonical `name`, and whether any member contributes its derivative
additively (`additive`, true as soon as one edge in the class declares `op = +`).
"""
struct ShareClass
    owner::Symbol
    owner_state::Symbol
    name::Symbol
    additive::Bool
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

    # `op` marks the edge's NON-OWNER endpoint. An endpoint that is some other edge's owner
    # endpoint therefore can never be additive: that would need either a second owner in the class
    # (rejected just below) or one component twice in one class (rejected by
    # `_check_one_state_per_class`). "The owner's equation is always kept" is true by construction
    # rather than by check.
    additive = Set{StateKey}()
    hard = Set{StateKey}()
    for sh in shares
        other = sh.owner === sh.a ? (sh.b, sh.b_state) : (sh.a, sh.a_state)
        push!(sh.op === (+) ? additive : hard, other)
        (other in additive && other in hard) && throw(
            ArgumentError(
                "share endpoint :$(other[1]).:$(other[2]) is declared both contributory (op = +) " *
                    "and hard-discard (op = overwrite); an endpoint either adds its derivative to " *
                    "the shared slot or has it discarded, not both. Give every edge naming it the " *
                    "same op."
            ),
        )
    end

    classes = Dict{StateKey, ShareClass}()
    for sh in shares
        r = _find_root(parent, (sh.a, sh.a_state))
        owner_state = sh.owner === sh.a ? sh.a_state : sh.b_state
        is_additive = sh.op === (+)
        prev = get(classes, r, nothing)
        if prev === nothing
            classes[r] = ShareClass(sh.owner, owner_state, sh.name, is_additive)
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
        else
            # `op` is per EDGE, so one contributor makes the whole class additive: the owner and
            # every hard-discard member then ride the same save/restore path in `_run!`.
            classes[r] = ShareClass(prev.owner, prev.owner_state, prev.name, prev.additive | is_additive)
        end
    end

    _check_one_state_per_class(shares, parent)
    return parent, classes, additive
end

# A component may contribute at most one state to a class. Two of its states in one class collapse
# onto a single global slot with one of them silently frozen — degenerate under hard-discard, and
# an outright double-count under `op = +`, where both local positions would save and restore the
# same slot. Iterating `shares` (not `parent`) keeps the message deterministic.
function _check_one_state_per_class(shares::Tuple, parent)
    seen = Dict{Tuple{StateKey, Symbol}, Symbol}()   # (class root, component) -> its state in that class
    for sh in shares, (comp, st) in ((sh.a, sh.a_state), (sh.b, sh.b_state))
        r = _find_root(parent, (comp, st))
        prev = get!(seen, (r, comp), st)
        prev === st || throw(
            ArgumentError(
                "states :$comp.:$prev and :$comp.:$st are transitively the same shared slot; a " *
                    "share merges states of two DIFFERENT components, so one component cannot " *
                    "contribute two states to one class. Check the edges naming :$comp."
            ),
        )
    end
    return nothing
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

# Operator-application order: every non-owner member of a hard-discard share class steps before
# that class's owner, so the owner's write into the shared slot is final. An ACCUMULATING class
# imposes no precedence at all — every member saves and restores the slot around its own write, so
# the total is order-free — which is what lets two components contribute to each other's slots, a
# topology hard-discard shares reject as cyclic. Stable topological sort over component keys;
# iterating components/states (not the `parent` dict) keeps it deterministic.
function _operator_order(comp_keys::NTuple{N, Symbol}, components::NamedTuple, parent, classes) where {N}
    edges = Tuple{Symbol, Symbol}[]   # (before, after)
    for ck in comp_keys
        for sname in state_names(components[ck])
            key = (ck, sname)
            haskey(parent, key) || continue
            cls = classes[_find_root(parent, key)]
            cls.additive && continue
            ck === cls.owner || push!(edges, (ck, cls.owner))
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
            throw(
                ArgumentError(
                    "share owners imply a cyclic operator order; restructure with nested " *
                        "couple(...), or declare the contested shares `op = +` — an accumulating " *
                        "class imposes no ordering"
                ),
            )
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
# and hard-discarded shared-slot derivatives are zeroed so the owner's write (last, by operator
# order) governs. A slot shared with `op = +` instead accumulates: each member saves the running
# value across its own write and adds it back, which sums the contributions and makes the class
# order-free. The plan is built once at couple() time and stored concretely-typed on the
# CoupledModel so the functor specializes and stays allocation-free.

"""
    CompEntry

One component's pre-resolved execution entry in a [`CoupledModel`](@ref)'s monolithic plan.
`block` is the component's slice of the global state (a `UnitRange` when contiguous, else an
index vector); `frozen` are local indices of shared states it neither owns nor contributes to
(zeroed after its functor, so the owner's write wins); `accumulated` are local indices belonging
to an accumulating class — saved before the functor and added back after, so contributions sum
whatever the operator order; `overwrites`/`adds` are state-sourced connect edges resolved to
`(src_global_index, dst_param_index)`, and `monitor_overwrites`/`monitor_adds` the
monitor-sourced ones resolved to `(monitor_scratch_index, dst_param_index)`; `params` is the
receiver's `writable_parameters` vector (a private deepcopy made by `couple`, so staging never
mutates the caller's model), or `nothing` when it has no incoming connect edges.
"""
struct CompEntry{M, B, P, F, AC, OW, AD, MOW, MAD}
    model::M
    block::B
    params::P
    frozen::F
    accumulated::AC
    overwrites::OW
    adds::AD
    monitor_overwrites::MOW
    monitor_adds::MAD
end

"""
    MonEntry

One monitor-sourcing component's entry in a [`CoupledModel`](@ref)'s monitor pre-pass. `block` is
the component's slice of the global state (as in [`CompEntry`](@ref)) and `range` its slice of
the coupling's flat monitor scratch vector.
"""
struct MonEntry{M, B}
    model::M
    block::B
    range::UnitRange{Int}
end

# Identity in base — zero overhead on the Float64 / explicit-solver path. `ext/ForwardDiffExt.jl`
# overrides this for `Dual`s, extracting the primal so a connect input is never stored as a Dual
# in a component's Float64 parameter slot under an implicit solver (it is frozen to its current
# value within the Newton step: correct fixed point, approximate Jacobian).
_connect_value(x) = x

# Build a concretely-typed tuple of entries in operator order (recursion keeps element types
# concrete → the functor specializes and stays allocation-free on the contiguous-block case).
_build_plan(components, parent, classes, additive, connects, layout, mon_offsets) =
    _entries(components, parent, classes, additive, connects, layout, mon_offsets, layout.operator_order, 1)

_entries(components, parent, classes, additive, connects, layout, mon_offsets, order, i) =
    i > length(order) ? () :
    (
        _entry(components, parent, classes, additive, connects, layout, mon_offsets, order[i]),
        _entries(components, parent, classes, additive, connects, layout, mon_offsets, order, i + 1)...,
    )

# The component's slice of the global state: a UnitRange when contiguous (fast view), else the
# index vector. The Union is build-time only — the concrete type is captured by the entry's `B`
# type parameter, so the functor stays specialized.
_block(idxs) = (idxs == first(idxs):last(idxs)) ? (first(idxs):last(idxs)) : copy(idxs)

function _entry(components, parent, classes, additive, connects, layout, mon_offsets, ck)
    si = layout.solution_indices
    frozen = _frozen_indices(components, parent, classes, additive, ck)
    accumulated = _accumulated_indices(components, parent, classes, ck)
    ow_l, ad_l, mow_l, mad_l = _connect_plan(components, connects, ck)
    toglobal(e) = (si[e[1]][e[2]], e[3])        # (src, src_local, dst_param) -> (src_global, dst_param)
    tomonitor(e) = (mon_offsets[e[1]] + e[2], e[3])   # -> (monitor_scratch_index, dst_param)
    overwrites = map(toglobal, ow_l) |> Tuple
    adds = map(toglobal, ad_l) |> Tuple
    monitor_overwrites = map(tomonitor, mow_l) |> Tuple
    monitor_adds = map(tomonitor, mad_l) |> Tuple
    has_edges = !(isempty(ow_l) && isempty(ad_l) && isempty(mow_l) && isempty(mad_l))
    params = has_edges ? writable_parameters(components[ck]) : nothing
    return CompEntry(
        components[ck], _block(si[ck]), params, frozen, accumulated,
        overwrites, adds, monitor_overwrites, monitor_adds,
    )
end

# Every global slot an accumulating class owns, collected from the very plan the walk will use, so
# the pre-walk reset and the per-component save/restore cannot disagree about which slots
# accumulate. Deduplicated — several components save and restore one slot, but it is reset once.
# `()` unless some share declares `op = +`, which is what makes `_reset_accumulators!` vanish.
function _accumulating_slots(plan::Tuple)
    slots = Int[]
    for e in plan, i in e.accumulated
        g = e.block[i]
        g in slots || push!(slots, g)
    end
    return Tuple(slots)
end

# Components sourcing at least one monitor, each with its slice of the coupling's flat monitor
# scratch. Order follows first appearance in the edge list; `_entry` resolves monitor edges
# against the same `offsets` map, so the two stay consistent. Empty when no edge sources a
# monitor — then the pre-pass is a `Tuple{}` method and compiles away entirely.
function _build_monitor_plan(components::NamedTuple, connects::Tuple, layout)
    order = Symbol[]
    for cn in connects
        _resolve_source(components[cn.src], cn.src, cn.src_state)[1] === :monitor || continue
        cn.src in order || push!(order, cn.src)
    end
    offsets = Dict{Symbol, Int}()
    total = 0
    for ck in order
        offsets[ck] = total
        total += num_monitors(components[ck])
    end
    return _mon_entries(components, layout, order, offsets, 1), offsets, total
end

_mon_entries(components, layout, order, offsets, i) =
    i > length(order) ? () :
    (
        _mon_entry(components, layout, offsets, order[i]),
        _mon_entries(components, layout, order, offsets, i + 1)...,
    )

function _mon_entry(components, layout, offsets, ck)
    off = offsets[ck]
    return MonEntry(components[ck], _block(layout.solution_indices[ck]), (off + 1):(off + num_monitors(components[ck])))
end

# Save the running value of every accumulating slot this component is about to overwrite. An
# `NTuple` in `eltype(dU)` — a stack temporary, so a `Dual` stays a `Dual`. A share must not lose
# derivative information the way a `connect` input does at `_connect_value`: flowing through `U`
# untouched by any Float64 scratch is the whole reason to prefer `share` under AD.
_save_accumulated(dU, block, ::Tuple{}) = ()
@inline _save_accumulated(dU, block, accumulated::Tuple) =
    map(i -> (@inbounds dU[block[i]]), accumulated)

# Add the saved running value back onto what the component just wrote, turning its wholesale write
# of the `du` view into an accumulation. Recursive on `Base.tail` so the index/value tuple pair
# unrolls with no dynamic tuple indexing.
_add_back!(dU, block, ::Tuple{}, ::Tuple{}) = nothing
@inline function _add_back!(dU, block, accumulated::Tuple, saved::Tuple)
    @inbounds dU[block[first(accumulated)]] += first(saved)
    return _add_back!(dU, block, Base.tail(accumulated), Base.tail(saved))
end

# Evaluate the assembled RHS: walk the plan (fully unrolled), each component staging its connect
# inputs, saving any accumulating slot it is about to overwrite, writing into its view of dU/U,
# zeroing its hard-discarded shared-slot derivatives, then adding the saved running sum back. The
# freeze must stay BETWEEN the model call and the add-back: reversed, a hard-discard member of an
# accumulating class would zero the running sum instead of just its own contribution.
_run!(dU, U, p, t, ::Tuple{}, mon) = nothing
function _run!(dU, U, p, t, plan, mon)
    e = first(plan)
    _connect!(U, mon, e.params, e.overwrites, e.adds, e.monitor_overwrites, e.monitor_adds)
    saved = _save_accumulated(dU, e.block, e.accumulated)
    e.model(view(dU, e.block), view(U, e.block), p, t)
    @inbounds for i in e.frozen
        dU[e.block[i]] = zero(eltype(dU))
    end
    _add_back!(dU, e.block, e.accumulated, saved)
    return _run!(dU, U, p, t, Base.tail(plan), mon)
end

# `dU` arrives dirty — the solver does not zero it — so every accumulating slot is reset once,
# before the walk. Same zero-then-accumulate discipline `_connect!` uses for an `op = +` parameter
# slot: the result is a cross-sectional sum of one evaluation's contributions, never a running
# total across evaluations.
_reset_accumulators!(dU, ::Tuple{}) = nothing
@inline function _reset_accumulators!(dU, slots::Tuple)
    @inbounds dU[first(slots)] = zero(eltype(dU))
    return _reset_accumulators!(dU, Base.tail(slots))
end

# Monitor pre-pass: fill the flat monitor scratch from each monitor-sourcing component's own
# `monitor_values!`. Monitors are algebraic in (U, t) and U does not change during an evaluation,
# so computing them once up front is exactly equivalent to recomputing per receiver, and cheaper.
# The source's state slice is copied through `_connect_value` first, which is what keeps a derived
# input frozen to its primal inside an implicit solver's Newton step (see ext/ForwardDiffExt.jl)
# and what keeps the scratch a plain real vector under ForwardDiff.
_monitors!(U, t, ::Tuple{}, mon, scratch) = nothing
function _monitors!(U, t, mplan, mon, scratch)
    e = first(mplan)
    n = length(e.block)
    @inbounds for i in 1:n
        scratch[i] = _connect_value(U[e.block[i]])
    end
    monitor_values!(view(mon, e.range), view(scratch, 1:n), t, e.model)
    return _monitors!(U, t, Base.tail(mplan), mon, scratch)
end

# Stage connect inputs into the receiver's parameter slots: state sources read live from U,
# monitor sources read the pre-pass scratch (already primal). `overwrite` = copy, `+` = a
# cross-sectional sum, reset to zero then accumulated each eval. Both source kinds resolve to
# homogeneous `(index, dst_param)` lists, so every loop below is a single concrete code path.
_connect!(U, mon, ::Nothing, ::Tuple{}, ::Tuple{}, ::Tuple{}, ::Tuple{}) = nothing
function _connect!(U, mon, params, overwrites, adds, monitor_overwrites, monitor_adds)
    @inbounds begin
        for (_, d) in adds
            params[d] = zero(eltype(params))
        end
        for (_, d) in monitor_adds
            params[d] = zero(eltype(params))
        end
        for (s, d) in overwrites
            params[d] = _connect_value(U[s])
        end
        for (s, d) in monitor_overwrites
            params[d] = mon[s]
        end
        for (s, d) in adds
            params[d] += _connect_value(U[s])
        end
        for (s, d) in monitor_adds
            params[d] += mon[s]
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

# Monolithic single-RHS: resolve every monitor-sourced connect input, then evaluate every
# component into the shared dU/U in operator order. This is the default solve path —
# `ODEProblem(cm, tspan)` + `solve` (see ext/SciMLBaseExt.jl).
function (cm::CoupledModel)(dU, U, p, t)
    _monitors!(U, t, cm.monitor_plan, cm.monitor_scratch, cm.source_scratch)
    _reset_accumulators!(dU, cm.additive_slots)
    _run!(dU, U, p, t, cm.plan, cm.monitor_scratch)
    return nothing
end

# `num_parameters` and `parameter_names` are deliberately left UNDEFINED for `CoupledModel`: a
# composite has no single flat parameter vector (parameters live on each component). Defining a
# method that only throws would make `hasmethod(num_parameters, Tuple{CoupledModel})` report a
# capability the type does not have — and this codebase feature-detects with `hasmethod` (see
# `_validate_specs`). Leaving it undefined yields a truthful `MethodError` instead.

# --- coupling-semantics helpers (consumed by the monolithic plan builder) ---

# Local state indices of `ck` that belong to a share class `ck` neither owns nor contributes to
# (the monolithic plan zeroes these so the class owner's write into the shared slot wins). Driven
# by the resolved equivalence classes, not by individual edges, so every hard-discard member of a
# multi-way share is frozen exactly once. `additive` carries the per-ENDPOINT ops, since one class
# may mix contributors and hard-discard members.
function _frozen_indices(components::NamedTuple, parent, classes, additive, ck::Symbol)
    snames = state_names(components[ck])
    frozen = Int[]
    for (li, sname) in enumerate(snames)
        key = (ck, sname)
        haskey(parent, key) || continue
        cls = classes[_find_root(parent, key)]
        (cls.owner === ck && cls.owner_state === sname) && continue   # the owner's write governs
        key in additive && continue                                   # a contributor ADDS its write
        push!(frozen, li)
    end
    return Tuple(frozen)
end

# Local state indices of `ck` belonging to an ACCUMULATING share class — whatever `ck`'s role in
# it. Owner, hard-discard member and contributor alike must save the slot's running value across
# their own write, because a component functor fills its whole `du` view and would otherwise
# clobber it. That uniformity is what makes an accumulating class independent of operator order.
# Empty unless some share declares `op = +`, so the save/restore compiles away on a hard-discard
# coupling.
function _accumulated_indices(components::NamedTuple, parent, classes, ck::Symbol)
    accumulated = Int[]
    for (li, sname) in enumerate(state_names(components[ck]))
        key = (ck, sname)
        haskey(parent, key) || continue
        classes[_find_root(parent, key)].additive || continue
        push!(accumulated, li)
    end
    return Tuple(accumulated)
end

# Resolve the connect edges targeting `ck` into four vectors of
# (src::Symbol, src_local::Int, dst_param::Int) tuples, partitioned by source kind and then by op.
# The source is named as (component, local index); `_entry` maps a state local index to a
# global-state index via `solution_indices` and a monitor local index to a monitor-scratch index
# via the offsets map. Partitioning up front keeps every per-eval write homogeneous (no dispatch).
function _connect_plan(components::NamedTuple, connects::Tuple, ck::Symbol)
    overwrites = Tuple{Symbol, Int, Int}[]
    adds = Tuple{Symbol, Int, Int}[]
    monitor_overwrites = Tuple{Symbol, Int, Int}[]
    monitor_adds = Tuple{Symbol, Int, Int}[]
    for cn in connects
        cn.dst === ck || continue
        kind, local_index = _resolve_source(components[cn.src], cn.src, cn.src_state)
        entry = (cn.src, local_index, parameter_index(components[ck], cn.dst_slot))
        if kind === :monitor
            cn.op === (+) ? push!(monitor_adds, entry) : push!(monitor_overwrites, entry)
        else
            cn.op === (+) ? push!(adds, entry) : push!(overwrites, entry)
        end
    end
    return overwrites, adds, monitor_overwrites, monitor_adds
end
