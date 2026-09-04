# ---------------------------------------------------------------------------
# State clamps — hold selected states fixed while the rest of the model evolves
# ---------------------------------------------------------------------------
#
# A clamp is a wrapper model, not a model feature: it calls the base RHS and then zeroes
# `du` at the held indices. That needs nothing but the functor plus the core state-side
# methods every `AbstractCellModel` already implements, so one implementation covers every
# model in the zoo — including a `CoupledModel`.
#
# The wrapper carries the held VALUES as well as the indices, so `default_initial_state`
# returns a seeded vector and the level a state is held at cannot disagree with the level it
# starts from.

"""
    ClampedCell(base, indices, values)
    ClampedCell(base, indices)

Wrap a cell model and hold the states at `indices` fixed at `values`, letting every other
state evolve normally: a voltage clamp, an ion clamp (`[Na]i`, `[Ca]i`), an `[ATP]i` clamp.
The two-argument form takes the held values from `default_initial_state(base)`, which holds
those states at baseline.

A clamp is **freeze plus seed**. The functor calls the base RHS and then writes `du[i] = 0`
at each held index, so a held state keeps whatever value it started with;
[`default_initial_state`](@ref) supplies that value from `values`. Both halves come from this
one object, so they cannot drift apart. To clamp by state name — and to seed a state vector
you already have, mid-protocol — use [`clamp_states`](@ref).

A clamped state is an algebraic input, not a conserved species: the clamp implicitly injects
or removes whatever flux is needed to hold it, exactly as a voltage clamp sources current.
Mass balance for that species no longer holds, by construction.

The base model's interface is forwarded unchanged — states, parameters, named lookup,
monitors, [`writable_parameters`](@ref), and Rush-Larsen. Model-specific accessors are not
(nothing generic can know about them); reach the wrapped model with [`base_model`](@ref).

# Examples

```julia
c = ClampedCell(ToRORd(), (1,))              # voltage clamp at the resting potential
c = ClampedCell(ToRORd(), (1,), (-20.0,))    # voltage clamp at -20 mV
solve(ODEProblem(c, default_initial_state(c), (0.0, 1000.0), nothing), FBDF())
```

!!! warning "Under a contributory `share`, clamp the coupling, not the component"
    Clamping a component of a [`couple`](@ref) zeroes only that component's own write. A
    contributory share (`op = +`) adds the other members' derivatives back afterwards, so the
    global slot still moves. Wrap the [`CoupledModel`](@ref) instead to hold it.

See also [`clamp_states`](@ref), [`base_model`](@ref).
"""
struct ClampedCell{M <: AbstractCardiacCellModel, N, T} <: AbstractCardiacCellModel
    base::M
    indices::NTuple{N, Int}
    values::NTuple{N, T}

    function ClampedCell(
        base::M,
        indices::NTuple{N, Int},
        values::NTuple{N, T},
    ) where {M <: AbstractCardiacCellModel, N, T}
        _check_clamp_indices(indices, length(default_initial_state(base)), M)
        return new{M, N, T}(base, indices, values)
    end
end

function _check_clamp_indices(indices, n::Int, M::Type)
    allunique(indices) ||
        throw(ArgumentError("clamped state indices must be unique, got $indices"))
    for i in indices
        1 <= i <= n || throw(
            ArgumentError("clamped state index $i is out of range for $M ($n states)"),
        )
    end
    return nothing
end

function ClampedCell(base::AbstractCardiacCellModel, indices::NTuple{N, Int}) where {N}
    u0 = default_initial_state(base)
    # Checked before the read, so an out-of-range index reports itself rather than surfacing as
    # a `BoundsError` from seeding the held values.
    _check_clamp_indices(indices, length(u0), typeof(base))
    return ClampedCell(base, indices, ntuple(k -> u0[indices[k]], N))
end

ClampedCell(base::AbstractCardiacCellModel, index::Integer) = ClampedCell(base, (Int(index),))
ClampedCell(base::AbstractCardiacCellModel, indices) =
    ClampedCell(base, Tuple(Int(i) for i in indices))
ClampedCell(base::AbstractCardiacCellModel, indices, values) =
    ClampedCell(base, Tuple(Int(i) for i in indices), Tuple(values))

# ---------------------------------------------------------------------------
# Functor — base RHS, then hold
# ---------------------------------------------------------------------------

function (c::ClampedCell)(du, u, p, t)
    c.base(du, u, p, t)
    @inbounds for i in c.indices
        du[i] = zero(eltype(du))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Interface forwarding — a clamp changes `du`, nothing else
# ---------------------------------------------------------------------------
#
# `num_parameters` / `parameter_names` are forwarded like the rest, which is right for the
# atomic models that have them and costs one thing for the composites that do not: a clamped
# `CoupledModel` answers `hasmethod` with `true` and then raises the base model's `MethodError`
# when called. The error names the wrapped type, so it still reads truthfully; capability
# detection on a wrapper should ask `base_model(c)`.

for f in (
    :num_states,
    :num_parameters,
    :transmembrane_potential_index,
    :state_names,
    :parameter_names,
    :writable_parameters,
    :num_monitors,
    :monitor_names,
)
    @eval $f(c::ClampedCell) = $f(c.base)
end

state_index(c::ClampedCell, name::Symbol) = state_index(c.base, name)
parameter_index(c::ClampedCell, name::Symbol) = parameter_index(c.base, name)
monitor_values!(mon, u, t, c::ClampedCell) = monitor_values!(mon, u, t, c.base)

# The one method a clamp does not forward verbatim: the held states start at `values`, so the
# seed and the hold are set from the same object.
function default_initial_state(c::ClampedCell)
    u = default_initial_state(c.base)
    @inbounds for (k, i) in enumerate(c.indices)
        u[i] = c.values[k]
    end
    return u
end

# Rush-Larsen writes `u_new` directly rather than through `du`, so zeroing derivatives would
# never reach it — restore each held state from `u` after the base step instead.
has_rush_larsen(c::ClampedCell) = has_rush_larsen(c.base)

function rush_larsen_step!(u_new, u, p, t, dt, c::ClampedCell)
    rush_larsen_step!(u_new, u, p, t, dt, c.base)
    @inbounds for i in c.indices
        u_new[i] = u[i]
    end
    return nothing
end

"""
    base_model(model::AbstractCellModel) -> AbstractCellModel

The model underneath any wrappers — the wrapped model for a [`ClampedCell`](@ref)
(recursively, for nested clamps), and `model` itself for everything else.

A wrapper forwards the CytoZoo interface, but it cannot forward accessors it has never heard
of. Unwrap to reach a model's own API:

```julia
c = ClampedCell(ToRORd(), (1,))
base_model(c).celltype
```
"""
base_model(model::AbstractCellModel) = model
base_model(c::ClampedCell) = base_model(c.base)

"""
    clamp_states(model, u0 = default_initial_state(model); state = value, ...) -> (clamped, u)

Clamp states **by name**, returning the [`ClampedCell`](@ref) and a copy of `u0` seeded to the
clamp values. Returning both is the point: a value in `u` without the matching hold decays
back to equilibrium, and a hold without the seed pins the state at baseline. Set from one
call, they cannot disagree.

Pass `u0` explicitly to continue a multi-segment protocol from the previous segment's final
state; the seed is re-applied to the state you hand in.

Names resolve against [`state_names`](@ref), so an unknown name raises an `ArgumentError`
naming the states that do exist. Values are converted to `eltype(u0)`, which keeps a `Float32`
model in `Float32`.

# Examples

```julia
# hold [Na]i at 20 mM while the rest of the cell runs
c, u = clamp_states(model; nai = 20.0)
sol = solve(ODEProblem(c, u, (0.0, 60_000.0), nothing), FBDF())

# second segment: same hold, new level, continuing from where the first ended
c2, u2 = clamp_states(model, sol.u[end]; nai = 7.5)
```

A qualified name reaches a component's state on a [`CoupledModel`](@ref) —
`clamp_states(cm; mito_cai = 0.1)` — which is also how to hold a coupled state that several
components write.

See also [`ClampedCell`](@ref), [`base_model`](@ref).
"""
function clamp_states(
    model::AbstractCardiacCellModel,
    u0::AbstractVector = default_initial_state(model);
    clamps...,
)
    isempty(clamps) &&
        throw(ArgumentError("clamp_states needs at least one `state = value` pair"))
    u = collect(u0)
    names = state_names(model)
    indices = Int[]
    values = eltype(u)[]
    for (name, value) in pairs(clamps)
        # Resolved against `state_names` rather than `state_index`: names are the core-tier
        # method every model implements, and a model whose `state_index` throws on an unknown
        # name (ToRORd's Dict lookup does) would bury this message under a `KeyError`.
        i = findfirst(==(name), names)
        i === nothing && throw(
            ArgumentError("$(typeof(model)) has no state :$name (states: $names)"),
        )
        u[i] = value
        push!(indices, i)
        push!(values, u[i])
    end
    return ClampedCell(model, Tuple(indices), Tuple(values)), u
end
