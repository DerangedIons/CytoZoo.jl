# Coupling Internals

For contributors, and for anyone who needs to reason about why a coupled model behaves the
way it does. Nothing here is public API — names and structures may change.

The guiding constraint: **all coupling work happens at `couple()` time**. The per-evaluation
path walks a precomputed, concretely-typed plan, with no dynamic dispatch, no allocation, and
no name lookups.

## What `couple` Does

`couple(nodes, edges)` runs a fixed sequence:

1. **Build the node collection.** Nodes become a `NamedTuple` keyed by node name, so
   component lookup is a compile-time field access rather than a dictionary probe.
2. **Partition the edges by type** into shares and connects. This is why the two edge kinds
   can be freely interleaved in one list.
3. **Validate every name.** Each `share` endpoint must resolve through
   [`state_index`](@ref); each `connect` source must resolve through `state_index` and each
   receiver through [`parameter_index`](@ref) plus [`writable_parameters`](@ref). Validation
   depends on those returning `nothing` for unknown names, which is how a typo becomes an
   actionable error instead of a `KeyError` from deep inside layout construction.
4. **Deepcopy the connect receivers**, so staging scratches a private copy rather than the
   caller's model instance.
5. **Resolve share equivalence classes.**
6. **Compute the global state layout.**
7. **Build the per-component execution plan.**

Note that `couple` never calls `num_states` or `num_parameters` on a component. The layout is
derived entirely from `state_names` and `default_initial_state`, which is why participating in
a coupling does not require the atomic tier of the interface.

## Share Resolution: Union-Find

Shares are **not** treated as pairs. Each `(component, state)` key is a node in a union-find
structure, and every `share` edge unions two keys. The result is a set of equivalence classes,
each collapsing to exactly one global slot.

This is what makes multi-way sharing work regardless of declaration order: fanning three edges
out from one owner puts all four endpoints in a single class, and the class becomes one slot.

Each class must declare exactly one owner. Because *every* edge in a class contributes its
`owner=`, a class that ends up with two different owners is contradictory, and `couple`
rejects it at construction. That is the mechanism behind the "fan out, do not chain" rule —
chaining `:A`–`:B` then `:B`–`:C` produces one class with two declared owners.

## The Global Layout

`CouplingLayout` is computed once and stored on the model:

| Field | Contents |
|---|---|
| `num_states` | number of distinct global slots |
| `solution_indices` | per component, the global indices it reads and writes, in that component's own local state order |
| `names` | canonical name of each global slot |
| `name_to_index` | reverse lookup for [`state_index`](@ref) |
| `u0` | assembled default initial global state |
| `operator_order` | component keys in the order their equations run |
| `vm_index` | global slot holding the primary component's transmembrane potential |

`solution_indices` is the crucial one. Each component gets an index vector mapping its
*local* state positions onto *global* slots. Under a share, two components' index vectors
point at the same global slot — which is precisely what "one variable" means.

Naming rules: the first node is the primary component and keeps bare state names; every other
component's states are prefixed `:<component>_<state>`; a shared slot takes the owner's name
unless `share` was given an explicit `name=`. Initial values for a shared slot come from the
owner.

## The Execution Plan

Rather than interpreting the edge lists on every evaluation, `couple` compiles them into a
tuple of per-component entries — a concretely-typed `Tuple`, so the whole walk specialises and
inlines. Each entry carries:

- the component model itself,
- its `solution_indices` view into the global vectors,
- its **frozen indices** — the local positions whose derivative must be zeroed because the
  component is a non-owner of that shared slot,
- its **connect plan**, partitioned into homogeneous overwrite and add lists, once for
  state-sourced edges and again for monitor-sourced ones.

Partitioning the connects by source kind and operation at construction is what keeps the staging
write free of dynamic dispatch: each list is a homogeneous tuple of `(index, dst_param)` pairs
processed by one concrete code path. State indices point into `U`, monitor indices into the
coupling's flat monitor scratch.

Alongside the plan, `couple` builds a **monitor pre-pass**: one entry per component that sources
at least one monitor, carrying that component's state block and its slice of the scratch. A
coupling with no monitor sources gets an empty tuple, so the pre-pass resolves to a
`Tuple{}` method and vanishes.

## Evaluating the Coupled Right-Hand Side

Each call to `(cm)(dU, U, p, t)` first runs the **monitor pre-pass**: for every monitor-sourcing
component, its state slice is copied out of `U` and its `monitor_values!` fills that component's
slice of the scratch. Monitors are algebraic in `(U, t)` and `U` does not change during an
evaluation, so computing them once up front is exactly equivalent to recomputing them at each
receiver, and cheaper when several receivers read the same source.

Then, for every component in `operator_order`:

1. **Stage connect inputs.** State sources are read **live from `U`**, monitor sources from the
   pre-pass scratch, and written into the receiver's parameter slots — overwrites first, then
   additive edges, whose target slots are reset to zero before summing. Because the read is from
   `U` rather than a cached buffer, the receiver always sees the current value.
2. **Run the component**, with each submodel writing into a `view` of the shared `dU` and
   reading a `view` of the shared `U`. A submodel therefore needs no awareness that it is
   coupled — it sees ordinary contiguous state and derivative vectors.
3. **Zero the frozen entries** of `dU`, discarding the non-owner's contribution to any shared
   slot.

The pre-pass running before all staging is why a monitor source may not also *receive* a connect
edge: its monitors would be computed from parameters staged on the previous evaluation. `couple`
rejects that overlap.

`operator_order` places the **owner last** for every shared slot, so the owner's write is the
final one and therefore the one that survives.

That combination — zero the non-owners, order the owner last — is how "which equation governs
this slot" is decided without any splitting schedule. It is a direct consequence of the
monolithic architecture; see [Design Notes](design.md).

## The ForwardDiff Seam

Connect staging passes each value through an internal `_connect_value` hook, as does the monitor
pre-pass when it copies a source's state slice. In the base package it is the identity function.

When ForwardDiff is loaded, the extension adds a method for `Dual` that extracts the primal
and recurses, so nested duals from higher-order differentiation also collapse. Without it,
storing a `Dual` in a receiver's `Float64` parameter vector would error outright.

This single seam is the origin of the automatic-differentiation limitation: the coupling term
is not carried through the dual, so derivatives taken across a `connect` lose it. `share` is
untouched, because it never leaves `U`. See
[Coupling Limitations](../guides/coupling/limitations.md).

## Monitor Aggregation

A `CoupledModel` sums [`num_monitors`](@ref) over components, concatenates
[`monitor_names`](@ref) with the same non-primary prefixing used for states, and slices each
component's own state portion via `solution_indices` before calling that component's
[`monitor_values!`](@ref).

Both the name concatenation and the value assembly walk components in **declaration order**,
not `operator_order`, so names and values stay aligned. This is a real distinction: operator
order is permuted by share ownership, and using it in one place but not the other would
silently misalign the output.

## Why `num_parameters` Is Undefined

A composite has no single flat parameter vector — parameters remain on the components. Rather
than defining a method that only throws, [`num_parameters`](@ref) and
[`parameter_names`](@ref) are simply left undefined for `CoupledModel`. Calling one raises a
`MethodError`, and `hasmethod` reports `false`, so capability detection gets a truthful
answer.

## Where the Code Lives

`src/coupling.jl` is deliberately free of any solver dependency. Layout, naming, validation,
interface queries, and the assembled right-hand side are all testable without a solver; only
the actual `ODEProblem` construction lives in `ext/SciMLBaseExt.jl`.

## See Also

- [Design Notes](design.md) — why this architecture was chosen.
- [Coupling Limitations](../guides/coupling/limitations.md) — the user-visible consequences.
