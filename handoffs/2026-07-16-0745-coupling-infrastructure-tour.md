---
slug: coupling-infrastructure-tour
created: 2026-07-16-0745
status: open
---

# Handoff: A guided tour of CytoZoo's coupling infrastructure

## Goal / why this matters

This is a **guided tour** of CytoZoo's model-coupling subsystem — the machinery that
composes two or more cardiac cell models into one combined model solved by a single ODE
solver. It's the densest part of the codebase, and its rationale has been scattered across
CLAUDE.md, the README, and two design docs that have since been deleted from `handoffs/`.
Read this cold and you should be able to reason about — and extend — `src/coupling.jl`
without re-reverse-engineering the API from the source.

The tour is code-focused: visitor-center overview first, then station-by-station through
each piece, with real `file:line` anchors and verbatim call sites. Everything cited is
current as of this writing (`src/coupling.jl` is 581 lines).

---

## Visitor center — the big idea

**The problem.** You have several standalone cell models (an electrophysiology model, a
mitochondrial model, a redox module…). Each is a self-contained ODE right-hand side. You
want them to run as *one* coupled system — sharing variables, feeding signals to each
other — and you want to solve that system **monolithically**: one solver over the whole
combined state, so there is no operator-splitting error and the whole thing stays
stiff-capable.

**The model.** Coupling is expressed as a **graph**:

- **Nodes** are `Subsystem`s — a cell model plus a `name`.
- **Edges** are directed and come in exactly two kinds:
  - **`share`** — two models' states are declared to be *the same physical variable*, living
    in one global slot. One model **owns** it (its equation governs); the other reads the
    value but its own derivative for that slot is discarded.
  - **`connect`** — a directed dataflow wire: before model B runs, model A's state value is
    written into one of B's *parameter* slots, which B reads live each evaluation.

`couple(nodes, edges)` compiles this graph — validates it, works out the global-state
layout, and pre-bakes an **execution plan** — into a `CoupledModel`. A `CoupledModel` is
itself a real functor `(cm)(dU, U, p, t)`, so it plugs straight into `ODEProblem` +
`solve`, *and* couplings nest (a `CoupledModel` is just another `AbstractCellModel`).

**The one principle to remember.** There is **no switch primitive** anywhere. Every
optional capability — a whole module, a single wire — is turned on by *composing it in* and
off by *leaving it out*. "Switching is composition." We'll come back to this at Station 8,
but keep it in mind the whole way: nothing in this subsystem toggles behavior through a
runtime flag or a `p` argument.

The tour route: the interface it all rests on (S1) → the building blocks you author (S2) →
how `couple()` assembles them (S3) → the two data structures it produces, layout (S4) and
plan (S5) → the semantic heart, share-vs-connect (S6) → how you actually solve and observe
(S7) → the composition principle (S8) → worked examples you can run (S9) → and the wing
that's still under construction.

---

## Station 1 — The type hierarchy & the interface contract

Everything coupling does is built on the plain-model interface in `src/interface.jl`. Two
abstract types form the spine:

```
AbstractCellModel                (interface.jl:36)
└── AbstractCardiacCellModel     (interface.jl:41)
    └── CoupledModel             (coupling.jl:155)   ← a coupling IS a cardiac cell model
```

Note that last line: `CoupledModel <: AbstractCardiacCellModel`. That's what lets couplings
nest and lets a `CoupledModel` be solved by exactly the same path as a leaf model.

**The required interface** — bare `function … end` stubs with no default, so a
non-implementing model errors:

| Method | Where | Meaning |
|---|---|---|
| `num_states` | `interface.jl:52` | total ODE state count |
| `num_parameters` | `interface.jl:61` | atomic models only (see below) |
| `transmembrane_potential_index` | `interface.jl:68` | index of Vm in the state vector |
| `default_initial_state` | `interface.jl:75` | length-`num_states` IC vector |

**The optional interface** — name-based access and traits, most with sensible defaults:

- `state_index` / `parameter_index` (`interface.jl:86,93`) — Symbol → index. No default.
- `state_names` / `parameter_names` (`interface.jl:100,107`).
- `writable_parameters(model) = model.parameters` (`interface.jl:121`) — **the one with a
  default**, and load-bearing for `connect` (below).
- Rush-Larsen: `has_rush_larsen` (default `false`, `interface.jl:132`), `rush_larsen_step!`.
- Monitors (derived observables): `num_monitors` (default `0`, `interface.jl:159`),
  `monitor_names` (default `()`, `:167`), `monitor_values!` (stub, `:175`).

**The atomic-vs-composite distinction** (interface docstring, `interface.jl:11-15`) is the
key idea that makes `CoupledModel` legal. A composite merges its components' states into one
global vector, so it *can* and *does* implement the state-side methods
(`coupling.jl:496-500`). But it has no single flat parameter vector — parameters stay on the
components — so `num_parameters` is required of *atomic* models only. `CoupledModel`
deliberately overrides it to **throw**:

```julia
# coupling.jl:543-545
num_parameters(::CoupledModel) = throw(
    ArgumentError("CoupledModel has no single parameter vector; parameters live on each component")
)
```

**What coupling requires of others.** A `share` needs nothing special — it works purely
through merged state. A `connect` *receiver* must implement exactly two things, and
`couple()` checks both at build time (Station 3):

1. `parameter_index(model, slot)` — to name the target slot.
2. `writable_parameters(model)` — the vector coupling writes into and the model reads.

That contract is stated verbatim in the `writable_parameters` docstring (`interface.jl:114-119`).
The default `writable_parameters` returns `model.parameters`, so any model with a
`parameters` field satisfies half of it for free.

*Questions on the interface before we move to the building blocks?*

---

## Station 2 — Building blocks: nodes & the two edges

These are the three things you author by hand. All live at the top of `src/coupling.jl`.

### `Subsystem` — a graph node

```julia
# coupling.jl:21-25
struct Subsystem{M}
    name::Symbol
    model::M
end
Subsystem(model; name::Symbol = gensym(:subsystem)) = Subsystem(name, model)
```

A node is just a model plus a `name` that edges reference. **Caveat:** the `gensym` default
name is only ergonomic for a *single-node* graph — the moment an edge references a node, you
must pass an explicit `name =` so the symbol is stable and predictable. The **first** node
passed to `couple` is the *primary* component (it keeps bare state names and supplies
`vm_index`; see Station 4).

### `share` — merge two states into one slot

```julia
# coupling.jl:35-42
struct ShareSpec
    a::Symbol; a_state::Symbol
    b::Symbol; b_state::Symbol
    owner::Symbol
    name::Symbol
end

# coupling.jl:52-60 (constructor)
share(:A => :d, :B => :x; owner = :A)          # A.d and B.x are one variable; A governs it
```

`a_state` of component `a` and `b_state` of `b` become the same physical quantity in one
global slot. `owner` (which must be `a` or `b`, validated at `coupling.jl:58`) wins the
equation; its state name is the default canonical `name` for the slot, and its initial value
seeds it. This is **hard-discard**: the non-owner's contribution to that slot's derivative is
zeroed (more at Station 5).

### `connect` — a directed dataflow wire

```julia
# coupling.jl:78-84
struct ConnectSpec{OP}
    src::Symbol; src_state::Symbol
    dst::Symbol; dst_slot::Symbol
    op::OP
end

overwrite(_old, new) = new                     # coupling.jl:68 — the default op (plain copy)

# coupling.jl:107-114 (constructor)
connect(:A => :Vm, :B => :Vm_ext)              # A.Vm → B's parameter slot :Vm_ext, each eval
connect(:R => :J,  :D => :w_in; op = +)        # additive fan-in into one slot
```

Before `dst` steps, `src`'s state value (read live from the global `U`) is combined into
`dst`'s parameter slot `dst_slot` via `op`. `op` must be `overwrite` (default) or `+`
— anything else is rejected at construction (`coupling.jl:108`). `OP` is a type parameter so
the op is captured in the type and the per-eval write stays dispatch-free.

*The rest of the tour is about what `couple()` does with these three. Ready?*

---

## Station 3 — `couple()`: from graph to `CoupledModel`

This is the orchestration hub. The whole body is nine lines, and reading it top-to-bottom is
the fastest way to understand the build pipeline:

```julia
# coupling.jl:172-182
function couple(nodes, edges = ())
    isempty(nodes) && throw(ArgumentError("couple requires at least one Subsystem node"))
    subsystems = _subsystems_namedtuple(nodes)          # keyed by node name
    components = map(s -> s.model, subsystems)
    shares, connects = _split_edges(edges)              # partition by type
    _validate_specs(components, shares, connects)        # fail fast, actionable errors
    components = _copy_connect_receivers(components, connects)   # deepcopy write targets
    layout = _compute_layout(components, shares)         # → Station 4
    plan = _build_plan(components, shares, connects, layout)     # → Station 5
    return CoupledModel(components, shares, connects, layout, plan)
end
```

The helpers, in order:

- **`_subsystems_namedtuple`** (`coupling.jl:198-209`) — builds a `NamedTuple` keyed by node
  name, erroring on duplicate names. Declaration order is preserved and matters (primary =
  first).
- **`_split_edges`** (`coupling.jl:212-229`) — partitions the edge list into a
  `Tuple{ShareSpec…}` and a `Tuple{ConnectSpec…}` by type; any other element type errors.
- **`_validate_specs`** (`coupling.jl:233-266`) — the gatekeeper. For each `share`: both
  components exist, no self-share (`:238`), both states exist. For each `connect`: src
  component/state exist; the dst implements `parameter_index` (`:251`); dst provides
  `writable_parameters` (`:256`, via `_provides_writable_parameters`, `:270-272`); the named
  slot exists (`:261`). Then `_check_connect_op_conflicts` (`:276-291`) rejects mixing
  `overwrite`/`+` into one `(dst, slot)`, or more than one `overwrite` into it (silent
  last-wins). Homogeneous `+` fan-in is allowed.
- **`_copy_connect_receivers`** (`coupling.jl:189-195`) — **deepcopies every component that
  receives a connect edge**. This is why the RHS can scratch a private parameter vector
  without ever mutating the caller's model instance. Read-only components are shared by
  reference (no copy). This is also the root of the thread-safety caveat (see Gotchas).

Then the two products — `layout` and `plan` — are the next two stations.

---

## Station 4 — The global-state layout

`_compute_layout` (`coupling.jl:304-365`) works out, once at build time, how every
component's local states map into the single shared global vector. The result:

```julia
# coupling.jl:129-137
struct CouplingLayout{SI <: NamedTuple, U <: AbstractVector}
    num_states::Int
    solution_indices::SI     # per-component global index vectors, in that component's local order
    names::Vector{Symbol}    # canonical name of each global slot
    name_to_index::Dict{Symbol, Int}
    u0::U                    # merged default initial state
    operator_order::Vector{Symbol}   # component keys in execution order (owner last)
    vm_index::Int
end
```

How the slots get assigned (`coupling.jl:320-343`): iterate components in declaration order,
and each component's states in local order. For each state:

- **Non-shared** → push a new global slot. The **prefixing rule** (`coupling.jl:328`):
  ```julia
  gname = ck === primary ? sname : Symbol(ck, :_, sname)
  ```
  The **primary** component (first node) keeps bare names (`:d`); every other component's
  names are prefixed with the node name (`:B_x`). Its IC value is pushed into `u0`.
- **Shared** → if the partner endpoint already made a slot, reuse it
  (`slot_of[key] = slot_of[partner]`, `:335`); otherwise make one new slot named `sh.name`
  seeded with the *owner's* initial value (`_share_initial_value`, `:367-370`).

After that: `u0`'s element type is promoted across all components' IC eltypes (`:318`);
`name_to_index` is built with a **collision check** that throws
`"coupled state name collision on :$nm; rename or share to resolve"` (`:345-350`);
`solution_indices` is a `NamedTuple` mapping each component to its global indices in local
order (`:352-357`); and `vm_index` is the global slot of the *primary's* transmembrane
potential (`:362`).

The subtle one is **`operator_order`** (`_operator_order`, `coupling.jl:374-400`): a stable
topological sort where every `share` adds a precedence edge `(non-owner → owner)`, so the
**owner always steps last** among the components sharing a slot. With no shares it's just
declaration order (`:380`); a cycle throws (`:391`). This ordering is what makes hard-discard
work — we'll see why at the next station.

---

## Station 5 — The execution plan & the functor walk

The second product of `couple()` is the **plan**: a pre-resolved, concretely-typed tuple that
the functor walks with zero allocation. Each component gets one entry:

```julia
# coupling.jl:421-428
struct CompEntry{M, B, P, F, OW, AD}
    model::M         # the component (a deepcopy if it's a connect receiver)
    block::B         # its slice of the global vector: a UnitRange if contiguous, else an index vector
    params::P        # its writable_parameters vector (private copy), or nothing if no incoming edges
    frozen::F        # local indices of shared states it does NOT own → zeroed after it runs
    overwrites::OW   # connect edges resolved to (src_global_index, dst_param_index)
    adds::AD         # ditto, for op = +
end
```

`_build_plan`/`_entries` (`coupling.jl:438-461`) build the tuple **by recursion, in operator
order**. Recursion (rather than a comprehension) is deliberate: it keeps every element's
concrete type in the tuple type `plan::PL`, so the functor specializes and stays
allocation-free. `_entry` computes the `block` (a fast `UnitRange` when the component's global
indices happen to be contiguous, else an index vector, `:453`), the `frozen` indices, the
`overwrites`/`adds` with sources already mapped local→global, and grabs
`writable_parameters` only if the component has incoming connect edges (else `nothing`, `:459`).

**The functor** is the payoff. `CoupledModel` is callable:

```julia
# coupling.jl:540
(cm::CoupledModel)(dU, U, p, t) = (_run!(dU, U, p, t, cm.plan); nothing)
```

`_run!` (`coupling.jl:466-474`) is fully unrolled via `Base.tail` recursion over the plan
tuple. For each entry, in order:

1. **Stage connect inputs** — `_connect!` writes this component's wired inputs into its
   private `params` (Station 6).
2. **Evaluate the component** into its own view of the shared vectors:
   `e.model(view(dU, e.block), view(U, e.block), p, t)`. It reads and writes *only its own
   slice*.
3. **Zero its frozen derivatives** — `dU[e.block[i]] = zero(eltype(dU))` for each frozen
   local index.

`frozen` (`_frozen_indices`, `coupling.jl:551-563`) is the set of local state positions a
component participates in via a `share` but does **not** own. Because a shared slot is one
global slot written by both components' views, the non-owner's derivative into it must be
discarded. Zeroing it *after* the non-owner runs — combined with owner-last operator order so
the owner runs afterward and writes the final value — is exactly what makes the owner's
equation govern. That's the whole hard-discard mechanism: **topological order + post-hoc
zeroing**.

---

## Station 6 — `share` vs `connect`, at the code level

This is the conceptual heart of the whole subsystem, so let's be precise about the
difference. They are *not* two flavors of the same thing — they touch different parts of the
system.

**`share` flows through `U` (the global state).** In `_compute_layout`, the two shared states
are assigned the *same global slot* (`slot_of[key] = slot_of[partner]`, `coupling.jl:335`), so
both components' `view`s of `dU`/`U` alias one slot. There is no parameter staging — the
non-owner's derivative is simply zeroed (`frozen`) and the owner's write wins. Because the
shared quantity is a genuine integrated **state**, an implicit solver differentiates through
it correctly.

**`connect` writes into a parameter slot.** The source state is read *live from `U`* each
evaluation and written into the receiver's private `params` vector before the receiver steps.
It never occupies a global state slot. Here's the staging routine:

```julia
# coupling.jl:478-492
_connect!(U, ::Nothing, ::Tuple{}, ::Tuple{}) = nothing   # fast no-op: component has no incoming edges
function _connect!(U, params, overwrites, adds)
    @inbounds begin
        for (_, d) in adds          # 1. zero every +-target slot
            params[d] = zero(eltype(params))
        end
        for (s, d) in overwrites    # 2. apply overwrites (plain copy)
            params[d] = _connect_value(U[s])
        end
        for (s, d) in adds          # 3. accumulate additive fan-in
            params[d] += _connect_value(U[s])
        end
    end
    return nothing
end
```

Two things to note. First, `+` is a **cross-sectional sum reset every eval**, not a running
total over time — the slot is zeroed then summed each call. Second, the ops are partitioned
into homogeneous `overwrites`/`adds` lists at build time by `_connect_plan`
(`coupling.jl:569-581`), so each loop is dispatch-free. `_connect_value` is the identity on
the plain path (`coupling.jl:434`) — its `Dual` override matters under implicit solvers, which
is Station 7.

The one-line summary to carry away: **`share` merges states in `U`; `connect` pipes a state
into a parameter.** A `share` changes what the *state vector* contains; a `connect` changes
what a model *reads as input*.

---

## Station 7 — Solving & observing (the extensions)

The base `src/coupling.jl` is pure Julia with no solver dependency. Two package extensions add
the solve path and AD support.

**`ext/SciMLBaseExt.jl` — the solve entry point.**

```julia
# SciMLBaseExt.jl:6-9
function SciMLBase.ODEProblem(model::CytoZoo.AbstractCellModel, tspan::Tuple;
                              u0 = CytoZoo.default_initial_state(model), p = nothing, kwargs...)
    return SciMLBase.ODEProblem{true}(model, u0, tspan, p; kwargs...)
end
```

This dispatches on *any* `AbstractCellModel`, so a `CoupledModel` hits the very same method —
there is **no** CoupledModel-specific `ODEProblem`. The model itself is passed as the in-place
RHS (`{true}` = isinplace), because every model is a `(model)(du,u,p,t)` functor. For a
coupling, that functor is the plan-walking `_run!` from Station 5. All the coupling machinery
was baked into `cm.plan` at `couple()` time; problem construction adds nothing. Default `u0`
is the merged `copy(cm.layout.u0)`; default `p = nothing` (the non-spatial path).

So the full solve path is just:

```julia
cm  = couple(nodes, edges)
sol = solve(ODEProblem(cm, (0.0, T)), Tsit5())
```

**`monitor_history(sol, model)`** (`SciMLBaseExt.jl:22-34`) surfaces DERIVED observables —
algebraic functions of state like conservation laws `ATPm = C_A − ADPm` — **post-solve**. It
allocates an `nmon × length(sol.t)` matrix (rows = monitors, cols = time), and for each saved
`sol.u[j]` calls `monitor_values!` into column `j`. A zero-monitor model returns an empty
`0×N` matrix without error. `CoupledModel` aggregates its components' monitors
(`coupling.jl:510-535`), walking them in *declaration* order (not operator order) so names and
values stay aligned, with the same `:comp_name` prefixing as states.

Why post-solve? The saved `sol.u` holds plain primal floats — no `Dual`s — so recomputing
monitors afterward never touches solver internals. Which is the exact concern of the last
extension.

**`ext/ForwardDiffExt.jl` — connect under implicit solvers.** One method:

```julia
# ForwardDiffExt.jl:12
CytoZoo._connect_value(x::ForwardDiff.Dual) = CytoZoo._connect_value(ForwardDiff.value(x))
```

Under an implicit solver, ForwardDiff threads `Dual`s through `U` to build the Jacobian. A
`connect` reads `U[s]` (a `Dual`) and tries to write it into the receiver's *`Float64`*
`params` slot — a type error, and conceptually wrong. `_connect_value` extracts the primal
(recursing so nested/higher-order Duals collapse to the underlying eltype), which **freezes the
connect input to its current value within the Newton step**. Consequence: the fixed point is
correct, but the Jacobian is *approximate* (the connect coupling's sensitivity is dropped from
the AD graph) — so a tightly-coupled stiff `connect` may converge poorly. `share` is
unaffected: it flows through `U` and never passes through a parameter slot. The base identity
`_connect_value(x) = x` keeps the `Float64`/GPU path allocation-free and ForwardDiff-free.

---

## Station 8 — Switching is composition

Now the principle from the visitor center, with teeth. **There is no switch primitive, no
construction-time switch kwarg, and nothing routes through `p`.** To turn an optional
capability off, you *compose without it* — omit a whole `Subsystem`, or omit a single edge —
and the rest of the system recovers its baseline. Omitting the element *is* the OFF-invariant.

Concretely, the same mechanism appears at three granularities (these are sockets 6–8 of the
canonical example, Station 9):

- **Module on/off** — include or omit a whole `Subsystem(ToyRedox())`.
- **Edge on/off** — include or omit a `connect`/`share` edge at the `couple()` call.
- **State↔param role flip** — include a `Subsystem(ToyH())` + a `connect` and `h` becomes a
  live integrated state; omit them and the receiver holds `h` as a fixed parameter.

Composition resolves entirely at `couple()`/construction time — never through the
time/space-varying `p`. This is why the architecture needs no new primitive to express "module
X is off": *off is just the graph without X's node or edge.*

**One thing that is deliberately NOT baked in as a permanent invariant:** "a shared state has
exactly one governing equation." Today's `share` is hard-discard (one owner governs), but some
feedback couplings need two models to each contribute a term to the *same* shared derivative.
That's the central gap — see the under-construction wing below.

---

## Station 9 — The worked examples (where to run it)

Two runnable files. Start with the first.

**`examples/coupling_toy.jl`** — the minimal 2-pattern intro (share + connect), ~90 lines.
Defines `ModelA` (owns shared `d`), `ModelB` (its `x` ≡ A's `d`), and a `Reader`. The two
call sites, verbatim:

```julia
# share: A.d ≡ B.x, owner A governs the slot (coupling_toy.jl:48-54)
coupled = couple(
    [Subsystem(ModelA(); name = :A), Subsystem(ModelB(); name = :B)],
    [share(:A => :d, :B => :x; owner = :A)],
)

# connect: Reader integrates A's d, read through a parameter slot (coupling_toy.jl:78-84)
coupled2 = couple(
    [Subsystem(ModelA(); name = :A), Subsystem(Reader(); name = :R)],
    [connect(:A => :d, :R => :d_ext)],
)
```

Run: `julia --project examples/coupling_toy.jl`. It prints `d(2) ≈ exp(-2)` (owner A wins,
B's `du(x)=3` discarded) and `acc(2) ≈ 1-exp(-2)` (the reader integrating A's `d`).

**`examples/coupling_mwe.jl` (+ `.md`)** — the **canonical API driver**. Two tiny models
(~8 states) that between them exercise the *entire* coupling taxonomy from the real
ECCMitoRedox architecture, small enough to hold in your head. The models:

- **`ToyDriver` D** (≈ Gauthier ECC) — states `u,v,a,w`; a monitor `b = Ca − a`.
- **`ToyResponder` R** (≈ Kembro mito) — states `y,m,e`; WIRE receiver slots `p_u,p_v`; a
  held-as-param `h`; a monitor `n = Cm − m`.
- **`ToyRedox`** — optional subsystem (state `z`), the redox module.
- **`ToyH`** — optional subsystem (state `h`) that flips `h` from param to live state.

The socket map (from `coupling_mwe.md`), ✅ = expressible today, ❌ = drives new API:

| # | Socket | Pattern | Toy expression | Status |
|---|--------|---------|----------------|--------|
| 1,2 | `u`,`v` | feedforward WIRE (overwrite): state → param | `connect(:D=>:u, :R=>:p_u)` | ✅ |
| 3 | `b` | DERIVED-source WIRE | *want* `connect(:D=>:b, :R=>:p_b)` | ❌ |
| 4 | `a`/`e` | adopt-native / drop receiver state | `share(:D=>:a, :R=>:e; owner=:D)` | ✅ |
| 5 | `w` | feedback additive contributed-flux into a shared derivative | *want* Redox's `+J` into D's `dw` | ❌ |
| 6 | redox | module on/off = compose w/without a subsystem | include/omit `Subsystem(ToyRedox())` | ✅ |
| 7 | edges | edge on/off = compose w/without an edge | include/omit the WIRE edges | ✅ |
| 8 | `h` | state↔param flip = compose w/without a subsystem | include/omit `Subsystem(ToyH())` + connect | ✅ |
| 9 | `b`,`n` | DERIVED / monitor (conservation law) | `monitor_names` / `monitor_values!` | ✅ |

The `.jl` runs each ✅ socket as a live demonstration with **guards** that fail loudly on
regression: **G1 isolation** (omit the WIRE edges ⇒ R's sub-trajectory reproduces standalone R
*exactly*, `coupling_mwe.jl:214-223`), **G2 one-way** (a pulse into `u` lifts R's `y`, which
peaks early then returns, `:225-241`), and **conservation closures** `a+b=Ca`, `m+n=Cm` checked
via `monitor_history` (`:314-327`). Run:
`julia --project=examples examples/coupling_mwe.jl`.

---

## The road not yet built (the under-construction wing)

Two patterns the taxonomy needs but CytoZoo **cannot express yet** — written as executable
target-API specs at the bottom of `coupling_mwe.jl` (`:329-357`):

1. **Additive contributed-flux into a shared derivative (socket 5) — the central feedback
   gap.** Some feedback couplings need model A's core equation *plus* model B's extra flux to
   drive the *same* shared state's derivative (e.g. `dw = Pw − Lw·w + J`, where `J` comes from
   the redox module). Today's `share` is hard-discard (one owner, non-owner zeroed) and
   `connect` writes a *parameter*, not a *derivative* — so neither can sum a non-owner flux into
   an owner's `dU`. Candidate spellings:
   ```julia
   share(:D => :w, :Redox => :w_flux; owner = :D, op = +)   # (a) additive share
   inject(:Redox => :J, :D => :w)                            # (b) a new flux-injection edge kind
   ```
   Implementation: an accumulate-into-owner's-`dU` path in `_run!` — a *contributory*
   alternative to today's `frozen`-index zeroing.

2. **DERIVED-source `connect` (socket 3).** Let a `connect` source resolve to a monitor/derived
   value (computed each RHS), not only a `state_index`:
   ```julia
   connect(:D => :b, :Redox => :p_b)   # b is a DERIVED monitor (Ca − a), not a state
   ```
   `connect` currently resolves its source only through `state_index` (`_connect_plan`,
   `coupling.jl:569-581`), so a monitor name fails. Needed: fall back to a monitor index and
   compute the derived value each eval.

Design each against the toy first, then apply the finalized API to the real 101-state
ECCMitoRedox model.

---

## Gotchas / constraints

- **Not thread-safe to solve concurrently.** `couple` deepcopies each connect receiver so the
  RHS scratches a *private* parameter vector (`coupling.jl:189-195`) — but a single
  `CoupledModel` still holds one such copy. For an `EnsembleProblem`, deepcopy the `cm` per
  trajectory in the `prob_func`. Full in-RHS reentrancy (routing connect inputs through the
  receiver's per-eval `p`) is a tracked follow-up.
- **`num_parameters(::CoupledModel)` throws by design** (`coupling.jl:543-545`) — a composite
  has no single flat parameter vector. Don't rely on it; use per-component access.
- **Implicit solver + tightly-stiff `connect` may converge poorly** — the connect input is
  frozen within the Newton step (approximate Jacobian; Station 7). `share` is unaffected.
- **No spatial resolution inside coupling.** `src/coupling.jl` is purely state/parameter-slot
  based — there is no `_resolve_spatial` here. `SpatialContext`/`overrides` live in the leaf
  models' RHS, not in the coupling layer.
- **First node is always primary** — it keeps bare state names and supplies `vm_index`. Order
  your `couple` node list intentionally.
- **Name collisions throw at `couple()` time** — two components with the same unprefixed name
  (or a bad share) error with `"coupled state name collision on :$nm; rename or share to
  resolve"` (`coupling.jl:345-350`).

---

## Where the design rationale lives now

The two prior design docs (`handoffs/2026-06-25-1515-coupling-monolithic-rhs.md` and
`handoffs/2026-07-13-1549-cytozoo-switches-are-subsystems.md`) referenced in CLAUDE.md have
been **deleted** — the `handoffs/` dir is otherwise empty. The surviving canonical sources are:

- **CLAUDE.md** → the "Coupling", "Variable roles", and "Derived observables" sections.
- **README** → the "Variable roles" table (the canonical version).
- **`examples/coupling_mwe.md`** → the socket map and the D↔Gauthier / R↔Kembro analogue mapping.
- **This doc** → the code-level tour tying them together.
