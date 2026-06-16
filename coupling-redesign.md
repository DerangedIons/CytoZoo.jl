# CytoZoo: Drop MTK, Add Operator-Splitting Coupling Interface

## Context

CytoZoo is being repositioned from a "registry of cardiac cell models" toward a **pure interface package for cell models AND for coupling them together**. Two motivating decisions, both made by the team:

1. **Move away from ModelingToolkit.jl.** The MTK-backed `BeelerReuter` model and the `MTKCardiacCellModelsExt` extension are removed. The `symbolic_system` / `has_symbolic_system` hooks lose their only consumer and are removed from the interface.

2. **Add a model-coupling interface.** Users need to compose two or more `AbstractCellModel`s into a single combined model. The driving example: model A has states `(b, c, d)`, model B has states `(x, y)`, where physically `d ≡ x`. The user wants to keep the equation from one side (their choice) and discard the other; the result must be solvable by OrdinaryDiffEq with the standard `f(du, u, p, t)` form. Read-only cross-references (e.g., model B reads model A's `Vm`) are also in scope.

Coupling is built on **OrdinaryDiffEqOperatorSplitting.jl** (OS) — operator splitting (Lie–Trotter, Strang) gives well-established convergence theory, lets each operator use its own inner solver, and keeps shared state in a single global vector (so aliased variables ARE the same memory — no callback dance to keep two copies consistent). PartitionedIntegrators.jl was considered as an alternative; OS is the better fit for intra-cell sub-model coupling where shared variables are literally the same physical quantity.

Top priorities, in order: **(1) high performance** (target: within ~1% of uncoupled component performance, no allocations in hot path); **(2) simple authoring API** for someone writing a single component model — they should write today's `f(du, u, p, t)` and not learn anything new for the alias case; **(3) simple coupling API** for the user composing models.

## Amendment (2026-06-16): coupling ships as an extension; OS is a weakdep

**Supersedes the "OS becomes a hard dependency, coupling lives in core" choice below** (see Architecture and the `Project.toml`/`src/coupling.jl` notes in "Files to modify / create"). Decided while implementing, after finding that `OrdinaryDiffEqOperatorSplitting` (OS) is **not in the General registry** and pulls in the full OrdinaryDiffEq solver tree.

A hard `[deps]` on OS would end CytoZoo's "zero runtime dependencies" property — every user, even one integrating a single ToRORd cell, would load the solver stack, and the base package would hard-depend on an unregistered package. Instead:

- **OS is a `[weakdeps]`** behind a new `CouplingExt`, alongside `SciMLBaseExt`/`ThunderboltExt`. The base stays zero-dep; coupling is opt-in via `using CytoZoo, OrdinaryDiffEqOperatorSplitting`.
- **Split point:** the pure-Julia parts — `CoupledModel`, `couple(...)`, layout computation, and the `AbstractCellModel` interface methods on `CoupledModel` — live in **base `src/coupling.jl`** (no OS import; fully testable without solving). The OS-dependent parts — `build_split_function(::CoupledModel)`, `OperatorSplittingProblem(::CoupledModel, tspan)`, and synchronizer wiring — live in **`ext/CouplingExt.jl`**.
- `couple` and `CoupledModel` are exported from base; the solving entry points become available once OS is loaded.

Everything else below (layout rules, alias-via-operator-order semantics, synchronizer-based cross-refs, the verification plan) stands unchanged.

## Architecture

```
┌──────────────────────────────────────────────────────────────────┐
│ CytoZoo  (deps: OrdinaryDiffEqOperatorSplitting + SciMLBase)     │
│                                                                  │
│   src/interface.jl   AbstractCellModel + interface methods       │
│   src/spatial.jl     SpatialContext (existing, unchanged)        │
│   src/coupling.jl    CoupledModel + couple(...)                  │
│                      build_split_function(::CoupledModel)        │
│                      OperatorSplittingProblem(::CoupledModel)    │
│                                                                  │
│   ext/               Thunderbolt, TWorld, SciMLBase (existing)   │
│                      *no MTKCardiacCellModelsExt*                │
└──────────────────────────────────────────────────────────────────┘
```

`CoupledModel <: AbstractCellModel` — a `CoupledModel` is itself a cell model, so couplings nest. Component models stay as plain `AbstractCellModel`s. **No authoring changes for the alias case.** Cross-refs require the receiving model to expose a parameter slot — the only place authoring is touched.

OS becomes a hard dependency, so coupling lives in core (no extension activation).

## Component model contract

A model that participates in coupling must provide (mostly already required):

- Functor `f(du, u, p, t) -> Nothing` with positional indexing into `u` and `du`.
- `num_states(model) -> Int`
- `default_initial_state(model) -> AbstractVector`
- `state_names(model) -> NTuple{N, Symbol}` (already optional, becomes required for participants in coupling)
- `state_index(model, name::Symbol) -> Int` (already optional, becomes required for participants in coupling)
- For cross-references: a parameter object (typically `model.params` or accessed via `parameter_index`) must contain a slot the receiver's `f` reads from. The slot's name is what the user references in the `refs=` declaration.

CytoZoo's interface methods get a small docstring update to indicate the new "required for coupling" status of `state_names` / `state_index`.

## CoupledModel data structure

```julia
struct CoupledModel{Cs<:NamedTuple, AS, RS, L} <: AbstractCardiacCellModel
    components::Cs        # NamedTuple, e.g., (A=model_A, B=model_B)
    aliases::AS           # Tuple of alias specs, see below
    refs::RS              # Tuple of cross-ref specs, see below
    layout::L             # Pre-computed global state layout (NamedTuple of NTuples)
end
```

Built by:

```julia
couple(components::NamedTuple;
       aliases = (),     # Tuple of alias specs
       refs    = ())     # Tuple of cross-ref specs
       -> CoupledModel
```

Spec types (concrete shapes are implementation details — pick whichever feels best at write time, candidates):

- Alias spec: `(comp_a => :state_a, comp_b => :state_b; owner = :comp_a)` — the symbol is the *component name* (NamedTuple key), not the model itself.
- Ref spec: `(source_comp => :state_or_param) => (target_comp => :param_slot)`.

The `couple(...)` constructor:
1. Validates names: each alias references states that exist on the named components; each ref's target slot exists on the target's parameters.
2. Computes the **global layout** — see next section.
3. Stores everything in the resulting `CoupledModel`.

## Global state layout & solution_indices

For example: A states `(b, c, d)`, B states `(x, y)`, alias `A.d ≡ B.x` with `owner=:A`.

| Global idx | Variable | Source | A's local idx | B's local idx |
|-----------:|----------|--------|---------------|---------------|
| 1 | b | A | 1 | — |
| 2 | c | A | 2 | — |
| 3 | d (≡x) | A (owner) | 3 | 1 |
| 4 | y | B | — | 2 |

**Layout rule**:
- Owner's full local-state vector is laid out first, in the owner's own order, contiguously starting at global index 1.
- For each non-owner component, its non-aliased states are appended in the component's local order; its aliased states do NOT take new global slots — they map to the existing global slot owned by the alias's owner.
- Each component's `solution_indices` is a tuple of global indices arranged so that, when `u[solution_indices_k]` is built by the OS solver, the result has the component's expected local order.

For the example:
- A's `solution_indices = (1, 2, 3)` — A's view `u[[1,2,3]] = (b, c, d)`. A's `f` is unchanged.
- B's `solution_indices = (3, 4)` — B's view `u[[3,4]] = (x, y)`. B's `f` is unchanged.

`state_names(coupled)` returns the global tuple, e.g. `(:b, :c, :d, :B_y)`. Naming rules:
- Owner's states keep their original names verbatim.
- Aliased states take a single canonical name (defaults to the owner's name; can be overridden in the alias spec).
- Non-owner components' non-aliased states are **always** prefixed by the component name (e.g., `:B_y`) — predictable and avoids surprising the user when a third component is added later that would have triggered a collision.
- Any residual collision after the above rules (e.g., owner has a state named `:B_y`) errors at `couple(...)` construction time.

`default_initial_state(coupled)` produces a vector where aliased slots get the owner's initial value and non-owner non-aliased slots come from the respective component's `default_initial_state`.

`state_index(coupled, name)` does a Dict lookup over the canonical names.

When CoupledModels nest: layout is computed flat at each `couple(...)` call; nested `couple(...)` calls produce a fresh flat layout from the nested components' (already flat) layouts. No recursion at solve time.

## Aliasing mechanic — how the owner is enforced

`GenericSplitFunction` applies operators in the order they appear in the operator tuple. For Lie–Trotter–Godunov over `[t, t+dt]`:

1. The non-owner B steps first: writes B's value of `x` into global slot 3.
2. The owner A steps next: overwrites global slot 3 with A's value of `d`.

A's write wins → A "owns" the equation. This is exactly the operator-splitting interpretation; no separate "skip B's equation" logic is required — both equations execute, and the splitting solver determines whose value is final at the timestep boundary. Both equations doing work is standard splitting cost; with reasonable `dt` this is mathematically correct (LT-G is 1st-order, Strang is 2nd-order).

`couple(...; aliases=[(A => :d, B => :x; owner=:A)])` orders the operators tuple as `(B_op, A_op)`. With Strang (`(B/2, A, B/2)`-style symmetric), owner semantics are preserved with 2nd-order accuracy. The CytoZoo coupler emits the ordered tuple based on the user's `owner` declaration.

If multiple aliases conflict on owner ordering (e.g., alias 1 wants A last, alias 2 wants B last), `couple(...)` errors at construction time with a clear message; users must restructure (e.g., split into nested CoupledModels) or pick consistent owners.

## Read-only cross-references via synchronizers

OS's `synchronizers` mechanism runs a callback before each operator step to update parameters from the current state vector. CytoZoo uses this for read-only cross-refs.

For `refs = [(A => :Vm) => (B => :Vm)]` (B reads A's Vm):

1. `couple(...)` validates that B's parameters have a `Vm` slot (via `parameter_index(B, :Vm)` or by introspecting `B.params`).
2. CytoZoo generates a synchronizer: before B's operator step, copy `u[global_idx_for_Vm]` into the appropriate slot of B's parameters.
3. B's `f` reads from `p.Vm` (or whatever the parameter slot is) — natural style for ion-channel modules already written this way.

The synchronizer is a small isbits closure-free struct (so it stays GPU-friendly when we eventually add GPU support). It captures: the global index of the source state, the parameter slot to write into, and the target component's parameter object reference.

## User-facing API surface

```julia
using CytoZoo
using OrdinaryDiffEq            # via OS reexport
using OrdinaryDiffEqOperatorSplitting  # hard dep, but importable

A = ToRORd()                    # AbstractCellModel
B = SomeChannelModel()          # AbstractCellModel — has B.params.Vm slot

coupled = couple(
    (A = A, B = B);
    aliases = [(A => :d, B => :x; owner = :A)],
    refs    = [(A => :Vm) => (B => :Vm)],
)

state_names(coupled)            # (:b, :c, :d, :B_y)
state_index(coupled, :d)        # 3
default_initial_state(coupled)  # owner's IC for aliased slots; concat otherwise
num_states(coupled)             # 4

prob = OperatorSplittingProblem(coupled, (0.0, 1000.0))
sol  = solve(prob, LieTrotterGodunov((Tsit5(), Tsit5())); dt = 0.05)
sol[:d]                         # trajectory of the aliased variable
```

For the simple "no alias" case (just running multiple operators sharing a global vector), the API still works — `aliases = ()` is the default and the layout becomes a straightforward concatenation.

## MTK removal scope

- Delete `ext/MTKCardiacCellModelsExt/` directory and `ext/MTKCardiacCellModelsExt.jl`.
- Remove `MTKCardiacCellModels`, `ModelingToolkit` from `[weakdeps]` and `[extensions]` in `Project.toml`.
- Remove `function BeelerReuter end` stub and `export BeelerReuter` from `src/CytoZoo.jl`.
- Remove `symbolic_system`, `has_symbolic_system` from `src/interface.jl` and from any exports.
- Delete `test/test_mtk_ext.jl`; remove inclusion from `test/runtests.jl`.
- Update `CLAUDE.md` to drop MTK references and document the new coupling layer.
- Update memory file `project_cytozoo.md` to mark MTK as removed (post-implementation).

## Files to modify / create

**Modify:**
- `Project.toml` — drop MTK weakdeps + extension; add OS as direct `[deps]` with appropriate `[compat]`; add SciMLBase to `[deps]` if not already (currently weakdep).
- `src/CytoZoo.jl` — drop `BeelerReuter` stub/export; add `include("coupling.jl")`; export `CoupledModel`, `couple`, `build_split_function`.
- `src/interface.jl` — remove `symbolic_system` / `has_symbolic_system`; small docstring updates noting `state_names` / `state_index` are required for coupling participants.
- `test/runtests.jl` — drop MTK include; add `test_coupling.jl` include.
- `CLAUDE.md` — update architecture description to drop MTK section, add coupling section.

**Create:**
- `src/coupling.jl` — `CoupledModel`, `couple(...)`, layout logic, `build_split_function`, `OperatorSplittingProblem(::CoupledModel, tspan; ...)`, synchronizer struct + generation.
- `test/test_coupling.jl` — alias correctness, cross-ref correctness, owner semantics, convergence under LT-G refinement.
- `examples/coupling_toy.jl` — the b/c/d + x/y minimal example.
- `examples/coupling_cardiac.jl` — one realistic cardiac case (ToRORd + small signaling sub-model, when an appropriate sub-model is available).

**Delete:**
- `ext/MTKCardiacCellModelsExt.jl`, `ext/MTKCardiacCellModelsExt/`
- `test/test_mtk_ext.jl`

## Verification

End-to-end checks (in `test/test_coupling.jl` and via REPL):

1. **Trivial single-component coupling**: `couple((A=ToRORd(),))` produces a `CoupledModel` whose RHS, run via `OperatorSplittingProblem`, matches the uncoupled `ToRORd` trajectory at the same `dt` to within solver tolerance.

2. **Alias correctness with hand-coded reference**: build the b/c/d + x/y toy with a known analytic combined ODE; confirm `couple(...; aliases=[(A=>:d, B=>:x; owner=:A)])` reproduces the analytic solution to a refinement-correct error as `dt → 0`.

3. **Owner semantics**: same toy, swap `owner=:A` vs `owner=:B`; verify operator order is reversed and the aliased trajectory matches the chosen owner's equation.

4. **Cross-ref correctness**: a synthetic case where B's RHS depends linearly on A.Vm; confirm the synchronizer delivers the correct value (compare against a hand-rolled monolithic version).

5. **Convergence**: LT-G with `dt = 1, 0.1, 0.01` on the toy; verify ~1st-order error decay. Strang once available: 2nd-order.

6. **Performance**: `@btime` the coupled RHS for ToRORd-only case vs uncoupled ToRORd RHS; confirm <1% overhead and zero allocations.

7. **MTK fully gone**: `using CytoZoo` doesn't reference MTK; `Pkg.test()` passes without MTK installed.

REPL check:
```bash
julia --project=. -e '
using CytoZoo, OrdinaryDiffEq, OrdinaryDiffEqOperatorSplitting
A = ToRORd()
coupled = couple((A=A,))
prob = OperatorSplittingProblem(coupled, (0.0, 100.0))
sol = solve(prob, LieTrotterGodunov((Tsit5(),)); dt=0.1)
println(sol[end][CytoZoo.transmembrane_potential_index(coupled)])
'
```

## Implementation order

1. **Drop MTK** (clean baseline; orthogonal). Verify tests pass without MTK.
2. **Add OS as a hard dep** in `Project.toml` (and `SciMLBase` if needed). Confirm `using CytoZoo` still works and existing tests pass.
3. **Add `CoupledModel`** to `src/coupling.jl`: data structure, `couple(...)` constructor, layout computation, and the `AbstractCellModel` interface methods (`num_states`, `state_names`, `state_index`, `default_initial_state`). This is pure Julia; testable without OS solving.
4. **Add `build_split_function(::CoupledModel)`** and `OperatorSplittingProblem(::CoupledModel, tspan; ...)`. Wraps each component's `f` so it gets a slice-aware view in operator form. Generates synchronizers from the `refs` declaration.
5. **Tests** in `test/test_coupling.jl` covering layout, alias correctness, owner semantics, cross-refs, convergence, and zero-allocation in the hot path.
6. **Examples**: the toy in `examples/coupling_toy.jl`; one realistic cardiac case once a candidate sub-model is identified.
7. **Documentation**: update `CLAUDE.md`; add a coupling section to README if/when the package gets a README pass.

## Open questions / future work

- **Strang and higher-order splitting**: LT-G ships first; Strang requires emitting `(B/2, A, B/2)`-style ordered tuples. Add once the LT-G path is solid.
- **GPU compatibility**: synchronizer structs are kept isbits-friendly so a future GPU pass needs only to wire OS's GPU support (when it lands). Current design does not preclude GPU but doesn't ship with it.
- **Hierarchical/nested coupling**: `couple(c1, couple(c2, c3))` — the design supports it (CoupledModel is AbstractCellModel) but the layout collapsing logic should be reviewed once we hit a real use case.
- **Alias ownership conflict resolution**: if multiple aliases imply contradictory operator orderings, error at `couple(...)` time. Future work: a dedicated topological resolver.
- **Non-ODE components**: the user noted "almost always ODEs (but not 100%)". DAE / algebraic-constraint coupling is out of scope for this design; revisit if a concrete need emerges.
