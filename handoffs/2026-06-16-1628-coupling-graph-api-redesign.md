---
slug: coupling-graph-api-redesign
created: 2026-06-16-1628
status: done
---

> **Resolved 2026-06-17 — implemented with refinements from a follow-up brainstorm.** The graph
> front-end shipped, but Kyle's feedback changed several names and added an edge-operation concept
> vs. the draft below. Final decisions (authoritative):
> - **Node type is `Subsystem`** (not `Domain`); constructor `Subsystem(model, alg = nothing; name)`.
> - **`alias → share`**, **`crossref → connect`** (both kept; renamed). Internal: `AliasSpec→ShareSpec`,
>   `CrossRef→ConnectSpec`, `CrossRefSync→ConnectSync`, `CoupledModel.aliases→shares`/`refs→connects`.
> - **`connect` carries an `op`** (`op(current_slot, source_value)`), default `overwrite` (copy);
>   non-default `op` errors for now (one-line generalization point in `forward_sync_external!`).
> - `coupled_algorithm(cm; scheme)` reads solvers off the nodes — `inner` arg removed; `_inner_solver→_node_alg`.
> - Exports: `couple, CoupledModel, Subsystem, share, connect, build_split_function, coupled_algorithm`.
>
> See `coupling-redesign.md` (2026-06-17 amendment) and the approved plan
> `~/.claude/plans/read-handoffs-2026-06-16-1628-coupling-g-fluttering-acorn.md`. The rest of this
> doc is the original draft; read names through the rename map above.

# Handoff: Reshape CytoZoo coupling API to PartitionedIntegrators-style graph construction

## Goal / why this matters
Make CytoZoo's model-coupling **front-end** feel like Kyle's own `PartitionedIntegrators.jl`
(`~/pro/dev/PartitionedIntegrators`): explicit `Domain` *nodes* (each a model + its solver) and a
list of *edges*, assembled positionally (`couple([n1, n2], [edge])`) instead of today's keyword-bag
`couple((A = ModelA(), B = ModelB()); aliases = [...])`. **The engine underneath does not change** —
this is purely an API-shape change. Kyle was explicit: *"i really just liked the way PI felt when
constructing the model / its API"* and *"use as much OS as possible, don't want too much custom
stuff."*

Work happens on branch `feat/coupling-interface`. The coupling feature is brand-new and unreleased,
so there is **no back-compat burden** — replace the old constructor outright.

## Background & current state
Current coupling lives in `src/coupling.jl` (pure Julia, no OS dep) + `ext/CouplingExt.jl`
(OrdinaryDiffEqOperatorSplitting solving). Two mechanisms, both kept as-is:
- **alias** `alias(:A => :d, :B => :x; owner = :A)` — two states share one global slot; owner's
  equation governs (hard-discard; non-owner derivative zeroed).
- **cross-ref** `crossref(:A => :Vm, :B => :Vm_ext)` — OS synchronizer copies a source state into a
  receiver's parameter slot before the receiver steps.

`CoupledModel <: AbstractCardiacCellModel` (couplings nest; usable as a Thunderbolt ionic model).
The full approved plan is at `/Users/kylebeggs/.claude/plans/i-want-the-api-keen-parrot.md` — read it
for the complete file-by-file detail. This handoff is the standalone brief.

## Key files / locations
- `src/coupling.jl` — `AliasSpec`/`alias`, `CrossRef`/`crossref`, `CouplingLayout`, `CoupledModel`,
  `couple`, `_compute_layout`, `_operator_order`, `_validate_specs`, interface methods. **Primary edits here.**
- `ext/CouplingExt.jl` — `ComponentOperator`, `build_split_function`, `_frozen_indices`,
  `OperatorSplittingProblem`, `coupled_algorithm`, `CrossRefSync` + sync hooks. **Only `coupled_algorithm` changes.**
- `src/CytoZoo.jl` — coupling `export` line (add `Domain`).
- `examples/coupling_toy.jl` — runnable demo to rewrite.
- `test/test_coupling.jl` (pure-Julia layout) and `test/test_coupling_ext.jl` (OS solving).
- `CLAUDE.md` (coupling section), `README.md` (coupling section), `coupling-redesign.md` (design log).
- Reference: `~/pro/dev/PartitionedIntegrators/src/domain.jl`, `partitioned_problem.jl`, `coupling.jl`.

## Decisions & conclusions (do not relitigate)
1. **Engine unchanged.** Keep shared-state OS backend + `CoupledModel <: AbstractCardiacCellModel`.
2. **Replace** `couple(::NamedTuple; aliases, refs)` with `couple(nodes, edges)`. No deprecation.
3. **Solver lives on the node** (PI style). `coupled_algorithm` reads solvers off the nodes; its
   `inner` argument is **removed**.
4. **Minimal custom machinery.** *No* custom `solve` wrapper, *no* custom solution type, *no*
   per-domain `sol[:A]` indexing, *no* manual step loop. Solve path stays pure OS
   (`OperatorSplittingProblem` + `coupled_algorithm` + `init`/`solve!`, read results from `integ.u`).
5. **`Domain` has no `kwargs` field.** OS 0.3 has no per-operator solve-option path, so PI-style
   `Domain(...; reltol=...)` kwargs would be inert — drop them.
6. **`alias`/`crossref` constructors stay unchanged** — they are already directed edge specs; do NOT
   adopt PI's generic `Coupling(from, to, affect!)` edge type (different semantics).

### Engine facts verified against OS v0.3.1 (`~/.julia/packages/OrdinaryDiffEqOperatorSplitting/3KKv4/`)
- `scheme(algs::Tuple)` with `scheme ∈ {LieTrotterGodunov, StrangMarchuk}` is the algorithm
  constructor (`src/solver.jl:10-12`, `:75-77`). `coupled_algorithm` keeps building `scheme(algs)`;
  only the *source* of `algs` changes.
- The OS solution object is **built empty and never populated** (no `savevalues!`/`push!` in the
  step loop; `fix_solution_buffer_sizes!` is dead code). Results are read off live `integ.u` — which
  is exactly what today's examples/tests already do. This is why there's no `sol[:A]` trajectory.

## What's left / next steps (suggested order)
1. **`src/coupling.jl`**: add the `Domain` node + constructor; change `CoupledModel` fields; rewrite
   `couple(nodes, edges)` + helpers; update docstrings.
   ```julia
   struct Domain{M, A}
       name::Symbol
       model::M     # an AbstractCellModel (a functor) — NB: PI's Domain.prob is an ODEProblem
       alg::A       # inner solver, e.g. Tsit5(); `nothing` allowed (base/layout work needs no solver)
   end
   Domain(model, alg = nothing; name::Symbol = gensym(:domain)) = Domain(name, model, alg)

   struct CoupledModel{Cs<:NamedTuple, AL<:NamedTuple, AS<:Tuple, RS<:Tuple, L<:CouplingLayout} <: AbstractCardiacCellModel
       components::Cs   # name => model   (KEEP — ext reads cm.components[ck])
       algs::AL         # name => inner solver (or nothing)  (NEW)
       aliases::AS
       refs::RS
       layout::L
   end

   function couple(nodes, edges = ())
       isempty(nodes) && throw(ArgumentError("couple requires at least one Domain node"))
       domains    = _domains_namedtuple(nodes)            # NamedTuple keyed by node.name; errors on dup names
       components = map(d -> d.model, domains)
       algs       = map(d -> d.alg,   domains)
       alias_tuple, ref_tuple = _split_edges(edges)        # partition by type: AliasSpec vs CrossRef
       _validate_specs(components, alias_tuple, ref_tuple) # REUSED unchanged
       layout = _compute_layout(components, alias_tuple)   # REUSED unchanged
       return CoupledModel(components, algs, alias_tuple, ref_tuple, layout)
   end
   ```
   - `_domains_namedtuple` pre-checks duplicate names (clear error before NamedTuple build).
   - `_split_edges` partitions `Tuple(edges)` into `(aliases, refs)` by `isa AliasSpec`/`isa CrossRef`;
     errors naming the offending type on anything else.
   - Node order is significant: first node = primary component (bare state names + `vm_index`).
   - **Reuse unchanged:** `alias`, `crossref`, `_compute_layout`, `_operator_order`, `_validate_specs`,
     `_check_component`, `_check_state`, `_coupled_initial_eltype`, `_alias_initial_value`, and all
     interface methods.
2. **`ext/CouplingExt.jl`**: rewrite `coupled_algorithm`, delete `_inner_solver`, add `_node_alg`.
   ```julia
   function CytoZoo.coupled_algorithm(cm::CoupledModel; scheme = OS.LieTrotterGodunov)
       algs = Tuple(_node_alg(cm, ck) for ck in cm.layout.operator_order)
       return scheme(algs)
   end
   function _node_alg(cm::CoupledModel, ck::Symbol)
       alg = cm.algs[ck]
       alg === nothing && throw(ArgumentError(
           "node :$ck has no inner solver; pass one to Domain(model, alg), e.g. Tsit5()"))
       return alg
   end
   ```
   Everything else in the ext is untouched.
3. **`src/CytoZoo.jl`**: add `Domain` to the coupling export line.
4. **`examples/coupling_toy.jl`**: rewrite both `couple(...)` calls + `coupled_algorithm(...)` (sketch below).
5. **Tests**: `test/test_coupling.jl` — rewrite each `couple(...)` to graph form (alg can be `nothing`);
   keep layout assertions; add tests for duplicate-name error, unknown-edge-type error, single-node
   graph, `gensym` name defaulting. `test/test_coupling_ext.jl` — graph-form construction with solvers
   on nodes, `coupled_algorithm(cm)` (no `inner`); add a nil-`alg` → `ArgumentError` test; zero-alloc
   test unchanged.
6. **Docs**: `CLAUDE.md` coupling section, `README.md` coupling block, `coupling-redesign.md` dated
   amendment (record graph front-end + the OS limitations that scoped it).

## API sketch (new toy — mirrors `examples/coupling_toy.jl`)
```julia
using CytoZoo
using OrdinaryDiffEqOperatorSplitting
using OrdinaryDiffEq

# (ModelA / ModelB / Reader struct + interface defs are identical to the current toy)

# 1. Alias: A.d ≡ B.x, owner = A  — nodes carry model + solver; edges are alias/crossref
coupled = couple(
    [Domain(ModelA(), Tsit5(); name = :A),
     Domain(ModelB(), Tsit5(); name = :B)],
    [alias(:A => :d, :B => :x; owner = :A)],
)

integ = init(
    OperatorSplittingProblem(coupled, (0.0, 2.0)),
    coupled_algorithm(coupled);          # solvers read off the nodes — no `inner` arg
    dt = 0.01, adaptive = false,
)
solve!(integ)
integ.u[state_index(coupled, :d)]        # owner A governs the shared slot (≈ exp(-2))

# 2. Cross-ref: a reader integrates A's d through a parameter slot
coupled2 = couple(
    [Domain(ModelA(), Tsit5(); name = :A),
     Domain(Reader(), Tsit5(); name = :R)],
    [crossref(:A => :d, :R => :d_ext)],
)
integ2 = init(OperatorSplittingProblem(coupled2, (0.0, 2.0)), coupled_algorithm(coupled2); dt = 0.01, adaptive = false)
solve!(integ2)
integ2.u[state_index(coupled2, :R_acc)]  # ≈ 1 - exp(-2)
```
Before → after at a glance:
- `couple((A = ModelA(), B = ModelB()); aliases = [...])`
  → `couple([Domain(ModelA(), Tsit5(); name=:A), Domain(ModelB(), Tsit5(); name=:B)], [...])`
- `coupled_algorithm(coupled, Tsit5())` → `coupled_algorithm(coupled)`

## Gotchas / constraints
- **Keep `CoupledModel.components`** (don't replace it with on-the-fly derivation) — the whole
  extension reads `cm.components[ck]`; keeping it means the ext machinery stays untouched.
- **Node order matters** — first node is the primary component (bare names + `vm_index`). Preserve
  iteration order when building the NamedTuple.
- **`gensym` default names** are non-deterministic; multi-node graphs that reference nodes in edges
  **must** pass explicit `name=`. The `gensym` default is only ergonomic for single-node graphs —
  document this.
- **`Domain.alg === nothing`** is fine for pure layout/base tests, but `coupled_algorithm` must throw
  a clear error if a solver is missing at solve time.
- **No `sol[:A]` / no custom solve** — results are read from `integ.u` after `solve!`, sliced by
  `cm.layout.solution_indices[name]` if a per-domain slice is needed (final-state only; OS saves no
  trajectory).
- No secrets involved in this task.

## Verification
1. `julia --project=. examples/coupling_toy.jl` — alias + cross-ref demos print the same numbers as
   today (`d(2) ≈ exp(-2)`; `acc(2) ≈ 1 - exp(-2)`).
2. `julia --project=. -e "using Pkg; Pkg.test()"` — coupling base + OS-ext tests pass, plus new
   edge-case tests.
