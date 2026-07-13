---
slug: coupling-inject-derived-connect
created: 2026-07-12-1542
status: open
---

# Handoff: Implement `inject` + DERIVED-source `connect` so the full coupling MWE runs live

## Goal / why this matters

`examples/coupling_mwe.jl` is the canonical toy that drives CytoZoo's coupling API — two tiny models
(`ToyDriver` D, `ToyResponder` R, ~8 states) covering the full coupling taxonomy from the
ECCMitoRedox architecture. Most sockets run **live and green today**, but **three patterns exist only
as commented "TARGET API" specs** because CytoZoo can't express them. Close these three so the MWE has
**no commented sockets** — the whole taxonomy runs live. This simultaneously unblocks ECCMitoRedox's
feedback stages (Stage 2+). The design is **decided** (see Decisions); this brief is the
implementation plan. No code has been written yet.

## Background & current state

The MWE was built and verified (`julia --project=examples examples/coupling_mwe.jl` — all live
assertions pass). Socket map (full table in `examples/coupling_mwe.md`):

- **Live today (✅):** feedforward WIRE overwrite (`connect(:D=>:u,:R=>:p_u)`, `:v→:p_v`); edge-gating
  switch (`cyto_ions_dynamic` = include/omit edges); hard-discard `share` / adopt-native
  (`share(:D=>:a,:R=>:e; owner=:D)`); DERIVED monitors + conservation closures; intra-model
  construction switches (`ToyResponder(; redox_on, dynamic_h)`) incl. the state-held-as-param inverse
  case.
- **TARGET API only (❌ — the gap):** three patterns, written as executable specs in the `TARGET API`
  block at the bottom of `examples/coupling_mwe.jl` (`:238-272`):
  1. **Additive feedback flux** (socket `w`): R contributes `J = 0.5·z` into D's `w` equation, so
     `dw = Pw − Lw·w + J`. Today `share` is hard-discard (owner wins, non-owner's derivative zeroed);
     nothing sums a contribution.
  2. **DERIVED-source `connect`** (socket `b`): wire D's derived `b = Ca − a` into R's `p_b` param
     slot. Today `connect` resolves its source only via `state_index`, so a derived/monitor name fails.
  3. **Module-switch protocol** (sockets `redox_on`, `dynamic_h`): already work as `ToyResponder`
     constructor kwargs; what's missing is a *blessed convention*.

The coupling core is a **monolithic single-RHS** design: `couple()` compiles the graph once into a
concretely-typed tuple of `CompEntry`; the functor `(cm)(dU,U,p,t)` walks it via `_run!`. It is
**allocation-free and that is a hard invariant** (`@allocated cm(dU,U,nothing,0.0) == 0` is asserted
for share/connect/scatter in `test/test_coupling.jl`).

## Decisions & conclusions (all confirmed with Kyle this session)

Supersedes the earlier open question of whether additive contribution should be spelled
`share(...; op=+)` (a phantom state) vs. a new edge kind — the new edge kind won (item 1 below).

**The unifying insight:** features 1 and 2 both need the *same* new capability — compute a named
derived scalar from a source component's local state each RHS, **allocation-free and Dual-safe**.
The existing `monitor_values!` is buffer-based (fills a vector) → allocates in-RHS and fights
ForwardDiff. The primitive is a **homogeneous tuple accessor** `monitor_values(model,u,t)->NTuple{N}`:
indexing a homogeneous tuple by a runtime `Int` is type-stable, needs no buffer, carries `Dual`s.
Both features read `monitor_values(src_model, view(U, src_block), t)[mi]`.

1. **Additive flux → new `inject` edge kind** (chosen over the `share(...; op=+)` "phantom state"
   spelling). `inject(:R=>:J, :D=>:w)`: R exposes `J` as a monitor; the value is **added** into D's
   `w` derivative after D writes its base. Rationale: J is genuinely algebraic (not an integrated
   state), so no phantom slot is needed; the state vector/Jacobian aren't inflated; and it keeps
   `Dual`s so the coupling is **fully implicit** (better Jacobian than `connect`).
2. **DERIVED sources ARE monitors** (chosen over a separate `coupled_output` hook). Reuses one
   concept for all deriveds; matches the existing "monitors = the DERIVED role" framing. Tradeoff
   accepted: an inject flux like `:J` becomes a monitor, so it also appears in `monitor_history`
   (socket-9 monitor-name assertions get updated — see Gotchas).
3. **Module-switch protocol → documentation only** (no new interface). YAGNI: sockets 6/8 already
   work as constructor kwargs; no consumer needs a declarative hook yet. Just bless the convention:
   construction-time, resolve at construction, **never** through `p`/`SpatialContext`; OFF recovers
   baseline.

ForwardDiff behaviour falls out correctly:
- **connect** (feature 2) stages into a `Float64` param slot → wrapped in `_connect_value`, so under
  an implicit solver the primal is frozen (approximate Jacobian — consistent with today's connect).
- **inject** (feature 1) adds into `dU` (carries `Dual`s) → **no** `_connect_value`, full Jacobian.

## Key files / locations

- `src/coupling.jl` — the machinery to extend. Load-bearing sites (verified this session):
  - `_run!` (`:465-474`) — the plan-walker. Per entry: `_connect!` stages inputs → `e.model(...)`
    writes its `dU` view → `frozen` loop zeroes non-owned shared slots. **inject's `+=` step slots
    in right after the model write.**
  - `CompEntry{M,B,P,F,OW,AD}` (`:421-428`) — per-component plan entry. Gains 3 new partitions.
  - `_entry` (`:448-461`) / `_connect_plan` (`:569-581`) — plan build; `_connect_plan:575` is the
    state-only src resolution (`state_index(components[cn.src], cn.src_state)`) to branch for monitor
    sources. Block logic at `:453` to factor out.
  - `_connect!` (`:478-492`) — stages connect inputs; widens to compute monitor-source values.
  - `_frozen_indices` (`:551-563`) — the frozen-zeroing plan helper.
  - `share`/`ShareSpec` (`:35-60`), `connect`/`ConnectSpec` (`:78-114`), `_split_edges` (`:212-229`),
    `_validate_specs`/`_check_state` (`:233-300`), `_operator_order` (`:374-400`), `couple` (`:172`),
    `CoupledModel` struct (`:155-161`), `_copy_connect_receivers` (`:189-195`).
  - CoupledModel monitor aggregation `monitor_values!` (`:524-536`) — delegates to per-component
    `monitor_values!`; unchanged (works via the new derived default).
- `src/interface.jl` — monitor hooks `num_monitors`/`monitor_names`/`monitor_values!` (`:153-175`);
  `writable_parameters` (`:109-121`). `state_index`/`parameter_index` are bare `function … end`.
- `examples/coupling_mwe.jl` — `ToyDriver` (`:33-75`), `ToyResponder` (`:82-128`),
  `responder_flux_J(m,u_R)=0.5*u_R[2]` (`:130-132`, the `J = gJ·z` flux R injects, defined and
  currently unused-by-design), live sockets (`:138-236`), TARGET API block (`:238-272`).
  `examples/coupling_mwe.md` — socket table + API backlog.
- `test/test_coupling.jl` — share/connect/zero-alloc test patterns; mocks `_MonoA`, `_MonoReader`,
  `_MonoP`/`_MonoQ`, `_ScatterA`/`_ScatterB`. `test/test_monitors.jl` — models to switch to the new
  `monitor_values` accessor (READ IT FIRST). `ext/ForwardDiffExt.jl` (`_connect_value(::Dual)`),
  `ext/SciMLBaseExt.jl` (`monitor_history` calls `monitor_values!` post-solve).

## What's left / next steps

### Step 1 — `src/interface.jl`: tuple monitor accessor
- Add primitive `monitor_values(::AbstractCellModel, u, t) = ()` (in-RHS, allocation-free, Dual-safe).
- Make `monitor_values!` a **derived default** filling the buffer from the tuple:
  ```julia
  function monitor_values!(mon, u, t, model::AbstractCellModel)
      vals = monitor_values(model, u, t)
      @inbounds for i in eachindex(vals); mon[i] = vals[i]; end
      return nothing
  end
  ```
- `num_monitors`/`monitor_names` unchanged. Update docstrings (name `monitor_values` the in-RHS
  primitive, `monitor_values!` the post-solve buffer variant). Models now implement `monitor_values`.

### Step 2 — `src/coupling.jl`: `inject` edge + monitor-sourced connect
- **New edge:**
  ```julia
  struct InjectSpec; src::Symbol; src_flux::Symbol; dst::Symbol; dst_state::Symbol; end
  inject(s::Pair{Symbol,Symbol}, d::Pair{Symbol,Symbol}) = InjectSpec(s.first,s.second,d.first,d.second)
  ```
  `_split_edges` partitions a third `injects::Tuple` (update the unknown-edge error). Thread
  `injects` through `couple`, `_validate_specs`, `_build_plan`, and the `CoupledModel` struct
  (new `injects::IJ` field + type param — store for introspection like `shares`/`connects`).
- **`CompEntry` grows 3 concretely-typed partitions:**
  ```julia
  struct CompEntry{M,B,P,F,OW,AD,MOW,MAD,IJ}
      model; block; params; frozen; overwrites; adds
      mon_overwrites::MOW   # (src_model, src_block, mi::Int, dst_param::Int) — receiver's entry
      mon_adds::MAD         #   same, op=+ monitor-sourced connects
      injects::IJ           # (src_model, src_block, mi::Int, owner_global::Int) — owner's entry
  end
  ```
- **Plan build:** factor block logic (`:453`) into `_as_block(idxs)`, reuse for source blocks.
  `_connect_plan`: `state_index(src,name)` non-`nothing` ⇒ state partition (existing); else
  `findfirst(==(name), monitor_names(src))` ⇒ monitor partition, split by op. `_entry` builds
  `mon_overwrites`/`mon_adds` = `(components[src], _as_block(si[src]), mi, dst_param)`, and `injects`
  (owner `ck`) = `(components[src], _as_block(si[src]), mi, si[ck][owner_local])` where
  `owner_local = state_index(components[ck], dst_state)`. Set `params` when any connect partition
  (state OR monitor) is non-empty; injects don't use `params`.
- **Eval (allocation-free):**
  ```
  _run!: _connect!(U, t, e.params, e.overwrites, e.adds, e.mon_overwrites, e.mon_adds)
         e.model(view(dU,e.block), view(U,e.block), p, t)
         <frozen zeroing, unchanged>
         _inject!(dU, U, t, e.injects)     # dU[g] += monitor_values(m, view(U,b), t)[mi]
  ```
  `_connect!` widens with `t` + the two monitor partitions (add a `::Nothing` base case): zero all
  add-target slots, apply overwrites (state `for`-loop + monitor recursion), apply adds (state
  `for` + monitor recursion); monitor values go through `_connect_value`. `_inject!` is a recursive
  helper with a `::Tuple{}` base case, added after the owner's model write.
- **Validation:** connect source check (`_check_state` at `:249`) widens to accept a state **or**
  monitor name (new `_check_connect_source`). New inject validation: `src`/`dst` exist,
  `src_flux ∈ monitor_names(src)`, `dst_state ∈ state_names(dst)`. Op-conflict check (`:276`) is
  source-agnostic — no change. `_copy_connect_receivers` unchanged (sources/inject-targets are
  read-only w.r.t. component state).

### Step 3 — `examples/coupling_mwe.jl`: go live
- `ToyDriver`: `monitor_values(m,u,t) = (m.parameters[6]-u[3],)`  (was `monitor_values!`).
- `ToyResponder`: monitors `(:n,:J)` → `num_monitors=2`, `monitor_names=(:n,:J)`,
  `monitor_values(m,u,t) = (m.parameters[11]-u[3], responder_flux_J(m,u))` (keep `responder_flux_J`).
- Socket 3 live: `connect(:D=>:b, :R=>:p_b)`; assert R's staged `p_b` tracks `Ca − a` live from `U`
  and R's redox derivative reflects it.
- Socket 5 live: `inject(:R=>:J, :D=>:w)`; assert `dU[w] ≈ (Pw−Lw·w)+0.5·z` structurally and that
  `w`'s steady state rises above `Pw/Lw=1` when `z>0`.
- Delete the socket-3/5 TARGET API block; keep a one-line pointer to the documented switch convention.

### Step 4 — Tests (`test/test_coupling.jl`)
Mirror existing style (`cm(dU,U,nothing,0.0)`, hand-computed `dU[state_index(cm,·)]`, warm-up +
`@allocated … == 0`): DERIVED-source connect (receiver derivative tracks the monitor live, zero-alloc);
inject (`dU[owner] == base + flux`, zero-alloc, `operator_order` unchanged); validation throws
(unknown connect state/monitor; non-monitor inject `src_flux`; unknown inject `dst_state`). Add a
monitor mock (extend a `_Mono*` model or add one). Switch `test/test_monitors.jl` models to
`monitor_values`.

### Step 5 — Docs
- `examples/coupling_mwe.md`: flip sockets 3 & 5 ❌→✅; update API backlog (1–2 done, 3 documented).
- `CLAUDE.md`: Coupling section (document `inject` as the **third edge kind** + DERIVED-source
  connect); Derived observables (`monitor_values` tuple primitive vs derived `monitor_values!`);
  Variable roles (the "feedback foot-gun (not built)" note is now BUILT via `inject`); bless the
  module-switch convention (construction-time, resolve at construction, **never** through
  `p`/`SpatialContext`; OFF recovers baseline).
- `runic --inplace src test` before finishing.

## Gotchas / constraints

- **Allocation-free is a hard invariant.** State partitions (`(Int,Int)` tuples) stay `for`-loops
  (homogeneous → stable). **Monitor/inject partitions can be heterogeneous** (different source model
  types), so a `for`-loop over them is type-unstable and may box/allocate — walk them with
  `first`/`Base.tail` **recursion** (same trick `_run!` uses), never a `for`-loop. Verify
  `@allocated cm(dU,U,nothing,0.0) == 0` for the new DERIVED-connect and inject paths.
- **`inject` needs NO `_operator_order` change.** It reads the source's *state* (`U`, fixed during
  the eval), not the source's *derivative* — so it doesn't matter whether the source's entry has run.
  The only ordering that matters is local (the `+=` after the owner writes its base), guaranteed
  because the inject step lives on the owner's own entry. (Contrast `share`, which forces owner-last.)
- **Don't break hard-discard `share`.** The `frozen`-zeroing path stays for plain shares; `inject`
  is additive and separate. Additive is a *new* mode, not a replacement.
- **ForwardDiff:** connect monitor-source → `_connect_value` (strip primal to a `Float64` slot);
  inject → keep the `Dual` and add into `dU` (full Jacobian). `monitor_values` is just arithmetic on
  `view(U,block)` + `model.parameters`, so a `Dual` state promotes cleanly.
- **Socket-9 test/assertion change:** R gaining a `:J` monitor means the plain-couple monitor history
  names become `(:b, :R_n, :R_J)` (was `(:b, :R_n)`). Update the socket-9 assertion in the MWE and
  any `test/test_monitors.jl` expectations.
- **`view(U, block)` is zero-alloc even for a `Vector` (non-contiguous) block** — confirmed by the
  existing scatter zero-alloc test — so monitor-source staging over a scattered source block is safe.
- **Source model reference in the plan:** resolve monitor/inject sources from the post-deepcopy
  `components[src]` in the plan builder (runs after `_copy_connect_receivers`), so a source that is
  also a connect receiver refers to the same (copied) instance.
- **No secrets** anywhere in this work; nothing redacted.
