# Coupling MWE — socket map

`coupling_mwe.jl` is the canonical **driver for CytoZoo's coupling API**. It is a minimal working
example: two tiny models (~8 states total) that between them exercise the *entire* coupling taxonomy
from the ECCMitoRedox architecture — feedforward **and** feedback, including the patterns CytoZoo
cannot express yet.

**Why a toy.** We were reverse-engineering the API from the real coupled model (Gauthier ECC, 76
states + Kembro mito, 25 states). Only the *set of coupling patterns* drives API design; the real
model's *scale* (more variables/parameters) adds nothing. The toy keeps the whole taxonomy in view
at once and turns each missing capability into a one-line, executable spec.

## The two models

- **Driver `D`** (`ToyDriver`, ≈ Gauthier ECC) — states `u, v, a, w`. Produces signals: a pulsed
  `u`, a slow `v`, a conserved-pool energy `a` (with derived partner `b = Ca − a`), and a shared
  feedback target `w`.
- **Responder `R`** (`ToyResponder`, ≈ Kembro mito) — states `y, z, m, e` (+ `h` when integrated),
  and WIRE receiver parameter slots `p_u, p_v, p_b`. Two construction switches gate its own
  terms/states: `redox_on` (the `z`/feedback term) and `dynamic_h` (whether `h` is a state or a held
  parameter). Switches resolve **at construction**, never through `p`.

## Socket map

Each row is a *role change* a coupling imposes. ✅ = expressible today; ❌ = drives new API.

| # | Socket | Coupling pattern | Toy expression | ECCMitoRedox analogue | CytoZoo |
|---|--------|------------------|----------------|-----------------------|---------|
| 1 | `u` | Feedforward **WIRE (overwrite)** — STATE → PARAM slot | `connect(:D=>:u, :R=>:p_u)` | Cai → K.Cai | ✅ |
| 2 | `v` | Feedforward WIRE (overwrite) | `connect(:D=>:v, :R=>:p_v)` | Nai → K.VNai | ✅ |
| 3 | `b` | **DERIVED-source WIRE** — DERIVED → PARAM slot | *want* `connect(:D=>:b, :R=>:p_b)` | ADPi = C_A − ATPi | ❌ |
| 4 | `a`/`e` | **Adopt-native / drop receiver state** — owner wins, other discarded | `share(:D=>:a, :R=>:e; owner=:D)` | ATPi G-native | ✅ (hard-discard `share`) |
| 5 | `w` | **Feedback additive contributed-flux** into a shared derivative | *want* R's `+J` summed into D's `dw` | NADH += −V_THD; ΔΨm += +V_IMAC | ❌ |
| 6 | `redox_on` | **Intra-model term switch** (OFF-invariant) | `ToyResponder(; redox_on=false)` | redox_on / CII_dynamic | ✅* |
| 7 | `cyto_ions_dynamic` | **Edge-gating switch** (OFF-invariant) | include/omit edges 1,2 at `couple()` | cyto_ions_dynamic | ✅ |
| 8 | `h` | **Inverse: state-held-as-param** | `ToyResponder(; dynamic_h=false)` | Hm/Pim held as G params | ✅* |
| 9 | `b`(D), `n`(R) | **DERIVED / monitor** (conservation law) | `monitor_names` / `monitor_values!` | ATPm = C_A − ADPm | ✅ |

`*` sockets 6 & 8 work today as model-construction switches; what's missing is a *blessed CytoZoo
convention* for them (see backlog item 3).

## Guards (built into the `.jl`)

- **G1 — isolation.** `cyto_ions_dynamic=false` ⇒ edges 1,2 dropped ⇒ R's sub-trajectory reproduces
  standalone R *exactly*. `redox_on=false` ⇒ `z ≡ 0` (D reproduces its baseline once `J` is wired).
- **G2 — one-way.** One pulse into `u` lifts R's `y`, which peaks early and returns toward baseline
  (no drift).
- **Closures.** `a + b = Ca` and `m + n = Cm` hold every step (asserted via `monitor_history`).

## API backlog (the ❌ rows)

Each is a concrete design question the toy exists to settle:

1. **Additive / contributory share** (socket 5) — *the central feedback gap.* Let a non-owner add a
   flux into an owner's shared-slot derivative, vs. today's hard-discard. Candidate spellings:
   `share(:D=>:w, :R=>:w_flux; owner=:D, op=+)`, or a new `inject(:R=>:J, :D=>:w)` edge. Implementation:
   an accumulate-into-owner's-`dU` path in `_run!` (`src/coupling.jl:466`), a contributory alternative
   to the `frozen`-index zeroing.
2. **DERIVED-source `connect`** (socket 3) — let a `connect` source resolve to a monitor/derived value
   (compute it each RHS), not only a `state_index` (`src/coupling.jl:575`).
3. **Module-switch protocol** (sockets 6, 8) — a declare/accept convention so any participant exposes
   construction-time term/state switches uniformly, with a documented guarantee they never route
   through `p`/`SpatialContext`. Items 1–2 are CytoZoo-core; item 3 is a model-package convention
   CytoZoo documents.

Design each against this toy, then apply the finalized API to ECCMitoRedox (downstream).

## Run

```bash
julia --project=examples examples/coupling_mwe.jl
```

Prints each socket's live demonstration; the assertions fail loudly if a pattern regresses. The
TARGET API section (bottom of the `.jl`) holds the not-yet-expressible calls as executable specs.
