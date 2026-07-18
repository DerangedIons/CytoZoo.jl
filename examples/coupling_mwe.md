# Coupling MWE — socket map

`coupling_mwe.jl` is the canonical **driver for CytoZoo's coupling API**. It is a minimal working
example: two tiny models (~8 states total) that between them exercise the *entire* coupling taxonomy
from the ECCMitoRedox architecture — feedforward **and** feedback, including the patterns CytoZoo
cannot express yet.

**Why a toy.** We were reverse-engineering the API from the real coupled model (Gauthier ECC, 76
states + Kembro mito, 25 states). Only the *set of coupling patterns* drives API design; the real
model's *scale* (more variables/parameters) adds nothing. The toy keeps the whole taxonomy in view
at once and turns each missing capability into a one-line, executable spec.

## The models

- **Driver `D`** (`ToyDriver`, ≈ Gauthier ECC) — states `u, v, a, w`. Produces signals: a pulsed
  `u`, a slow `v`, a conserved-pool energy `a` (with derived partner `b = Ca − a`), and a shared
  feedback target `w`.
- **Responder `R`** (`ToyResponder`, ≈ Kembro mito) — the core downstream model, states `y, m, e`,
  with WIRE receiver parameter slots `p_u, p_v` and a held-as-param `h`.
- **`ToyRedox`** — an *optional* subsystem (state `z`, the redox module). Compose it in to switch
  redox on; leave it out to switch off.
- **`ToyH`** — an *optional* subsystem (state `h`). Compose it in (with a `connect` edge) and `R`
  reads `h` as a live state; leave it out and `R` holds `h` as a fixed parameter.

Optional capabilities are their **own subsystems**, composed in or left out — there is no
construction-time switch, and nothing routes through `p`.

## Socket map

Each row is a *role change* a coupling imposes. ✅ = expressible today; ❌ = drives new API.

| # | Socket | Coupling pattern | Toy expression | ECCMitoRedox analogue | CytoZoo |
|---|--------|------------------|----------------|-----------------------|---------|
| 1 | `u` | Feedforward **WIRE (overwrite)** — STATE → PARAM slot | `connect(:D=>:u, :R=>:p_u)` | Cai → K.Cai | ✅ |
| 2 | `v` | Feedforward WIRE (overwrite) | `connect(:D=>:v, :R=>:p_v)` | Nai → K.VNai | ✅ |
| 3 | `b` | **DERIVED-source WIRE** — DERIVED → PARAM slot | *want* `connect(:D=>:b, :R=>:p_b)` | ADPi = C_A − ATPi | ❌ |
| 4 | `a`/`e` | **Adopt-native / drop receiver state** — owner wins, other discarded | `share(:D=>:a, :R=>:e; owner=:D)` | ATPi G-native | ✅ (hard-discard `share`) |
| 5 | `w` | **Feedback additive contributed-flux** into a shared derivative | *want* Redox's `+J` summed into D's `dw` | NADH += −V_THD; ΔΨm += +V_IMAC | ❌ |
| 6 | redox | **Module on/off = compose with/without a subsystem** | include/omit `Subsystem(ToyRedox())` | redox_on / CII_dynamic | ✅ |
| 7 | edges | **Edge on/off = compose with/without an edge** | include/omit the WIRE edges at `couple()` | cyto_ions_dynamic | ✅ |
| 8 | `h` | **State↔param role flip = compose with/without a subsystem** | include/omit `Subsystem(ToyH())` + `connect` | Hm/Pim held as G params | ✅ |
| 9 | `b`(D), `n`(R) | **DERIVED / monitor** (conservation law) | `monitor_names` / `monitor_values!` | ATPm = C_A − ADPm | ✅ |

Rows 6–8 are the **same mechanism** at three granularities: switching *is* composition — omit an
element (a whole subsystem, or an edge) and the rest recovers its baseline (the OFF-invariant).
There is no switch primitive and no construction-time switch kwarg.

## Guards (built into the `.jl`)

- **G1 — isolation.** Omit the WIRE edges ⇒ R's sub-trajectory reproduces standalone R *exactly*;
  omit the `ToyRedox` subsystem ⇒ D and R recover their baseline (`z` is a leaf nothing reads).
- **G2 — one-way.** One pulse into `u` lifts R's `y`, which peaks early and returns toward baseline
  (no drift).
- **Closures.** `a + b = Ca` and `m + n = Cm` hold every step (asserted via `monitor_history`).

## API backlog (the ❌ rows)

Each is a concrete design question the toy exists to settle:

1. **Additive / contributory share** (socket 5) — *the central feedback gap.* Let a non-owner add a
   flux into an owner's shared-slot derivative, vs. today's hard-discard. Candidate spellings:
   `share(:D=>:w, :Redox=>:w_flux; owner=:D, op=+)`, or a new `inject(:Redox=>:J, :D=>:w)` edge.
   Implementation: an accumulate-into-owner's-`dU` path in `_run!` (`src/coupling.jl`), a contributory
   alternative to the `frozen`-index zeroing.
2. **DERIVED-source `connect`** (socket 3) — let a `connect` source resolve to a monitor/derived value
   (compute it each RHS), not only a `state_index` (`src/coupling.jl`).

Design each against this toy, then apply the finalized API to ECCMitoRedox (downstream).

## Run

```bash
julia --project=examples examples/coupling_mwe.jl
```

Prints each socket's live demonstration; the assertions fail loudly if a pattern regresses. The
TARGET API section (bottom of the `.jl`) holds the not-yet-expressible calls as executable specs.
