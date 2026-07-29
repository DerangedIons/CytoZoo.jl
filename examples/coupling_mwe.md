# Coupling MWE — socket map

`coupling_mwe.jl` is the canonical **driver for CytoZoo's coupling API**. It is a minimal working
example: two tiny models (~8 states total) that between them exercise the *entire* coupling taxonomy
from the ECCMitoRedox architecture — feedforward **and** feedback. Every socket in the map is now
live.

**Why a toy.** We were reverse-engineering the API from the real coupled model (Gauthier ECC, 76
states + Kembro mito, 25 states). Only the *set of coupling patterns* drives API design; the real
model's *scale* (more variables/parameters) adds nothing. The toy keeps the whole taxonomy in view
at once and turned each missing capability into a one-line, executable spec until it was built.

## The models

- **Driver `D`** (`ToyDriver`, ≈ Gauthier ECC) — states `u, v, a, w`. Produces signals: a pulsed
  `u`, a slow `v`, a conserved-pool energy `a` (with derived partner `b = Ca − a`), and a shared
  feedback target `w`.
- **Responder `R`** (`ToyResponder`, ≈ Kembro mito) — the core downstream model, states `y, m, e`,
  with WIRE receiver parameter slots `p_u, p_v` and a held-as-param `h`.
- **`ToyRedox`** — an *optional* subsystem (states `z` and `w_flux`, the redox module). Compose it
  in to switch redox on; leave it out to switch off. `w_flux` carries the term it contributes back
  to `D`: its derivative is the flux `J = gJ·z`, summed into `D`'s `dw` by a contributory share.
- **`ToyH`** — an *optional* subsystem (state `h`). Compose it in (with a `connect` edge) and `R`
  reads `h` as a live state; leave it out and `R` holds `h` as a fixed parameter.

Optional capabilities are their **own subsystems**, composed in or left out — there is no
construction-time switch, and nothing routes through `p`.

## Socket map

Each row is a *role change* a coupling imposes. ✅ = expressible today.

| # | Socket | Coupling pattern | Toy expression | ECCMitoRedox analogue | CytoZoo |
|---|--------|------------------|----------------|-----------------------|---------|
| 1 | `u` | Feedforward **WIRE (overwrite)** — STATE → PARAM slot | `connect(:D=>:u, :R=>:p_u)` | Cai → K.Cai | ✅ |
| 2 | `v` | Feedforward WIRE (overwrite) | `connect(:D=>:v, :R=>:p_v)` | Nai → K.VNai | ✅ |
| 3 | `b` | **DERIVED-source WIRE** — DERIVED → PARAM slot | `connect(:D=>:b, :Redox=>:p_b)` | ADPi = C_A − ATPi | ✅ (monitor source) |
| 4 | `a`/`e` | **Adopt-native / drop receiver state** — owner wins, other discarded | `share(:D=>:a, :R=>:e; owner=:D)` | ATPi G-native | ✅ (hard-discard `share`) |
| 5 | `w` | **Feedback additive contributed-flux** into a shared derivative | `share(:D=>:w, :Redox=>:w_flux; owner=:D, op=+)` | NADH += −V_THD; ΔΨm += +V_IMAC | ✅ (contributory `share`) |
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
- **G3 — contribution.** Socket 5's `dw` is checked against a hand-rolled `Pw − Lw·w + J`, and the
  *same graph* with the default `op` must give `Pw − Lw·w` — one keyword apart, so a leak in either
  direction fails loudly.
- **Closures.** `a + b = Ca` and `m + n = Cm` hold every step (asserted via `monitor_history`).

Socket 3 needed **no change to the toy models at all**: `ToyDriver` already declared `b = Ca − a`
as a monitor for socket 9, and a `connect` source resolves against states *then* monitors. Wiring a
derived quantity and observing one are the same declaration — the conservation law stays inside the
model that owns it instead of being restated at the edge.

Socket 5 needed **one new state on `ToyRedox`** and nothing else. A module contributes a term by
carrying it as an ordinary state whose derivative is the flux; `op = +` decides that the shared slot
sums rather than discards. Written standalone, `ToyRedox` is unaware it is coupled — the same
authoring property hard-discard `share` has.

## API backlog

Empty. Every socket is live.

What remains is an alternative *spelling* for socket 5, not a missing capability: a new
`inject(:Redox => :J, :D => :w)` edge sourcing a named DERIVED flux. It was not built, for two
reasons worth keeping. A monitor source is resolved in a pre-pass, so the sourcing component may not
itself receive a `connect` edge — which the real consumer does, to read its driver's `V_O2`,
`V_SDH`, `V_He_F`; and monitor values pass through `_connect_value`, which extracts the primal of a
`Dual`, so an injected flux would silently vanish from every derivative. A share never leaves `U`
and has neither problem. Revisit only for a module that genuinely cannot carry the term as a state.

## Run

```bash
julia --project=examples examples/coupling_mwe.jl
```

Prints each socket's live demonstration; the assertions fail loudly if a pattern regresses. The
TARGET API section (bottom of the `.jl`) holds the not-yet-expressible calls as executable specs.
