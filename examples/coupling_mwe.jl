# Coupling MWE — the canonical driver for CytoZoo's coupling API.
#
# Two tiny models (Driver D, Responder R; ~8 states total) that between them exercise the ENTIRE
# coupling taxonomy from the ECCMitoRedox architecture (feedforward + feedback), in a form small
# enough to reason about whole. It is deliberately abstract: only the *set of coupling patterns*
# drives API design — the real model's 101-state *scale* does not. See `coupling_mwe.md` for the
# socket map and the D↔Gauthier / R↔Kembro analogue mapping.
#
# New to CytoZoo coupling? Start with examples/coupling_toy.jl — a minimal 2-pattern (share +
# connect) intro. This file is the comprehensive driver, not a gentle first read.
#
# Socket map (✅ = expressible today, ❌ = drives new API):
#   1 u   feedforward WIRE (overwrite): D state → R param slot            ✅  connect(:D=>:u, :R=>:p_u)
#   2 v   feedforward WIRE (overwrite)                                    ✅  connect(:D=>:v, :R=>:p_v)
#   3 b   DERIVED-source WIRE: D derived → R param slot                   ❌  connect sources a STATE only
#   4 a/e adopt-native / drop receiver state (owner wins)                 ✅  share(:D=>:a, :R=>:e; owner=:D)
#   5 w   feedback additive contributed-flux into a shared derivative     ❌  share=hard-discard; connect writes a param
#   6 redox_on   intra-model term switch (OFF-invariant)                  ✅* model construction switch (protocol ❌)
#   7 cyto_ions  edge-gating switch (OFF-invariant)                       ✅  include/omit edges at couple()
#   8 h   inverse: state-held-as-param                                    ✅* model construction switch (protocol ❌)
#   9 b,n DERIVED / monitor (conservation law)                           ✅  monitor hooks
#
# Run (needs an OrdinaryDiffEq solver):
#   julia --project examples/coupling_mwe.jl

using Pkg
Pkg.activate(@__DIR__)

using CytoZoo
using OrdinaryDiffEq

# =====================================================================================
# Driver D — the upstream model (≈ Gauthier ECC). States u,v,a,w; nominal Vm = u.
# =====================================================================================

const DRIVER_STATES = (:u, :v, :a, :w)
const DRIVER_PARAMS = (:ku, :kv, :v0, :kf, :kr, :Ca, :Pw, :Lw)

struct ToyDriver{T, S} <: AbstractCardiacCellModel
    parameters::Vector{T}          # DRIVER_PARAMS, in order
    stim::S                        # Istim for u; a CytoZoo AbstractStimulus
end

# Defaults: a* = kf·Ca/(kf+kr) = 2, b* = Ca-a* = 2, w* = Pw/Lw = 1.
function ToyDriver(;
        stim = Stimulus(; amplitude = 20.0, period = 1.0, duration = 0.05, start = 0.2),
        elT::Type = Float64,
    )
    θ = elT[1.0, 0.5, 2.0, 1.0, 1.0, 4.0, 1.0, 1.0]   # ku,kv,v0,kf,kr,Ca,Pw,Lw
    return ToyDriver(θ, stim)
end

CytoZoo.num_states(::ToyDriver) = 4
CytoZoo.num_parameters(m::ToyDriver) = length(m.parameters)
CytoZoo.state_names(::ToyDriver) = DRIVER_STATES
CytoZoo.parameter_names(::ToyDriver) = DRIVER_PARAMS
CytoZoo.state_index(::ToyDriver, n::Symbol) = findfirst(==(n), DRIVER_STATES)
CytoZoo.parameter_index(::ToyDriver, n::Symbol) = findfirst(==(n), DRIVER_PARAMS)
CytoZoo.transmembrane_potential_index(::ToyDriver) = 1
CytoZoo.default_initial_state(::ToyDriver{T}) where {T} = T[0, 0.5, 1, 0]   # u,v,a,w

function (m::ToyDriver)(du, u, p, t)
    θ = m.parameters
    ku, kv, v0, kf, kr, Ca, Pw, Lw = θ[1], θ[2], θ[3], θ[4], θ[5], θ[6], θ[7], θ[8]
    du[1] = -ku * u[1] + m.stim(nothing, t)   # u — pulsed driver (≈ Cai)
    du[2] = kv * (v0 - u[2])                   # v — slow accumulator (≈ Nai)
    du[3] = kf * (Ca - u[3]) - kr * u[3]       # a — energy; conserved pool a + b = Ca (≈ ATPi)
    du[4] = Pw - Lw * u[4]                      # w — shared feedback target (≈ ΔΨm/NADH); R adds +J
    return nothing
end

# DERIVED monitor: b = Ca - a (conservation a + b = Ca).
CytoZoo.num_monitors(::ToyDriver) = 1
CytoZoo.monitor_names(::ToyDriver) = (:b,)
function CytoZoo.monitor_values!(mon, u, t, m::ToyDriver)
    mon[1] = m.parameters[6] - u[3]            # Ca - a
    return nothing
end

# =====================================================================================
# Responder R — the downstream model (≈ Kembro mito). States y,z,m,e (+h if dynamic_h).
# Two construction switches gate its OWN terms/states — resolved at construction, never through p.
# =====================================================================================

const RESPONDER_PARAMS =
    (:p_u, :p_v, :p_b, :kin, :ky, :shunt, :kz, :ke, :e0, :km, :Cm, :kh, :h0)

struct ToyResponder{T} <: AbstractCardiacCellModel
    parameters::Vector{T}          # RESPONDER_PARAMS, in order (p_u/p_v/p_b are WIRE receiver slots)
    redox_on::Bool                 # socket 6 — gate the z/J redox term (OFF ⇒ z ≡ 0)
    dynamic_h::Bool                # socket 8 — h integrated (true) or held as param h0 (false)
end

function ToyResponder(; redox_on::Bool = true, dynamic_h::Bool = false, elT::Type = Float64)
    #    p_u  p_v  p_b  kin  ky  shunt kz  ke  e0  km  Cm  kh  h0
    θ = elT[0.1, 1.0, 0.0, 1.0, 1.0, 0.2, 1.0, 0.5, 3.0, 1.0, 3.0, 0.5, 0.5]
    return ToyResponder(θ, redox_on, dynamic_h)
end

# num_states / state_names depend on the dynamic_h switch — the inverse-case (state ↔ param).
CytoZoo.state_names(m::ToyResponder) = m.dynamic_h ? (:y, :z, :m, :e, :h) : (:y, :z, :m, :e)
CytoZoo.num_states(m::ToyResponder) = m.dynamic_h ? 5 : 4
CytoZoo.num_parameters(m::ToyResponder) = length(m.parameters)
CytoZoo.parameter_names(::ToyResponder) = RESPONDER_PARAMS
CytoZoo.state_index(m::ToyResponder, n::Symbol) = findfirst(==(n), state_names(m))
CytoZoo.parameter_index(::ToyResponder, n::Symbol) = findfirst(==(n), RESPONDER_PARAMS)
CytoZoo.transmembrane_potential_index(::ToyResponder) = 1
CytoZoo.default_initial_state(m::ToyResponder{T}) where {T} =
    m.dynamic_h ? T[0, 0, 0, 0, 0.5] : T[0, 0, 0, 0]   # y,z,m,e[,h]

function (m::ToyResponder)(du, u, p, t)
    θ = m.parameters
    p_u, p_v, p_b, kin, ky = θ[1], θ[2], θ[3], θ[4], θ[5]
    shunt, kz, ke, e0, km = θ[6], θ[7], θ[8], θ[9], θ[10]
    du[1] = p_u * kin - p_v * ky * u[1]                       # y — driven by wired u (influx), v (efflux)
    du[2] = m.redox_on ? (shunt * (p_u + p_b) - kz * u[2]) : zero(eltype(du))  # z — redox, gated
    du[3] = km * (u[4] - u[3])                                # m — tracks energy e; pool m + n = Cm
    du[4] = ke * (e0 - u[4])                                  # e — own energy (shared with D.a ⇒ discarded)
    if m.dynamic_h
        du[5] = θ[12] * (θ[13] - u[5])                       # h — kh·(h0 - h), only when integrated
    end
    return nothing
end

# DERIVED monitor: n = Cm - m (conservation m + n = Cm).
CytoZoo.num_monitors(::ToyResponder) = 1
CytoZoo.monitor_names(::ToyResponder) = (:n,)
function CytoZoo.monitor_values!(mon, u, t, m::ToyResponder)
    mon[1] = m.parameters[11] - u[3]           # Cm - m
    return nothing
end

# The feedback flux R contributes upstream (socket 5). It is a function of R's own state, ready to
# be injected into D's `w` derivative — but no edge kind can carry it yet (see TARGET API below).
responder_flux_J(m::ToyResponder, u_R) = 0.5 * u_R[2]   # J = gJ·z

# =====================================================================================
# Live demonstrations (everything CytoZoo expresses today)
# =====================================================================================

approx(a, b; atol = 1e-8) = isapprox(a, b; atol = atol)

# --- sockets 1,2,7 — feedforward WIREs + cyto_ions_dynamic gate, with guards G1/G2 ------------
function build_wired(; cyto_ions_dynamic::Bool, stim)
    D = ToyDriver(; stim)
    R = ToyResponder()
    edges = cyto_ions_dynamic ?
        (connect(:D => :u, :R => :p_u), connect(:D => :v, :R => :p_v)) : ()
    return couple((Subsystem(D; name = :D), Subsystem(R; name = :R)), edges)
end

println("── sockets 1,2,7 — feedforward WIREs + edge-gating switch ──")
let stim = ToyDriver().stim, tspan = (0.0, 5.0), opts = (dt = 0.01, adaptive = false)
    # G1 isolation: edges OFF ⇒ R's sub-trajectory reproduces standalone R exactly.
    cm_off = build_wired(cyto_ions_dynamic = false, stim = stim)
    @assert length(cm_off.connects) == 0
    R = ToyResponder()
    solR = solve(ODEProblem(R, tspan), Tsit5(); opts...)
    solC = solve(ODEProblem(cm_off, tspan), Tsit5(); opts...)
    kidx = cm_off.layout.solution_indices.R
    g1 = maximum(maximum(abs, solC.u[j][kidx] .- solR.u[j]) for j in eachindex(solR.t))
    @assert approx(g1, 0.0; atol = 1e-9) "G1 isolation broken: max |coupled_R - standalone_R| = $g1"
    println("  G1 isolation (edges off ⇒ standalone R):  max |Δ| = ", g1)

    # G2 one-way: one pulse into u lifts R's y, which then returns (difference vs a never-paced run).
    onepulse = Stimulus(; amplitude = 20.0, period = 1e6, duration = 0.05, start = 0.5)
    norest = Stimulus(; amplitude = 20.0, period = 1e6, duration = 0.05, start = 1e9)
    cm_ap = build_wired(cyto_ions_dynamic = true, stim = onepulse)
    cm_rest = build_wired(cyto_ions_dynamic = true, stim = norest)
    @assert length(cm_ap.connects) == 2
    sap = solve(ODEProblem(cm_ap, tspan), Tsit5(); opts...)
    srest = solve(ODEProblem(cm_rest, tspan), Tsit5(); opts...)
    yidx = state_index(cm_ap, :R_y)
    dY = [sap.u[j][yidx] - srest.u[j][yidx] for j in eachindex(sap.t)]
    ipeak = argmax(dY)
    @assert dY[ipeak] > 0 "G2: the pulse did not lift y"
    @assert sap.t[ipeak] < (tspan[2] / 2) "G2: lift is not an early transient"
    @assert dY[end] < dY[ipeak] "G2: y is not returning toward baseline"
    println("  G2 one-way (pulse ⇒ y lifts then returns):  peak Δy = ",
            round(dY[ipeak]; digits = 4), " at t = ", sap.t[ipeak], ",  Δy(end) = ",
            round(dY[end]; digits = 4))
end

# --- socket 4 — share / adopt-native (owner D wins; R's own equation discarded) --------------
println("── socket 4 — share (adopt-native / drop receiver state) ──")
let tspan = (0.0, 8.0)
    cm = couple(
        (Subsystem(ToyDriver(); name = :D), Subsystem(ToyResponder(); name = :R)),
        (share(:D => :a, :R => :e; owner = :D),),
    )
    # Structural: the shared slot carries D's da; R's de is discarded (frozen).
    U = default_initial_state(cm)
    dU = similar(U)
    cm(dU, U, nothing, 0.0)
    a_slot = state_index(cm, :a)
    D = ToyDriver()
    dU_D = similar(default_initial_state(D)); D(dU_D, default_initial_state(D), nothing, 0.0)
    @assert approx(dU[a_slot], dU_D[3]) "share: owner D's da did not win the shared slot"

    # Behavioural: R's m tracks the shared energy ⇒ converges to D's a* (=2), NOT R's own e0 (=3).
    sol = solve(ODEProblem(cm, tspan), Tsit5(); dt = 0.01, adaptive = false)
    m_final = sol.u[end][state_index(cm, :R_m)]
    @assert approx(m_final, 2.0; atol = 1e-2) "share: R.m did not adopt D's energy a*"
    println("  shared slot = D's a (R's de discarded);  R.m → ", round(m_final; digits = 4),
            "  (D's a* = 2.0, not R's e0 = 3.0)")
end

# --- sockets 6,8 — intra-model construction switches + OFF-invariants ------------------------
println("── sockets 6,8 — module switches (construction-time, not through p) ──")
let tspan = (0.0, 5.0), opts = (dt = 0.01, adaptive = false)
    # socket 6: redox_on = false ⇒ z ≡ 0 (the OFF-invariant / regression baseline).
    Roff = ToyResponder(redox_on = false)
    solz = solve(ODEProblem(Roff, tspan), Tsit5(); opts...)
    zmax = maximum(abs(u[2]) for u in solz.u)
    @assert approx(zmax, 0.0) "redox OFF-invariant broken: z drifted from 0 (max |z| = $zmax)"
    println("  redox_on=false ⇒ z ≡ 0:  max |z| = ", zmax)

    # socket 8: dynamic_h flips h between an integrated STATE and a held PARAM — the state count changes.
    Rheld = ToyResponder(dynamic_h = false)
    Rdyn = ToyResponder(dynamic_h = true)
    @assert num_states(Rheld) == 4 && state_names(Rheld) == (:y, :z, :m, :e)
    @assert num_states(Rdyn) == 5 && state_names(Rdyn) == (:y, :z, :m, :e, :h)
    println("  dynamic_h=false ⇒ h is PARAM (4 states); dynamic_h=true ⇒ h is STATE (5 states)")
end

# --- socket 9 — DERIVED monitors / conservation-law closures --------------------------------
println("── socket 9 — DERIVED monitors (conservation closures) ──")
let tspan = (0.0, 5.0), opts = (dt = 0.05, adaptive = false)
    cm = couple((Subsystem(ToyDriver(); name = :D), Subsystem(ToyResponder(); name = :R)))
    sol = solve(ODEProblem(cm, tspan), Tsit5(); opts...)
    h = monitor_history(sol, cm)
    @assert h.names == (:b, :R_n) "coupled monitor names wrong: $(h.names)"
    a_i, m_i = state_index(cm, :a), state_index(cm, :R_m)
    Ca, Cm = ToyDriver().parameters[6], ToyResponder().parameters[11]
    okb = all(approx(h.values[1, j], Ca - sol.u[j][a_i]) for j in eachindex(sol.t))   # a + b = Ca
    okn = all(approx(h.values[2, j], Cm - sol.u[j][m_i]) for j in eachindex(sol.t))   # m + n = Cm
    @assert okb && okn "monitor conservation closure violated"
    println("  monitors ", h.names, ";  closures a+b=Ca and m+n=Cm hold every step")
end

# =====================================================================================
# TARGET API — patterns the taxonomy needs but CytoZoo cannot express yet.
# These are executable specs: the exact calls we WANT, plus what CytoZoo must add.
# =====================================================================================
#
# ── socket 5 — feedback additive contributed-flux into a shared derivative ───────────────────
#   R's flux J = gJ·z (responder_flux_J) must be ADDED into D's w equation:  dw = Pw - Lw·w + J.
#   `share` is hard-discard (owner's equation wins entirely; the non-owner's is zeroed via `frozen`)
#   and `connect` writes into a *parameter* slot, not a *derivative* — so neither can sum a
#   non-owner flux into an owner's dU. This is THE feedback gap (ECCMitoRedox §6: NADH += −V_THD,
#   ΔΨm += +V_IMAC). Two candidate spellings to decide between:
#
#       # (a) additive share — R exposes a local `w_flux` state whose derivative sums into the slot
#       share(:D => :w, :R => :w_flux; owner = :D, op = +)
#
#       # (b) flux injection — a new edge kind adding a named derived flux into an owner's derivative
#       inject(:R => :J, :D => :w)
#
#   Needed in CytoZoo: an accumulate-into-owner's-dU path in `_run!` (src/coupling.jl:466), i.e. a
#   contributory alternative to today's `frozen`-index zeroing.
#
# ── socket 3 — DERIVED-source WIRE ───────────────────────────────────────────────────────────
#   Wire D's DERIVED b (= Ca − a) into R's p_b slot so R's redox reads it live:
#
#       connect(:D => :b, :R => :p_b)
#
#   `connect` resolves its source only through `state_index` (src/coupling.jl:575), so a monitor /
#   derived name fails today. Needed: fall back to a monitor index (compute the derived value each
#   RHS) when the source name is not a state.
#
# ── sockets 6,8 — module-switch protocol ─────────────────────────────────────────────────────
#   `redox_on` / `dynamic_h` above WORK as model-construction switches. What is missing is a BLESSED
#   CytoZoo convention (a declare/accept protocol) so any participant can expose such switches
#   uniformly — and a documented guarantee they resolve at construction, never through `p` /
#   `SpatialContext` (the time/space-varying payload). Until then each model spells them ad hoc.

println("\nAll live demonstrations passed. See the TARGET API section for the API backlog.")
