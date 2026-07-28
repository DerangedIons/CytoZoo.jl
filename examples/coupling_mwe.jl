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
#   3 b   DERIVED-source WIRE: D derived → Redox param slot               ✅  connect(:D=>:b, :Redox=>:p_b)
#   4 a/e adopt-native / drop receiver state (owner wins)                 ✅  share(:D=>:a, :R=>:e; owner=:D)
#   5 w   feedback additive contributed-flux into a shared derivative     ❌  share=hard-discard; connect writes a param
#   6 redox    module on/off = compose with/without the ToyRedox subsystem  ✅  include/omit Subsystem(ToyRedox)
#   7 cyto_ions edge on/off = compose with/without the WIRE edges           ✅  include/omit edges at couple()
#   8 h        state↔param role flip = compose with/without the ToyH subsystem  ✅  include/omit Subsystem(ToyH) + connect
#   9 b,n DERIVED / monitor (conservation law)                           ✅  monitor hooks
#
# Sockets 6–8 are ONE mechanism at three granularities: switching is composition. To turn an
# optional capability off you compose without it (a whole subsystem, or an edge) and the rest
# recovers its baseline (the OFF-invariant). There is no switch primitive and no construction-time
# switch kwarg — anything switchable is a subsystem you include or leave out.
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
    du[4] = Pw - Lw * u[4]                      # w — shared feedback target (≈ ΔΨm/NADH); Redox adds +J
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
# Responder R — the downstream CORE model (≈ Kembro mito). States y, m, e.
# Optional capabilities are their OWN subsystems (ToyRedox, ToyH), composed in or left out —
# there is no construction-time switch. Here `h` is a held parameter; the ToyH subsystem (below)
# turns it into a live state via a connect edge (socket 8, the state↔param role flip).
# =====================================================================================

const RESPONDER_PARAMS = (:p_u, :p_v, :kin, :ky, :ke, :e0, :km, :Cm, :h)

struct ToyResponder{T} <: AbstractCardiacCellModel
    parameters::Vector{T}          # RESPONDER_PARAMS (p_u/p_v are WIRE receiver slots; h is held-as-param)
end

function ToyResponder(; elT::Type = Float64)
    #    p_u  p_v  kin  ky  ke  e0  km  Cm  h
    θ = elT[0.1, 1.0, 1.0, 1.0, 0.5, 3.0, 1.0, 3.0, 1.0]
    return ToyResponder(θ)
end

CytoZoo.num_states(::ToyResponder) = 3
CytoZoo.num_parameters(m::ToyResponder) = length(m.parameters)
CytoZoo.state_names(::ToyResponder) = (:y, :m, :e)
CytoZoo.parameter_names(::ToyResponder) = RESPONDER_PARAMS
CytoZoo.state_index(::ToyResponder, n::Symbol) = findfirst(==(n), (:y, :m, :e))
CytoZoo.parameter_index(::ToyResponder, n::Symbol) = findfirst(==(n), RESPONDER_PARAMS)
CytoZoo.transmembrane_potential_index(::ToyResponder) = 1
CytoZoo.default_initial_state(::ToyResponder{T}) where {T} = T[0, 0, 0]   # y,m,e

function (m::ToyResponder)(du, u, p, t)
    θ = m.parameters
    p_u, p_v, kin, ky, ke = θ[1], θ[2], θ[3], θ[4], θ[5]
    e0, km, h = θ[6], θ[7], θ[9]
    du[1] = p_u * kin * h - p_v * ky * u[1]   # y — wired influx (scaled by held/dynamic h), efflux
    du[2] = km * (u[3] - u[2])                # m — tracks energy e; pool m + n = Cm
    du[3] = ke * (e0 - u[3])                  # e — own energy (shared with D.a ⇒ discarded)
    return nothing
end

# DERIVED monitor: n = Cm - m (conservation m + n = Cm).
CytoZoo.num_monitors(::ToyResponder) = 1
CytoZoo.monitor_names(::ToyResponder) = (:n,)
function CytoZoo.monitor_values!(mon, u, t, m::ToyResponder)
    mon[1] = m.parameters[8] - u[2]            # Cm - m
    return nothing
end

# =====================================================================================
# ToyRedox — an OPTIONAL subsystem (≈ the redox module). State z. Compose it in to switch redox
# ON; leave it out to switch OFF (the rest recovers baseline, since z is a leaf nothing reads).
# It reads a wired driver signal through p_u; its flux J = gJ·z is the feedback R would contribute
# upstream (socket 5, not yet expressible).
# =====================================================================================

const REDOX_PARAMS = (:p_u, :p_b, :shunt, :kz)

struct ToyRedox{T} <: AbstractCardiacCellModel
    parameters::Vector{T}          # REDOX_PARAMS (p_u/p_b are WIRE receiver slots)
end

function ToyRedox(; elT::Type = Float64)
    θ = elT[0.0, 0.0, 0.2, 1.0]   # p_u, p_b, shunt, kz
    return ToyRedox(θ)
end

CytoZoo.num_states(::ToyRedox) = 1
CytoZoo.num_parameters(m::ToyRedox) = length(m.parameters)
CytoZoo.state_names(::ToyRedox) = (:z,)
CytoZoo.parameter_names(::ToyRedox) = REDOX_PARAMS
CytoZoo.state_index(::ToyRedox, n::Symbol) = findfirst(==(n), (:z,))
CytoZoo.parameter_index(::ToyRedox, n::Symbol) = findfirst(==(n), REDOX_PARAMS)
CytoZoo.transmembrane_potential_index(::ToyRedox) = 1
CytoZoo.default_initial_state(::ToyRedox{T}) where {T} = T[0]   # z

function (m::ToyRedox)(du, u, p, t)
    p_u, p_b, shunt, kz = m.parameters[1], m.parameters[2], m.parameters[3], m.parameters[4]
    du[1] = shunt * (p_u + p_b) - kz * u[1]    # z — redox accumulator driven by wired inputs
    return nothing
end

# The feedback flux the redox module contributes upstream (socket 5): J = gJ·z, ready to be added
# into D's `w` derivative — but no edge kind can carry it yet (see TARGET API below).
redox_flux_J(m::ToyRedox, u) = 0.5 * u[1]   # J = gJ·z

# =====================================================================================
# ToyH — an OPTIONAL subsystem carrying the `h` dynamics. State h. Compose it in (with the connect
# edge below) and R reads `h` as a live STATE; leave it out and R holds `h` as a fixed parameter.
# This is the state↔param role flip (socket 8), expressed purely by composition.
# =====================================================================================

const H_PARAMS = (:kh, :h0)

struct ToyH{T} <: AbstractCardiacCellModel
    parameters::Vector{T}          # H_PARAMS, in order
end

ToyH(; elT::Type = Float64) = ToyH(elT[0.5, 1.0])   # kh, h0

CytoZoo.num_states(::ToyH) = 1
CytoZoo.num_parameters(m::ToyH) = length(m.parameters)
CytoZoo.state_names(::ToyH) = (:h,)
CytoZoo.parameter_names(::ToyH) = H_PARAMS
CytoZoo.state_index(::ToyH, n::Symbol) = findfirst(==(n), (:h,))
CytoZoo.parameter_index(::ToyH, n::Symbol) = findfirst(==(n), H_PARAMS)
CytoZoo.transmembrane_potential_index(::ToyH) = 1
CytoZoo.default_initial_state(::ToyH{T}) where {T} = T[0.2]   # h starts below h0 ⇒ visibly ramps

function (m::ToyH)(du, u, p, t)
    kh, h0 = m.parameters[1], m.parameters[2]
    du[1] = kh * (h0 - u[1])   # h — relaxes toward h0
    return nothing
end

# =====================================================================================
# Live demonstrations (everything CytoZoo expresses today)
# =====================================================================================

approx(a, b; atol = 1e-8) = isapprox(a, b; atol = atol)

# --- sockets 1,2,7 — feedforward WIREs + compose-with/without the WIRE edges, guards G1/G2 ------
function build_wired(; wire_edges::Bool, stim)
    D = ToyDriver(; stim)
    R = ToyResponder()
    edges = wire_edges ?
        (connect(:D => :u, :R => :p_u), connect(:D => :v, :R => :p_v)) : ()
    return couple((Subsystem(D; name = :D), Subsystem(R; name = :R)), edges)
end

println("── sockets 1,2,7 — feedforward WIREs + edge on/off by composition ──")
let stim = ToyDriver().stim, tspan = (0.0, 5.0), opts = (dt = 0.01, adaptive = false)
    # G1 isolation: edges OFF ⇒ R's sub-trajectory reproduces standalone R exactly.
    cm_off = build_wired(wire_edges = false, stim = stim)
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
    cm_ap = build_wired(wire_edges = true, stim = onepulse)
    cm_rest = build_wired(wire_edges = true, stim = norest)
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

# --- socket 3 — DERIVED-source WIRE (D's monitor b = Ca − a drives Redox's p_b slot) ---------
println("── socket 3 — DERIVED-source WIRE (monitor → param slot) ──")
let tspan = (0.0, 5.0), opts = (dt = 0.01, adaptive = false)
    D = ToyDriver()
    cm = couple(
        (Subsystem(D; name = :D), Subsystem(ToyRedox(); name = :Redox)),
        (connect(:D => :b, :Redox => :p_b),),
    )
    # `b` is DERIVED, not integrated: wiring it adds no state, and the conservation law a + b = Ca
    # stays inside ToyDriver rather than being restated at the edge.
    @assert state_index(cm, :b) === nothing "the monitor b leaked into the global state"
    @assert num_states(cm) == 5 "expected D's 4 states + Redox's z, got $(num_states(cm))"

    # Structural: Redox's dz reads b = Ca − a live, not its standalone p_b default of 0.
    U = default_initial_state(cm)
    dU = similar(U)
    cm(dU, U, nothing, 0.0)
    Ca = D.parameters[6]
    b0 = Ca - U[state_index(cm, :a)]
    shunt, kz = ToyRedox().parameters[3], ToyRedox().parameters[4]
    dz_want = shunt * (0.0 + b0) - kz * U[state_index(cm, :Redox_z)]
    @assert approx(dU[state_index(cm, :Redox_z)], dz_want) "socket 3: dz did not read the derived b"

    # Behavioural + OFF-invariant: omit the edge and p_b holds its default, so z never moves.
    off = couple((Subsystem(ToyDriver(); name = :D), Subsystem(ToyRedox(); name = :Redox)))
    son = solve(ODEProblem(cm, tspan), Tsit5(); opts...)
    soff = solve(ODEProblem(off, tspan), Tsit5(); opts...)
    zon = son.u[end][state_index(cm, :Redox_z)]
    zoff = soff.u[end][state_index(off, :Redox_z)]
    @assert zon > 1.0e-3 "socket 3: wiring the derived b made no difference to z"
    @assert approx(zoff, 0.0; atol = 1.0e-9) "socket 3 OFF-invariant broken: z moved without the edge"
    println("  z driven by DERIVED b = Ca − a:  z(end) = ", round(zon; digits = 4),
            ";  omit the edge ⇒ z(end) = ", zoff)
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

# --- socket 6 — redox module on/off by composition (include vs. omit the ToyRedox subsystem) --
println("── socket 6 — module on/off by composition (ToyRedox subsystem) ──")
let tspan = (0.0, 5.0), opts = (dt = 0.01, adaptive = false)
    D, R = ToyDriver(), ToyResponder()
    redox_edge = connect(:D => :u, :Redox => :p_u)   # drive the redox module from D's u
    on = couple(
        (Subsystem(D; name = :D), Subsystem(R; name = :R), Subsystem(ToyRedox(); name = :Redox)),
        (redox_edge,),
    )
    off = couple((Subsystem(D; name = :D), Subsystem(R; name = :R)))
    son = solve(ODEProblem(on, tspan), Tsit5(); opts...)
    soff = solve(ODEProblem(off, tspan), Tsit5(); opts...)

    # ON: the redox module contributes a live z state driven by the wired input.
    zmax = maximum(abs(u[state_index(on, :Redox_z)]) for u in son.u)
    @assert zmax > 0 "redox composed in but z never moved"

    # OFF-invariant: omitting ToyRedox recovers D & R exactly — z is a leaf, nothing downstream reads it.
    dr, rr = on.layout.solution_indices.D, on.layout.solution_indices.R
    dr0, rr0 = off.layout.solution_indices.D, off.layout.solution_indices.R
    dmax = maximum(maximum(abs, son.u[j][dr] .- soff.u[j][dr0]) for j in eachindex(soff.t))
    rmax = maximum(maximum(abs, son.u[j][rr] .- soff.u[j][rr0]) for j in eachindex(soff.t))
    @assert approx(dmax, 0.0; atol = 1e-9) && approx(rmax, 0.0; atol = 1e-9) "redox OFF-invariant broken: D/R drifted (D $dmax, R $rmax)"
    println("  compose in ⇒ live z (max |z| = ", round(zmax; digits = 4),
            ");  omit ⇒ D,R recover baseline (max |Δ| = ", max(dmax, rmax), ")")
end

# --- socket 8 — h as a live state (compose ToyH) vs held param (omit it) ----------------------
println("── socket 8 — state↔param role flip by composition (ToyH subsystem) ──")
let tspan = (0.0, 5.0), opts = (dt = 0.01, adaptive = false)
    R = ToyResponder()
    # OFF: no ToyH ⇒ R holds h as its fixed parameter (default h0 = 1.0) — the baseline.
    soff = solve(ODEProblem(R, tspan), Tsit5(); opts...)
    # ON: compose ToyH and wire its live state into R's h slot ⇒ h is now an integrated STATE.
    on = couple((Subsystem(R; name = :R), Subsystem(ToyH(); name = :H)), (connect(:H => :h, :R => :h),))
    son = solve(ODEProblem(on, tspan), Tsit5(); opts...)

    @assert state_index(on, :H_h) !== nothing "ToyH's h did not become a global state"
    yidx = state_index(on, :y)            # R is primary here ⇒ bare :y
    # With ToyH, h ramps 0.2→1.0, so R's influx (∝ h) starts throttled: y differs from the held-h baseline.
    dy = maximum(abs(son.u[j][yidx] - soff.u[j][1]) for j in eachindex(soff.t))
    @assert dy > 1e-3 "held-vs-dynamic h produced no observable difference in y"
    println("  compose in ⇒ h is a live state (:H_h); omit ⇒ h held as param;  max |Δy| = ",
            round(dy; digits = 4))
end

# --- socket 9 — DERIVED monitors / conservation-law closures --------------------------------
println("── socket 9 — DERIVED monitors (conservation closures) ──")
let tspan = (0.0, 5.0), opts = (dt = 0.05, adaptive = false)
    cm = couple((Subsystem(ToyDriver(); name = :D), Subsystem(ToyResponder(); name = :R)))
    sol = solve(ODEProblem(cm, tspan), Tsit5(); opts...)
    h = monitor_history(sol, cm)
    @assert h.names == (:b, :R_n) "coupled monitor names wrong: $(h.names)"
    a_i, m_i = state_index(cm, :a), state_index(cm, :R_m)
    Ca, Cm = ToyDriver().parameters[6], ToyResponder().parameters[8]
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
#   The redox module's flux J = gJ·z (redox_flux_J) must be ADDED into D's w equation:
#   dw = Pw - Lw·w + J.  `share` is hard-discard (owner's equation wins entirely; the non-owner's
#   is zeroed via `frozen`) and `connect` writes into a *parameter* slot, not a *derivative* — so
#   neither can sum a non-owner flux into an owner's dU. This is THE feedback gap (ECCMitoRedox §6:
#   NADH += −V_THD, ΔΨm += +V_IMAC). Two candidate spellings to decide between:
#
#       # (a) additive share — Redox exposes a local `w_flux` state whose derivative sums into the slot
#       share(:D => :w, :Redox => :w_flux; owner = :D, op = +)
#
#       # (b) flux injection — a new edge kind adding a named derived flux into an owner's derivative
#       inject(:Redox => :J, :D => :w)
#
#   Needed in CytoZoo: an accumulate-into-owner's-dU path in `_run!` (src/coupling.jl), a
#   contributory alternative to today's `frozen`-index zeroing.

println("\nAll live demonstrations passed. See the TARGET API section for the API backlog.")
