# Coupling tests: layout, naming, sharing, owner semantics, validation, and the monolithic
# single-RHS functor. Solver-driven coupling solves live in test_scimlbase_ext.jl.

using ForwardDiff   # loads ext/ForwardDiffExt.jl — the additive-share Jacobian test needs it

# Minimal mock components implementing the participant interface.
struct _CplMockA <: CytoZoo.AbstractCellModel end
CytoZoo.num_states(::_CplMockA) = 3
CytoZoo.state_names(::_CplMockA) = (:b, :c, :d)
CytoZoo.default_initial_state(::_CplMockA) = [1.0, 2.0, 3.0]
CytoZoo.state_index(::_CplMockA, n::Symbol) = findfirst(==(n), (:b, :c, :d))
CytoZoo.transmembrane_potential_index(::_CplMockA) = 1

struct _CplMockB <: CytoZoo.AbstractCellModel end
CytoZoo.num_states(::_CplMockB) = 2
CytoZoo.state_names(::_CplMockB) = (:x, :y)
CytoZoo.default_initial_state(::_CplMockB) = [10.0, 20.0]
CytoZoo.state_index(::_CplMockB, n::Symbol) = findfirst(==(n), (:x, :y))
CytoZoo.transmembrane_potential_index(::_CplMockB) = 1

# Mock with a parameter slot :in, for connect validation. A connect receiver must expose a
# writable `parameters` vector (the monolithic plan binds it), so the mock carries one.
struct _CplMockP <: CytoZoo.AbstractCellModel
    parameters::Vector{Float64}
end
_CplMockP() = _CplMockP([0.0])
CytoZoo.num_states(::_CplMockP) = 1
CytoZoo.state_names(::_CplMockP) = (:z,)
CytoZoo.default_initial_state(::_CplMockP) = [0.0]
CytoZoo.state_index(::_CplMockP, n::Symbol) = findfirst(==(n), (:z,))
CytoZoo.transmembrane_potential_index(::_CplMockP) = 1
CytoZoo.parameter_index(::_CplMockP, n::Symbol) = n === :in ? 1 : nothing

# Functor-bearing toy models for the monolithic-RHS functor tests (the mocks above have no
# functor). _MonoA is a connect/share source; _MonoReader consumes a connect edge through its
# parameter slot; _MonoP owns a (stiff) shared state v that _MonoQ reads but does not drive.
struct _MonoA <: CytoZoo.AbstractCardiacCellModel end
CytoZoo.num_states(::_MonoA) = 2
CytoZoo.state_names(::_MonoA) = (:c, :d)
CytoZoo.default_initial_state(::_MonoA) = [5.0, 1.0]
CytoZoo.state_index(::_MonoA, n::Symbol) = findfirst(==(n), (:c, :d))
CytoZoo.transmembrane_potential_index(::_MonoA) = 2
function (::_MonoA)(du, u, p, t)
    du[1] = -0.2 * u[1]
    du[2] = -u[2]                      # d(t) = exp(-t)
    return nothing
end

struct _MonoReader <: CytoZoo.AbstractCardiacCellModel
    parameters::Vector{Float64}
end
_MonoReader() = _MonoReader([0.0])
CytoZoo.num_states(::_MonoReader) = 1
CytoZoo.state_names(::_MonoReader) = (:acc,)
CytoZoo.default_initial_state(::_MonoReader) = [0.0]
CytoZoo.state_index(::_MonoReader, n::Symbol) = findfirst(==(n), (:acc,))
CytoZoo.transmembrane_potential_index(::_MonoReader) = 1
CytoZoo.parameter_index(::_MonoReader, n::Symbol) = n === :d_ext ? 1 : nothing
(m::_MonoReader)(du, u, p, t) = (du[1] = m.parameters[1]; nothing)

const _MONO_K = 500.0
struct _MonoP <: CytoZoo.AbstractCardiacCellModel end     # owner of stiff shared v
CytoZoo.num_states(::_MonoP) = 2
CytoZoo.state_names(::_MonoP) = (:a, :v)
CytoZoo.default_initial_state(::_MonoP) = [1.0, 10.0]
CytoZoo.state_index(::_MonoP, n::Symbol) = findfirst(==(n), (:a, :v))
CytoZoo.transmembrane_potential_index(::_MonoP) = 2
function (::_MonoP)(du, u, p, t)
    du[1] = -0.1 * u[1]
    du[2] = -_MONO_K * (u[2] - u[1])  # v relaxes (stiffly) to a
    return nothing
end

struct _MonoQ <: CytoZoo.AbstractCardiacCellModel end     # non-owner reads shared v
CytoZoo.num_states(::_MonoQ) = 2
CytoZoo.state_names(::_MonoQ) = (:v, :w)
CytoZoo.default_initial_state(::_MonoQ) = [10.0, 0.0]
CytoZoo.state_index(::_MonoQ, n::Symbol) = findfirst(==(n), (:v, :w))
CytoZoo.transmembrane_potential_index(::_MonoQ) = 1
function (::_MonoQ)(du, u, p, t)
    du[1] = 12345.0                   # v — discarded (P owns the shared slot)
    du[2] = u[1] - u[2]               # w — driven by the shared v
    return nothing
end

struct _MonoR <: CytoZoo.AbstractCardiacCellModel end     # third member of a 3-way share on v
CytoZoo.num_states(::_MonoR) = 2
CytoZoo.state_names(::_MonoR) = (:v, :s)
CytoZoo.default_initial_state(::_MonoR) = [10.0, 0.0]
CytoZoo.state_index(::_MonoR, n::Symbol) = findfirst(==(n), (:v, :s))
CytoZoo.transmembrane_potential_index(::_MonoR) = 1
function (::_MonoR)(du, u, p, t)
    du[1] = 777.0                     # v — discarded (P owns the shared slot)
    du[2] = 2 * u[1]                  # s — driven by the shared v
    return nothing
end

# Contributory-share toys (`op = +`). Each writes its whole `du` view like any component; the
# difference is what the coupling does with the shared slot afterwards — added, not discarded. Two
# of them so one slot can be shown to SUM several contributions.
struct _MonoJ <: CytoZoo.AbstractCardiacCellModel end
CytoZoo.num_states(::_MonoJ) = 2
CytoZoo.state_names(::_MonoJ) = (:v, :j)
CytoZoo.default_initial_state(::_MonoJ) = [10.0, 0.0]
CytoZoo.state_index(::_MonoJ, n::Symbol) = findfirst(==(n), (:v, :j))
CytoZoo.transmembrane_potential_index(::_MonoJ) = 1
function (::_MonoJ)(du, u, p, t)
    du[1] = 20.0                      # v — CONTRIBUTED into the shared slot
    du[2] = u[1]                      # j — also reads the shared v
    return nothing
end

struct _MonoJ2 <: CytoZoo.AbstractCardiacCellModel end
CytoZoo.num_states(::_MonoJ2) = 1
CytoZoo.state_names(::_MonoJ2) = (:v,)
CytoZoo.default_initial_state(::_MonoJ2) = [10.0]
CytoZoo.state_index(::_MonoJ2, n::Symbol) = findfirst(==(n), (:v,))
CytoZoo.transmembrane_potential_index(::_MonoJ2) = 1
(::_MonoJ2)(du, u, p, t) = (du[1] = 300.0; nothing)   # v — a second contribution into the same slot

# A contribution that DEPENDS on the shared state. This is what distinguishes an additive share
# from a `connect`: the term stays in the Jacobian instead of being frozen to its primal.
struct _MonoJState <: CytoZoo.AbstractCardiacCellModel end
CytoZoo.num_states(::_MonoJState) = 1
CytoZoo.state_names(::_MonoJState) = (:v,)
CytoZoo.default_initial_state(::_MonoJState) = [10.0]
CytoZoo.state_index(::_MonoJState, n::Symbol) = findfirst(==(n), (:v,))
CytoZoo.transmembrane_potential_index(::_MonoJState) = 1
(::_MonoJState)(du, u, p, t) = (du[1] = -3.0 * u[1]; nothing)   # contributes -3v into the shared slot

# Mutual-contribution toys: each owns one slot and contributes to the other's. Under hard-discard
# this topology is a cyclic operator order; an accumulating class imposes no ordering, so it is
# legal here.
struct _MonoMutX <: CytoZoo.AbstractCardiacCellModel end
CytoZoo.num_states(::_MonoMutX) = 2
CytoZoo.state_names(::_MonoMutX) = (:f, :g)
CytoZoo.default_initial_state(::_MonoMutX) = [0.0, 0.0]
CytoZoo.state_index(::_MonoMutX, n::Symbol) = findfirst(==(n), (:f, :g))
CytoZoo.transmembrane_potential_index(::_MonoMutX) = 1
(::_MonoMutX)(du, u, p, t) = (du[1] = 1.0; du[2] = 3.0; nothing)

struct _MonoMutY <: CytoZoo.AbstractCardiacCellModel end
CytoZoo.num_states(::_MonoMutY) = 2
CytoZoo.state_names(::_MonoMutY) = (:f, :g)
CytoZoo.default_initial_state(::_MonoMutY) = [0.0, 0.0]
CytoZoo.state_index(::_MonoMutY, n::Symbol) = findfirst(==(n), (:f, :g))
CytoZoo.transmembrane_potential_index(::_MonoMutY) = 1
(::_MonoMutY)(du, u, p, t) = (du[1] = 10.0; du[2] = 30.0; nothing)

# Monitor-source toy: owns a conserved pool a + b = _MONO_C but integrates only `a`, surfacing
# the complement `b` as a DERIVED monitor. This is the shape a real monitor source has (Gauthier's
# ADPi = C_A - ATPi). The second monitor `ag` depends on TWO states, which no per-edge scalar
# transform could express. `monitor_calls` counts evaluations so the pre-pass can be shown to run
# once per eval rather than once per receiving edge.
const _MONO_C = 8.0
mutable struct _MonoDerived <: CytoZoo.AbstractCardiacCellModel
    monitor_calls::Int
end
_MonoDerived() = _MonoDerived(0)
CytoZoo.num_states(::_MonoDerived) = 2
CytoZoo.state_names(::_MonoDerived) = (:a, :g)
CytoZoo.default_initial_state(::_MonoDerived) = [3.0, 1.0]
CytoZoo.state_index(::_MonoDerived, n::Symbol) = findfirst(==(n), (:a, :g))
CytoZoo.transmembrane_potential_index(::_MonoDerived) = 1
CytoZoo.num_monitors(::_MonoDerived) = 2
CytoZoo.monitor_names(::_MonoDerived) = (:b, :ag)
function CytoZoo.monitor_values!(mon, u, t, m::_MonoDerived)
    m.monitor_calls += 1
    mon[1] = _MONO_C - u[1]            # conserved-pool complement
    mon[2] = u[1] + u[2]               # depends on two states
    return nothing
end
function (::_MonoDerived)(du, u, p, t)
    du[1] = -u[1]                      # a(t) = 3exp(-t)
    du[2] = 0.0
    return nothing
end

# Declares :s as BOTH a state and a monitor, so a connect source naming it is ambiguous.
struct _AmbiguousName <: CytoZoo.AbstractCardiacCellModel end
CytoZoo.num_states(::_AmbiguousName) = 1
CytoZoo.state_names(::_AmbiguousName) = (:s,)
CytoZoo.default_initial_state(::_AmbiguousName) = [1.0]
CytoZoo.state_index(::_AmbiguousName, n::Symbol) = findfirst(==(n), (:s,))
CytoZoo.transmembrane_potential_index(::_AmbiguousName) = 1
CytoZoo.num_monitors(::_AmbiguousName) = 1
CytoZoo.monitor_names(::_AmbiguousName) = (:s,)
CytoZoo.monitor_values!(mon, u, t, ::_AmbiguousName) = (mon[1] = 0.0; nothing)
(::_AmbiguousName)(du, u, p, t) = (du[1] = 0.0; nothing)

# Sources a monitor AND exposes a receiver slot whose staged value its monitor reads — the
# one-eval lag hazard `couple` rejects, since the pre-pass runs before any staging.
struct _DerivedReceiver <: CytoZoo.AbstractCardiacCellModel
    parameters::Vector{Float64}
end
_DerivedReceiver() = _DerivedReceiver([0.0])
CytoZoo.num_states(::_DerivedReceiver) = 1
CytoZoo.state_names(::_DerivedReceiver) = (:a,)
CytoZoo.default_initial_state(::_DerivedReceiver) = [1.0]
CytoZoo.state_index(::_DerivedReceiver, n::Symbol) = findfirst(==(n), (:a,))
CytoZoo.transmembrane_potential_index(::_DerivedReceiver) = 1
CytoZoo.parameter_index(::_DerivedReceiver, n::Symbol) = n === :in ? 1 : nothing
CytoZoo.num_monitors(::_DerivedReceiver) = 1
CytoZoo.monitor_names(::_DerivedReceiver) = (:b,)
CytoZoo.monitor_values!(mon, u, t, m::_DerivedReceiver) = (mon[1] = m.parameters[1] - u[1]; nothing)
(m::_DerivedReceiver)(du, u, p, t) = (du[1] = m.parameters[1]; nothing)

# Reports more state names than its initial-state vector has entries — `couple` must reject
# this at construction rather than emitting a BoundsError from inside the layout loop.
struct _BadLength <: CytoZoo.AbstractCardiacCellModel end
CytoZoo.num_states(::_BadLength) = 3
CytoZoo.state_names(::_BadLength) = (:m, :n, :o)
CytoZoo.default_initial_state(::_BadLength) = [1.0, 2.0]
CytoZoo.state_index(::_BadLength, n::Symbol) = findfirst(==(n), (:m, :n, :o))
CytoZoo.transmembrane_potential_index(::_BadLength) = 1

# Non-contiguous-block toy models: sharing B's *last* state with A's *first* scatters B's
# global indices (x appended after A, y folded onto A's early slot) → B.solution_indices = [3,1],
# which forces the `copy(idxs)` / Vector-indexed `view` branch in the plan.
struct _ScatterA <: CytoZoo.AbstractCardiacCellModel end
CytoZoo.num_states(::_ScatterA) = 2
CytoZoo.state_names(::_ScatterA) = (:p, :q)
CytoZoo.default_initial_state(::_ScatterA) = [2.0, 3.0]
CytoZoo.state_index(::_ScatterA, n::Symbol) = findfirst(==(n), (:p, :q))
CytoZoo.transmembrane_potential_index(::_ScatterA) = 1
function (::_ScatterA)(du, u, p, t)
    du[1] = -u[1]                     # p — owns the shared slot
    du[2] = 2u[2]                     # q
    return nothing
end

struct _ScatterB <: CytoZoo.AbstractCardiacCellModel end
CytoZoo.num_states(::_ScatterB) = 2
CytoZoo.state_names(::_ScatterB) = (:x, :y)
CytoZoo.default_initial_state(::_ScatterB) = [7.0, 9.0]
CytoZoo.state_index(::_ScatterB, n::Symbol) = findfirst(==(n), (:x, :y))
CytoZoo.transmembrane_potential_index(::_ScatterB) = 1
function (::_ScatterB)(du, u, p, t)
    du[1] = 5.0                       # x
    du[2] = 999.0                     # y — discarded (A owns the shared slot)
    return nothing
end

# Receiver that implements parameter_index but exposes no writable_parameters (no `parameters`
# field, no override) — for the connect-participation contract error.
struct _NoParamsField <: CytoZoo.AbstractCardiacCellModel end
CytoZoo.num_states(::_NoParamsField) = 1
CytoZoo.state_names(::_NoParamsField) = (:z,)
CytoZoo.default_initial_state(::_NoParamsField) = [0.0]
CytoZoo.state_index(::_NoParamsField, n::Symbol) = findfirst(==(n), (:z,))
CytoZoo.transmembrane_potential_index(::_NoParamsField) = 1
CytoZoo.parameter_index(::_NoParamsField, n::Symbol) = n === :in ? 1 : nothing

@testset "couple — single node" begin
    cm = couple([Subsystem(_CplMockA(); name = :A)])
    @test cm isa CoupledModel
    @test num_states(cm) == 3
    @test state_names(cm) == (:b, :c, :d)
    @test cm.layout.solution_indices.A == [1, 2, 3]
    @test default_initial_state(cm) == [1.0, 2.0, 3.0]
    @test cm.layout.operator_order == [:A]
    @test transmembrane_potential_index(cm) == 1
end

@testset "couple — gensym default name (single node)" begin
    cm = couple([Subsystem(_CplMockA())])      # name auto-generated; fine for a single node
    @test num_states(cm) == 3
    @test state_names(cm) == (:b, :c, :d)
    @test length(cm.components) == 1
end

@testset "couple — share, owner = first component" begin
    cm = couple(
        [Subsystem(_CplMockA(); name = :A), Subsystem(_CplMockB(); name = :B)],
        [share(:A => :d, :B => :x; owner = :A)],
    )
    @test num_states(cm) == 4
    @test state_names(cm) == (:b, :c, :d, :B_y)          # owner A keeps :d; B's y prefixed
    @test cm.layout.solution_indices.A == [1, 2, 3]
    @test cm.layout.solution_indices.B == [3, 4]         # x shares slot 3, y is slot 4
    @test default_initial_state(cm) == [1.0, 2.0, 3.0, 20.0]  # slot 3 = owner A's d IC
    @test cm.layout.operator_order == [:B, :A]           # owner A steps last, wins the slot
    @test state_index(cm, :d) == 3
    @test state_index(cm, :B_y) == 4
    @test state_index(cm, :nope) === nothing
    @test transmembrane_potential_index(cm) == 1
end

@testset "couple — share, owner = second component" begin
    cm = couple(
        [Subsystem(_CplMockA(); name = :A), Subsystem(_CplMockB(); name = :B)],
        [share(:A => :d, :B => :x; owner = :B)],
    )
    @test state_names(cm) == (:b, :c, :x, :B_y)          # canonical name from owner B
    @test cm.layout.solution_indices.A == [1, 2, 3]
    @test cm.layout.solution_indices.B == [3, 4]
    @test default_initial_state(cm) == [1.0, 2.0, 10.0, 20.0]  # slot 3 = owner B's x IC
    @test cm.layout.operator_order == [:A, :B]           # owner B steps last
end

@testset "couple — explicit canonical share name" begin
    cm = couple(
        [Subsystem(_CplMockA(); name = :A), Subsystem(_CplMockB(); name = :B)],
        [share(:A => :d, :B => :x; owner = :A, name = :shared)],
    )
    @test state_names(cm) == (:b, :c, :shared, :B_y)
    @test state_index(cm, :shared) == 3
end

@testset "couple — validation and parameter contract" begin
    @test_throws ArgumentError share(:A => :d, :B => :x; owner = :C)                       # owner not a participant
    @test_throws ArgumentError couple(
        [Subsystem(_CplMockA(); name = :A), Subsystem(_CplMockB(); name = :B)],
        [share(:A => :nope, :B => :x; owner = :A)]
    )                                        # unknown state
    @test_throws ArgumentError couple(
        [Subsystem(_CplMockA(); name = :A), Subsystem(_CplMockB(); name = :B)],
        [share(:Z => :d, :B => :x; owner = :Z)]
    )                                           # unknown component
    @test_throws ArgumentError couple(
        [
            Subsystem(_CplMockA(); name = :A),
            Subsystem(_CplMockB(); name = :A),
        ]
    )                                                # duplicate node name
    @test_throws ArgumentError couple([Subsystem(_CplMockA(); name = :A)], [:not_an_edge]) # unknown edge type
    cm = couple([Subsystem(_CplMockA(); name = :A)])
    @test_throws MethodError num_parameters(cm)                                            # no single parameter vector
    @test !hasmethod(num_parameters, Tuple{typeof(cm)})                                    # and does not claim to have one
end

@testset "couple — connect validation" begin
    cm = couple(
        [Subsystem(_CplMockA(); name = :A), Subsystem(_CplMockP(); name = :P)],
        [connect(:A => :b, :P => :in)],
    )
    @test cm.connects[1] isa CytoZoo.ConnectSpec
    @test connect(:A => :b, :P => :in; op = +) isa CytoZoo.ConnectSpec                     # + (sum) is supported
    @test_throws ArgumentError connect(:A => :b, :P => :in; op = max)                      # unsupported op
    @test_throws ArgumentError couple(
        [Subsystem(_CplMockA(); name = :A), Subsystem(_CplMockP(); name = :P)],
        [connect(:A => :b, :P => :nope)]
    )                                                  # no such param slot
    @test_throws ArgumentError couple(
        [Subsystem(_CplMockA(); name = :A), Subsystem(_CplMockP(); name = :P)],
        [connect(:A => :nostate, :P => :in)]
    )                                              # no such source state
end

@testset "couple — monolithic functor" begin
    # connect: receiver reads a source state live from U through its parameter slot.
    cm = couple(
        [Subsystem(_MonoA(); name = :A), Subsystem(_MonoReader(); name = :R)],
        [connect(:A => :d, :R => :d_ext)],
    )
    U = default_initial_state(cm)                       # [c=5, d=1, acc=0]
    dU = similar(U)
    cm(dU, U, nothing, 0.0)
    @test dU[state_index(cm, :c)] == -0.2 * 5.0         # A's own equations untouched
    @test dU[state_index(cm, :d)] == -1.0
    @test dU[state_index(cm, :R_acc)] == 1.0            # acc' = d_ext = d(0) = 1

    # connect is read live from U each eval (not a stale parameter): change U, derivative tracks.
    U2 = copy(U); U2[state_index(cm, :d)] = 7.0
    cm(dU, U2, nothing, 0.0)
    @test dU[state_index(cm, :R_acc)] == 7.0

    # zero allocation on the assembled RHS (contiguous blocks).
    cm(dU, U, nothing, 0.0)                             # warm up
    @test (@allocated cm(dU, U, nothing, 0.0)) == 0

    # connect op = +: two edges sum into one slot (reset to zero then accumulated).
    cm_sum = couple(
        [Subsystem(_MonoA(); name = :A), Subsystem(_MonoReader(); name = :R)],
        [connect(:A => :c, :R => :d_ext; op = +), connect(:A => :d, :R => :d_ext; op = +)],
    )
    Us = default_initial_state(cm_sum)
    dUs = similar(Us)
    cm_sum(dUs, Us, nothing, 0.0)
    @test dUs[state_index(cm_sum, :R_acc)] == 5.0 + 1.0          # c + d
    # Cross-sectional sum, NOT a running total: the slot is reset to zero before each eval, so
    # repeated evaluation gives the same answer. Without the reset this drifts to 12, 18, ...
    cm_sum(dUs, Us, nothing, 0.0)
    cm_sum(dUs, Us, nothing, 0.0)
    @test dUs[state_index(cm_sum, :R_acc)] == 5.0 + 1.0

    # share owner-last: only the owner's derivative reaches the shared slot.
    cm_sh = couple(
        [Subsystem(_MonoP(); name = :P), Subsystem(_MonoQ(); name = :Q)],
        [share(:P => :v, :Q => :v; owner = :P)],
    )
    @test cm_sh.layout.operator_order == [:Q, :P]               # owner P steps last
    Ush = default_initial_state(cm_sh)                          # [a=1, v=10, Q_w=0]
    dUsh = similar(Ush)
    cm_sh(dUsh, Ush, nothing, 0.0)
    @test dUsh[state_index(cm_sh, :a)] == -0.1 * 1.0
    @test dUsh[state_index(cm_sh, :v)] == -_MONO_K * (10.0 - 1.0)  # P's value, NOT Q's 12345
    @test dUsh[state_index(cm_sh, :Q_w)] == 10.0 - 0.0            # Q's w driven by shared v
    cm_sh(dUsh, Ush, nothing, 0.0)                               # warm up
    @test (@allocated cm_sh(dUsh, Ush, nothing, 0.0)) == 0        # share path (frozen zeroing) allocation-free
end

# Shares are resolved as equivalence classes, not pairwise, so a state may take part in several
# share edges and the whole group still collapses to ONE global slot. Before this was union-find,
# a second edge touching an already-shared state silently split the group and left the surviving
# slot with an identically-zero derivative.
@testset "couple — multi-way share topologies" begin
    # Fan-out from the owner: P.v ≡ Q.v ≡ R.v is a single slot governed by P.
    nodes = [Subsystem(_MonoP(); name = :P), Subsystem(_MonoQ(); name = :Q), Subsystem(_MonoR(); name = :R)]
    edges = [share(:P => :v, :Q => :v; owner = :P), share(:P => :v, :R => :v; owner = :P)]
    cm = couple(nodes, edges)
    @test num_states(cm) == 4                                   # a, v, Q_w, R_s — one shared v
    @test state_names(cm) == (:a, :v, :Q_w, :R_s)
    vslot = state_index(cm, :v)
    @test cm.layout.solution_indices.Q[1] == vslot               # every member maps to the same slot
    @test cm.layout.solution_indices.R[1] == vslot
    @test cm.layout.operator_order[end] == :P                    # owner steps last

    U = default_initial_state(cm)                                # [a=1, v=10, Q_w=0, R_s=0]
    @test U[vslot] == 10.0                                       # IC comes from the owner
    dU = similar(U)
    cm(dU, U, nothing, 0.0)
    @test dU[vslot] == -_MONO_K * (10.0 - 1.0)                   # P governs; 12345 and 777 discarded
    @test dU[state_index(cm, :Q_w)] == 10.0 - 0.0                # non-owners still READ the shared v
    @test dU[state_index(cm, :R_s)] == 2 * 10.0
    cm(dU, U, nothing, 0.0)                                      # warm up
    @test (@allocated cm(dU, U, nothing, 0.0)) == 0              # multi-way share stays allocation-free

    # Declaration order of the nodes must not change the physics: the owner may be listed last.
    cm_reordered = couple(
        [Subsystem(_MonoQ(); name = :Q), Subsystem(_MonoR(); name = :R), Subsystem(_MonoP(); name = :P)],
        edges,
    )
    @test num_states(cm_reordered) == 4
    @test state_names(cm_reordered) == (:v, :w, :R_s, :P_a)      # shared slot keeps the owner's name
    Ur = default_initial_state(cm_reordered)
    dUr = similar(Ur)
    cm_reordered(dUr, Ur, nothing, 0.0)
    @test dUr[state_index(cm_reordered, :v)] == -_MONO_K * (10.0 - 1.0)   # still P's equation
    @test cm_reordered.layout.operator_order[end] == :P

    # A chained share names one group through two edges, so the two `owner=`s contradict each
    # other. That must be an actionable error, never a silently split layout.
    @test_throws ArgumentError couple(
        nodes, [share(:P => :v, :Q => :v; owner = :P), share(:Q => :v, :R => :v; owner = :Q)]
    )

    # A component whose state_names and default_initial_state disagree is caught at couple time.
    @test_throws ArgumentError couple([Subsystem(_BadLength(); name = :X)])
end

@testset "couple — non-contiguous block" begin
    cm = couple(
        [Subsystem(_ScatterA(); name = :A), Subsystem(_ScatterB(); name = :B)],
        [share(:A => :p, :B => :y; owner = :A)],
    )
    @test cm.layout.solution_indices.B == [3, 1]        # x→3 (appended), y→1 (shared) — non-contiguous
    U = default_initial_state(cm)
    @test U == [2.0, 3.0, 7.0]                          # slot1=A.p IC, slot2=A.q, slot3=B.x
    dU = similar(U)
    cm(dU, U, nothing, 0.0)
    @test dU[1] == -2.0                                 # A.p owns shared slot 1 (B's 999 discarded)
    @test dU[2] == 2 * 3.0                              # A.q
    @test dU[3] == 5.0                                  # B.x, written through the scattered view

    cm(dU, U, nothing, 0.0)                             # warm up
    @test (@allocated cm(dU, U, nothing, 0.0)) == 0     # scattered (Vector-indexed) view stays allocation-free
end

@testset "couple — contributory share (op = +)" begin
    nodes = [Subsystem(_MonoP(); name = :P), Subsystem(_MonoJ(); name = :J)]
    cm = couple(nodes, [share(:P => :v, :J => :v; owner = :P, op = +)])

    vslot = state_index(cm, :v)
    @test num_states(cm) == 3                                    # a, v, J_j — still one shared slot
    @test state_names(cm) == (:a, :v, :J_j)
    @test cm.layout.solution_indices.J[1] == vslot                # the contributor aliases the slot
    @test default_initial_state(cm)[vslot] == 10.0                # IC still comes from the OWNER

    # The contributor is accumulated, never frozen; the plan carries the slot for the pre-walk reset.
    entry = cm.plan[findfirst(e -> e.model isa _MonoJ, cm.plan)]
    @test entry.frozen == ()
    @test entry.accumulated == (1,)
    @test cm.additive_slots == (vslot,)

    # dU arrives dirty from the solver — the accumulate path must not inherit whatever was there.
    U = default_initial_state(cm)
    dU = similar(U)
    fill!(dU, 1.0e6)
    cm(dU, U, nothing, 0.0)
    @test dU[vslot] == -_MONO_K * (10.0 - 1.0) + 20.0             # owner's equation PLUS the contribution
    @test dU[state_index(cm, :J_j)] == 10.0                       # the contributor still reads the shared v

    # Cross-sectional, NOT a running total: every accumulating slot is reset before the walk, so
    # repeated evaluation gives the same answer. Without the reset this drifts by +20 each call.
    cm(dU, U, nothing, 0.0)
    cm(dU, U, nothing, 0.0)
    @test dU[vslot] == -_MONO_K * (10.0 - 1.0) + 20.0

    cm(dU, U, nothing, 0.0)                                       # warm up
    @test (@allocated cm(dU, U, nothing, 0.0)) == 0               # save/restore stays allocation-free
end

@testset "couple — contributory share sums, and mixes with hard-discard" begin
    # One class, three roles at once: P owns, Q is hard-discarded, J and J2 both contribute.
    nodes = [
        Subsystem(_MonoP(); name = :P), Subsystem(_MonoQ(); name = :Q),
        Subsystem(_MonoJ(); name = :J), Subsystem(_MonoJ2(); name = :J2),
    ]
    edges = [
        share(:P => :v, :Q => :v; owner = :P),                    # hard discard — 12345 dropped
        share(:P => :v, :J => :v; owner = :P, op = +),            # +20
        share(:P => :v, :J2 => :v; owner = :P, op = +),           # +300
    ]
    cm = couple(nodes, edges)
    vslot = state_index(cm, :v)
    want = -_MONO_K * (10.0 - 1.0) + 20.0 + 300.0

    U = default_initial_state(cm)
    dU = similar(U)
    fill!(dU, -7.0)
    cm(dU, U, nothing, 0.0)
    @test dU[vslot] == want                                       # contributions sum; the poison is gone
    @test dU[state_index(cm, :Q_w)] == 10.0                       # the frozen member still READS the slot

    # An accumulating class is order-free, so every declaration order must give the same physics.
    # This is the property that lets a mixed class exist at all: hard-discard wants the owner last,
    # accumulation wants it first, and save/restore removes the question.
    for perm in ([2, 1, 4, 3], [4, 3, 2, 1], [3, 4, 1, 2], [2, 4, 1, 3])
        cmp_ = couple(nodes[perm], edges)
        Up = default_initial_state(cmp_)
        dUp = similar(Up)
        fill!(dUp, 1.0e6)
        cmp_(dUp, Up, nothing, 0.0)
        @test dUp[state_index(cmp_, :v)] == want
    end

    cm(dU, U, nothing, 0.0)                                       # warm up
    @test (@allocated cm(dU, U, nothing, 0.0)) == 0
end

@testset "couple — contributory share on a non-contiguous block" begin
    # Same scattered topology as above (B.solution_indices == [3, 1]), but B now CONTRIBUTES its
    # 999 instead of having it discarded. `block[i]` compiles differently for a Vector-indexed
    # block than for a UnitRange, so the save/restore has to be allocation-free on both.
    cm = couple(
        [Subsystem(_ScatterA(); name = :A), Subsystem(_ScatterB(); name = :B)],
        [share(:A => :p, :B => :y; owner = :A, op = +)],
    )
    @test cm.layout.solution_indices.B == [3, 1]
    @test cm.additive_slots == (1,)

    U = default_initial_state(cm)
    dU = similar(U)
    fill!(dU, 1.0e6)
    cm(dU, U, nothing, 0.0)
    @test dU[1] == -2.0 + 999.0                       # A.p's equation plus B.y's contribution
    @test dU[2] == 2 * 3.0                            # A.q untouched
    @test dU[3] == 5.0                                # B.x, written through the scattered view

    cm(dU, U, nothing, 0.0)                           # warm up
    @test (@allocated cm(dU, U, nothing, 0.0)) == 0
end

@testset "couple — contributory shares may be mutual" begin
    # X owns f and contributes to g; Y owns g and contributes to f. Under hard-discard the two
    # precedence edges (X before Y, Y before X) are a cycle; an accumulating class emits none.
    cm = couple(
        [Subsystem(_MonoMutX(); name = :X), Subsystem(_MonoMutY(); name = :Y)],
        [share(:X => :f, :Y => :f; owner = :X, op = +),
         share(:Y => :g, :X => :g; owner = :Y, op = +)],
    )
    U = default_initial_state(cm)
    dU = similar(U)
    fill!(dU, 1.0e6)
    cm(dU, U, nothing, 0.0)
    @test dU[state_index(cm, :f)] == 1.0 + 10.0       # X's equation + Y's contribution
    @test dU[state_index(cm, :g)] == 30.0 + 3.0       # Y's equation + X's contribution

    # The same pair under hard-discard is still, correctly, a cyclic operator order.
    @test_throws ArgumentError couple(
        [Subsystem(_MonoMutX(); name = :X), Subsystem(_MonoMutY(); name = :Y)],
        [share(:X => :f, :Y => :f; owner = :X), share(:Y => :g, :X => :g; owner = :Y)],
    )
end

@testset "couple — contributory share validation" begin
    nodes = [Subsystem(_MonoP(); name = :P), Subsystem(_MonoJ(); name = :J), Subsystem(_MonoJ2(); name = :J2)]

    # R1 — an unsupported op is rejected where it is written, not at couple time.
    @test_throws ArgumentError share(:P => :v, :J => :v; owner = :P, op = *)
    @test share(:P => :v, :J => :v; owner = :P, op = +).op === (+)          # ...but `+` is fine
    @test share(:P => :v, :J => :v; owner = :P).op === overwrite            # ...and `overwrite` is the default

    # R3 — one endpoint cannot both contribute and be discarded.
    @test_throws ArgumentError couple(
        nodes,
        [share(:P => :v, :J => :v; owner = :P, op = +), share(:P => :v, :J => :v; owner = :P)],
    )

    # R4 — one component may contribute at most one state to a class. Two of its states in one
    # class collapse onto a single slot, which double-counts under `op = +` and silently freezes
    # one of them under hard-discard.
    @test_throws ArgumentError couple(
        [Subsystem(_MonoP(); name = :P), Subsystem(_MonoQ(); name = :Q)],
        [share(:P => :v, :Q => :v; owner = :P), share(:P => :v, :Q => :w; owner = :P)],
    )

    # R2 — `op` does not relax the one-owner rule: chaining still names two owners for one class.
    @test_throws ArgumentError couple(
        nodes,
        [share(:P => :v, :J => :v; owner = :P, op = +), share(:J => :v, :J2 => :v; owner = :J, op = +)],
    )
end

@testset "couple — a contributory share survives ForwardDiff" begin
    # The sharp difference between an additive share and a `connect`. A `connect` input is frozen
    # to its primal (ext/ForwardDiffExt.jl), so a derivative taken through one silently loses the
    # coupling term. A share never leaves `U`, and the accumulate path holds its saved values in
    # `eltype(dU)` — a Dual stays a Dual — so a STATE-DEPENDENT contribution stays in the Jacobian.
    cm = couple(
        [Subsystem(_MonoP(); name = :P), Subsystem(_MonoJState(); name = :J)],
        [share(:P => :v, :J => :v; owner = :P, op = +)],
    )
    U = default_initial_state(cm)
    J = ForwardDiff.jacobian(u -> (d = similar(u); cm(d, u, nothing, 0.0); d), U)

    a, v = state_index(cm, :a), state_index(cm, :v)
    # dv/dt = -K(v - a) + (-3v)  =>  ∂/∂v = -K - 3, and the "-3" is the contributed term.
    @test J[v, v] == -_MONO_K - 3.0
    @test J[v, a] == _MONO_K
    @test J[a, a] == -0.1
    @test all(isfinite, J)
end

@testset "couple — hard-discard couplings are structurally inert" begin
    # The accumulate path must be invisible unless some share asks for it: with no `op = +` edge
    # every entry's `accumulated` is an empty tuple, so `_save_accumulated` / `_add_back!` /
    # `_reset_accumulators!` resolve to their `::Tuple{}` methods and compile away. Asserting the
    # structure (not just the numbers) is what makes "the numbers cannot move" a property rather
    # than an observation.
    couplings = [
        couple([Subsystem(_MonoP(); name = :P), Subsystem(_MonoQ(); name = :Q)],
               [share(:P => :v, :Q => :v; owner = :P)]),
        couple([Subsystem(_ScatterA(); name = :A), Subsystem(_ScatterB(); name = :B)],
               [share(:A => :p, :B => :y; owner = :A)]),
        couple([Subsystem(_MonoA(); name = :A), Subsystem(_MonoReader(); name = :R)],
               [connect(:A => :d, :R => :d_ext)]),
    ]
    for cm in couplings
        @test cm.additive_slots === ()
        @test all(e -> e.accumulated === (), cm.plan)
    end
end

@testset "couple — connect does not mutate the caller's model" begin
    reader = _MonoReader()
    reader.parameters[1] = 42.0                         # sentinel the coupling must not clobber
    cm = couple(
        [Subsystem(_MonoA(); name = :A), Subsystem(reader; name = :R)],
        [connect(:A => :d, :R => :d_ext)],
    )
    U = default_initial_state(cm)
    dU = similar(U)
    cm(dU, U, nothing, 0.0)                             # stages d(0)=1 into the deepcopied receiver's slot
    @test dU[state_index(cm, :R_acc)] == 1.0            # coupling still works on the copy
    @test reader.parameters[1] == 42.0                  # caller's original model untouched
end

@testset "couple — connect op conflicts rejected" begin
    @test_throws ArgumentError couple(                  # two overwrites into one slot (last-wins)
        [Subsystem(_CplMockA(); name = :A), Subsystem(_CplMockP(); name = :P)],
        [connect(:A => :b, :P => :in), connect(:A => :c, :P => :in)],
    )
    @test_throws ArgumentError couple(                  # mixed overwrite/+ into one slot
        [Subsystem(_CplMockA(); name = :A), Subsystem(_CplMockP(); name = :P)],
        [connect(:A => :b, :P => :in), connect(:A => :c, :P => :in; op = +)],
    )
end

# A connect source may name a STATE or a MONITOR. A monitor source lets a quantity an algebraic
# law determines (a conserved pool's complement) be wired without promoting it to an integrated
# state — the law stays inside the model that owns it.
@testset "couple — monitor-sourced connect edges" begin
    src = _MonoDerived()
    cm = couple(
        [Subsystem(src; name = :D), Subsystem(_MonoReader(); name = :R)],
        [connect(:D => :b, :R => :d_ext)],
    )
    @test num_states(cm) == 3
    @test state_names(cm) == (:a, :g, :R_acc)      # the monitor did NOT become a state
    @test state_index(cm, :b) === nothing

    U = default_initial_state(cm)                   # [a=3, g=1, acc=0]
    dU = similar(U)
    cm(dU, U, nothing, 0.0)
    @test dU[state_index(cm, :R_acc)] == _MONO_C - 3.0     # acc' = b = C - a
    @test dU[state_index(cm, :a)] == -3.0                  # source's own equations untouched

    # Recomputed live from U each eval, exactly like a state source.
    U2 = copy(U)
    U2[state_index(cm, :a)] = 5.0
    cm(dU, U2, nothing, 0.0)
    @test dU[state_index(cm, :R_acc)] == _MONO_C - 5.0

    cm(dU, U, nothing, 0.0)                                # warm up
    @test (@allocated cm(dU, U, nothing, 0.0)) == 0        # pre-pass + scratch stay allocation-free

    # A monitor over several states — the case no per-edge scalar transform could express.
    cm_multi = couple(
        [Subsystem(_MonoDerived(); name = :D), Subsystem(_MonoReader(); name = :R)],
        [connect(:D => :ag, :R => :d_ext)],
    )
    Um = default_initial_state(cm_multi)
    dUm = similar(Um)
    cm_multi(dUm, Um, nothing, 0.0)
    @test dUm[state_index(cm_multi, :R_acc)] == 3.0 + 1.0

    # State and monitor sources sum into one slot, and the sum stays cross-sectional (the slot is
    # reset each eval) — the same invariant the all-state `+` path holds.
    cm_sum = couple(
        [Subsystem(_MonoDerived(); name = :D), Subsystem(_MonoReader(); name = :R)],
        [connect(:D => :a, :R => :d_ext; op = +), connect(:D => :b, :R => :d_ext; op = +)],
    )
    Us = default_initial_state(cm_sum)
    dUs = similar(Us)
    cm_sum(dUs, Us, nothing, 0.0)
    @test dUs[state_index(cm_sum, :R_acc)] == 3.0 + (_MONO_C - 3.0)
    cm_sum(dUs, Us, nothing, 0.0)
    cm_sum(dUs, Us, nothing, 0.0)
    @test dUs[state_index(cm_sum, :R_acc)] == 3.0 + (_MONO_C - 3.0)
    @test (@allocated cm_sum(dUs, Us, nothing, 0.0)) == 0

    # The pre-pass evaluates each source ONCE per eval, not once per receiving edge.
    counted = _MonoDerived()
    cm_two = couple(
        [
            Subsystem(counted; name = :D),
            Subsystem(_MonoReader(); name = :R1),
            Subsystem(_MonoReader(); name = :R2),
        ],
        [connect(:D => :b, :R1 => :d_ext), connect(:D => :ag, :R2 => :d_ext)],
    )
    Ut = default_initial_state(cm_two)
    dUt = similar(Ut)
    counted.monitor_calls = 0
    cm_two(dUt, Ut, nothing, 0.0)
    @test counted.monitor_calls == 1
    @test dUt[state_index(cm_two, :R1_acc)] == _MONO_C - 3.0
    @test dUt[state_index(cm_two, :R2_acc)] == 3.0 + 1.0

    # A coupling with no monitor source carries no pre-pass at all, so the existing path is
    # untouched (the empty-tuple method compiles away).
    plain = couple(
        [Subsystem(_MonoA(); name = :A), Subsystem(_MonoReader(); name = :R)],
        [connect(:A => :d, :R => :d_ext)],
    )
    @test plain.monitor_plan === ()
    @test isempty(plain.monitor_scratch)
end

@testset "couple — monitor-source validation" begin
    # Unknown name: neither a state nor a monitor.
    @test_throws ArgumentError couple(
        [Subsystem(_MonoDerived(); name = :D), Subsystem(_MonoReader(); name = :R)],
        [connect(:D => :nope, :R => :d_ext)],
    )
    # A monitor has no derivative to own, so it cannot be a `share` endpoint.
    @test_throws ArgumentError couple(
        [Subsystem(_MonoDerived(); name = :D), Subsystem(_MonoQ(); name = :Q)],
        [share(:D => :b, :Q => :w; owner = :D)],
    )
    # A name declared as both a state and a monitor is ambiguous as a connect source.
    @test_throws ArgumentError couple(
        [Subsystem(_AmbiguousName(); name = :X), Subsystem(_MonoReader(); name = :R)],
        [connect(:X => :s, :R => :d_ext)],
    )
    # A monitor source that also receives an edge would compute its monitors from last eval's
    # staged parameters — rejected rather than silently lagged.
    @test_throws ArgumentError couple(
        [Subsystem(_DerivedReceiver(); name = :D), Subsystem(_MonoReader(); name = :R)],
        [connect(:D => :b, :R => :d_ext), connect(:R => :acc, :D => :in)],
    )
    # ...but receiving an edge is fine when the component sources no monitor.
    @test couple(
        [Subsystem(_DerivedReceiver(); name = :D), Subsystem(_MonoReader(); name = :R)],
        [connect(:R => :acc, :D => :in)],
    ) isa CytoZoo.CoupledModel
end

@testset "couple — self-share rejected" begin
    @test_throws ArgumentError couple(
        [Subsystem(_CplMockA(); name = :A)],
        [share(:A => :b, :A => :c; owner = :A)],
    )
end

@testset "couple — connect-participation contract" begin
    @test_throws ArgumentError couple(                  # receiver implements no parameter_index
        [Subsystem(_MonoA(); name = :S), Subsystem(_CplMockA(); name = :A)],
        [connect(:S => :c, :A => :in)],
    )
    @test_throws ArgumentError couple(                  # receiver has parameter_index but no writable_parameters
        [Subsystem(_MonoA(); name = :S), Subsystem(_NoParamsField(); name = :N)],
        [connect(:S => :c, :N => :in)],
    )
end
