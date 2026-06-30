# Coupling tests: layout, naming, sharing, owner semantics, validation, and the monolithic
# single-RHS functor. Solver-driven coupling solves live in test_scimlbase_ext.jl.

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
    @test_throws ArgumentError num_parameters(cm)                                          # no single parameter vector
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
end
