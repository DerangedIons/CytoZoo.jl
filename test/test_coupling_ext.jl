# CouplingExt solving tests. Included only when OrdinaryDiffEqOperatorSplitting is
# resolvable (gated in runtests.jl, like the TWorld tests); skipped in CI where OS is absent.

using OrdinaryDiffEqOperatorSplitting
using OrdinaryDiffEq

# Toy component A: state `d` decays at rate 1, `v` is constant. Analytic: d(t) = d0·exp(-t).
struct _ExtToyA <: CytoZoo.AbstractCellModel end
CytoZoo.num_states(::_ExtToyA) = 2
CytoZoo.state_names(::_ExtToyA) = (:v, :d)
CytoZoo.default_initial_state(::_ExtToyA) = [0.0, 1.0]
CytoZoo.state_index(::_ExtToyA, n::Symbol) = findfirst(==(n), (:v, :d))
CytoZoo.transmembrane_potential_index(::_ExtToyA) = 1
function (::_ExtToyA)(du, u, p, t)
    du[1] = 0.0
    du[2] = -u[2]
    return nothing
end

# Toy component B with a *competing* equation for the shared state (`du(x) = +10`) that
# hard-discard must throw away, plus `y` that reads the shared value (`du(y) = -x`).
struct _ExtToyB <: CytoZoo.AbstractCellModel end
CytoZoo.num_states(::_ExtToyB) = 2
CytoZoo.state_names(::_ExtToyB) = (:x, :y)
CytoZoo.default_initial_state(::_ExtToyB) = [99.0, 0.0]   # x IC ignored — owner supplies the shared slot
CytoZoo.state_index(::_ExtToyB, n::Symbol) = findfirst(==(n), (:x, :y))
CytoZoo.transmembrane_potential_index(::_ExtToyB) = 1
function (::_ExtToyB)(du, u, p, t)
    du[1] = 10.0      # discarded for the shared slot (owner A governs it)
    du[2] = -u[1]     # y reads the shared value
    return nothing
end

@testset "CouplingExt — single node matches analytic" begin
    cm = couple([Subsystem(_ExtToyA(), Tsit5(); name = :A)])
    prob = OperatorSplittingProblem(cm, (0.0, 1.0))
    integ = init(prob, coupled_algorithm(cm); dt = 0.01, adaptive = false)
    solve!(integ)
    @test integ.u[state_index(cm, :v)] ≈ 0.0 atol = 1e-8
    @test integ.u[state_index(cm, :d)] ≈ exp(-1.0) rtol = 1e-4
end

@testset "CouplingExt — missing inner solver errors" begin
    cm = couple([Subsystem(_ExtToyA(); name = :A)])      # no alg passed to the node
    @test_throws ArgumentError coupled_algorithm(cm)
end

@testset "CouplingExt — share is hard-discard (owner equation only)" begin
    cm = couple(
        [Subsystem(_ExtToyA(), Tsit5(); name = :A), Subsystem(_ExtToyB(), Tsit5(); name = :B)],
        [share(:A => :d, :B => :x; owner = :A)],
    )
    @test num_states(cm) == 3
    @test default_initial_state(cm) == [0.0, 1.0, 0.0]          # shared d from owner A; B's y appended
    prob = OperatorSplittingProblem(cm, (0.0, 1.0))
    integ = init(prob, coupled_algorithm(cm); dt = 0.005, adaptive = false)
    solve!(integ)
    # shared slot follows ONLY owner A's decay; B's du(x) = +10 is discarded
    @test integ.u[state_index(cm, :d)] ≈ exp(-1.0) rtol = 1e-3
    # B read the owner-governed shared value (x = d): dy = -d ⇒ y(1) = exp(-1) - 1
    @test integ.u[state_index(cm, :B_y)] ≈ (exp(-1.0) - 1.0) atol = 2e-2
end

# Connect source: constant Vm = 2.0. Receiver reads it via a parameter slot.
struct _RefSrc <: CytoZoo.AbstractCellModel end
CytoZoo.num_states(::_RefSrc) = 2
CytoZoo.state_names(::_RefSrc) = (:Vm, :w)
CytoZoo.default_initial_state(::_RefSrc) = [2.0, 0.0]
CytoZoo.state_index(::_RefSrc, n::Symbol) = findfirst(==(n), (:Vm, :w))
CytoZoo.transmembrane_potential_index(::_RefSrc) = 1
(::_RefSrc)(du, u, p, t) = (du[1] = 0.0; du[2] = 0.0; nothing)

struct _RefRecv <: CytoZoo.AbstractCellModel
    parameters::Vector{Float64}     # one slot :Vm_ext, fed by the connect edge
end
_RefRecv() = _RefRecv([0.0])
CytoZoo.num_states(::_RefRecv) = 1
CytoZoo.state_names(::_RefRecv) = (:z,)
CytoZoo.default_initial_state(::_RefRecv) = [0.0]
CytoZoo.state_index(::_RefRecv, n::Symbol) = findfirst(==(n), (:z,))
CytoZoo.transmembrane_potential_index(::_RefRecv) = 1
CytoZoo.parameter_index(::_RefRecv, n::Symbol) = n === :Vm_ext ? 1 : nothing
(m::_RefRecv)(du, u, p, t) = (du[1] = m.parameters[1]; nothing)   # z integrates the synced Vm

@testset "CouplingExt — connect feeds a source state into a receiver param" begin
    recv = _RefRecv()
    cm = couple(
        [Subsystem(_RefSrc(), Tsit5(); name = :A), Subsystem(recv, Tsit5(); name = :B)],
        [connect(:A => :Vm, :B => :Vm_ext)],
    )
    @test num_states(cm) == 3                              # Vm, w, B_z
    prob = OperatorSplittingProblem(cm, (0.0, 1.0))
    integ = init(prob, coupled_algorithm(cm); dt = 0.1, adaptive = false)
    solve!(integ)
    @test integ.u[state_index(cm, :B_z)] ≈ 2.0 rtol = 1e-6  # z = ∫ Vm dt = 2·1, Vm constant 2
    @test recv.parameters[1] ≈ 2.0                          # connect edge wrote Vm into the receiver param
end

# Sum (op = +): two source states summed into one receiver parameter slot each step.
struct _SumSrc <: CytoZoo.AbstractCellModel end
CytoZoo.num_states(::_SumSrc) = 2
CytoZoo.state_names(::_SumSrc) = (:p, :q)
CytoZoo.default_initial_state(::_SumSrc) = [2.0, 3.0]
CytoZoo.state_index(::_SumSrc, n::Symbol) = findfirst(==(n), (:p, :q))
CytoZoo.transmembrane_potential_index(::_SumSrc) = 1
(::_SumSrc)(du, u, p, t) = (du[1] = 0.0; du[2] = 0.0; nothing)   # p, q held constant at 2, 3

struct _SumRecv <: CytoZoo.AbstractCellModel
    parameters::Vector{Float64}     # one slot :s, fed by two summed connect edges
end
_SumRecv() = _SumRecv([0.0])
CytoZoo.num_states(::_SumRecv) = 1
CytoZoo.state_names(::_SumRecv) = (:acc,)
CytoZoo.default_initial_state(::_SumRecv) = [0.0]
CytoZoo.state_index(::_SumRecv, n::Symbol) = findfirst(==(n), (:acc,))
CytoZoo.transmembrane_potential_index(::_SumRecv) = 1
CytoZoo.parameter_index(::_SumRecv, n::Symbol) = n === :s ? 1 : nothing
(m::_SumRecv)(du, u, p, t) = (du[1] = m.parameters[1]; nothing)   # acc integrates the summed slot

@testset "CouplingExt — connect op = + sums sources into one slot" begin
    recv = _SumRecv()
    cm = couple(
        [Subsystem(_SumSrc(), Tsit5(); name = :A), Subsystem(recv, Tsit5(); name = :B)],
        [connect(:A => :p, :B => :s; op = +), connect(:A => :q, :B => :s; op = +)],
    )
    integ = init(OperatorSplittingProblem(cm, (0.0, 1.0)), coupled_algorithm(cm); dt = 0.1, adaptive = false)
    solve!(integ)
    @test recv.parameters[1] ≈ 5.0                              # 2 + 3, reset-and-summed each step
    @test integ.u[state_index(cm, :B_acc)] ≈ 5.0 rtol = 1e-6   # ∫₀¹ (p + q) dt = 5
end

# Convergence toy: A's `a` decays (a = exp(-t)); B integrates it via a connect edge. The exact
# b(1) = ∫₀¹ exp(-t) dt = 1 - exp(-1). Splitting freezes `a` per step ⇒ O(dt) (1st-order) error.
struct _ConvA <: CytoZoo.AbstractCellModel end
CytoZoo.num_states(::_ConvA) = 1
CytoZoo.state_names(::_ConvA) = (:a,)
CytoZoo.default_initial_state(::_ConvA) = [1.0]
CytoZoo.state_index(::_ConvA, n::Symbol) = findfirst(==(n), (:a,))
CytoZoo.transmembrane_potential_index(::_ConvA) = 1
(::_ConvA)(du, u, p, t) = (du[1] = -u[1]; nothing)

struct _ConvB <: CytoZoo.AbstractCellModel
    parameters::Vector{Float64}
end
_ConvB() = _ConvB([0.0])
CytoZoo.num_states(::_ConvB) = 1
CytoZoo.state_names(::_ConvB) = (:b,)
CytoZoo.default_initial_state(::_ConvB) = [0.0]
CytoZoo.state_index(::_ConvB, n::Symbol) = findfirst(==(n), (:b,))
CytoZoo.transmembrane_potential_index(::_ConvB) = 1
CytoZoo.parameter_index(::_ConvB, n::Symbol) = n === :a_in ? 1 : nothing
(m::_ConvB)(du, u, p, t) = (du[1] = m.parameters[1]; nothing)

@testset "CouplingExt — LT-G shows ~1st-order convergence" begin
    function solve_b(dt)
        cm = couple(
            [Subsystem(_ConvA(), Tsit5(); name = :A), Subsystem(_ConvB(), Tsit5(); name = :B)],
            [connect(:A => :a, :B => :a_in)],
        )
        integ = init(OperatorSplittingProblem(cm, (0.0, 1.0)), coupled_algorithm(cm); dt = dt, adaptive = false)
        solve!(integ)
        return integ.u[state_index(cm, :B_b)]
    end
    exact = 1.0 - exp(-1.0)
    e1, e2, e3 = abs(solve_b(0.1) - exact), abs(solve_b(0.05) - exact), abs(solve_b(0.025) - exact)
    @test e1 > e2 > e3                       # error shrinks with dt
    @test 1.7 < e1 / e2 < 2.3                # halving dt ≈ halves error (1st order)
    @test 1.7 < e2 / e3 < 2.3
end

@testset "CouplingExt — operators allocate nothing in the hot path" begin
    cm = couple(
        [Subsystem(_ExtToyA(), Tsit5(); name = :A), Subsystem(_ExtToyB(), Tsit5(); name = :B)],
        [share(:A => :d, :B => :x; owner = :A)],
    )
    gsf = build_split_function(cm)
    for (i, ck) in enumerate(cm.layout.operator_order)
        op = gsf.functions[i].f                  # the ComponentOperator inside the ODEFunction
        n = num_states(cm.components[ck])
        du, u = zeros(n), ones(n)
        op(du, u, nothing, 0.0)                  # warmup
        @test @allocated(op(du, u, nothing, 0.0)) == 0
    end
end
