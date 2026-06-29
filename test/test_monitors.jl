# DERIVED-observable tests: the monitor hook (num_monitors / monitor_names / monitor_values!),
# the post-solve monitor_history helper (SciMLBase ext), and CoupledModel aggregation
# (prefixed names + per-component state slicing). Needs a solver for monitor_history.

using OrdinaryDiffEq

# Mock model exercising a conservation-law-style monitor `mon = C - s`. States decay (s(t) =
# s0·exp(-t)) so the trajectory is analytic; the functor ignores p so it solves standalone.
struct _MonMock <: CytoZoo.AbstractCardiacCellModel
    C::Float64
end
CytoZoo.num_states(::_MonMock) = 2
CytoZoo.state_names(::_MonMock) = (:vm, :s)
CytoZoo.default_initial_state(::_MonMock) = [0.5, 1.0]
CytoZoo.state_index(::_MonMock, n::Symbol) = findfirst(==(n), (:vm, :s))
CytoZoo.transmembrane_potential_index(::_MonMock) = 1
function (::_MonMock)(du, u, p, t)
    du[1] = -u[1]
    du[2] = -u[2]
    return nothing
end
CytoZoo.num_monitors(::_MonMock) = 1
CytoZoo.monitor_names(::_MonMock) = (:deriv,)
function CytoZoo.monitor_values!(mon, u, t, m::_MonMock)
    mon[1] = m.C - u[2]
    return nothing
end

@testset "Monitors" begin
    @testset "single model — interface hooks" begin
        m = _MonMock(3.0)
        @test num_monitors(m) == 1
        @test monitor_names(m) == (:deriv,)
        mon = zeros(1)
        monitor_values!(mon, [0.5, 0.25], 0.0, m)
        @test mon[1] == 3.0 - 0.25
    end

    @testset "single model — defaults" begin
        # A model that opts out reports zero monitors / no names (base interface defaults).
        @test num_monitors(ToRORd()) == 0
        @test monitor_names(ToRORd()) == ()
    end

    @testset "single model — monitor_history recovers the law" begin
        m = _MonMock(3.0)
        sol = solve(ODEProblem(m, (0.0, 2.0)), Tsit5(); dt = 0.1, adaptive = false)
        h = monitor_history(sol, m)
        @test h.names == (:deriv,)
        @test h.t == sol.t
        @test size(h.values) == (1, length(sol.t))
        # Each column recomputed from that step's state: mon = C - s.
        @test all(h.values[1, j] == m.C - sol.u[j][2] for j in eachindex(sol.t))
        # Conservation law tracks the analytic trajectory s(t) = exp(-t).
        @test all(isapprox(h.values[1, j], m.C - exp(-sol.t[j]); atol = 1.0e-4) for j in eachindex(sol.t))
    end

    @testset "coupled — prefixed names + per-component slicing" begin
        cm = couple([Subsystem(_MonMock(2.0); name = :A), Subsystem(_MonMock(5.0); name = :B)])
        @test num_monitors(cm) == 2
        @test monitor_names(cm) == (:deriv, :B_deriv)          # primary bare, non-primary prefixed

        # Global state [vm_A, s_A, vm_B, s_B]; each component's monitor reads its OWN slice.
        U = [0.5, 0.25, 0.7, 0.6]
        mon = zeros(2)
        monitor_values!(mon, U, 0.0, cm)
        @test mon == [2.0 - 0.25, 5.0 - 0.6]
    end

    @testset "coupled — monitor_history end-to-end" begin
        cm = couple([Subsystem(_MonMock(2.0); name = :A), Subsystem(_MonMock(5.0); name = :B)])
        s_A = state_index(cm, :s)
        s_B = state_index(cm, :B_s)
        sol = solve(ODEProblem(cm, (0.0, 1.0)), Tsit5(); dt = 0.1, adaptive = false)
        h = monitor_history(sol, cm)
        @test h.names == (:deriv, :B_deriv)
        @test size(h.values) == (2, length(sol.t))
        @test all(h.values[1, j] == 2.0 - sol.u[j][s_A] for j in eachindex(sol.t))
        @test all(h.values[2, j] == 5.0 - sol.u[j][s_B] for j in eachindex(sol.t))
    end

    @testset "zero monitors — empty history, no error" begin
        model = ToRORd()
        sol = solve(ODEProblem(model, (0.0, 1.0)), Tsit5(); adaptive = true)
        h = monitor_history(sol, model)
        @test h.names == ()
        @test size(h.values) == (0, length(sol.t))
    end
end
