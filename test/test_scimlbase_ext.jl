using OrdinaryDiffEq

@testset "SciMLBase Extension" begin
    model = ToRORd()
    u0 = default_initial_state(model)

    @testset "ODEProblem(model, tspan)" begin
        prob = ODEProblem(model, (0.0, 1.0))
        @test prob.u0 == u0
        @test prob.tspan == (0.0, 1.0)
    end

    @testset "ODEProblem(model, tspan; u0=...)" begin
        custom_u0 = copy(u0)
        custom_u0[1] = -80.0
        prob = ODEProblem(model, (0.0, 1.0); u0 = custom_u0)
        @test prob.u0 == custom_u0
    end

    @testset "solve produces finite results" begin
        prob = ODEProblem(model, (0.0, 1.0))
        sol = solve(prob, Tsit5(); adaptive = true)
        @test sol.retcode == ReturnCode.Success
        @test all(isfinite, sol.u[end])
    end
end

# Monolithic single-RHS coupling solved through ODEProblem(cm, tspan). The toy models
# (_MonoA/_MonoReader/_MonoP/_MonoQ) are defined in test_coupling.jl, included earlier.
@testset "Monolithic coupling solve" begin
    analytic = 1 - exp(-2.0)              # acc(2) for acc' = d, d = exp(-t)
    cm = couple(
        [Subsystem(_MonoA(); name = :A), Subsystem(_MonoReader(); name = :R)],
        [connect(:A => :d, :R => :d_ext)],
    )
    acc_idx = state_index(cm, :R_acc)

    @testset "accuracy vs analytic (explicit, fixed dt)" begin
        sol = solve(ODEProblem(cm, (0.0, 2.0)), Tsit5(); dt = 0.05, adaptive = false)
        @test sol.retcode == ReturnCode.Success
        @test isapprox(sol.u[end][acc_idx], analytic; atol = 1.0e-8)
    end

    @testset "connect under an implicit solver (ForwardDiffExt freeze)" begin
        # ForwardDiff threads Duals through U; the connect write must not store a Dual in the
        # receiver's Float64 parameter slot. ForwardDiffExt freezes the connect input to its primal.
        sol = solve(ODEProblem(cm, (0.0, 2.0)), Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-10)
        @test sol.retcode == ReturnCode.Success
        @test isapprox(sol.u[end][acc_idx], analytic; atol = 1.0e-6)
    end

    @testset "stiff share coupling is stable (implicit)" begin
        cm_s = couple(
            [Subsystem(_MonoP(); name = :P), Subsystem(_MonoQ(); name = :Q)],
            [share(:P => :v, :Q => :v; owner = :P)],
        )
        v_idx = state_index(cm_s, :v)
        sol = solve(ODEProblem(cm_s, (0.0, 5.0)), Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-10)
        @test sol.retcode == ReturnCode.Success
        @test all(isfinite, sol.u[end])
        @test isapprox(sol.u[end][v_idx], exp(-0.1 * 5.0); atol = 1.0e-2)   # v relaxed to a
    end
end
