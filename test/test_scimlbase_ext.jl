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

    @testset "monitor-sourced connect solves, explicit and implicit" begin
        # acc' = b = C - a with a(t) = 3exp(-t)  =>  acc(2) = 2C - 3(1 - exp(-2)).
        cm_m = couple(
            [Subsystem(_MonoDerived(); name = :D), Subsystem(_MonoReader(); name = :R)],
            [connect(:D => :b, :R => :d_ext)],
        )
        idx = state_index(cm_m, :R_acc)
        exact = 2 * _MONO_C - 3 * (1 - exp(-2.0))

        sol = solve(ODEProblem(cm_m, (0.0, 2.0)), Tsit5(); dt = 0.01, adaptive = false)
        @test sol.retcode == ReturnCode.Success
        @test isapprox(sol.u[end][idx], exact; atol = 1.0e-6)

        # Under Rodas5P the state slice reaching monitor_values! holds Duals; the pre-pass must
        # extract primals before computing, so neither the scratch nor the receiver's Float64
        # parameter slot ever sees one.
        sol_i = solve(ODEProblem(cm_m, (0.0, 2.0)), Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-10)
        @test sol_i.retcode == ReturnCode.Success
        @test isapprox(sol_i.u[end][idx], exact; atol = 1.0e-6)
    end

    @testset "contributory share integrates to the analytic solution" begin
        # A.p owns the slot with p' = -p; B.y contributes a constant 999 into it. Together
        # p' = -p + 999 with p(0) = 2, so p(t) = 999 - 997exp(-t) — a closed form that is wrong
        # under hard-discard (which would give 2exp(-t)) and wrong again if the contribution were
        # accumulated across evaluations instead of reset each one.
        cm_c = couple(
            [Subsystem(_ScatterA(); name = :A), Subsystem(_ScatterB(); name = :B)],
            [share(:A => :p, :B => :y; owner = :A, op = +)],
        )
        p_idx = state_index(cm_c, :p)
        exact = 999 - 997 * exp(-2.0)

        sol = solve(ODEProblem(cm_c, (0.0, 2.0)), Tsit5(); dt = 0.001, adaptive = false)
        @test sol.retcode == ReturnCode.Success
        @test isapprox(sol.u[end][p_idx], exact; rtol = 1.0e-6)

        # Under an implicit solver the accumulate path carries Duals. Unlike a `connect` input,
        # nothing extracts a primal, so the contributed term stays in the Jacobian.
        sol_i = solve(ODEProblem(cm_c, (0.0, 2.0)), Rodas5P(); reltol = 1.0e-10, abstol = 1.0e-12)
        @test sol_i.retcode == ReturnCode.Success
        @test isapprox(sol_i.u[end][p_idx], exact; rtol = 1.0e-6)
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
