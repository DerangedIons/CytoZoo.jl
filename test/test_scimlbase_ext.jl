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
        prob = ODEProblem(model, (0.0, 1.0); u0=custom_u0)
        @test prob.u0 == custom_u0
    end

    @testset "solve produces finite results" begin
        prob = ODEProblem(model, (0.0, 1.0))
        sol = solve(prob, Tsit5(); adaptive=true)
        @test sol.retcode == ReturnCode.Success
        @test all(isfinite, sol.u[end])
    end
end
