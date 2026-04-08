using MTKCardiacCellModels
using ModelingToolkit
using OrdinaryDiffEq

@testset "MTKCardiacCellModels Extension" begin
    @testset "BeelerReuter model creation" begin
        model = BeelerReuter()
        @test num_states(model) == 8
        @test num_parameters(model) > 0
        @test has_symbolic_system(model)
        @test symbolic_system(model) !== nothing
        @test transmembrane_potential_index(model) > 0
    end

    @testset "BeelerReuter interface compliance" begin
        model = BeelerReuter()
        u0 = default_initial_state(model)
        @test length(u0) == num_states(model)
        @test length(state_names(model)) == num_states(model)
        @test length(parameter_names(model)) == num_parameters(model)
    end

    @testset "BeelerReuter functor" begin
        model = BeelerReuter()
        u0 = default_initial_state(model)
        du = similar(u0)
        model(du, u0, nothing, 0.0)
        @test all(isfinite, du)
    end

    @testset "BeelerReuter ODEProblem + solve" begin
        model = BeelerReuter()
        prob = ODEProblem(model, (0.0, 100.0))
        sol = solve(prob, Tsit5())
        @test sol.retcode == ReturnCode.Success
        @test all(isfinite, sol.u[end])
    end

    @testset "BeelerReuter caching" begin
        m1 = BeelerReuter()
        m2 = BeelerReuter()
        @test m1 === m2
    end
end
