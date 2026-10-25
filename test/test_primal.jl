using ForwardDiff
using Test

@testset "primal" begin
    @testset "plain numbers pass through" begin
        @test CytoZoo.primal(1.5) === 1.5
        @test CytoZoo.primal(1.5f0) === 1.5f0
        @test CytoZoo.primal(2) === 2
    end

    @testset "Duals are unwrapped" begin
        @test CytoZoo.primal(ForwardDiff.Dual(0.0, 1.0)) === 0.0
        nested = ForwardDiff.Dual(ForwardDiff.Dual(0.5, 1.0), ForwardDiff.Dual(1.0, 0.0))
        @test CytoZoo.primal(nested) === 0.5
    end

    @testset "guarded comparison differentiates cleanly at the guard" begin
        f(x) = (CytoZoo.primal(x) > 0) ? x / x : zero(x)
        g(x) = (x > 0) ? x / x : zero(x)
        @test ForwardDiff.derivative(f, 0.0) == 0.0
        if pkgversion(ForwardDiff) >= v"1"
            @test isnan(ForwardDiff.derivative(g, 0.0))
        end
    end

    @testset "type stability and allocations" begin
        d = ForwardDiff.Dual(1.0, 1.0)
        @test (@inferred CytoZoo.primal(d)) === 1.0
        alloc_primal(x) = @allocated CytoZoo.primal(x)
        alloc_primal(d)
        @test alloc_primal(d) == 0
    end
end
