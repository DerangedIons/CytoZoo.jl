using CytoZoo
using Test
using JET

@testset "CytoZoo.jl" begin
    @testset "Code linting (JET.jl)" begin
        JET.test_package(CytoZoo; target_defined_modules = true)
    end
    # Write your tests here.
end
