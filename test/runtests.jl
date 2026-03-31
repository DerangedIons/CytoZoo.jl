using CytoZoo
using Test

@testset "CytoZoo.jl" begin
    include("test_interface.jl")
    include("test_torord.jl")
    include("test_torord_correctness.jl")
end
