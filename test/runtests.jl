using CytoZoo
using Test

@testset "CytoZoo.jl" begin
    include("test_interface.jl")
    include("test_torord.jl")
    include("test_torord_correctness.jl")
    include("test_scimlbase_ext.jl")
    include("test_mtk_ext.jl")
end
