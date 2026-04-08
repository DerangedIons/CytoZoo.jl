using CytoZoo
using Test

@testset "CytoZoo.jl" begin
    include("test_interface.jl")
    include("test_torord.jl")
    include("test_torord_correctness.jl")
    include("test_scimlbase_ext.jl")
    if Base.find_package("MTKCardiacCellModels") !== nothing
        include("test_mtk_ext.jl")
    else
        @warn "MTKCardiacCellModels not available, skipping MTK extension tests"
    end
end
