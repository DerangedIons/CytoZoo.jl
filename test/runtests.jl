using CytoZoo
using Test

@testset "CytoZoo.jl" begin
    include("test_interface.jl")
    include("test_stimulus.jl")
    include("test_coupling.jl")
    if Base.find_package("OrdinaryDiffEqOperatorSplitting") !== nothing
        include("test_coupling_ext.jl")
    else
        @warn "OrdinaryDiffEqOperatorSplitting not available, skipping coupling extension tests"
    end
    include("test_torord.jl")
    include("test_torord_correctness.jl")
    include("test_scimlbase_ext.jl")
    if Base.find_package("TWorld") !== nothing
        include("test_tworld.jl")
    else
        @warn "TWorld not available, skipping TWorld extension tests"
    end
end
