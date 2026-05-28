using TWorld

# `TWorldCellModel` is now exported by TWorld itself (native adherence to the
# CytoZoo interface). No CytoZoo wrapper / ext involved.

@testset "TWorldCellModel — construction" begin
    for celltype in (0.0, 1.0, 2.0), sex in (0.0, 1.0, 2.0)
        model = TWorldCellModel(; celltype, sex)
        @test model.params.celltype == celltype
        @test model.params.sex == sex
    end
end

@testset "TWorldCellModel — ODE evaluation" begin
    model = TWorldCellModel()
    u = default_initial_state(model)
    du = similar(u)

    model(du, u, nothing, 0.0)

    @test !any(isnan, du)
    @test !any(isinf, du)
    @test du[1] != 0.0  # dVm/dt should be nonzero (stimulus at t=0)
end

@testset "TWorldCellModel — cross-validation vs direct TWorld" begin
    model = TWorldCellModel()
    p = TWorld.TWorldParameters()
    u = TWorld.tworld_initial_conditions()

    du_cytozoo = similar(u)
    du_direct = similar(u)

    model(du_cytozoo, u, nothing, 0.0)
    TWorld.tworld_ode!(du_direct, u, p, 0.0)

    @test du_cytozoo == du_direct  # bitwise identical (pure delegation)
end

@testset "TWorldCellModel — interface completeness" begin
    model = TWorldCellModel()

    @test num_states(model) == 92
    @test num_parameters(model) > 0
    @test transmembrane_potential_index(model) == 1
    # TWorld does not expose a single-step Rush-Larsen primitive; users go
    # through `solve(ODEProblem(model, tspan), TWorld.RushLarsen(); dt=...)`.
    @test has_rush_larsen(model) == false

    u0 = default_initial_state(model)
    @test length(u0) == 92
    @test !any(isnan, u0)

    snames = state_names(model)
    @test length(snames) == 92
    @test snames[1] == :v
    @test state_index(model, :v) == 1
    @test state_index(model, :ca_junc) == TWorld.IDX_CA_JUNC

    pnames = parameter_names(model)
    @test length(pnames) == num_parameters(model)
    @test :celltype in pnames
    @test :nao in pnames
    @test :stim_fn ∉ pnames
    @test :x_coord ∉ pnames

    @test parameter_index(model, :nao) isa Int
end

@testset "TWorldCellModel — SpatialContext dispatch" begin
    model = TWorldCellModel()
    u = default_initial_state(model)
    du = similar(u)

    p = SpatialContext([0.5, 0.0, 0.0], nothing)
    model(du, u, p, 0.0)

    @test !any(isnan, du)
    @test !any(isinf, du)
end

@testset "TWorldCellModel — allocation measurement" begin
    model = TWorldCellModel()
    u = default_initial_state(model)
    du = similar(u)

    # Warmup
    model(du, u, nothing, 0.0)

    alloc = @allocated model(du, u, nothing, 0.0)
    @info "TWorldCellModel functor allocations: $alloc bytes"
end

@testset "Stimulus / StimulusParametric — exported by both packages" begin
    # Both CytoZoo and TWorld expose the same Stimulus and StimulusParametric
    # types (CytoZoo owns them; TWorld re-exports). Identity check confirms
    # no shadowing.
    @test TWorld.Stimulus === CytoZoo.Stimulus
    @test TWorld.StimulusParametric === CytoZoo.StimulusParametric
    @test TWorld.stim_eval === CytoZoo.stim_eval
end
