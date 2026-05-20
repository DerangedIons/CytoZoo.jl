@testset "Interface compliance — ToRORd" begin
    model = ToRORd()

    @test num_states(model) == 65
    @test num_parameters(model) == 181
    @test transmembrane_potential_index(model) == 1
    @test has_rush_larsen(model) == true

    u0 = default_initial_state(model)
    @test length(u0) == 65
    @test u0[1] ≈ -91.33918  # Vm at rest

    @test state_index(model, :v) == 1
    @test state_index(model, :hL) == 63
    @test parameter_index(model, :GNa) == 14
    @test parameter_index(model, :IKr_Multiplier) == 27

    @test length(state_names(model)) == 65
    @test length(parameter_names(model)) == 181
    @test state_names(model)[1] == :v
    @test parameter_names(model)[14] == :GNa

    @test model.parameters[14] ≈ 11.7802  # GNa default
    @test model.celltype == 0
end

@testset "SpatialContext" begin
    sc = SpatialContext([1.0, 0.0, 0.0], (IKr_Multiplier = (x, t) -> 0.5,))
    @test sc.x == [1.0, 0.0, 0.0]
    @test sc.spatial_funcs.IKr_Multiplier([0.0], 0.0) == 0.5

    sc_scalar = SpatialContext([1.0, 0.0, 0.0], (IKr_Multiplier = 0.5,))
    @test CytoZoo._resolve_spatial(sc_scalar.spatial_funcs.IKr_Multiplier, [0.0], 0.0) == 0.5

    sc_step = SpatialContext(
        (1.0, 0.0, 0.0),
        (IKr_Multiplier = SpatialStep(1, 1.5, 1.0, 0.5),),
    )
    @test isbitstype(typeof(sc_step))
end

@testset "Constructors" begin
    m1 = ToRORd()
    @test eltype(m1.parameters) == Float64
    @test m1.celltype == 0

    m2 = ToRORd(Float32; celltype = 1)
    @test eltype(m2.parameters) == Float32
    @test m2.celltype == 1

    m3 = ToRORd(; celltype = 2)
    @test m3.celltype == 2
end
