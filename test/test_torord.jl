@testset "ToRORd — ODE evaluation" begin
    model = ToRORd()
    u = default_initial_state(model)
    du = similar(u)

    model(du, u, nothing, 0.0)

    @test !any(isnan, du)
    @test !any(isinf, du)
    @test du[1] != 0.0  # dVm/dt should be nonzero (stimulus at t=0)
end

@testset "ToRORd — Rush-Larsen step" begin
    model = ToRORd()
    u = default_initial_state(model)
    u_new = similar(u)

    rush_larsen_step!(u_new, u, 0.0, 0.01, model)

    @test !any(isnan, u_new)
    @test !any(isinf, u_new)
    @test u_new[1] != u[1]  # Vm should change
end

@testset "ToRORd — Float32 support" begin
    model = ToRORd(Float32)
    u = default_initial_state(model)
    du = similar(u)

    @test eltype(u) == Float32
    model(du, u, nothing, Float32(0))

    @test !any(isnan, du)
    @test !any(isinf, du)
end

@testset "ToRORd — zero allocations (functor)" begin
    model = ToRORd()
    u = default_initial_state(model)
    du = similar(u)

    # Warmup
    model(du, u, nothing, 0.0)

    alloc = @allocated model(du, u, nothing, 0.0)
    @test alloc == 0
end

@testset "ToRORd — DiffEq functor signature" begin
    model = ToRORd()
    u = default_initial_state(model)
    du = similar(u)

    # Should work with any p (ignored)
    model(du, u, nothing, 0.0)
    model(du, u, :unused, 0.0)
    model(du, u, 42, 0.0)

    @test !any(isnan, du)
end

@testset "ToRORd — Spatial wrapper evaluation" begin
    base = ToRORd()
    spatial = Spatial(base, (
        IKr_Multiplier = (x, t) -> 0.5,
        isHypoxic = (x, t) -> 0.0,
        celltype = (x, t) -> 0.0,
    ))

    u = default_initial_state(spatial)
    du = similar(u)
    x = [0.5, 0.0, 0.0]

    spatial(du, u, x, 0.0)

    @test !any(isnan, du)
    @test !any(isinf, du)
end
