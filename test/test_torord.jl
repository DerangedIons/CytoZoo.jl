@testset "ToRORd — ODE evaluation" begin
    model = ToRORd()
    u = default_initial_state(model)
    du = similar(u)

    model(du, u, nothing, 0.0)

    @test !any(isnan, du)
    @test !any(isinf, du)
    @test du[1] != 0.0  # dVm/dt should be nonzero (stimulus at t=0)
end

@testset "ToRORd — custom stimulus" begin
    # Constant zero stim → dVm/dt at rest matches the default-stim "off" phase
    model_off = ToRORd(; stim = FunctionStimulus((x, t) -> 0.0))
    model_default = ToRORd()
    u = default_initial_state(model_default)

    du_off = similar(u)
    du_at_pulse_end = similar(u)

    model_off(du_off, u, nothing, 0.0)
    # t = 2.0 → default stim is OFF (pulse 0–1 ms, period 1000); off-stim model
    # at any t should be identical to default at t = 2.0.
    model_default(du_at_pulse_end, u, nothing, 2.0)

    @test du_off[1] ≈ du_at_pulse_end[1] atol = 1e-12
end

@testset "ToRORd — Rush-Larsen step" begin
    model = ToRORd()
    u = default_initial_state(model)
    u_new = similar(u)

    rush_larsen_step!(u_new, u, nothing, 0.0, 0.01, model)

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

@testset "ToRORd — SpatialContext evaluation" begin
    model = ToRORd()
    u = default_initial_state(model)
    du = similar(u)
    p = SpatialContext([0.5, 0.0, 0.0], (
        IKr_Multiplier = (x, t) -> 0.5,
        isHypoxic = (x, t) -> 0.0,
        celltype = (x, t) -> 0.0,
    ))

    model(du, u, p, 0.0)

    @test !any(isnan, du)
    @test !any(isinf, du)
end

@testset "ToRORd — SpatialContext with scalar overrides" begin
    model = ToRORd()
    u = default_initial_state(model)
    du = similar(u)
    p = SpatialContext([0.5, 0.0, 0.0], (IKr_Multiplier = 0.5,))

    model(du, u, p, 0.0)

    @test !any(isnan, du)
    @test !any(isinf, du)
end

@testset "ToRORd — Rush-Larsen with SpatialContext" begin
    model = ToRORd()
    u = default_initial_state(model)
    u_new = similar(u)
    p = SpatialContext([0.5, 0.0, 0.0], (
        IKr_Multiplier = (x, t) -> 0.5,
    ))

    rush_larsen_step!(u_new, u, p, 0.0, 0.01, model)

    @test !any(isnan, u_new)
    @test !any(isinf, u_new)
end
