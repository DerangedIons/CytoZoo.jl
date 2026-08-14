using OrdinaryDiffEq

@testset "FHN — interface compliance" begin
    model = ParametrizedFHNModel()

    @test num_states(model) == 2
    @test num_parameters(model) == 5
    @test transmembrane_potential_index(model) == 1
    @test has_rush_larsen(model) == false
    @test num_monitors(model) == 0

    @test state_names(model) == (:v, :s)
    @test parameter_names(model) == (:a, :b, :c, :d, :e)

    @test state_index(model, :v) == 1
    @test state_index(model, :s) == 2
    @test parameter_index(model, :a) == 1
    @test parameter_index(model, :e) == 5

    # The interface convention: unknown names return `nothing`, never throw. `couple`'s
    # connect validation tests `=== nothing`, so a throwing lookup breaks its error message.
    @test state_index(model, :not_a_state) === nothing
    @test parameter_index(model, :not_a_parameter) === nothing

    u0 = default_initial_state(model)
    @test length(u0) == 2
    @test eltype(u0) == Float64
    @test u0 == [0.0, 0.0]
end

@testset "FHN — constructors and defaults" begin
    # `ParametrizedFHNModel` is Thunderbolt's spelling of the same type, not a second one.
    # Every construction below therefore exercises both names at once.
    @test ParametrizedFHNModel === FHNModel

    m = ParametrizedFHNModel()
    @test (m.a, m.b, m.c, m.d, m.e) == (0.1, 0.5, 1.0, 0.0, 0.01)
    @test m isa FHNModel{Float64}

    m32 = ParametrizedFHNModel(Float32)
    @test m32.a isa Float32
    @test m32.e isa Float32
    @test eltype(default_initial_state(m32)) == Float32
    @test m32 isa FHNModel                 # the unparameterized name matches any eltype...
    @test !(m32 isa FHNModel{Float64})     # ...so the eltype claim needs the parameter

    # Per-parameter overrides keep the rest at their defaults and promote to T.
    mk = ParametrizedFHNModel(; a = 0.25, e = 1 // 200)
    @test mk.a == 0.25
    @test mk.e == 0.005
    @test mk.b == 0.5

    # The whole model is isbits with the default (isbits) stimulus — this is what lets it
    # be captured by value inside a GPU kernel.
    @test isbitstype(typeof(m))
    @test isbitstype(typeof(m32))

    # `FunctionStimulus` is isbits *iff* its function is, which is the contract its
    # docstring states. A non-capturing closure is a singleton type, so it stays isbits
    # and can still ride into a GPU kernel...
    @test isbitstype(typeof(ParametrizedFHNModel(; stim = FunctionStimulus((x, t) -> 0.0))))

    # ...while one that closes over heap data is the CPU-only escape hatch.
    waveform = [0.0, -0.5]
    @test !isbitstype(
        typeof(ParametrizedFHNModel(; stim = FunctionStimulus((x, t) -> waveform[1]))),
    )
end

@testset "FHN — resting state is a genuine fixed point" begin
    model = ParametrizedFHNModel()
    u = default_initial_state(model)
    du = fill(NaN, 2)

    model(du, u, nothing, 0.0)

    # Unlike ToRORd's un-paced default, (0, 0) really is the steady state, so a solve
    # started there stays there until something stimulates it.
    @test du == [0.0, 0.0]
end

@testset "FHN — cubic threshold behaviour" begin
    model = ParametrizedFHNModel()          # a = 0.1
    du = fill(NaN, 2)

    # Below threshold the cubic pulls v back down; above it, v runs away to the upper
    # branch. This sign flip across v = a IS the FitzHugh–Nagumo model — if it is wrong,
    # nothing built on top of it will propagate.
    for (v, want_positive) in ((0.05, false), (0.15, true))
        model(du, [v, 0.0], nothing, 0.0)
        @test (du[1] > 0) == want_positive
    end

    # v = a and v = 0 and v = 1 are the three roots of the cubic, so with s = 0 the
    # voltage derivative vanishes at each.
    for v in (0.0, 0.1, 1.0)
        model(du, [v, 0.0], nothing, 0.0)
        @test du[1] == 0.0
    end
end

@testset "FHN — recovery equation" begin
    model = ParametrizedFHNModel(; b = 0.5, c = 1.0, d = 0.2, e = 0.01)
    du = fill(NaN, 2)
    v, s = 0.8, 0.3

    model(du, [v, s], nothing, 0.0)

    @test du[2] ≈ model.e * (model.b * v - model.c * s - model.d)
    @test du[1] ≈ v * (1 - v) * (v - model.a) - s
end

@testset "FHN — stimulus sign convention" begin
    # Negative amplitude depolarizes, matching ToRORd and the rest of the zoo.
    stim = Stimulus(; amplitude = -0.5, period = 100.0, duration = 1.0, start = 0.0)
    model = ParametrizedFHNModel(; stim)
    u = default_initial_state(model)
    du = fill(NaN, 2)

    model(du, u, nothing, 0.0)          # inside the pulse
    @test du[1] ≈ 0.5

    model(du, u, nothing, 5.0)          # outside the pulse
    @test du[1] == 0.0

    # And the default model is quiet at every phase of the default period — a tissue
    # framework owns stimulation, so a nonzero default would double-stimulate.
    quiet = ParametrizedFHNModel()
    for t in (0.0, 0.5, 1.5, 500.0, 1000.0)
        quiet(du, u, nothing, t)
        @test du[1] == 0.0
    end
end

@testset "FHN — SpatialContext overrides" begin
    model = ParametrizedFHNModel()
    du_plain = fill(NaN, 2)
    du_spatial = fill(NaN, 2)
    u = [0.15, 0.0]

    # Raising the threshold above v turns a firing cell into a decaying one — a change
    # of sign, not just of magnitude, so this cannot pass by accident.
    p = SpatialContext((0.0, 0.0, 0.0), (a = Constant(0.3),))
    model(du_plain, u, nothing, 0.0)
    model(du_spatial, u, p, 0.0)
    @test du_plain[1] > 0
    @test du_spatial[1] < 0

    # A position-dependent override: SpatialStep in x picks a different `a` on each side.
    step = SpatialStep(1, 0.5, 0.05, 0.3)      # a = 0.05 left of x = 0.5, 0.3 right of it
    du_left, du_right = fill(NaN, 2), fill(NaN, 2)
    model(du_left, u, SpatialContext((0.25, 0.0, 0.0), (a = step,)), 0.0)
    model(du_right, u, SpatialContext((0.75, 0.0, 0.0), (a = step,)), 0.0)
    @test du_left[1] > 0
    @test du_right[1] < 0

    # Every parameter is overridable, and an override of a parameter the context does not
    # mention falls back to the model field.
    p_all = SpatialContext(
        (0.0,),
        (a = 0.2, b = 0.6, c = 1.1, d = 0.05, e = 0.02),
    )
    model(du_spatial, [0.8, 0.3], p_all, 0.0)
    @test du_spatial[1] ≈ 0.8 * (1 - 0.8) * (0.8 - 0.2) - 0.3
    @test du_spatial[2] ≈ 0.02 * (0.6 * 0.8 - 1.1 * 0.3 - 0.05)

    @test isbitstype(typeof(SpatialContext((0.0, 0.0), (a = step,))))
end

@testset "FHN — Float32 genericity" begin
    model = ParametrizedFHNModel(Float32; stim = Stimulus(Float32; amplitude = -0.5f0))
    u = default_initial_state(model)
    du = fill(NaN32, 2)

    model(du, u, nothing, Float32(0))

    @test eltype(du) == Float32
    @test all(isfinite, du)
    @test du[1] ≈ 0.5f0

    # Spatial path too — `_resolve_spatial` results are wrapped in T.
    p = SpatialContext((0.0f0,), (a = Constant(0.2f0),))
    model(du, Float32[0.8, 0.3], p, Float32(5))
    @test eltype(du) == Float32
    @test all(isfinite, du)
end

@testset "FHN — zero allocations (functor)" begin
    model = ParametrizedFHNModel(; stim = Stimulus(; amplitude = -0.5))
    u = default_initial_state(model)
    du = similar(u)
    p = SpatialContext((0.0, 0.0, 0.0), (a = Constant(0.2),))

    model(du, u, nothing, 0.0)          # warm up both dispatches
    model(du, u, p, 0.0)

    @test @allocated(model(du, u, nothing, 0.0)) == 0
    @test @allocated(model(du, u, p, 0.0)) == 0
end

@testset "FHN — full excitation cycle under a solver" begin
    # Real usage: one suprathreshold pulse, then watch the pulse rise, plateau on the
    # upper branch, and recover. `e = 0.01` sets the recovery timescale, so a few hundred
    # time units covers a whole cycle.
    model = ParametrizedFHNModel(;
        stim = Stimulus(; amplitude = -0.5, period = 1.0e6, duration = 1.0, start = 0.0),
    )
    prob = ODEProblem(model, (0.0, 400.0))
    sol = solve(prob, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-8, saveat = 1.0)

    v = [u[1] for u in sol.u]
    s = [u[2] for u in sol.u]
    vmax, imax = findmax(v)

    @info "FHN excitation cycle" vmax t_peak = sol.t[imax] v_end = v[end] s_max = maximum(s) s_end = s[end]

    @test sol.retcode == ReturnCode.Success
    @test vmax > 0.8                       # depolarizes onto the upper branch
    @test sol.t[imax] < 100.0              # upstroke is fast relative to recovery
    # `s` has to charge past the knee of the v-nullcline `s = v(1-v)(v-a)` — that is what
    # tips `v` off the upper branch and forces the recovery. With `a = 0.1` the cubic peaks
    # at s ≈ 0.1262 (v ≈ 0.6846), so the knee, not some round number, is the bound worth
    # asserting; the run reaches ≈ 0.155.
    @test maximum(s) > 0.1262
    @test v[end] < 0.05                    # and the cell returns to rest
    @test s[end] < 0.1

    # No stimulus at all → the cell never leaves the fixed point.
    quiet = ParametrizedFHNModel()
    solq = solve(ODEProblem(quiet, (0.0, 400.0)), Tsit5(); abstol = 1.0e-10, reltol = 1.0e-8)
    @test maximum(abs, reduce(vcat, solq.u)) < 1.0e-8
end

@testset "FHN — connect edges are rejected with an actionable message" begin
    # The parameters are immutable struct fields, so the model cannot receive a `connect`.
    # `couple` must say so rather than failing later inside `setindex!`.
    src = ParametrizedFHNModel()
    dst = ParametrizedFHNModel()
    nodes = (Subsystem(src; name = :src), Subsystem(dst; name = :dst))
    err = try
        couple(nodes, (connect(:src => :v, :dst => :a),))
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("parameters", err.msg)
end
