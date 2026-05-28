@testset "Stimulus" begin
    @testset "default periodic pulse" begin
        s = Stimulus()
        # Default: amp=-53, period=1000, duration=1, start=0 → fires at t∈[0,1), [1000,1001), …
        @test s(0.0) == -53.0
        @test s(0.5) == -53.0
        @test s(1.0) == 0.0
        @test s(999.9) == 0.0
        @test s(1000.0) == -53.0
        @test stim_eval(s, 0.0, 0.0, 0.0, 0.0) == -53.0
        @test stim_eval(s, 0.0, 0.0, 0.0, 500.0) == 0.0
    end

    @testset "custom amplitude / timing kwargs" begin
        s = Stimulus(; amplitude = -40.0, period = 500.0, duration = 2.0, start = 10.0)
        @test s(9.999) == 0.0
        @test s(10.0) == -40.0
        @test s(11.5) == -40.0
        @test s(12.0) == 0.0
        @test s(510.0) == -40.0
    end

    @testset "constant amplitude convenience" begin
        s = Stimulus(-25.0)
        @test s(0.0) == -25.0
        @test s(0.999) == -25.0
        @test s(1.0) == 0.0
    end

    @testset "custom callable" begin
        s = Stimulus(t -> ifelse(t < 5.0, -10.0, 5.0))
        @test s(0.0) == -10.0
        @test s(10.0) == 5.0
    end

    @testset "idempotent on Stimulus" begin
        s = Stimulus(; amplitude = -53.0)
        @test Stimulus(s) === s
    end

    @testset "typed constructor preserves element type" begin
        s64 = Stimulus(Float64)
        @test s64(0.0) === -53.0          # default amplitude, Float64

        s32 = Stimulus(Float32)
        v = s32(Float32(0.0))
        @test v === Float32(-53)          # native Float32 — no widening
        @test typeof(v) === Float32

        # promote() based kwarg constructor: mixing types yields the wider type
        s_mixed = Stimulus(; amplitude = Float32(-53), period = 1000.0,
                            duration = 1.0, start = 0.0)
        @test typeof(s_mixed(0.0)) === Float64  # period/duration/start are Float64
    end
end

@testset "StimulusParametric" begin
    sp = StimulusParametric(; amplitude = -53.0, period = 1000.0, duration = 1.0, start = 0.0)
    @test sp(0.0) == -53.0
    @test sp(0.5) == -53.0
    @test sp(1.0) == 0.0
    @test sp(1000.0) == -53.0
    @test stim_eval(sp, 0.0, 0.0, 0.0, 500.0) == 0.0

    sp32 = StimulusParametric(; amplitude = Float32(-53), period = Float32(1000),
                              duration = Float32(1), start = Float32(0))
    @test sp32(Float32(0.5)) === Float32(-53)
    @test sp32 isa StimulusParametric{Float32}

    # isbits is the GPU/Rush-Larsen prerequisite
    @test isbits(sp)
    @test isbits(sp32)
end
