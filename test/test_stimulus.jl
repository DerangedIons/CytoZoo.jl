@testset "Stimulus" begin
    @testset "default periodic pulse" begin
        s = Stimulus()
        # Default: amp=-53, period=1000, duration=1, start=0 → fires at t∈[0,1), [1000,1001), …
        # Position-independent: `x` is ignored (here `nothing`, the non-spatial path).
        @test s(nothing, 0.0) == -53.0
        @test s(nothing, 0.5) == -53.0
        @test s(nothing, 1.0) == 0.0
        @test s(nothing, 999.9) == 0.0
        @test s(nothing, 1000.0) == -53.0
    end

    @testset "custom amplitude / timing kwargs" begin
        s = Stimulus(; amplitude = -40.0, period = 500.0, duration = 2.0, start = 10.0)
        @test s(nothing, 9.999) == 0.0
        @test s(nothing, 10.0) == -40.0
        @test s(nothing, 11.5) == -40.0
        @test s(nothing, 12.0) == 0.0
        @test s(nothing, 510.0) == -40.0
    end

    @testset "isbits / typed constructor preserves element type" begin
        @test isbits(Stimulus())

        s64 = Stimulus(Float64)
        @test s64(nothing, 0.0) === -53.0          # default amplitude, Float64

        s32 = Stimulus(Float32)
        @test s32 isa Stimulus{Float32}
        @test isbits(s32)
        v = s32(nothing, Float32(0.0))
        @test v === Float32(-53)                   # native Float32 — no widening

        # promote()-based kwarg constructor: mixing types yields the wider type
        s_mixed = Stimulus(; amplitude = Float32(-53), period = 1000.0,
                           duration = 1.0, start = 0.0)
        @test typeof(s_mixed(nothing, 0.0)) === Float64
    end
end

@testset "FunctionStimulus" begin
    @testset "arbitrary time-only waveform" begin
        s = FunctionStimulus((x, t) -> ifelse(t < 5.0, -10.0, 5.0))
        @test s(nothing, 0.0) == -10.0
        @test s(nothing, 10.0) == 5.0
    end

    @testset "first-class spatial — fires only in-region" begin
        # Localized periodic pulse: active only where x[1] > 0.5
        s = FunctionStimulus((x, t) -> (mod(t, 1000.0) < 1.0 && x[1] > 0.5) ? -53.0 : 0.0)
        @test s([1.0, 0.0, 0.0], 0.0) == -53.0    # in-region, in-pulse
        @test s([0.0, 0.0, 0.0], 0.0) == 0.0      # out-of-region
        @test s([1.0, 0.0, 0.0], 500.0) == 0.0    # in-region, out-of-pulse

        # isbits when the wrapped function captures nothing boxed
        @test isbits(FunctionStimulus((x, t) -> zero(t)))
    end
end
