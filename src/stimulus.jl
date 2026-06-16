# Stimulus current models. See `AbstractStimulus` (interface.jl) for the
# `(x, t) -> current` contract. Spatial dependence is first-class: a stimulus
# may index the position `x`; a position-independent one ignores it.

"""
    Stimulus{T}

Periodic rectangular pulse: returns `amplitude` for `duration` time units every
`period`, starting at `start`, and zero otherwise. Position-independent (the
built-in default), isbits, and GPU/Rush-Larsen safe.

# Constructors
```julia
Stimulus()                                            # default: -53, period 1000, 1 ms
Stimulus(; amplitude=-53.0, period=1000.0, duration=1.0, start=0.0)
Stimulus(Float32)                                     # element type T (Float32/GPU paths)
```
"""
struct Stimulus{T} <: AbstractStimulus
    amplitude::T
    period::T
    duration::T
    start::T
end

@inline (s::Stimulus)(x, t) =
    ifelse(mod(t - s.start, s.period) < s.duration, s.amplitude, zero(s.amplitude))

function Stimulus(;
        amplitude::Real = -53.0, period::Real = 1000.0,
        duration::Real = 1.0, start::Real = 0.0,
    )
    amp, per, dur, st = promote(amplitude, period, duration, start)
    return Stimulus(amp, per, dur, st)
end

# Typed convenience — all parameters in element type `T`. Useful for Float32 /
# GPU paths where Float64 defaults would force per-step conversion.
function Stimulus(::Type{T};
        amplitude::Real = -53, period::Real = 1000,
        duration::Real = 1, start::Real = 0,
    ) where {T <: Real}
    return Stimulus(T(amplitude), T(period), T(duration), T(start))
end

"""
    FunctionStimulus{F}

Stimulus wrapping an arbitrary `(x, t) -> current` function, for waveforms the
parametric `Stimulus` can't express (biphasic, S1–S2, ramps). Time-only
functions just ignore `x`. Isbits iff `F` is — closures capturing boxed values
are CPU-only.
"""
struct FunctionStimulus{F} <: AbstractStimulus
    f::F
end

@inline (s::FunctionStimulus)(x, t) = s.f(x, t)
