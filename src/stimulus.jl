# Time-only stimulus types shared across cell models.
#
# Spatial modulation is layered on top via `SpatialContext` — the stimulus
# itself takes only `t`, and the consumer composes it with position when
# needed.

# GPU-safe `mod` — Float64/Float32 use `Base.mod`; the Real method covers
# ForwardDiff.Dual; the generic fallback handles exotic float types.
@inline _safe_mod(a::Float64, b::Float64) = mod(a, b)
@inline _safe_mod(a::Float32, b::Float32) = mod(a, b)
@inline _safe_mod(a::Real, b::Real) = mod(a, b)
@inline _safe_mod(a, b) = a - b * oftype(a, floor(BigFloat(a / b)))

"""
    Stimulus{F}

Callable wrapper for a time-only stimulus current function. `s(t)` ≡ `s.f(t)`.

Spatial modulation is not part of this type — it belongs to a separate layer
(e.g. `SpatialContext` composition).

# Constructors
```julia
Stimulus()                                            # default periodic pulse
Stimulus(; amplitude=-53.0, period=1000.0, duration=1.0, start=0.0)
Stimulus(-40.0)                                       # constant amplitude
Stimulus(t -> ifelse(mod(t, 500) < 1, -53.0, 0.0))    # custom function
```
"""
struct Stimulus{F}
    f::F
end

Stimulus(f::Function) = Stimulus{typeof(f)}(f)
(s::Stimulus)(t) = s.f(t)

function Stimulus(;
        amplitude::Real = -53.0, period::Real = 1000.0,
        duration::Real = 1.0, start::Real = 0.0,
    )
    amp, per, dur, st = promote(amplitude, period, duration, start)
    z = zero(amp)
    return Stimulus(t -> ifelse(mod(t - st, per) < dur, amp, z))
end

# Typed convenience — constructs a default periodic pulse with all parameters
# in element type `T`. Useful for Float32 / GPU paths where Float64-captured
# closures would force per-step conversion.
function Stimulus(::Type{T};
        amplitude::Real = -53, period::Real = 1000,
        duration::Real = 1, start::Real = 0,
    ) where {T <: Real}
    return Stimulus(;
        amplitude = T(amplitude), period = T(period),
        duration = T(duration), start = T(start),
    )
end

Stimulus(amp::Real) = Stimulus(; amplitude = amp)
Stimulus(s::Stimulus) = s

"""
    StimulusParametric{T}

Closure-free, isbits stimulus for GPU and integrators that require scalar
amplitude/period/duration/start (e.g. Rush-Larsen / gotranx-style).
"""
struct StimulusParametric{T}
    amplitude::T
    period::T
    duration::T
    start::T
end

function StimulusParametric(;
        amplitude::T = -53.0, period::T = 1000.0,
        duration::T = 1.0, start::T = 0.0,
    ) where {T <: Real}
    return StimulusParametric{T}(amplitude, period, duration, start)
end

@inline function (s::StimulusParametric)(t)
    return ifelse(_safe_mod(t - s.start, s.period) < s.duration, s.amplitude, zero(s.amplitude))
end

"""
    stim_eval(s, x, y, z, t)

Evaluate a stimulus at position `(x, y, z)` and time `t`. The base types
(`Stimulus`, `StimulusParametric`) ignore position; the signature exists so
that downstream packages (e.g. CytoZoo's MTK extension, were one to return)
can register a single symbolic function.
"""
@inline stim_eval(s, _x, _y, _z, t) = s(t)
