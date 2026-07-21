# Stimulus

A cell model needs an external stimulus current to fire. In CytoZoo the stimulus is a
first-class object attached to the model, not a parameter buried in the parameter vector and
not something you patch in through a callback.

## The Contract

[`AbstractStimulus`](@ref) is the supertype. A subtype must be callable as

```julia
(s::AbstractStimulus)(x, t) -> current
```

and return the **full** stimulus current — not a multiplier, not an amplitude to be scaled
later. `x` is a position vector (the same one [`SpatialFunction`](@ref)s receive) and `t` is
time.

Two rules follow from that signature:

1. **Spatial dependence is first-class.** A stimulus may index `x`, so "stimulate only this
   region of tissue" needs no special support.
2. **A stimulus used on the non-spatial path must ignore `x`.** When a model is evaluated as
   `model(du, u, nothing, t)` there is no position to pass, so the stimulus is called as
   `s(nothing, t)`. A stimulus that dereferences `x` is only valid under a
   [`SpatialContext`](@ref).

## The Built-In Periodic Pulse

[`Stimulus`](@ref CytoZoo.Stimulus) is a rectangular periodic pulse — the common case, and the default on
`ToRORd`:

```@example stim
using CytoZoo

s = Stimulus()
```

It returns `amplitude` for `duration` time units every `period`, starting at `start`, and
zero otherwise:

```@example stim
[s(nothing, t) for t in (0.0, 0.5, 1.5, 1000.0, 1000.5)]
```

Configure it with keywords:

```@example stim
Stimulus(; amplitude = -40.0, period = 500.0, duration = 2.0, start = 50.0)
```

`Stimulus` holds no closure, so it is isbits and safe on a GPU and in a Rush-Larsen step:

```@example stim
isbitstype(typeof(s))
```

Passing an element type builds all four fields in that type, which avoids per-step
conversions on `Float32` and GPU paths:

```@example stim
typeof(Stimulus(Float32))
```

Attach one at construction:

```@example stim
model = ToRORd(; stim = Stimulus(; amplitude = -80.0, period = 500.0))
model.stim.period
```

## Arbitrary Waveforms

For anything the parametric pulse cannot express — biphasic pulses, S1–S2 protocols, ramps —
wrap a function in [`FunctionStimulus`](@ref):

```@example stim
biphasic = FunctionStimulus() do x, t
    τ = mod(t, 1000.0)
    τ < 1.0  ? -53.0 :
    τ < 2.0  ?  20.0 : 0.0
end

[biphasic(nothing, t) for t in (0.5, 1.5, 3.0)]
```

An S1–S2 protocol — a train of conditioning beats followed by one premature beat — is just
as direct:

```@example stim
function s1s2(x, t; s1_period = 1000.0, n_s1 = 8, s2_delay = 300.0)
    t_s2 = n_s1 * s1_period + s2_delay
    in_s1 = t < n_s1 * s1_period && mod(t, s1_period) < 1.0
    in_s2 = t_s2 <= t < t_s2 + 1.0
    return (in_s1 || in_s2) ? -53.0 : 0.0
end

model_s1s2 = ToRORd(; stim = FunctionStimulus(s1s2))
nothing # hide
```

`FunctionStimulus` is isbits if and only if the function it wraps is. A plain top-level
function or a closure over isbits values stays GPU-safe; a closure capturing a boxed or
heap-allocated value is CPU-only.

## Position-Dependent Stimulus

Because `x` is passed through, restricting a stimulus to part of a domain is a normal
predicate. This one fires only where the first coordinate exceeds 1.5:

```@example stim
local_stim = FunctionStimulus((x, t) -> (mod(t, 1000.0) < 1.0 && x[1] > 1.5) ? -53.0 : 0.0)

tissue_model = ToRORd(; stim = local_stim)

p_inside  = SpatialContext([2.0, 0.0, 0.0], nothing)
p_outside = SpatialContext([1.0, 0.0, 0.0], nothing)

(local_stim(p_inside.x, 0.5), local_stim(p_outside.x, 0.5))
```

This stimulus indexes `x`, so it must not be used on the non-spatial path — `model(du, u,
nothing, t)` would pass `nothing` as `x` and fail.

## Writing a Custom Stimulus

Subtype [`AbstractStimulus`](@ref) when you want an isbits stimulus with its own parameters —
the reason to prefer this over `FunctionStimulus` is keeping GPU compatibility while still
carrying configuration:

```@example stim
struct RampStimulus{T} <: AbstractStimulus
    amplitude::T
    period::T
    rise::T
end

@inline function (s::RampStimulus)(x, t)
    τ = mod(t, s.period)
    return τ < s.rise ? s.amplitude * (τ / s.rise) : zero(s.amplitude)
end

r = RampStimulus(-53.0, 1000.0, 2.0)
([r(nothing, t) for t in (0.0, 1.0, 2.0, 3.0)], isbitstype(typeof(r)))
```

Stimulus types are owned by CytoZoo and re-exported by model packages that adhere to the
interface natively, so the same stimulus objects work across every model you load.

## See Also

- [Spatial Heterogeneity](spatial.md) — the `x` that stimuli receive, and per-cell parameters.
- [Rush-Larsen Integration](rush_larsen.md) — why an isbits stimulus matters there.
