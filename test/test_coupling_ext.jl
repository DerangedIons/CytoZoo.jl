# CouplingExt solving tests. Included only when OrdinaryDiffEqOperatorSplitting is
# resolvable (gated in runtests.jl, like the TWorld tests); skipped in CI where OS is absent.

using OrdinaryDiffEqOperatorSplitting
using OrdinaryDiffEq

# Toy component A: state `d` decays at rate 1, `v` is constant. Analytic: d(t) = d0·exp(-t).
struct _ExtToyA <: CytoZoo.AbstractCellModel end
CytoZoo.num_states(::_ExtToyA) = 2
CytoZoo.state_names(::_ExtToyA) = (:v, :d)
CytoZoo.default_initial_state(::_ExtToyA) = [0.0, 1.0]
CytoZoo.state_index(::_ExtToyA, n::Symbol) = findfirst(==(n), (:v, :d))
CytoZoo.transmembrane_potential_index(::_ExtToyA) = 1
function (::_ExtToyA)(du, u, p, t)
    du[1] = 0.0
    du[2] = -u[2]
    return nothing
end

# Toy component B: inert (no dynamics). Used to check the shared slot is driven by the
# alias owner's equation, not the non-owner's.
struct _ExtToyB <: CytoZoo.AbstractCellModel end
CytoZoo.num_states(::_ExtToyB) = 2
CytoZoo.state_names(::_ExtToyB) = (:x, :y)
CytoZoo.default_initial_state(::_ExtToyB) = [5.0, 7.0]
CytoZoo.state_index(::_ExtToyB, n::Symbol) = findfirst(==(n), (:x, :y))
CytoZoo.transmembrane_potential_index(::_ExtToyB) = 1
function (::_ExtToyB)(du, u, p, t)
    du[1] = 0.0
    du[2] = 0.0
    return nothing
end

@testset "CouplingExt — single component matches analytic" begin
    cm = couple((A = _ExtToyA(),))
    prob = OperatorSplittingProblem(cm, (0.0, 1.0))
    integ = init(prob, coupled_algorithm(cm, Tsit5()); dt = 0.01, adaptive = false)
    solve!(integ)
    @test integ.u[state_index(cm, :v)] ≈ 0.0 atol = 1e-8
    @test integ.u[state_index(cm, :d)] ≈ exp(-1.0) rtol = 1e-4
end

@testset "CouplingExt — alias shares a slot, driven by owner dynamics" begin
    cm = couple((A = _ExtToyA(), B = _ExtToyB()); aliases = [alias(:A => :d, :B => :x; owner = :A)])
    @test num_states(cm) == 3
    @test default_initial_state(cm) == [0.0, 1.0, 7.0]      # d from owner A; B's y appended
    prob = OperatorSplittingProblem(cm, (0.0, 1.0))
    integ = init(prob, coupled_algorithm(cm, Tsit5()); dt = 0.01, adaptive = false)
    solve!(integ)
    @test integ.u[state_index(cm, :d)] ≈ exp(-1.0) rtol = 1e-4   # owner A decays the shared slot
    @test integ.u[state_index(cm, :B_y)] ≈ 7.0 atol = 1e-8       # B inert
end
