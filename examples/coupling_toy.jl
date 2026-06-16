# Minimal demonstration of CytoZoo's operator-splitting model coupling:
#   1. an alias  — two models share one state; the owner's equation governs it (hard-discard)
#   2. a cross-ref — a model reads another model's state through a parameter slot
#
# Run (needs OrdinaryDiffEqOperatorSplitting + an OrdinaryDiffEq solver):
#   julia --project examples/coupling_toy.jl

using CytoZoo
using OrdinaryDiffEqOperatorSplitting
using OrdinaryDiffEq

# --- two toy component models ---

# Model A owns the shared variable `d` (exponential decay) plus its own state `c`.
struct ModelA <: AbstractCardiacCellModel end
CytoZoo.num_states(::ModelA) = 2
CytoZoo.state_names(::ModelA) = (:c, :d)
CytoZoo.default_initial_state(::ModelA) = [5.0, 1.0]
CytoZoo.state_index(::ModelA, n::Symbol) = findfirst(==(n), (:c, :d))
CytoZoo.transmembrane_potential_index(::ModelA) = 2
function (::ModelA)(du, u, p, t)
    du[1] = -0.2 * u[1]        # c
    du[2] = -u[2]              # d — the shared variable, owned by A
    return nothing
end

# Model B has state `x` (≡ A's `d` when aliased) and `y` driven by it.
struct ModelB <: AbstractCardiacCellModel end
CytoZoo.num_states(::ModelB) = 2
CytoZoo.state_names(::ModelB) = (:x, :y)
CytoZoo.default_initial_state(::ModelB) = [0.0, 0.0]
CytoZoo.state_index(::ModelB, n::Symbol) = findfirst(==(n), (:x, :y))
CytoZoo.transmembrane_potential_index(::ModelB) = 1
function (::ModelB)(du, u, p, t)
    du[1] = 3.0               # x — DISCARDED under the alias (A owns the shared slot)
    du[2] = u[1] - u[2]       # y — driven by the shared value x
    return nothing
end

# --- 1. Alias: A.d ≡ B.x, owner = A ---------------------------------------------------------
coupled = couple((A = ModelA(), B = ModelB()); aliases = [alias(:A => :d, :B => :x; owner = :A)])

println("alias  → states: ", state_names(coupled), "  num_states: ", num_states(coupled))

integ = init(
    OperatorSplittingProblem(coupled, (0.0, 2.0)),
    coupled_algorithm(coupled, Tsit5());
    dt = 0.01, adaptive = false,
)
solve!(integ)
println("  d(2) = ", round(integ.u[state_index(coupled, :d)]; digits = 5),
    "   (owner A: exp(-2) = ", round(exp(-2.0); digits = 5), "; B's du(x)=3 discarded)")
println("  y(2) = ", round(integ.u[state_index(coupled, :B_y)]; digits = 5))

# --- 2. Cross-ref: a reader integrates A's `d`, read through a parameter slot ---------------
struct Reader <: AbstractCardiacCellModel
    parameters::Vector{Float64}     # slot 1 = :d_ext, written by the synchronizer each step
end
Reader() = Reader([0.0])
CytoZoo.num_states(::Reader) = 1
CytoZoo.state_names(::Reader) = (:acc,)
CytoZoo.default_initial_state(::Reader) = [0.0]
CytoZoo.state_index(::Reader, n::Symbol) = findfirst(==(n), (:acc,))
CytoZoo.transmembrane_potential_index(::Reader) = 1
CytoZoo.parameter_index(::Reader, n::Symbol) = n === :d_ext ? 1 : nothing
(m::Reader)(du, u, p, t) = (du[1] = m.parameters[1]; nothing)   # acc integrates A's d

coupled2 = couple((A = ModelA(), R = Reader()); refs = [crossref(:A => :d, :R => :d_ext)])
integ2 = init(
    OperatorSplittingProblem(coupled2, (0.0, 2.0)),
    coupled_algorithm(coupled2, Tsit5());
    dt = 0.01, adaptive = false,
)
solve!(integ2)
println("cross-ref → acc(2) = ", round(integ2.u[state_index(coupled2, :R_acc)]; digits = 5),
    "   (∫₀² exp(-t) dt = ", round(1 - exp(-2.0); digits = 5), ")")
