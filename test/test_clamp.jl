# State-clamp tests: the ClampedCell wrapper (freeze `du`, seed `u0`), the name-keyed
# constructor and `seed!`, interface forwarding, the Rush-Larsen path, and how a clamp composes with
# `couple`. Needs a solver for the "a held state does not move" trajectory test.

using OrdinaryDiffEq

# Mock with a parameter slot (:k, a connect target), a monitor, and a state whose derivative
# depends on the state below it, so freezing one state is visible in another's trajectory.
struct _ClampMock <: CytoZoo.AbstractCardiacCellModel
    parameters::Vector{Float64}
end
_ClampMock() = _ClampMock([2.0])
CytoZoo.num_states(::_ClampMock) = 3
CytoZoo.num_parameters(::_ClampMock) = 1
CytoZoo.state_names(::_ClampMock) = (:v, :na, :ca)
CytoZoo.parameter_names(::_ClampMock) = (:k,)
CytoZoo.default_initial_state(::_ClampMock) = [-80.0, 10.0, 0.1]
CytoZoo.state_index(::_ClampMock, n::Symbol) = findfirst(==(n), (:v, :na, :ca))
CytoZoo.parameter_index(::_ClampMock, n::Symbol) = n === :k ? 1 : nothing
CytoZoo.transmembrane_potential_index(::_ClampMock) = 1
function (m::_ClampMock)(du, u, p, t)
    du[1] = -0.5 * (u[1] + 80.0)
    du[2] = -u[2]
    du[3] = m.parameters[1] * u[2] - u[3]
    return nothing
end
CytoZoo.num_monitors(::_ClampMock) = 1
CytoZoo.monitor_names(::_ClampMock) = (:total,)
function CytoZoo.monitor_values!(mon, u, t, ::_ClampMock)
    mon[1] = u[2] + u[3]
    return nothing
end

@testset "ClampedCell" begin
    @testset "the functor zeroes only the held derivatives" begin
        m = _ClampMock()
        u = [-70.0, 10.0, 0.1]
        du_base = similar(u)
        m(du_base, u, nothing, 0.0)

        c = ClampedCell(m, (2,))
        du = fill(1.0e6, 3)             # du arrives dirty
        c(du, u, nothing, 0.0)

        @test du[2] == 0.0
        @test du[1] == du_base[1]
        @test du[3] == du_base[3]       # still sees the held state's VALUE, only its rate is gone
    end

    @testset "several states at once" begin
        c = ClampedCell(_ClampMock(), (1, 3))
        du = fill(1.0e6, 3)
        c(du, [-70.0, 10.0, 0.1], nothing, 0.0)
        @test du[1] == 0.0
        @test du[3] == 0.0
        @test du[2] == -10.0
    end

    @testset "zero allocations (functor)" begin
        c = ClampedCell(_ClampMock(), (2,))
        u = default_initial_state(c)
        du = similar(u)
        c(du, u, nothing, 0.0)          # warmup
        @test (@allocated c(du, u, nothing, 0.0)) == 0
    end

    @testset "the seed and the hold come from one object" begin
        m = _ClampMock()

        # Two-argument form: held at the base model's own initial value.
        @test default_initial_state(ClampedCell(m, (2,))) == default_initial_state(m)

        # Explicit values: the wrapper seeds them, so the level held is the level started from.
        c = ClampedCell(m, (2,), (20.0,))
        @test default_initial_state(c) == [-80.0, 20.0, 0.1]
        @test c.indices == (2,)
        @test c.values == (20.0,)
    end

    @testset "the keyword constructor resolves names and seeds the state vector" begin
        m = _ClampMock()
        c = ClampedCell(m; na = 20.0)

        @test c isa ClampedCell
        @test c.indices == (2,)
        @test c.values == (20.0,)
        u = default_initial_state(c)
        @test u == [-80.0, 20.0, 0.1]                # the wrapper alone supplies the seed

        du = similar(u)
        c(du, u, nothing, 0.0)
        @test du[2] == 0.0

        # Several names, order independent.
        c2 = ClampedCell(m; ca = 0.5, v = -20.0)
        @test sort(collect(c2.indices)) == [1, 3]
        @test default_initial_state(c2) == [-20.0, 10.0, 0.5]
    end

    @testset "seed! continues a protocol from a handed-in state" begin
        m = _ClampMock()
        u_prev = [-12.0, 20.0, 7.0]
        c = ClampedCell(m; na = 7.5)
        u = seed!(copy(u_prev), c)

        @test u == [-12.0, 7.5, 7.0]                 # only the clamped slot is re-seeded
        @test u_prev == [-12.0, 20.0, 7.0]           # copied first, so the caller's vector is untouched
        @test default_initial_state(c) == [-80.0, 7.5, 0.1]  # base IC, same hold

        # In place, returning its argument.
        w = zeros(3)
        @test seed!(w, c) === w
        @test w == [0.0, 7.5, 0.0]

        # Nested clamps seed every level; a bare model seeds nothing.
        nested = ClampedCell(c, (1,), (-20.0,))
        @test seed!(zeros(3), nested) == [-20.0, 7.5, 0.0]
        @test seed!([1.0, 2.0, 3.0], m) == [1.0, 2.0, 3.0]
    end

    @testset "the keyword constructor rejects unknown names and empty clamps" begin
        m = _ClampMock()
        @test_throws ArgumentError ClampedCell(m; nope = 1.0)
        @test_throws ArgumentError ClampedCell(m)
        err = try
            ClampedCell(m; nope = 1.0)
        catch e
            e
        end
        @test occursin("nope", err.msg) && occursin("(:v, :na, :ca)", err.msg)
    end

    @testset "the model's element type survives the clamp" begin
        c = ClampedCell(FHNModel(Float32); v = 0.5)
        @test c.values === (0.5f0,)
        @test eltype(default_initial_state(c)) === Float32
    end

    @testset "invalid indices are rejected at construction" begin
        m = _ClampMock()
        @test_throws ArgumentError ClampedCell(m, (2, 2))
        @test_throws ArgumentError ClampedCell(m, (4,))
        @test_throws ArgumentError ClampedCell(m, (0,))
    end

    @testset "index containers" begin
        m = _ClampMock()
        @test ClampedCell(m, 2).indices == (2,)
        @test ClampedCell(m, [1, 3]).indices == (1, 3)
        @test ClampedCell(m, [2], [20.0]).values == (20.0,)
    end

    @testset "the base interface is forwarded" begin
        m = _ClampMock()
        c = ClampedCell(m, (2,))

        @test num_states(c) == 3
        @test num_parameters(c) == 1
        @test transmembrane_potential_index(c) == 1
        @test state_names(c) == (:v, :na, :ca)
        @test parameter_names(c) == (:k,)
        @test state_index(c, :ca) == 3
        @test state_index(c, :nope) === nothing
        @test parameter_index(c, :k) == 1
        @test writable_parameters(c) === m.parameters

        @test num_monitors(c) == 1
        @test monitor_names(c) == (:total,)
        mon = zeros(1)
        monitor_values!(mon, [-80.0, 10.0, 0.1], 0.0, c)
        @test mon[1] == 10.1
    end

    @testset "base_model unwraps, however deep" begin
        m = _ClampMock()
        @test base_model(m) === m
        @test base_model(ClampedCell(m, (2,))) === m
        @test base_model(ClampedCell(ClampedCell(m, (2,)), (1,))) === m
    end

    @testset "nested clamps hold both states" begin
        m = _ClampMock()
        c = ClampedCell(ClampedCell(m, (2,), (20.0,)), (1,), (-20.0,))
        @test default_initial_state(c) == [-20.0, 20.0, 0.1]

        du = fill(1.0e6, 3)
        c(du, default_initial_state(c), nothing, 0.0)
        @test du[1] == 0.0
        @test du[2] == 0.0
        @test du[3] != 0.0
    end

    @testset "a held state does not move over a solve" begin
        m = _ClampMock()
        c = ClampedCell(m; na = 20.0)
        u = seed!([-60.0, 10.0, 0.1], c)
        sol = solve(ODEProblem(c, u, (0.0, 5.0), nothing), Tsit5(); reltol = 1.0e-10, abstol = 1.0e-10)
        vend, naend, caend = sol.u[end]

        @test all(s[2] == 20.0 for s in sol.u)          # held exactly, at every saved step
        @test naend == 20.0
        @test vend ≈ -80.0 + 20.0 * exp(-2.5) atol = 1.0e-6      # v relaxes, unaffected
        # ca is driven by the held na: dca/dt = k·na − ca, k = 2 ⇒ ca(∞) = 40.
        @test caend ≈ 40.0 - 39.9 * exp(-5.0) atol = 1.0e-6

        # Unclamped, na decays and ca follows it back down — the contrast the clamp buys.
        free = solve(ODEProblem(m, u, (0.0, 5.0), nothing), Tsit5(); reltol = 1.0e-10, abstol = 1.0e-10)
        @test free.u[end][2] ≈ 20.0 * exp(-5.0) atol = 1.0e-6
        @test free.u[end][3] ≈ exp(-5.0) * (0.1 + 40.0 * 5.0) atol = 1.0e-6
    end

    @testset "monitors survive the wrapper post-solve" begin
        c = ClampedCell(_ClampMock(); na = 20.0)
        sol = solve(ODEProblem(c, default_initial_state(c), (0.0, 5.0), nothing), Tsit5())
        h = monitor_history(sol, c)
        @test h.names == (:total,)
        @test h.values[1, end] ≈ sol.u[end][2] + sol.u[end][3]
    end

    @testset "Rush-Larsen holds the clamped state too" begin
        base = ToRORd()
        c = ClampedCell(base, (1,))                             # voltage clamp
        @test has_rush_larsen(c) == has_rush_larsen(base) == true

        u = default_initial_state(base)
        u_new = similar(u)
        u_new_base = similar(u)
        rush_larsen_step!(u_new, u, nothing, 0.0, 0.01, c)
        rush_larsen_step!(u_new_base, u, nothing, 0.0, 0.01, base)

        @test u_new[1] == u[1]                                  # Vm held, not stepped
        @test u_new_base[1] != u[1]                             # the base step does move it
        @test u_new[2:end] == u_new_base[2:end]                 # nothing else changed
    end
end

@testset "ClampedCell — composition with couple" begin
    @testset "a clamped component clamps its own states, in component-local indices" begin
        clamped_a = ClampedCell(_ClampMock(); na = 20.0)
        cm = couple([Subsystem(clamped_a; name = :A), Subsystem(_ClampMock(); name = :B)])

        # The component's seeded IC flows into the coupled layout — no hand-patching of U0.
        U = default_initial_state(cm)
        @test U[state_index(cm, :na)] == 20.0
        @test U[state_index(cm, :B_na)] == 10.0

        dU = fill(1.0e6, num_states(cm))
        cm(dU, U, nothing, 0.0)
        @test dU[state_index(cm, :na)] == 0.0
        @test dU[state_index(cm, :B_na)] == -10.0               # the other component is untouched
    end

    @testset "under a contributory share, clamp the coupling and not the component" begin
        # J contributes +20 into P's shared v slot; P's own equation drives it as well.
        clamped_j = ClampedCell(_MonoJ(), (1,))
        cm = couple(
            [Subsystem(_MonoP(); name = :P), Subsystem(clamped_j; name = :J)],
            [share(:P => :v, :J => :v; owner = :P, op = +)],
        )
        vslot = state_index(cm, :v)
        U = default_initial_state(cm)
        dU = similar(U)
        cm(dU, U, nothing, 0.0)

        # Clamping the contributor removes ITS contribution (+20) and nothing else — the owner's
        # term still drives the slot, so the state is not held.
        @test dU[vslot] != 0.0
        @test dU[vslot] == -_MONO_K * (U[vslot] - U[state_index(cm, :a)])

        # Wrapping the coupling holds it.
        outer = ClampedCell(cm; v = 10.0)
        U_outer = default_initial_state(outer)
        dU_outer = similar(U_outer)
        outer(dU_outer, U_outer, nothing, 0.0)
        @test dU_outer[vslot] == 0.0
    end

    @testset "a clamped receiver still takes a connect edge" begin
        clamped_p = ClampedCell(_ClampMock(); v = -20.0)
        cm = couple(
            [Subsystem(_MonoA(); name = :S), Subsystem(clamped_p; name = :R)],
            [connect(:S => :d, :R => :k)],
        )
        U = default_initial_state(cm)
        dU = similar(U)
        cm(dU, U, nothing, 0.0)

        # The staged input reached the wrapped model's parameter vector: dca/dt = d·na − ca.
        d = U[state_index(cm, :d)]
        @test dU[state_index(cm, :R_ca)] ≈ d * U[state_index(cm, :R_na)] - U[state_index(cm, :R_ca)]
        @test dU[state_index(cm, :R_v)] == 0.0
    end

    @testset "a clamp does not fake a parameter vector the base lacks" begin
        # FHNModel keeps its parameters in struct fields, so it cannot receive a connect edge —
        # and wrapping it in a clamp must not make `couple` believe otherwise.
        clamped_fhn = ClampedCell(FHNModel(), (1,))
        @test_throws ArgumentError couple(
            [Subsystem(_MonoA(); name = :S), Subsystem(clamped_fhn; name = :F)],
            [connect(:S => :d, :F => :a)],
        )
    end
end
