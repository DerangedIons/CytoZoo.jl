# A model that gotranx splits in two, coupled back together with `connect`, must be
# the model it was split from. Subtracting a component from a gotranx ODE turns
# every name the two halves share into a *missing variable* of the half that reads
# it; the generated adapter keeps those as trailing writable parameters, NaN until
# something stages them, so a `connect` edge supplies them by name.
#
# split_toy.ode has the shape of an excitation-contraction split: the troponin half
# reads the membrane's calcium `c` (a state), and the membrane reads troponin's
# binding flux `J` (a monitor that depends on the `c` it was just sent, so it needs
# the feedthrough ordering of the monitor pre-pass). The fixtures are committed
# because this repo's CI has no Python. Regenerate all three (the whole model and
# its two halves), from the repository root, with gotranx 2.5.0 or newer:
#
#   gotranx ode2cytozoo test/fixtures/split_toy.ode \
#       -o test/fixtures/split_toy_cytozoo --model-name SplitToy
#   gotranx ode2cytozoo test/fixtures/split_toy.ode \
#       -o test/fixtures/split_toy_membrane_cytozoo --model-name SplitToyMembrane \
#       --remove-component troponin
#   gotranx ode2cytozoo test/fixtures/split_toy.ode \
#       -o test/fixtures/split_toy_troponin_cytozoo --model-name SplitToyTroponin \
#       --keep-component troponin --v-name b
#
# Only a coupling's first component supplies its transmembrane potential, and the
# troponin half has none, so `--v-name` points at its only state.

using OrdinaryDiffEq

module GeneratedSplitToy
include("fixtures/split_toy_cytozoo.jl")
end

module GeneratedSplitToyMembrane
include("fixtures/split_toy_membrane_cytozoo.jl")
end

module GeneratedSplitToyTroponin
include("fixtures/split_toy_troponin_cytozoo.jl")
end

function _split_toy_coupled()
    return couple(
        [
            Subsystem(GeneratedSplitToyMembrane.SplitToyMembrane(); name = :M),
            Subsystem(GeneratedSplitToyTroponin.SplitToyTroponin(); name = :Tn),
        ],
        [connect(:M => :c, :Tn => :c), connect(:Tn => :J, :M => :J)],
    )
end

@testset "Generated split adapters" begin
    membrane = GeneratedSplitToyMembrane.SplitToyMembrane()
    troponin = GeneratedSplitToyTroponin.SplitToyTroponin()
    full = GeneratedSplitToy.SplitToy()

    @testset "Inputs are trailing parameters, by name" begin
        @test parameter_names(membrane)[end] === :J
        @test parameter_index(membrane, :J) == num_parameters(membrane)
        @test parameter_names(troponin)[end] === :c
        @test parameter_index(troponin, :c) == num_parameters(troponin)
        # The inputs are counted, so the vector `connect` writes into is all of them.
        @test length(writable_parameters(membrane)) == num_parameters(membrane)
        # ...and an input is not a parameter of the unsplit model.
        @test num_parameters(membrane) + num_parameters(troponin) == num_parameters(full) + 2
    end

    @testset "An unconnected input poisons the RHS" begin
        for T in (Float64, Float32)
            m = GeneratedSplitToyTroponin.SplitToyTroponin(T)
            @test isnan(m.parameters[parameter_index(m, :c)])
            u = default_initial_state(m)
            du = zeros(T, num_states(m))
            m(du, u, nothing, zero(T))
            @test eltype(du) == T
            @test all(isnan, du)
        end
    end

    @testset "An input computes in the state's element type, like a parameter" begin
        # A half built in Float64 and evaluated on Float32 states converts its
        # parameters to Float32; an input it left in Float64 would widen every term
        # that reads it, and the result would no longer be the Float32 model's.
        for (M, input) in (
                (GeneratedSplitToyMembrane.SplitToyMembrane, :J),
                (GeneratedSplitToyTroponin.SplitToyTroponin, :c),
            )
            m64 = M()
            m32 = M(Float32)
            m64.parameters[parameter_index(m64, input)] = 0.1
            m32.parameters[parameter_index(m32, input)] = 0.1
            u = default_initial_state(m32)
            du64 = zeros(Float32, num_states(m32))
            du32 = zeros(Float32, num_states(m32))
            m64(du64, u, nothing, 0.0f0)
            m32(du32, u, nothing, 0.0f0)
            @test du64 == du32
        end
    end

    check_interface_conformance(
        membrane; name = "split toy membrane (generated)", expect_monitors = true,
        inputs = (J = 0.1,),
    )
    check_interface_conformance(
        troponin; name = "split toy troponin (generated)", expect_monitors = true,
        inputs = (c = 0.4,),
    )

    @testset "Coupled halves are the unsplit model" begin
        cm = _split_toy_coupled()
        # Each of the unsplit model's states at its slot in the coupled vector. The
        # first component's states keep bare names; the second's are prefixed.
        slot(s) = s in state_names(membrane) ? state_index(cm, s) : state_index(cm, Symbol(:Tn_, s))
        perm = [slot(s) for s in state_names(full)]
        @test sort(perm) == 1:num_states(full)
        @test default_initial_state(cm)[perm] == default_initial_state(full)

        # Bitwise, not approximately: each half prints the expressions of the whole,
        # with an input standing where the whole computes the value inline. At the
        # initial state and away from it, so the test cannot pass on the initial
        # condition alone.
        for u_full in (default_initial_state(full), [0.9, -0.3, 1.7])
            U = similar(u_full)
            U[perm] = u_full
            dU = fill(NaN, length(U))
            cm(dU, U, nothing, 0.0)
            du_full = fill(NaN, length(u_full))
            full(du_full, u_full, nothing, 0.0)
            @test dU[perm] == du_full
        end

        # And through a solve, where every evaluation stages afresh.
        tspan = (0.0, 5.0)
        sol = solve(ODEProblem(cm, tspan), Tsit5(); dt = 0.01, adaptive = false)
        sol_full = solve(ODEProblem(full, tspan), Tsit5(); dt = 0.01, adaptive = false)
        @test sol.u[end][perm] == sol_full.u[end]
        @test sol.u[end] != default_initial_state(cm)

        # A stiff solver evaluates the halves on Dual states, to which the inputs are
        # converted like the parameters. Its Jacobian omits the coupling terms, so the
        # trajectory agrees only to the tolerances. FBDF, not a Rosenbrock method,
        # which takes hundreds of times the steps here (see the `connect`
        # limitations).
        sol_ad = solve(ODEProblem(cm, tspan), FBDF(); reltol = 1.0e-10, abstol = 1.0e-12)
        sol_ad_full = solve(ODEProblem(full, tspan), FBDF(); reltol = 1.0e-10, abstol = 1.0e-12)
        @test sol_ad.retcode == ReturnCode.Success
        @test sol_ad.u[end][perm] ≈ sol_ad_full.u[end] rtol = 1.0e-7
    end
end
