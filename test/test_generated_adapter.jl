# The generated adapter, exercised against the real interface. The fixtures are
# committed because this repo's CI has no Python.
#
# They come to roughly 5,200 lines between them, which is not carelessness: about
# 2,000 of those are the plain fixture's string-keyed `parameter_index`/
# `state_index`/`monitor_index` chains -- one `if name == "..."` branch per name,
# which `ode2julia` always emits and no test here calls. They are kept because the
# plain fixture has to stay byte-for-byte what `ode2julia` produces, which is the
# whole point of comparing against it.
#
# Regenerate both, from this directory:
#
#   /workspaces/main_repo/.venv/bin/gotranx ode2cytozoo \
#       /workspaces/main_repo/third-party/gotranx/tests/odefiles/ORdmm_Land.ode \
#       -o test/fixtures/ordmm_land_cytozoo --model-name ORdmmLand
#   /workspaces/main_repo/.venv/bin/gotranx ode2julia \
#       /workspaces/main_repo/third-party/gotranx/tests/odefiles/ORdmm_Land.ode \
#       -o test/fixtures/ordmm_land_plain
#
# The CytoZoo fixture carries its own generator-emitted banner; the plain one
# carries none, because it must match `ode2julia` output exactly.

module GeneratedCytoZoo
include("fixtures/ordmm_land_cytozoo.jl")
end

module GeneratedPlain
include("fixtures/ordmm_land_plain.jl")
end

@testset "Generated CytoZoo adapter" begin
    model = GeneratedCytoZoo.ORdmmLand()
    check_interface_conformance(model; name = "ORdmmLand (generated)", expect_monitors = true)

    @testset "Unknown names raise rather than returning -1" begin
        # gotranx's own state_index returns -1 for an unknown name, which would
        # silently become u[-1]. The adapter must raise at the point of the mistake.
        @test_throws KeyError state_index(model, :definitely_not_a_state)
        @test_throws KeyError parameter_index(model, :definitely_not_a_parameter)
    end

    @testset "Float32 end to end" begin
        m32 = GeneratedCytoZoo.ORdmmLand(Float32)
        @test eltype(m32.parameters) == Float32
        u32 = default_initial_state(m32)
        @test eltype(u32) == Float32
        du32 = similar(u32)
        m32(du32, u32, nothing, 0.0f0)
        @test eltype(du32) == Float32
        @test all(isfinite, du32)
    end

    @testset "RHS agrees bitwise with ode2julia" begin
        n = num_states(model)
        u = default_initial_state(model)
        p = copy(model.parameters)
        du_cz = zeros(n)
        du_plain = zeros(n)
        model(du_cz, u, nothing, 0.0)
        GeneratedPlain.rhs(0.0, u, p, du_plain)
        @test du_cz == du_plain          # bitwise, not approximate
    end

    @testset "Spatial override takes effect" begin
        u = default_initial_state(model)
        base = zeros(num_states(model)); model(base, u, nothing, 0.0)

        # Not parameter_names(model)[1]: for ORdmm_Land that is :Aff, which only
        # reaches dv/dt through ICaL/ICaK/ICaNa, each gated by the L-type Ca
        # activation state `d` — exactly 0.0 at the default (resting) initial
        # state, so every term built on it is bitwise 0.0 regardless of Aff.
        # Overriding Aff there is a true no-op, not a bug in the override
        # machinery (verified by scanning all 139 parameters: most, including
        # :F used here, do change du). Faraday's constant feeds Nernst
        # potentials used throughout the model and is never zero-gated.
        name = :F
        i = parameter_index(model, name)
        sc = SpatialContext([0.0, 0.0, 0.0], NamedTuple{(name,)}((model.parameters[i] * 2 + 1,)))
        varied = zeros(num_states(model)); model(varied, u, sc, 0.0)
        @test varied != base

        # p = nothing must be untouched by the override machinery.
        again = zeros(num_states(model)); model(again, u, nothing, 0.0)
        @test again == base
    end
end
