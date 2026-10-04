# A model-agnostic conformance harness. Any AbstractCellModel should pass it;
# test_interface.jl keeps the ToRORd-specific value assertions.

function check_interface_conformance(model; name::String, expect_monitors::Bool = false)
    @testset "conformance — $name" begin
        n = num_states(model)
        @test n isa Int && n > 0

        np = num_parameters(model)
        @test np isa Int && np > 0

        vi = transmembrane_potential_index(model)
        @test vi isa Int && 1 <= vi <= n

        u0 = default_initial_state(model)
        @test length(u0) == n

        sn = state_names(model)
        @test length(sn) == n
        @test allunique(sn)
        for (i, s) in enumerate(sn)
            @test state_index(model, s) == i
        end

        pn = parameter_names(model)
        @test length(pn) == np
        @test allunique(pn)
        for (i, s) in enumerate(pn)
            @test parameter_index(model, s) == i
        end

        # `writable_parameters` is documented as *optional* interface (src/interface.jl):
        # its default reads `model.parameters`, and a model with no such field (e.g. FHN's
        # isbits struct fields) legitimately opts out, rejecting `connect` instead of
        # supporting it (see `_provides_writable_parameters` in src/coupling.jl, which the
        # connect-validation path itself uses to decide this). Only check the length when
        # the model actually provides it.
        if CytoZoo._provides_writable_parameters(model)
            @test length(writable_parameters(model)) == np
        end

        # The RHS must WRITE every du entry, never accumulate into an unassigned slot.
        # Poisoning du with NaN catches a slot the RHS forgets: it stays NaN.
        du = fill(NaN, n)
        model(du, copy(u0), nothing, 0.0)
        @test !any(isnan, du)

        if expect_monitors
            nm = num_monitors(model)
            @test nm > 0
            @test length(monitor_names(model)) == nm
            mon = zeros(eltype(u0), nm)
            monitor_values!(mon, copy(u0), 0.0, model)
            @test all(isfinite, mon)
        end
    end
end

@testset "Conformance harness on existing models" begin
    check_interface_conformance(ToRORd(); name = "ToRORd")
    check_interface_conformance(FHNModel(); name = "FHN")
end
