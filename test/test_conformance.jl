# A model-agnostic conformance harness. Any AbstractCellModel should pass it;
# test_interface.jl keeps the ToRORd-specific value assertions.

function check_interface_conformance(model; name::String, expect_monitors::Bool = false)
    @testset "conformance — $name" begin
        n = num_states(model)
        @test n isa Int && n > 0

        np = num_parameters(model)
        # Not `np > 0`: a model whose kinetics carry no tunable constant is a
        # legitimate `AbstractCellModel`, and nothing in the interface requires
        # one. Only the agreement between `np` and the name/value vectors below
        # is part of the contract.
        @test np isa Int && np >= 0

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

        # The element type the model was constructed with must survive to the
        # state vector and through the right-hand side. A model that silently
        # widens to Float64 defeats the whole point of a Float32 (or Dual, or
        # GPU) run: every cell in a tissue solve pays for it.
        #
        # `M` is the model's own UnionAll, so this stays model-agnostic -- the
        # `M(::Type{T})` constructor is part of the pattern every model here
        # follows (ToRORd, FHNModel and the gotranx-generated adapters alike).
        M = Base.typename(typeof(model)).wrapper
        @testset "element type round-trip" begin
            for T in (Float64, Float32)
                @test eltype(default_initial_state(M(T))) == T
            end

            m32 = M(Float32)
            u32 = default_initial_state(m32)
            @test eltype(u32) == Float32

            # Poisoned, as above: a slot the Float32 path forgets stays NaN.
            du32 = fill(Float32(NaN), num_states(m32))
            m32(du32, copy(u32), nothing, 0.0f0)
            @test eltype(du32) == Float32
            @test !any(isnan, du32)
            @test all(isfinite, du32)
        end

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
