# Base (pure-Julia) coupling tests: layout, naming, aliasing, owner semantics, validation.
# OS-based solving is exercised separately in the CouplingExt tests (require OS loaded).

# Minimal mock components implementing the participant interface.
struct _CplMockA <: CytoZoo.AbstractCellModel end
CytoZoo.num_states(::_CplMockA) = 3
CytoZoo.state_names(::_CplMockA) = (:b, :c, :d)
CytoZoo.default_initial_state(::_CplMockA) = [1.0, 2.0, 3.0]
CytoZoo.state_index(::_CplMockA, n::Symbol) = findfirst(==(n), (:b, :c, :d))
CytoZoo.transmembrane_potential_index(::_CplMockA) = 1

struct _CplMockB <: CytoZoo.AbstractCellModel end
CytoZoo.num_states(::_CplMockB) = 2
CytoZoo.state_names(::_CplMockB) = (:x, :y)
CytoZoo.default_initial_state(::_CplMockB) = [10.0, 20.0]
CytoZoo.state_index(::_CplMockB, n::Symbol) = findfirst(==(n), (:x, :y))
CytoZoo.transmembrane_potential_index(::_CplMockB) = 1

@testset "couple — single component" begin
    cm = couple((A = _CplMockA(),))
    @test cm isa CoupledModel
    @test num_states(cm) == 3
    @test state_names(cm) == (:b, :c, :d)
    @test cm.layout.solution_indices.A == [1, 2, 3]
    @test default_initial_state(cm) == [1.0, 2.0, 3.0]
    @test cm.layout.operator_order == [:A]
    @test transmembrane_potential_index(cm) == 1
end

@testset "couple — alias, owner = first component" begin
    cm = couple((A = _CplMockA(), B = _CplMockB()); aliases = [alias(:A => :d, :B => :x; owner = :A)])
    @test num_states(cm) == 4
    @test state_names(cm) == (:b, :c, :d, :B_y)          # owner A keeps :d; B's y prefixed
    @test cm.layout.solution_indices.A == [1, 2, 3]
    @test cm.layout.solution_indices.B == [3, 4]         # x aliases to slot 3, y is slot 4
    @test default_initial_state(cm) == [1.0, 2.0, 3.0, 20.0]  # slot 3 = owner A's d IC
    @test cm.layout.operator_order == [:B, :A]           # owner A steps last, wins the slot
    @test state_index(cm, :d) == 3
    @test state_index(cm, :B_y) == 4
    @test state_index(cm, :nope) === nothing
    @test transmembrane_potential_index(cm) == 1
end

@testset "couple — alias, owner = second component" begin
    cm = couple((A = _CplMockA(), B = _CplMockB()); aliases = [alias(:A => :d, :B => :x; owner = :B)])
    @test state_names(cm) == (:b, :c, :x, :B_y)          # canonical name from owner B
    @test cm.layout.solution_indices.A == [1, 2, 3]
    @test cm.layout.solution_indices.B == [3, 4]
    @test default_initial_state(cm) == [1.0, 2.0, 10.0, 20.0]  # slot 3 = owner B's x IC
    @test cm.layout.operator_order == [:A, :B]           # owner B steps last
end

@testset "couple — explicit canonical alias name" begin
    cm = couple((A = _CplMockA(), B = _CplMockB()); aliases = [alias(:A => :d, :B => :x; owner = :A, name = :shared)])
    @test state_names(cm) == (:b, :c, :shared, :B_y)
    @test state_index(cm, :shared) == 3
end

@testset "couple — validation and non-functor contract" begin
    @test_throws ArgumentError alias(:A => :d, :B => :x; owner = :C)                       # owner not a participant
    @test_throws ArgumentError couple((A = _CplMockA(), B = _CplMockB());
        aliases = [alias(:A => :nope, :B => :x; owner = :A)])                              # unknown state
    @test_throws ArgumentError couple((A = _CplMockA(), B = _CplMockB());
        aliases = [alias(:Z => :d, :B => :x; owner = :Z)])                                 # unknown component
    cm = couple((A = _CplMockA(),))
    @test_throws ArgumentError cm(nothing, nothing, nothing, 0.0)                          # not a single functor
    @test_throws ArgumentError num_parameters(cm)
end
