module ForwardDiffExt

using CytoZoo
import ForwardDiff

# Under an implicit solver, ForwardDiff threads Dual numbers through the global state U. A
# `connect` edge stages a source state into the receiver's (Float64) parameter slot during the
# monolithic RHS — storing a Dual there would error. Extract the primal so the connect input is
# frozen to its current value within the Newton step (semi-implicit: correct fixed point,
# approximate Jacobian). Recurses so nested Duals (higher-order AD) also collapse to the slot's
# eltype. `share` coupling is unaffected — it flows entirely through U.
CytoZoo._connect_value(x::ForwardDiff.Dual) = CytoZoo._connect_value(ForwardDiff.value(x))

# Comparisons in generated code decide on plain values: `Dual(0.0, 1.0) > 0` is `true` on
# ForwardDiff >= 1.0, which would send a guard like `(S > 0) ? A / S : 0` down the 0/0 branch.
CytoZoo.primal(x::ForwardDiff.Dual) = CytoZoo.primal(ForwardDiff.value(x))

end
