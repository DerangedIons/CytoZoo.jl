const FHN_PARAMETER_NAMES = (:a, :b, :c, :d, :e)

const FHN_NUM_PARAMS = length(FHN_PARAMETER_NAMES)

# Nagumo's original excitability threshold and the linear recovery coefficients, in the
# scaling Thunderbolt's `ParametrizedFHNModel` uses. `e` is the recovery/excitation
# timescale ratio: at 0.01 the recovery variable is a hundred times slower than the
# upstroke, which is what makes the pulse shape action-potential-like.
const FHN_DEFAULT_PARAMETERS = (a = 0.1, b = 0.5, c = 1.0, d = 0.0, e = 0.01)
