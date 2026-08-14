# `:v` for the transmembrane potential, matching `ToRORd`; Thunderbolt writes the same
# variable `φₘ`, and tissue frameworks that prefer that spelling alias it to
# `transmembrane_potential_index` rather than renaming the state here.
const FHN_STATE_NAMES = (:v, :s)

const FHN_NUM_STATES = length(FHN_STATE_NAMES)

# Both states rest at zero: v = 0 is the stable node of the cubic when 0 < a < 1, and
# s = 0 is its matching recovery value (`d = 0` by default). Unlike `ToRORd`, this is a
# genuine steady state — no pre-pacing needed.
const FHN_DEFAULT_INITIAL_STATE = (0.0, 0.0)
