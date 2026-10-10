"""Regenerate the split-toy fixtures of test/test_generated_split.jl.

    python test/fixtures/generate_split_toy.py

`gotranx ode2cytozoo` cannot split a model, so the three adapters are generated
here: the whole of split_toy.ode, and its two halves, split at "troponin".
Needs a gotranx whose ode2cytozoo supports missing variables (newer than 2.4.0).
"""

from pathlib import Path

import gotranx
from gotranx.cli.gotran2cytozoo import get_code

here = Path(__file__).parent
source = here / "split_toy.ode"
ode = gotranx.load_ode(source)
troponin = ode.get_component("troponin")
adapters = {
    "split_toy_cytozoo.jl": get_code(ode, model_name="SplitToy", source=source),
    "split_toy_membrane_cytozoo.jl": get_code(
        ode - troponin, model_name="SplitToyMembrane", source=source
    ),
    # Only a coupling's first component supplies its transmembrane potential, and
    # the troponin half has none, so it points at its only state.
    "split_toy_troponin_cytozoo.jl": get_code(
        troponin.to_ode(), model_name="SplitToyTroponin", v_name="b", source=source
    ),
}
for name, code in adapters.items():
    (here / name).write_text(code)
