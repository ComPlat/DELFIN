"""The TS -> GOAT plain-text constraints equal OPI's typed rendering.

``delfin.ts_goat.render_ts_goat_input`` writes the input without OPI, so
the dashboard works where the ``delfin[orca]`` extra is absent. This
pins the text to OPI's own model: each ``{ C n C }`` line must be what
``opi.input.blocks.block_geom.Constraint(mode="C", atom1=n)`` prints,
and the set must be accepted by ``BlockGeom``. Runs in the OPI CI job
(``tests/test_opi_*.py``, DELFIN_OPI_REQUIRED=1).

OPI's ``Constraint.from_string`` is not used: in orca-pi 2.0.0 it raises
"String not valid" for ``"{ C 12 C }"`` (its pattern requires a second
atom index), so the typed constructor is the only OPI route.
"""

import re

import pytest

pytest.importorskip("opi")

from opi.input.blocks.block_geom import (  # noqa: E402
    BlockGeom,
    Constraint,
    Constraints,
)

from delfin.ts_goat import render_ts_goat_input  # noqa: E402


def _norm(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def test_constraint_lines_match_opi():
    indices = [0, 3, 12]
    elements = ["C"] * 14
    coordinates = [(float(i), 0.0, 0.0) for i in range(14)]
    text = render_ts_goat_input(elements, coordinates, indices)
    ours = [_norm(line) for line in text.splitlines()
            if line.strip().startswith("{")]
    typed = [Constraint(mode="C", atom1=i) for i in indices]
    assert ours == [_norm(str(c)) for c in typed]
    block = BlockGeom(constraints=Constraints(constraints=typed))
    assert block.constraints is not None
    assert len(block.constraints.constraints) == len(indices)
