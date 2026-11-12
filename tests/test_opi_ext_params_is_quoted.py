"""ORCA reads ``Ext_Params`` only as a quoted string.

Manual 4.1.10 (structurereactivity/optimizations.txt, "Optimizing with
External Methods") writes ``Ext_Params "optional command line arguments"``.
Measured with ORCA 6.1.1 on an input written by the OPI path: the
unquoted line ``Ext_Params --device cpu`` is refused by the
parser ("Invalid assignment in METHOD block ... Last token : -"); the
quoted line is accepted and ORCA calls ``<ProgExt> <base>_EXT.extinp.tmp
--device cpu``.  The plain-text fallback always quoted it;
the OPI path must write the same line.

Input: water through both ExtOpt builders (orca_extopt, orca_mlips).
Output: exactly one ``Ext_Params "<args>"`` line.  Runs with OPI only.
"""

from __future__ import annotations

import re

import pytest

pytest.importorskip("opi")

from delfin.mlp_tools.orca_extopt import build_extopt_input  # noqa: E402
from delfin.mlp_tools.orca_mlips import build_orca_mlips_input  # noqa: E402

WATER = (["O", "H", "H"], [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]])
_LINE = re.compile(r"(?im)^\s*Ext_Params\s+(.*)$")


def _ext_params_values(text: str) -> list[str]:
    return [m.group(1).strip() for m in _LINE.finditer(text)]


def test_orca_mlips_writes_ext_params_quoted(tmp_path):
    inp = build_orca_mlips_input("opt", elements=WATER[0], coordinates=WATER[1],
                                 work_dir=tmp_path)
    assert _ext_params_values(inp.read_text()) == ['"--device cpu"']


def test_delfin_wrapper_writes_ext_params_quoted(tmp_path):
    inp = build_extopt_input("opt", elements=WATER[0], coordinates=WATER[1],
                             work_dir=tmp_path)
    assert _ext_params_values(inp.read_text()) == ['"--backend ani2x --device cpu"']
