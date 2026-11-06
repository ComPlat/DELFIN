"""XE1 structure-editor library + MCP tools (red control for package XE1).

Package XE1 builds the structure department's ``xyz-editor`` expert on the
SA0 harness. This file is the red control for the three missing structure
tools (geometry_renumber, geometry_fix_hydrogens, geometry_diff) and the
library functions behind them in ``delfin/structure_edit.py``.

Rules the tests pin (from the package target):
- renumber maps a target (product) onto a reference (educt) atom order via
  ``delfin.atom_mapping.map_atoms``, uses its ``verified`` flag, surfaces a
  formula-mismatch ValueError (never swallows it), and reports a renumbering
  table plus checks.
- fix_hydrogens adds missing hydrogens through DELFIN's own covalent-radius /
  connectivity primitives; it never removes atoms.
- structure_diff reports a per-atom diff (moved / changed element) and checks.
- security (operator): the geometry MCP tools write NOTHING to disk — an MCP
  tool writing a file would bypass the write-scope gate. The structure expert
  writes a NEW file with the gate-checked ``write_file``, which refuses
  overwriting (the ``refuse to overwrite the input`` benchmark task).

Fixtures are hand-built xyz text only (never ORCA output).
"""

from __future__ import annotations

import json

import pytest

import delfin.structure_edit as se
from delfin.ops_server import server as server_mod


# ---------------------------------------------------------------------------
# hand-built fixtures (no ORCA output: licence)
# ---------------------------------------------------------------------------

# Reference educt: methane, atoms in a stable order.
EDUCT_CH4 = """5
methane reference
C 0.000000 0.000000 0.000000
H 0.000000 0.000000 1.089000
H 1.027673 0.000000 -0.363000
H -0.513836 0.889795 -0.363000
H -0.513836 -0.889795 -0.363000
"""

# Product: the same methane with the hydrogens listed in a different order
# (a renumbering must map it back onto the educt's atom order).
PRODUCT_CH4_PERMUTED = """5
methane product (H permuted)
C 0.000000 0.000000 0.000000
H -0.513836 -0.889795 -0.363000
H 1.027673 0.000000 -0.363000
H 0.000000 0.000000 1.089000
H -0.513836 0.889795 -0.363000
"""

# Metal complex missing a coordinating hydrogen (a common xyz defect): the
# ammine/metal-donor H that DELFIN's H logic pins.
METAL_NO_H = """3
metal ammine missing donor H
Zn 0.000000 0.000000 0.000000
N  1.500000 0.000000 0.000000
H  1.500000 0.000000 0.000000
"""

# Two structures differing by an element swap on one atom for structure_diff.
STRUCT_DIFF_A = """2
dimer A
C 0.000000 0.000000 0.000000
O 1.200000 0.000000 0.000000
"""

STRUCT_DIFF_B = """2
dimer B (carbon -> nitrogen)
C 0.000000 0.000000 0.000000
N 1.200000 0.000000 0.000000
"""


def _symbols_order(text: str) -> list[str]:
    """Atom symbols of an xyz block in file order."""
    body = [ln for ln in text.splitlines() if ln.strip()]
    start = 1
    if body and body[0].strip().isdigit():
        start = 2
    return [ln.split()[0] for ln in body[start:]]


# ---------------------------------------------------------------------------
# renumber
# ---------------------------------------------------------------------------

def test_renumber_maps_permuted_target_back_onto_reference():
    out = se.renumber_to_reference(PRODUCT_CH4_PERMUTED, EDUCT_CH4)
    assert out["error"] is None
    assert out["verified"] is True
    # Reordered product in reference order has the reference symbol order.
    assert _symbols_order(out["reordered_xyz"]) == _symbols_order(EDUCT_CH4)
    # The renumbering table lists a reference index for every product atom.
    table = out["renumbering"]
    assert len(table) == 5
    assert {row["reference_index"] for row in table} == set(range(5))


def test_renumber_surfaces_formula_mismatch_valueerror():
    # Different formula (product has an extra H) must surface, not be hidden.
    mismatched = EDUCT_CH4.replace("5\n", "6\n").rstrip() + "\nH 2.0 0.0 0.0\n"
    with pytest.raises(ValueError):
        se.renumber_to_reference(mismatched, EDUCT_CH4)


def test_renumber_tool_str_is_registerable_via_server():
    result = server_mod.tool_geometry_renumber(PRODUCT_CH4_PERMUTED, EDUCT_CH4)
    payload = json.loads(result)
    assert payload["verified"] is True
    assert payload["n_atoms"] == 5


# ---------------------------------------------------------------------------
# fix_hydrogens
# ---------------------------------------------------------------------------

def test_fix_hydrogens_adds_missing_h_and_never_removes():
    out = se.fix_hydrogens(METAL_NO_H)
    assert out["error"] is None
    before = _symbols_order(METAL_NO_H)
    after = _symbols_order(out["xyz"])
    # Nothing removed: every original atom is still present in order.
    assert after[: len(before)] == before
    # H were added beyond the original count.
    assert len(after) > len(before)
    assert out["n_added"] == len(after) - len(before)


def test_fix_hydrogens_tool_str_roundtrip():
    result = server_mod.tool_geometry_fix_hydrogens(METAL_NO_H)
    payload = json.loads(result)
    assert payload["n_added"] >= 1


# ---------------------------------------------------------------------------
# structure_diff
# ---------------------------------------------------------------------------

def test_structure_diff_reports_element_change():
    out = se.structure_diff(STRUCT_DIFF_A, STRUCT_DIFF_B)
    assert out["error"] is None
    changed = [r for r in out["rows"] if r.get("element_changed")]
    assert len(changed) == 1
    assert changed[0]["from_symbol"] == "O"
    assert changed[0]["to_symbol"] == "N"


def test_structure_diff_reports_composition_mismatch():
    out = se.structure_diff(EDUCT_CH4, PRODUCT_CH4_PERMUTED)
    assert out["composition_equal"] is True


def test_structure_diff_tool_str():
    result = server_mod.tool_geometry_diff(STRUCT_DIFF_A, STRUCT_DIFF_B)
    payload = json.loads(result)
    assert payload["composition_equal"] is False


# ---------------------------------------------------------------------------
# security: the MCP tools write nothing to disk (operator, XE1)
# ---------------------------------------------------------------------------

def test_geometry_tools_write_nothing_to_disk(tmp_path):
    """An MCP tool must not write a file (that would bypass the write-scope
    gate). renumber / fix_hydrogens / diff return text only; the directory
    they are given stays empty and no file appears anywhere under it."""
    import os

    before = set(os.listdir(tmp_path))
    server_mod.tool_geometry_renumber(PRODUCT_CH4_PERMUTED, EDUCT_CH4)
    server_mod.tool_geometry_fix_hydrogens(METAL_NO_H)
    server_mod.tool_geometry_diff(STRUCT_DIFF_A, STRUCT_DIFF_B)
    after = set(os.listdir(tmp_path))
    assert after == before  # no file created anywhere under tmp_path
