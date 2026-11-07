"""XE1 adversarial tests (reviewer nacht-s7) — attack the structure_edit
library and its MCP tools on their real call paths.

Points attacked (operator ruling + QS directives):
- (a) renumber surfaces map_atoms' ``verified`` when it is False;
- (b) a FORMULA/element mismatch (same count, different elements) surfaces as
      an error through the REAL tool path, never as a bogus zero mapping;
- (c) identical structures give a REAL diff: rows present, zero moved, green
      checks — not missing or empty;
- (d) per-atom diff is COMPLETE: every atom of both sides accounted for when
      the two structures have different lengths;
- chemistry: fix_hydrogens must not place a hydrogen on a terminal
      double-bonded oxygen (CO2) — a valence-only count is wrong there.

Fixtures are hand-built xyz text only (never ORCA output).
"""

from __future__ import annotations

import json

import pytest

import delfin.structure_edit as se
from delfin.ops_server import server as server_mod


# Hand-built fixture: carbon dioxide, terminal double-bonded O atoms.
CO2 = """3
carbon dioxide
C 0.000000 0.000000 0.000000
O 1.160000 0.000000 0.000000
O -1.160000 0.000000 0.000000
"""


def _symbols(text: str) -> list[str]:
    body = [ln for ln in text.splitlines() if ln.strip()]
    start = 2 if body and body[0].strip().isdigit() else 1
    return [ln.split()[0] for ln in body[start:]]


# ---------------------------------------------------------------------------
# (a) verified=False must be visible to the caller
# ---------------------------------------------------------------------------

def test_renumber_surfaces_unverified_mapping():
    """(a) A same-count/element-mismatch pair must either raise (surfaced
    from map_atoms) or return verified=False — never a silent verified=True.
    Evidence (auto-check 2026-11-06): map_atoms raises ValueError
    'different molecular formulae' (delfin/atom_mapping.py:404) — that is
    proper surfacing, so the raise path is asserted here."""
    ref = (
        "5\nref\nC 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n"
        "H 0.0 0.0 1.0\nH -1.0 0.0 0.0\n"
    )
    tgt = (
        "5\ntgt (F instead of H)\nC 0.0 0.0 0.0\nF 1.0 0.0 0.0\n"
        "H 0.0 1.0 0.0\nH 0.0 0.0 1.0\nH -1.0 0.0 0.0\n"
    )
    with pytest.raises(ValueError) as excinfo:
        se.renumber_to_reference(tgt, ref)
    assert "formula" in str(excinfo.value)


# ---------------------------------------------------------------------------
# (b) element mismatch through the REAL tool path
# ---------------------------------------------------------------------------

def test_renumber_tool_surfaces_element_mismatch_not_bogus_mapping():
    """(b) Same atom count, different elements: the mismatch must be VISIBLE
    through the real tool path — the ValueError from map_atoms
    (atom_mapping.py:404) propagates through tool_geometry_renumber
    (server.py:2055). Evidence (gate 2026-11-06): it raises."""
    ref = (
        "3\nwater ref\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n"
    )
    tgt = (
        "3\nH2S tgt\nS 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n"
    )
    with pytest.raises(ValueError) as excinfo:
        server_mod.tool_geometry_renumber(tgt, ref)
    assert "formula" in str(excinfo.value)


# ---------------------------------------------------------------------------
# (c) identical structures -> a real, complete diff
# ---------------------------------------------------------------------------

def test_diff_identical_structures_is_real_and_green():
    """(c) Identical inputs must give a REAL diff: one row per atom, none
    moved, none changed, checks green — never an empty/missing payload."""
    out = se.structure_diff(CO2, CO2)
    assert out["error"] is None
    assert len(out["rows"]) == 3
    assert all(not r["moved"] for r in out["rows"])
    assert all(not r["element_changed"] for r in out["rows"])
    assert out["composition_equal"] is True
    assert out["charge_mult_match"] is True


# ---------------------------------------------------------------------------
# (d) unequal lengths -> every atom accounted for
# ---------------------------------------------------------------------------

def test_diff_unequal_lengths_accounts_for_every_atom():
    """(d) When one side has extra atoms, the diff must cover BOTH totals:
    len(a) + (len(b) - len(a)) rows, with the extras flagged as changed."""
    longer = CO2.replace("3\n", "4\n", 1) + "H 0.0 0.0 2.0\n"
    out = se.structure_diff(CO2, longer)
    assert out["error"] is None
    assert len(out["rows"]) == 4
    assert out["composition_equal"] is False
    extras = [r for r in out["rows"] if r["index"] >= 3]
    assert len(extras) == 1 and extras[0]["to_symbol"] == "H"
    # And the reverse direction too.
    out_rev = se.structure_diff(longer, CO2)
    assert len(out_rev["rows"]) == 4
    assert out_rev["rows"][3]["from_symbol"] == "H"
    assert out_rev["rows"][3]["to_symbol"] is None


# ---------------------------------------------------------------------------
# chemistry: no guessed H on a terminal double-bonded oxygen
# ---------------------------------------------------------------------------

def test_fix_hydrogens_does_not_hydroxylate_co2():
    """A valence-only count adds an H to each terminal O of CO2 (heavy
    degree 1, valence 2). That turns CO2 into HO-C=O(H) — a chemically
    wrong structure. fix_hydrogens must not guess H there."""
    out = se.fix_hydrogens(CO2)
    assert out["error"] is None
    added_symbols = _symbols(out["xyz"])
    assert added_symbols == ["C", "O", "O"], added_symbols
    assert out["n_added"] == 0, out["added_rows"]


def test_fix_hydrogens_adds_missing_hs_to_a_bare_carbon():
    """Positive direction (mirror of the CO2 case): a single isolated sp3
    carbon (heavy-bond order 0, target valence 4) MUST get four hydrogens
    added -- methane, not a bare C. Fires when it should."""
    bare_c = "1\nsingle carbon\nC 0.000000 0.000000 0.000000\n"
    out = se.fix_hydrogens(bare_c)
    assert out["error"] is None
    added_symbols = _symbols(out["xyz"])
    assert added_symbols == ["C", "H", "H", "H", "H"], added_symbols
    assert out["n_added"] == 4, out["added_rows"]
    # Every added H is bonded to the carbon (already placed, not stray).
    bonded_to = {r["bonded_to"] for r in out["added_rows"]}
    assert bonded_to == {0}, bonded_to


# ---------------------------------------------------------------------------
# robustness: None / empty inputs refuse cleanly, never crash deep
# ---------------------------------------------------------------------------

def test_tools_refuse_none_and_empty_inputs_cleanly():
    with pytest.raises((ValueError, TypeError)):
        server_mod.tool_geometry_renumber(None, None)
    with pytest.raises((ValueError, TypeError)):
        server_mod.tool_geometry_fix_hydrogens("")
    with pytest.raises((ValueError, TypeError)):
        server_mod.tool_geometry_diff("", CO2)


# ---------------------------------------------------------------------------
# security contract (operator ruling, binding): the tools never write files;
# an output/output_path attempt is refused; ".." / overwrite are impossible
# because no tool accepts a path argument at all.
# ---------------------------------------------------------------------------

def test_geometry_tools_never_accept_an_output_path():
    """(e) The tools take no path/output parameter (their schema is
    xyz-text-in only): passing one is a TypeError. There is therefore no
    ``..`` traversal and no input/overwrite possible -- a tool can never
    write a file, which would bypass the preset's write-scope gate."""
    import os

    tmp = os.getcwd()  # real cwd inside the gate; nothing must appear here
    before = set(os.listdir(tmp))
    with pytest.raises(TypeError):
        server_mod.tool_geometry_renumber(CO2, CO2, output_path="new.xyz")
    with pytest.raises(TypeError):
        server_mod.tool_geometry_fix_hydrogens(CO2, output="new.xyz")
    with pytest.raises(TypeError):
        server_mod.tool_geometry_diff(CO2, CO2, out_path="../esc.xyz")
    after = set(os.listdir(tmp))
    assert after == before  # no file, no directory, nothing written


# ---------------------------------------------------------------------------
# physical realism (package target #1): the diff must report whether the
# structure has two heavy atoms closer than a stated distance (a clash).
# A clashing pair must surface a red no_clash_below check -- today's
# structure_diff returns only two string checks and evaluates no clash.
# ---------------------------------------------------------------------------

def test_structure_diff_surfaces_a_clash_check():
    """A diff of two identical-composition structures where one heavy-atom
    pair sits 0.5 A apart (far below any covalent-radii sum) must evaluate a
    no-clash check to False. Composition equal and charge/mult equal are
    NOT enough: the pair is chemically impossible and must be flagged."""
    clash_a = (
        "2\nclash\nC 0.000000 0.000000 0.000000\n"
        "C 0.500000 0.000000 0.000000\n"
    )
    clash_b = (
        "2\nclash\nC 0.000000 0.000000 0.000000\n"
        "C 0.500000 0.000000 0.000000\n"
    )
    out = se.structure_diff(clash_a, clash_b)
    assert out["error"] is None
    assert out["composition_equal"] is True
    assert out["charge_mult_match"] is True
    # The package target requires a "no clashes below a stated distance"
    # check: a clashing pair must surface it. The builder's fix exposes
    # clash_detected (structure_edit.py:368) evaluated against
    # _CLASH_ANGSTROM = 1.0 A (structure_edit.py:55). My 0.5 A pair is
    # far below that, so it MUST be flagged -- a chemically impossible
    # structure reported as a clean diff is the defect this catches.
    assert out["clash_detected"] is True
    assert out["min_distance"] < out["clash_angstrom"]
