"""Adjust H: whether an edit fills the valences it opened, or leaves them.

Filling them is right most of the time -- a carbon that has just lost a bond
should get a hydrogen back.  It is wrong when the open valence is the point:
a radical, a coordination site waiting for a ligand, a fragment about to be
joined to something else.  So it is a switch, and it stands beside the element
that a click draws with.
"""

from __future__ import annotations

import pytest

from delfin.dashboard.molecule_builder import (
    delete_atoms, grow_from, place_atom, set_bond_order, set_element,
    structure_from_xyz,
)

ETHANE = (
    "8\nethane\n"
    "C 0 0 0\nC 1.53 0 0\n"
    "H -0.36 1.02 0\nH -0.36 -0.51 0.88\nH -0.36 -0.51 -0.88\n"
    "H 1.89 1.02 0\nH 1.89 -0.51 0.88\nH 1.89 -0.51 -0.88\n"
)


def hydrogens(structure):
    return sum(1 for symbol in structure.symbols if symbol == "H")


def test_a_bond_that_becomes_double_gives_two_hydrogens_back():
    """Ethane to ethene is not a label change: two hydrogens have to go."""
    on = structure_from_xyz(ETHANE, {})
    set_bond_order(on, 0, 1, 2, adjust_h=True)
    assert hydrogens(on) == 4

    off = structure_from_xyz(ETHANE, {})
    set_bond_order(off, 0, 1, 2, adjust_h=False)
    assert hydrogens(off) == 6, "switched off, nothing is taken away"


def test_changing_the_element_re_satisfies_it_only_when_asked():
    on = structure_from_xyz(ETHANE, {})
    set_element(on, 0, "N", adjust_h=True)
    assert hydrogens(on) == 5, "nitrogen takes one fewer than carbon"

    off = structure_from_xyz(ETHANE, {})
    set_element(off, 0, "N", adjust_h=False)
    assert hydrogens(off) == 6


def test_every_edit_takes_the_switch():
    """One switch, and no edit that quietly ignores it."""
    import inspect

    for edit in (place_atom, grow_from, set_element, delete_atoms,
                 set_bond_order):
        parameters = inspect.signature(edit).parameters
        assert "adjust_h" in parameters, f"{edit.__name__} ignores the switch"
        assert parameters["adjust_h"].default is True, (
            f"{edit.__name__} must fill valences unless it is told otherwise"
        )


def test_the_toolbar_hands_it_to_every_edit():
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def on_submit_cmd")[1].split("\n    def ")[0]
    assert "keep_h = not bool(submit_adjust_h_btn.value)" in handler
    assert handler.count("adjust_h=not keep_h") >= 6, (
        "place, grow, set element, delete and both bond edits"
    )


# ---------------------------------------------------------------------------
# the element it draws with
# ---------------------------------------------------------------------------
@pytest.fixture
def editor(tmp_path):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_submit
    from delfin.dashboard.context import DashboardContext

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = lambda _script: None
    _widget, refs = tab_submit.create_tab(ctx)
    return refs


def test_every_element_is_there(editor):
    choices = list(editor["submit_element_dd"].options)

    assert len(choices) == 118, f"{len(choices)} elements offered"
    for symbol in ("H", "C", "Fe", "Pt", "U", "Og", "Tc", "Xe", "La", "Am"):
        assert symbol in choices, symbol
    assert len(set(choices)) == len(choices), "no element twice"


def test_typing_a_symbol_reaches_the_element_it_names(editor):
    """A native select types ahead by label and matches in the order the
    options stand in.  Sorted by atomic number, pressing P would land on
    phosphorus but I would land on indium, which comes first.  The common ones
    are therefore at the top, and each of them is the first of its letter.
    """
    choices = list(editor["submit_element_dd"].options)

    for symbol in ("C", "H", "N", "O", "S", "P", "F", "B", "I"):
        first = next(c for c in choices if c.upper().startswith(symbol))
        assert first == symbol, (
            f"pressing {symbol} would land on {first}, not on {symbol}"
        )
    # and the rest are still reachable, by typing both letters
    assert "Pd" in choices and "In" in choices


def test_the_switch_stands_beside_the_element_it_draws_with(editor):
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    assert "submit_element_dd, submit_adjust_h_btn," in source, (
        "it belongs behind the element dropdown, where drawing happens"
    )
    assert editor["submit_adjust_h_btn"].value is True, (
        "filling the valences is what is wanted most of the time"
    )
