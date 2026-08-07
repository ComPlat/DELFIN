"""The list of held values: nothing preselected, and what is selected can be seen.

A held value used to be a line of text with a delete button.  It now answers the
two questions worth asking about it -- which atoms does it hold, and to what --
and the second one can be changed without setting the whole constraint again.
Nothing is selected to begin with, because selecting marks atoms in the picture
and a preselected entry would mark a set nobody asked for.
"""

from __future__ import annotations

import pathlib
import tempfile

import pytest

from delfin.dashboard.context import DashboardContext
from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

_WATER = "3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"


@pytest.fixture
def editor(tmp_path):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_submit

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    sent: list[str] = []
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = sent.append
    _widget, refs = tab_submit.create_tab(ctx)
    refs["coords_widget"].value = _WATER
    sent.clear()
    return refs, sent


def _hold(refs, atoms, value, kind="distance", mode="pull"):
    state = refs["editor_state"]
    state["constraints"] = list(state.get("constraints") or []) + [
        {"atoms": list(atoms), "value": float(value), "kind": kind, "mode": mode}
    ]
    refs["refresh_constraints"]()


def test_nothing_is_selected_when_a_value_is_held(editor):
    refs, _sent = editor

    _hold(refs, (0, 1), 1.20)

    dropdown = refs["submit_constraint_dd"]
    assert dropdown.value == "", "a preselected entry would mark atoms unasked"
    assert dropdown.options[0][1] == "", "the first entry must select nothing"
    assert len(dropdown.options) == 2, "the held value is missing from the list"


def test_selecting_a_held_value_shows_its_atoms_and_its_value(editor):
    refs, sent = editor
    _hold(refs, (0, 1), 1.20)
    sent.clear()

    refs["submit_constraint_dd"].value = "c0"

    script = "\n".join(sent)
    assert "setPicks" in script, "nothing told the picture which atoms these are"
    assert "[0, 1]" in script, f"the wrong atoms were marked: {script}"
    assert "1.200" in script, "the box was not given the value being held"


def test_deselecting_drops_the_marks_again(editor):
    refs, sent = editor
    _hold(refs, (0, 1), 1.20)
    refs["submit_constraint_dd"].value = "c0"
    sent.clear()

    refs["submit_constraint_dd"].value = ""

    assert "clearSelection" in "\n".join(sent)


def test_editing_the_value_retunes_the_selected_constraint(editor):
    refs, _sent = editor
    _hold(refs, (0, 1), 1.20)
    refs["submit_constraint_dd"].value = "c0"

    refs["submit_internal_value"].value = 1.55

    held = refs["editor_state"]["constraints"]
    assert held[0]["value"] == pytest.approx(1.55)
    assert held[0]["atoms"] == [0, 1], "retuning must not touch the atoms"
    assert held[0]["mode"] == "pull", "retuning must not touch the mode"


def test_the_value_box_alone_does_not_invent_a_constraint(editor):
    refs, _sent = editor
    _hold(refs, (0, 1), 1.20)
    # nothing selected in the list
    refs["submit_internal_value"].value = 9.99

    assert refs["editor_state"]["constraints"][0]["value"] == pytest.approx(1.20)


# ---------------------------------------------------------------------------
# the browser side
# ---------------------------------------------------------------------------
def test_the_editor_can_be_told_which_atoms_to_mark():
    editor_js = submit_manip_bootstrap_js()

    assert "setPicks: setPicks" in editor_js, "the API does not expose it"
    start = editor_js.index("function setPicks(")
    body = editor_js[start:start + 1600]
    assert "redrawHighlights" in body, "the marks would never be drawn"
    assert "pushPicksToPython" in body, "Python would not learn what is picked"
    # the held value has to be written after the readout, or the measured
    # value the readout just wrote would stand instead
    assert body.index("redrawHighlights") < body.index("submit-internal-value")


def test_reset_goes_back_to_the_loaded_structure(editor):
    refs, _sent = editor
    state = refs["editor_state"]
    _hold(refs, (0, 1), 1.20)
    state["bond_edits"] = {(0, 1): 2}
    state["poly_applied"] = "sqp_4"
    # Pulled apart in the viewer: that writes the coordinates back the way a
    # drag does, which is an edit of this structure and not a new one.
    state["manip_inflight"] = True
    refs["coords_widget"].value = (
        "3\nmangled\nO 0.0 0.0 0.0\nH 3.96 0.0 0.0\nH -0.24 0.93 0.0\n"
    )

    refs["submit_reset_btn"].click()

    assert refs["coords_widget"].value == _WATER
    assert state["constraints"] == []
    assert state["bond_edits"] == {}
    assert state["poly_applied"] is None
    assert state["structure_undo"] == []
