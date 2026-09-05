"""Set is a mode: while it is on, the value box turns the selection live.

Set and Hold answer two different questions.  Set is "put it there, and let me
turn it by hand while I watch"; Hold is "keep it there while the field runs,
with pull or fix".  Set used to be a single push that also dropped the picks,
so turning a dihedral by hand meant re-picking four atoms for every step.
"""

from __future__ import annotations

import pytest

from delfin.dashboard.context import DashboardContext

# Hydrogen peroxide: four atoms, so a dihedral is a real selection here.
_PEROXIDE = (
    "4\nhydrogen peroxide\n"
    "O -0.70 0.10 0.00\n"
    "O  0.70 -0.10 0.00\n"
    "H -0.95 -0.50 0.72\n"
    "H  0.95 0.50 0.72\n"
)


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
    refs["coords_widget"].value = _PEROXIDE
    sent.clear()
    return refs, sent


def _pick(refs, indices):
    refs["submit_pick_sync"].value = ",".join(str(i) for i in indices)


def test_set_is_a_mode_that_can_be_switched_off_again(editor):
    refs, _sent = editor
    assert refs["submit_internal_btn"].value is False

    refs["submit_internal_btn"].value = True
    refs["submit_internal_btn"].value = False

    assert refs["submit_internal_btn"].value is False


def test_switching_set_on_puts_the_selection_at_the_value(editor):
    refs, sent = editor
    _pick(refs, [0, 1])
    refs["submit_internal_value"].value = 1.20
    sent.clear()

    refs["submit_internal_btn"].value = True

    assert "setInternal" in "\n".join(sent)


def test_while_set_is_on_every_step_of_the_box_turns_it(editor):
    refs, sent = editor
    _pick(refs, [0, 1, 2])
    refs["submit_internal_btn"].value = True
    sent.clear()

    refs["submit_internal_value"].value = 104.5      # one press of an arrow key

    script = "\n".join(sent)
    assert "setInternal" in script
    assert "104.5" in script


def test_the_picks_survive_a_step_so_it_can_be_turned_again(editor):
    refs, sent = editor
    _pick(refs, [0, 1, 2])
    refs["submit_internal_btn"].value = True
    sent.clear()

    refs["submit_internal_value"].value = 104.6

    assert "clearSelection" not in "\n".join(sent), (
        "letting go of the picks is what made turning by hand impossible"
    )


def test_with_set_off_the_box_moves_nothing(editor):
    refs, sent = editor
    _pick(refs, [0, 1, 2])
    sent.clear()

    refs["submit_internal_value"].value = 120.0

    assert "setInternal" not in "\n".join(sent)


def test_each_quantity_gets_the_step_it_is_turned_in(editor):
    """Three quantities, three steps: a hundredth of an Angstrom for a bond,
    a tenth of a degree for an angle, half a degree for a dihedral -- which is
    the one that gets swept through a whole rotation by hand."""
    refs, _sent = editor

    _pick(refs, [0, 1])
    assert refs["submit_internal_value"].step == 0.01

    _pick(refs, [0, 1, 2])
    assert refs["submit_internal_value"].step == 0.1

    _pick(refs, [0, 1, 2, 3])
    assert refs["submit_internal_value"].step == 0.5


def test_a_held_value_being_retuned_keeps_the_box_to_itself(editor):
    """With an entry selected in the held list, the box is that entry's."""
    refs, sent = editor
    state = refs["editor_state"]
    state["constraints"] = [
        {"atoms": [0, 1], "value": 1.2, "kind": "distance", "mode": "pull"}
    ]
    refs["refresh_constraints"]()
    refs["submit_constraint_dd"].value = "c0"
    refs["submit_internal_btn"].value = True
    sent.clear()

    refs["submit_internal_value"].value = 1.55

    assert state["constraints"][0]["value"] == pytest.approx(1.55)
    assert "setInternal" not in "\n".join(sent), (
        "retuning a held value must not also shove the structure"
    )
