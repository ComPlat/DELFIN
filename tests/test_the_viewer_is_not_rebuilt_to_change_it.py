"""Browsing molecules and switching the numbers must not rebuild the viewer.

Both used to re-emit the whole viewer HTML: a new div, a new WebGL context, the
old one disposed.  A browser grants a handful of contexts, so doing that on
every click is both the slow path and the reason a viewer can end up black until
the tab is switched.  Changing what a living viewer shows is a model swap.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

TWO_MOLECULES = (
    "a.xyz;\n2\n\nH 0 0 0\nH 0 0 0.74\n*\n\n"
    "b.xyz;\n2\n\nO 0 0 0\nO 0 0 1.21\n*\n"
)


@pytest.fixture
def tab(tmp_path):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_orca_builder

    scripts: list[str] = []
    ctx = SimpleNamespace(
        backend=SimpleNamespace(),
        calc_dir=Path(tmp_path),
        run_js=scripts.append,
        add_init_js=scripts.append,
    )
    _widget, widgets_map = tab_orca_builder.create_tab(ctx)
    widgets_map["orca_coords"].value = TWO_MOLECULES
    scripts.clear()
    return widgets_map, scripts


def test_browsing_to_the_next_molecule_keeps_the_viewer(tab):
    widgets_map, scripts = tab

    widgets_map["orca_mol_next_btn"].click()

    assert len(scripts) == 1, "browsing should be one message to the page"
    script = scripts[0]
    assert "removeAllModels" in script, "the model was not swapped"
    assert "__delfinCreateViewer" not in script, "a new WebGL context was created"


def test_browsing_does_not_move_the_camera(tab):
    widgets_map, scripts = tab

    widgets_map["orca_mol_next_btn"].click()

    script = scripts[0]
    for camera in ("zoomTo", "setView", "center("):
        assert camera not in script, (
            f"{camera} would throw away the orientation the user set"
        )


def test_switching_the_numbers_off_keeps_the_viewer(tab):
    """And keeps the model, too.

    The switch is the editor's ``#`` now -- the Builder had a second pair of
    numbering controls beside the block stepper, and one tab wanting two of
    them was one too many.

    The numbers used to go away by swapping the model for the same
    coordinates. That is how the Submit tab's editor lost its bonds -- a model
    rebuilt from coordinates has them perceived from distances again -- and
    the Builder is getting that editor. So the sprites go and nothing else.
    """
    widgets_map, scripts = tab

    widgets_map["submit_labels_btn"].value = True
    scripts.clear()
    widgets_map["submit_labels_btn"].value = False

    assert len(scripts) == 1
    script = scripts[0]
    assert "__delfinAtomNumbers.set(" in script and ",false," in script
    assert "_submitMolViewerByScope" in script
    assert "removeAllModels" not in script
    assert "addModel" not in script
    assert "__delfinCreateViewer" not in script


def test_switching_the_numbers_back_on_keeps_the_viewer(tab):
    widgets_map, scripts = tab
    scripts.clear()

    widgets_map["submit_labels_btn"].value = True

    assert len(scripts) == 1
    assert "__delfinAtomNumbers.set(" in scripts[0] and ",true," in scripts[0]
    assert "_submitMolViewerByScope" in scripts[0]
    assert "removeAllModels" not in scripts[0]
    assert "__delfinCreateViewer" not in scripts[0]


def test_the_update_waits_for_a_viewer_that_is_still_starting(tab):
    widgets_map, scripts = tab

    widgets_map["orca_mol_next_btn"].click()

    # Clicking right after the viewer was displayed must not be a no-op.
    assert "setTimeout(go,50)" in scripts[0]
