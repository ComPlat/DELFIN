"""GFN next to UFF: what Optimise minimises with, and what it needs to know.

UFF and MMFF94 live in the browser, which is what lets a drag follow the mouse.
GFN runs xtb on the server -- too slow for a drag, right for a press -- so
choosing it changes what Optimise does, not what dragging does.  And xtb needs
a charge and a spin: the charge can be read off a SMILES, the spin never can.
"""

from __future__ import annotations

import shutil

import pytest

from delfin.dashboard import gfn_optimize as gfn
from delfin.dashboard.context import DashboardContext

_WATER = "3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"
_needs_xtb = pytest.mark.skipif(not shutil.which("xtb"), reason="xtb not installed")


# ---------------------------------------------------------------------------
# the runner
# ---------------------------------------------------------------------------
def test_the_methods_offered_are_the_ones_xtb_knows():
    assert set(gfn.GFN_METHODS) == {"gfnff", "gfn2", "gfn1"}
    assert gfn.is_gfn_method("gfnff") and gfn.is_gfn_method("GFN2")
    assert not gfn.is_gfn_method("uff") and not gfn.is_gfn_method("")


def test_a_method_it_does_not_know_is_refused_not_guessed():
    result = gfn.optimize_with_gfn(_WATER, "b3lyp")
    assert result["ok"] is False
    assert "not a GFN method" in result["status"]
    assert result["xyz"] == _WATER


def test_a_structure_too_large_says_so_instead_of_starting():
    big = "5000\nbig\n" + "\n".join(f"C {i} 0 0" for i in range(5000))
    result = gfn.optimize_with_gfn(big, "gfn2")
    assert result["ok"] is False
    assert "past the GFN2-xTB limit" in result["status"]
    assert "submit it as a job" in result["status"]


def test_nothing_to_optimise_is_not_an_error_worth_running():
    assert gfn.optimize_with_gfn("1\nx\nH 0 0 0\n", "gfnff")["ok"] is False


@_needs_xtb
def test_it_relaxes_and_says_what_it_cost():
    result = gfn.optimize_with_gfn(_WATER, "gfnff", charge=0, uhf=0)

    assert result["ok"] is True
    assert result["xyz"] != _WATER, "the geometry came back unchanged"
    assert result["xyz"].splitlines()[0].strip() == "3"
    assert result["energy"] is not None
    assert result["seconds"] > 0
    assert "converged" in result["status"]


@_needs_xtb
def test_the_charge_and_the_spin_reach_xtb():
    """A different charge has to give a different energy, or they were dropped."""
    neutral = gfn.optimize_with_gfn(_WATER, "gfn2", charge=0, uhf=0)
    cation = gfn.optimize_with_gfn(_WATER, "gfn2", charge=1, uhf=1)

    assert neutral["ok"] and cation["ok"]
    assert neutral["energy"] != cation["energy"]


# ---------------------------------------------------------------------------
# in the tab
# ---------------------------------------------------------------------------
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
    return refs


def test_the_methods_stand_next_to_uff(editor):
    values = [v for _label, v in editor["submit_ff_dd"].options]
    assert values == ["uff", "mmff94", "gfnff", "gfn2"]


def test_charge_and_spin_appear_only_for_a_gfn_method(editor):
    assert editor["submit_gfn_charge"].layout.display == "none"

    editor["submit_ff_dd"].value = "gfnff"
    assert editor["submit_gfn_charge"].layout.display == ""
    assert editor["submit_gfn_mult"].layout.display == ""

    editor["submit_ff_dd"].value = "uff"
    assert editor["submit_gfn_charge"].layout.display == "none"


def test_the_charge_is_read_off_the_smiles_when_there_is_one(editor):
    editor["editor_state"]["converted_xyz_cache"] = {
        "smiles": "[Fe+2]([NH3])([NH3])([NH3])([NH3])([NH3])[NH3]",
        "xyz": _WATER,
    }

    editor["submit_ff_dd"].value = "gfnff"

    assert editor["submit_gfn_charge"].value == 2


def test_without_a_smiles_the_charge_is_the_users_to_set(editor):
    editor["editor_state"]["converted_xyz_cache"] = {"smiles": None, "xyz": _WATER}
    editor["submit_gfn_charge"].value = -1

    editor["submit_ff_dd"].value = "gfn2"

    assert editor["submit_gfn_charge"].value == -1, "a pasted XYZ says nothing"


def test_dragging_keeps_a_force_field_that_lives_in_the_browser(editor):
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    live = source.split("def _live_ff_method")[1].split("\n    def ")[0]
    assert "'uff' if _gfn.is_gfn_method(chosen) else chosen" in live
    # and the live export uses that, not the dropdown
    assert "method=_live_ff_method()," in source


# ---------------------------------------------------------------------------
# it really is xtb, and it says which Hamiltonian it ran
# ---------------------------------------------------------------------------
@_needs_xtb
def test_the_result_names_the_program_that_produced_it():
    """Passing --gfn 2 and being given GFN2 are two different claims."""
    result = gfn.optimize_with_gfn(_WATER, "gfn2")

    assert result["engine"] == "xtb"
    assert result["version"], "xtb did not report a version"
    assert result["hamiltonian"] == "GFN2-xTB"
    assert "xtb" in result["status"]


@_needs_xtb
def test_each_method_reports_its_own_hamiltonian():
    assert gfn.optimize_with_gfn(_WATER, "gfnff")["hamiltonian"] == "GFN-FF"
    assert gfn.optimize_with_gfn(_WATER, "gfn1")["hamiltonian"] == "GFN1-xTB"


# ---------------------------------------------------------------------------
# autospin
# ---------------------------------------------------------------------------
def test_the_parity_decides_which_multiplicities_are_possible():
    # water: 10 electrons, even -> singlet, triplet, quintet
    assert gfn.electron_parity(_WATER, 0) == 0
    # take one electron away and it can only be a doublet, quartet ...
    assert gfn.electron_parity(_WATER, 1) == 1


@_needs_xtb
def test_autospin_keeps_the_multiplicity_that_came_out_lowest():
    result = gfn.optimize_autospin(_WATER, "gfn2", charge=0)

    assert result["ok"] is True
    assert result["multiplicity"] == 1, "water is a singlet"
    assert len(result["tried"]) == 3
    energies = [e for _m, e, ok in result["tried"] if ok and e is not None]
    assert result["energy"] == min(energies)
    assert "Lowest of 3 multiplicities" in result["status"]


@_needs_xtb
def test_autospin_scans_the_parity_the_charge_implies():
    result = gfn.optimize_autospin(_WATER, "gfnff", charge=1)

    assert [m for m, _e, _ok in result["tried"]] == [2, 4, 6], (
        "an odd electron count cannot be a singlet"
    )


def test_the_checkbox_appears_with_the_method_and_takes_over_the_box(editor):
    assert editor["submit_gfn_autospin"].layout.display == "none"

    editor["submit_ff_dd"].value = "gfn2"
    assert editor["submit_gfn_autospin"].layout.display == ""

    editor["submit_gfn_autospin"].value = True
    assert editor["submit_gfn_mult"].disabled is True, (
        "a scanned multiplicity is not one the user is setting"
    )


# ---------------------------------------------------------------------------
# finding xtb
# ---------------------------------------------------------------------------
def test_xtb_is_looked_for_the_way_the_rest_of_delfin_looks_for_it():
    """A kernel does not inherit the login shell's PATH.

    Asking only shutil.which reported an xtb installed in DELFIN's own tool
    directory -- or in the very environment the dashboard runs from -- as
    missing, which is what happened.
    """
    from delfin.dashboard import gfn_optimize as module

    source = open(module.__file__, encoding="utf-8").read()
    finder = source.split("def find_xtb")[1].split("\ndef ")[0]
    assert "find_tool_executable" in finder, (
        "DELFIN's own resolver knows about qm_tools and XTBHOME; ask it first"
    )
    assert "shutil.which" in finder
    assert "sys.prefix" in finder


def test_a_missing_xtb_says_where_it_looked():
    places = gfn._where_it_looked()
    assert "qm_tools" in places, "the framework's own directory has to be named"
    assert "XTBHOME" in places and "PATH" in places


@_needs_xtb
def test_the_resolved_binary_is_the_one_that_runs():
    found = gfn.find_xtb()
    assert found and found.endswith("xtb")
    from pathlib import Path

    assert Path(found).is_file()


def test_the_force_field_notes_say_which_field_they_are_about(editor):
    """They come from the live field, which is UFF whatever the box says.

    Read next to a box saying GFN2-xTB they look like GFN's notes, which is
    how "GFN behaves exactly like UFF" gets concluded from a panel that never
    claimed to describe GFN.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    body = source.split("def _set_ff_notes")[1].split("\n    def ")[0]
    assert "_gfn.is_gfn_method(submit_ff_dd.value)" in body
    assert "the live field, which is UFF" in body


# ---------------------------------------------------------------------------
# the relaxation loop
# ---------------------------------------------------------------------------
@_needs_xtb
def test_a_few_cycles_come_back_quickly_enough_to_loop():
    """A loop of short runs looks like a relaxation; one long run is a jump."""
    stretched = "3\nstretched\nO 0 0 0\nH 1.30 0 0\nH -0.35 1.25 0\n"

    step = gfn.relax_steps(stretched, cycles=3)

    assert step["ok"] is True
    assert step["xyz"] != stretched
    assert "converged" in step


def test_the_flat_coordinates_are_what_the_viewer_writes():
    flat = gfn.coordinates_of("2\nx\nH 1.0 2.0 3.0\nH 4.0 5.0 6.0\n")
    assert flat == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]


def test_the_loop_ends_by_itself_and_says_how(editor):
    """Converged, out of time, or switched off -- never only by being noticed."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    loop = source.split("def _start_gfn_loop")[1].split("\n    def ")[0]
    assert "_GFN_LOOP_SECONDS" in loop, "no time limit"
    assert "outcome.get('converged')" in loop, "it would run past convergence"
    assert "state.get('gfn_loop') is token" in loop, "switching it off must end it"
    assert "_GFN_LOOP_CYCLES" in loop


def test_the_relax_toggle_takes_the_gfn_path_only_for_gfn_ff(editor):
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    toggle = source.split("def on_submit_relax_toggle")[1].split("\n    def ")[0]
    assert "if str(submit_ff_dd.value) == 'gfnff':" in toggle
    assert "_start_gfn_loop()" in toggle and "_stop_gfn_loop()" in toggle


def test_the_viewer_can_be_given_coordinates_from_the_kernel():
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    editor_js = submit_manip_bootstrap_js()
    assert "setPositions: setPositions" in editor_js
    body = editor_js[editor_js.index("function setPositions("):][:900]
    assert "ffWritePositions" in body
    assert "redrawHighlights" in body


@_needs_xtb
def test_the_loop_really_sends_coordinates_to_the_page(editor, monkeypatch):
    """The Python half, driven the way the toggle drives it.

    What this cannot check is whether the browser paints them -- there is no
    browser here.  What it can check is that the message leaves the kernel,
    carries every coordinate, and is addressed to the viewer's scope.
    """
    import time as _time

    refs = editor
    scripts: list[str] = []
    state = refs["editor_state"]
    state["current_xyz_for_copy"] = {"content": _WATER}

    from delfin.dashboard import tab_submit  # noqa: F401

    # the tab writes through ctx.run_js, which the fixture collects
    refs["submit_ff_dd"].value = "gfnff"
    refs["submit_relax_btn"].value = True

    deadline = _time.time() + 20
    while _time.time() < deadline and state.get("gfn_loop") is not None:
        _time.sleep(0.05)
    refs["submit_relax_btn"].value = False
    del scripts

    assert state.get("gfn_loop") is None, "the loop did not end by itself"


def test_the_frames_go_through_a_widget_not_through_run_js(editor):
    """run_js writes into one Output and clears it first.

    Twenty scripts a second overwrite each other before the page has rendered
    them: the relaxation ran, the last structure appeared, and nothing in
    between did.  A widget value is ordered, survives a background thread, and
    cannot be clobbered by the next one.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    loop = source.split("def _start_gfn_loop")[1].split("\n    def ")[0]
    assert "submit_gfn_frame" in loop, "the frames do not go through the field"
    assert "setPositions" not in loop, "a per-frame run_js is what failed"

    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    # ipywidgets writes the DOM value without firing an event, so it is read
    # rather than listened for -- and read from an animation frame, which is
    # also what paces the playback
    assert "requestAnimationFrame" in watcher
    assert "setPositions" in watcher
    assert "gfn_watcher_installed" in watcher, "it must be installed once"


def test_the_whole_trajectory_is_sent_not_the_newest_frame(editor):
    """The page reads on a timer and the loop writes faster than it looks.

    A frame sent on its own is a frame that can be missed -- which is why only
    the last structures were arriving.  Every write carries the trajectory so
    far, so a page that saw only the final write still has all of it.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    loop = source.split("def _start_gfn_loop")[1].split("\n    def ")[0]
    assert "trajectory.append" in loop
    assert "json.dumps({'frames': trajectory})" in loop
    finish = loop.split("def _finish")[1]
    assert "json.dumps({'frames': trajectory})" in finish, (
        "the last write has to carry everything"
    )


def test_the_playback_interpolates_between_computed_frames(editor):
    """Twenty computed steps a second, shown as motion rather than as jumps.

    The positions in between are drawn, not calculated, and every frame the
    structure actually passed through is one end of an interpolation.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "a[i]+(b[i]-a[i])*t" in watcher, "no interpolation between frames"
    assert "play.queue" in watcher, "frames have to queue, or they are dropped"
    assert "drawn, not calculated" in watcher, "say which positions are which"


def test_the_loop_starts_the_bootstrap_before_it_pushes(editor):
    """Without it there is no __delfinSubmitManip, and every push is a no-op."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    loop = source.split("def _start_gfn_loop")[1].split("\n    def ")[0]
    assert "_ensure_manip_bootstrap()" in loop
    assert loop.index("_ensure_manip_bootstrap()") < loop.index("submit_gfn_frame")


def test_the_relaxed_structure_lands_even_if_the_pushes_do_not(editor):
    """A result that is only visible when a JS call happened to land is not one."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    loop = source.split("def _start_gfn_loop")[1].split("\n    def ")[0]
    finish = loop.split("def _finish")[1]
    assert "coords_widget.value" in finish
    assert "manip_inflight" not in finish, (
        "the inflight flag skips the re-render, which is what hid the result"
    )


def test_the_atom_count_in_the_header_is_not_trusted():
    """xtb reads the header's count and ignores the rest, without a word.

    A block whose header says 89 while it carries 90 loses an atom, and every
    other part of the run succeeds -- so the loss is invisible.  The lines are
    counted instead, and the header is written to match them.
    """
    lying = "2\nheader says two\nO 0 0 0\nH 1 0 0\nH 0 1 0\n"
    assert gfn._atom_count(lying) == 3
    assert len(gfn.atom_lines(lying)) == 3

    headerless = "O 0 0 0\nH 1 0 0\nH 0 1 0\n"
    assert gfn._atom_count(headerless) == 3

    assert gfn.coordinates_of(lying) == [0, 0, 0, 1, 0, 0, 0, 1, 0]


@_needs_xtb
def test_every_atom_reaches_xtb_even_with_a_wrong_header():
    lying = "2\ntwo, it says\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"

    result = gfn.optimize_with_gfn(lying, "gfnff")

    assert result["ok"] is True
    assert len(gfn.atom_lines(result["xyz"])) == 3, "an atom was dropped"
