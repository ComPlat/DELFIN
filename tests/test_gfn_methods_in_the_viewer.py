"""GFN next to UFF: what Optimise minimises with, and what it needs to know.

UFF and MMFF94 live in the browser, which is what lets a drag follow the mouse.
GFN runs xtb on the server -- too slow for a drag, right for a press -- so
choosing it changes what Optimise does, not what dragging does.  And xtb needs
a charge and a spin: the charge can be read off a SMILES, the spin never can.
"""

from __future__ import annotations

import pathlib
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


@pytest.fixture
def player_js(tmp_path):
    """The trajectory player as the page receives it.

    Reading it out of the Python source instead means reading string literals
    that are split across lines wherever they happened to be too long.
    """
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_submit

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = lambda _script: None
    tab_submit.create_tab(ctx)
    startup = "\n".join(ctx.init_js_parts)
    return startup[startup.index("window.__delfinGfnPlay"):]


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


def test_the_follow_ends_with_the_drag_and_with_the_method(editor):
    """A loop of processes on a shared machine that only a user can stop is a
    loop that will not be stopped.

    This one has no clock because it needs none: it runs while a hand is on an
    atom and stops when the hand comes off, when the switch goes up, or when a
    method that has no xtb behind it is chosen.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    follow = source.split("def _gfn_follow_step")[1].split("\n    def ")[0]
    assert "while state.get('gfn_follow')" in follow, (
        "it has to stop when the drag does"
    )

    for name in ("on_submit_relax_toggle", "on_submit_ff_changed",
                 "on_submit_cmd"):
        body = source.split(f"def {name}")[1].split("\n    def ")[0]
        assert "_end_gfn_follow()" in body, name


def test_switching_relax_on_starts_the_relaxation_at_once(editor):
    """Pressed, and nothing happening, is the switch claiming to be live.

    It used to wait for a drag, which is not what the same switch does under
    the browser's field: there the structure starts settling the moment it goes
    down.  It ends by itself, because a structure that has converged has
    nothing left to ask for -- that is what keeps a loop of processes on a
    shared machine from being a loop that will not be stopped.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    toggle = source.split("def on_submit_relax_toggle")[1].split("\n    def ")[0]
    gfn_branch = toggle.split("if _gfn.is_gfn_method(submit_ff_dd.value):")[1]
    assert "_arm_gfn_takeup()" in gfn_branch, "the press is what starts it"

    settle = source.split("def _gfn_settle_now")[1].split("\n    def ")[0]
    assert "not outcome.get('converged') and _gfn_live_is_on()" in settle, (
        "a round that did not converge has to be followed by another"
    )
    assert "state.get('gfn_follow')" in settle, (
        "and a hand on an atom takes over from a relaxation of the whole thing"
    )


def test_a_held_value_is_taken_up_without_pressing_optimise(editor):
    """Setting a value, watching nothing happen and having to press Optimise
    for it is the switch claiming to be live and not being it."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    for name in ("on_submit_hold", "on_submit_hold_mode", "_apply_internal_now"):
        body = source.split(f"def {name}")[1].split("\n    def ")[0]
        assert "_arm_gfn_takeup(" in body, name

    takeup = source.split("def _arm_gfn_takeup")[1].split("\n    def ")[0]
    assert "_gfn_live_is_on()" in takeup, "only while the switch is down"
    assert "state['gfn_settle_forced'] = True" in takeup


def test_the_runner_says_whether_it_converged(editor):
    """String-matching a status line to decide whether to go round again is a
    decision resting on a sentence someone may reword."""
    import inspect

    source = inspect.getsource(gfn.optimize_with_gfn)
    assert "'converged': True," in source and "'converged': False," in source


def test_the_viewer_can_be_given_coordinates_from_the_kernel():
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    editor_js = submit_manip_bootstrap_js()
    assert "setPositions: setPositions" in editor_js
    body = editor_js[editor_js.index("function setPositions("):][:900]
    assert "ffWritePositions" in body
    assert "redrawHighlights" in body


@_needs_xtb
def test_relaxing_ends_by_itself_and_then_costs_nothing(editor):
    """A loop of processes on a shared machine that only a user can stop is a
    loop that will not be stopped.  This one stops when the structure has
    converged, and nothing runs after that until something changes.
    """
    import time as _time

    refs = editor
    state = refs["editor_state"]
    propane = (
        "9\npropane\n"
        "C -1.26 0.00 0.00\nC 0.00 1.16 0.00\nC 1.26 0.00 0.00\n"
        "H -2.15 0.63 0.00\nH -1.30 -0.64 0.88\nH 0.00 1.50 0.89\n"
        "H 0.00 1.50 -0.89\nH 2.15 0.63 0.00\nH 1.30 -0.64 0.88\n"
    )
    refs["coords_widget"].value = propane
    state["current_xyz_for_copy"] = {"content": propane}
    refs["submit_ff_dd"].value = "gfnff"

    refs["submit_relax_btn"].value = True
    deadline = _time.time() + 60
    while _time.time() < deadline and "converged after" not in refs["mol_status"].value:
        _time.sleep(0.05)
    assert "converged after" in refs["mol_status"].value, refs["mol_status"].value
    assert "Settled with GFN-FF" in refs["coords_widget"].value

    # and then it is quiet: nothing armed, nothing running, nothing to run
    _time.sleep(1.0)
    assert state.get("gfn_settle_busy") is False
    assert state.get("gfn_settle_forced") is False
    assert state.get("gfn_settle_armed") is False

    refs["submit_relax_btn"].value = False


def test_the_frames_go_through_a_widget_not_through_run_js(editor):
    """run_js writes into one Output and clears it first.

    Twenty scripts a second overwrite each other before the page has rendered
    them: the relaxation ran, the last structure appeared, and nothing in
    between did.  A widget value is ordered, survives a background thread, and
    cannot be clobbered by the next one.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    for name in ("_gfn_follow_step", "on_submit_optimize(change=None, "
                 "every_frame=False)"):
        body = source.split(f"def {name}")[1].split("\n    def ")[0]
        assert "submit_gfn_frame" in body, f"{name} does not use the field"
        assert "setPositions" not in body, "a per-frame run_js is what failed"

    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    # ipywidgets writes the DOM value without firing an event, so it is read
    # rather than listened for -- and read from an animation frame, which is
    # also what paces the playback
    assert "requestAnimationFrame" in watcher
    assert "setPositions" in watcher
    assert "gfn_watcher_installed" in watcher, "it must be installed once"


def test_the_trail_is_sent_not_the_newest_frame_alone(editor):
    """The page reads on a timer and the kernel writes faster than it looks.

    A frame sent on its own is a frame that can be missed -- which is why only
    the last structures were arriving.  Every write carries a run of frames, so
    a page that looked once still has the ones it did not see.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    follow = source.split("def _gfn_follow_step")[1].split("\n    def ")[0]
    assert "frames.append(" in follow
    assert "frames[-40:]" in follow, "one frame per write can be missed"

    handler = source.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0]
    assert "walked[-400:]" in handler


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


def test_the_bootstrap_is_on_the_page_before_anything_pushes(editor):
    """Without it there is no __delfinSubmitManip, and every push is a no-op --
    which is exactly what "Relaxing..." and nothing moving looked like."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    toggle = source.split("def on_submit_relax_toggle")[1].split("\n    def ")[0]
    gfn_branch = toggle.split("if _gfn.is_gfn_method(submit_ff_dd.value):")[1]
    assert "_ensure_manip_bootstrap()" in gfn_branch
    assert "_install_gfn_frame_watcher()" in gfn_branch

    handler = source.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0]
    assert handler.index("_ensure_manip_bootstrap()") < handler.index(
        "def _push_frames")


def test_the_optimised_structure_lands_even_if_the_pushes_do_not(editor):
    """A result that is only visible when a JS call happened to land is not one.

    The playback is the nice part; the coordinate box is the part that has to
    be true, and Copy and Submit both read it.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0]
    apply_step = handler.split("def _apply()")[1]
    assert "coords_widget.value = (" in apply_step
    assert "if played[0]:" in apply_step, (
        "the re-render is skipped only when the picture is already right"
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


def test_what_is_wrong_with_the_structure_is_said_without_xtb(monkeypatch):
    """CI has no xtb, and neither does every machine a user sits at.

    Too large is too large on any machine; hearing about a missing program
    first leaves a caller unable to act on either problem.
    """
    monkeypatch.setattr(gfn, "find_xtb", lambda: None)
    big = "5000\nbig\n" + "\n".join(f"C {i} 0 0" for i in range(5000))

    assert "past the GFN2-xTB limit" in gfn.optimize_with_gfn(big, "gfn2")["status"]
    assert gfn.optimize_with_gfn("1\nx\nH 0 0 0\n", "gfnff")["status"] == (
        "There is nothing to optimise."
    )
    # a structure with nothing wrong with it does report the missing program
    fine = "3\nx\nO 0 0 0\nH 1 0 0\nH 0 1 0\n"
    assert "needs xtb" in gfn.optimize_with_gfn(fine, "gfnff")["status"]


def test_fullscreen_has_a_status_line_of_its_own():
    """Both views need it, and neither may take the other's.

    Fullscreen relocates its members by hand; ipywidgets knows nothing about
    that, so the line borrowed for the big view came back somewhere else and
    was lost from the small one.  Two widgets carrying the same text is the
    version where nothing is moved that the small view needs.
    """
    from delfin.dashboard import tab_submit
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    source = open(tab_submit.__file__, encoding="utf-8").read()
    assert "mol_status_fs.add_class('submit-fs-member-status')" in source
    setter = source.split("def _set_mol_status")[1].split("\n    def ")[0]
    assert "mol_status.value = rendered_html" in setter
    assert "mol_status_fs.value = '' if prompt else rendered_html" in setter, (
        "they must agree on everything except a prompt to load a structure"
    )
    assert "mol_status_fs.value = ''" in setter, "and agree when empty too"

    enter = submit_manip_bootstrap_js().split("function enterFullscreen")[1][:1100]
    assert "'.submit-fs-member-status'" in enter


def test_relax_means_the_molecule_follows_the_drag_under_gfn(editor):
    """There is no GFN engine in the browser to run once per frame.

    So under GFN this switch means the other half of the same idea: while an
    atom is being dragged the rest of the molecule follows it, one short xtb
    run per push and nothing at all when nothing is being dragged.  It used to
    be switched off and hidden, which was right while there was nothing for it
    to do.
    """
    assert editor["submit_relax_btn"].disabled is False

    editor["submit_ff_dd"].value = "gfnff"
    assert editor["submit_relax_btn"].disabled is False
    assert "follow" in editor["submit_relax_btn"].tooltip
    assert "GFN-FF on the server" in editor["submit_relax_btn"].tooltip

    editor["submit_ff_dd"].value = "uff"
    assert editor["submit_relax_btn"].disabled is False
    assert "continuously" in editor["submit_relax_btn"].tooltip


def test_optimise_sends_the_path_for_the_viewer_to_play(editor):
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "outcome.get('frames')" in handler
    assert "submit_gfn_frame" in handler
    assert "_install_gfn_frame_watcher" in handler
    # both ways of showing a result rebuild the viewer, and both would tear
    # the playback down
    assert "if not played[0]:" in handler, "the isomer path re-renders too"


@_needs_xtb
def test_an_xtb_error_is_a_failure_not_a_partial_result():
    """xtb writes xtbopt.log as it goes, so the files exist even when it dies.

    Two atoms on top of each other make it stop with "something is totally
    wrong" -- and the geometry it had reached was being handed back as though
    it were a result.  Coordinates no one should use must not arrive looking
    like coordinates someone can.
    """
    broken = "3\noverlapping\nO 0 0 0\nH 0.0001 0 0\nH 0 0.0001 0\n"

    result = gfn.optimize_with_gfn(broken, "gfnff")

    if result["ok"]:            # some builds cope with this one
        pytest.skip("this xtb build survived the overlap")
    assert result["xyz"] == broken, "the structure must be left as it was"
    assert result["frames"] == []
    assert "charge" in result["status"] and "overlap" in result["status"]


def test_each_run_is_told_apart_so_a_short_one_still_plays(editor):
    """The player counted the frames it had seen, and the count carried over.

    A run with fewer frames than the one before it therefore played nothing --
    which is why the playback looked like it worked only sometimes.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "run!==play.run" in watcher, "a new run has to reset the count"
    assert "play.seen=0" in watcher

    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "'run': r" in handler, "every payload has to name its run"
    assert "state['gfn_run']" in handler


def test_leaving_fullscreen_puts_every_member_back():
    """The status line is needed in both views, not only the big one.

    Fullscreen moves the widgets into an overlay and the overlay is removed on
    exit -- so a member that could not be put back where it came from would be
    carried out of the page with it and be missing from the small view.
    """
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    editor_js = submit_manip_bootstrap_js()
    exit_body = editor_js.split("function exitFullscreen")[1][:1600]
    assert "insertBefore" in exit_body and "appendChild" in exit_body
    assert "isConnected" in exit_body, "an orphaned member is a lost control"
    assert "root.appendChild(el)" in exit_body


def test_the_finished_geometry_does_not_tear_down_the_playback(editor):
    """Writing the coordinates the ordinary way rebuilds the viewer.

    Done while the trajectory was playing, that destroyed it milliseconds
    after it started -- which is why only the end of the optimisation was ever
    seen.  The playback's last frame is that geometry already.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "played = [False]" in handler
    assert "played[0] = True" in handler
    apply_body = handler.split("def _apply")[1]
    assert "if played[0]:" in apply_body
    assert "state['manip_inflight'] = True" in apply_body
    # the flag has to be set before the write that would re-render
    assert apply_body.index("if played[0]:") < apply_body.index("coords_widget.value = (")


def test_the_playback_finds_its_field_in_fullscreen_too(editor):
    """Fullscreen moves the viewer into an overlay carrying the same scope
    class, and the frame field is not one of the things it takes.  Looking
    only inside the first element with that class found the overlay and no
    field -- so the playback worked small and showed nothing big."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "querySelectorAll" in watcher, "one root is not enough in fullscreen"
    assert "if(!field) field=document.querySelector(" in watcher


def test_optimise_is_a_switch_that_can_be_turned_off(editor):
    """On starts it, off stops it, and it goes back up by itself.

    A push button that cannot be un-pushed leaves the only way out of a long
    optimisation being to wait for it.
    """
    import ipywidgets as w

    button = editor["submit_optimize_btn"]
    assert isinstance(button, w.ToggleButton)
    assert isinstance(editor["submit_optimize_all_btn"], w.ToggleButton)

    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "state['optimize_run'] = None" in handler, "off has to end the run"
    assert "should_stop=_stopped" in handler, "the run has to watch for it"
    apply_body = handler.split("def _apply")[1]
    assert "switch.value = False" in apply_body, (
        "the switches have to release themselves when the work is over"
    )


@_needs_xtb
def test_a_run_that_is_switched_off_ends_rather_than_being_waited_out():
    import threading
    import time

    stretched = "3\nstretched\nO 0 0 0\nH 1.40 0 0\nH -0.40 1.35 0\n"
    stop = [False]
    threading.Timer(0.05, lambda: stop.__setitem__(0, True)).start()

    started = time.perf_counter()
    result = gfn.optimize_with_gfn(
        stretched, "gfnff", should_stop=lambda: stop[0], timeout=60)
    elapsed = time.perf_counter() - started

    if result["ok"]:
        pytest.skip("it finished before the stop arrived")
    assert "was stopped" in result["status"]
    assert result["xyz"] == stretched, "a stopped run must change nothing"
    assert elapsed < 30


def test_the_page_says_what_the_playback_is_doing(editor):
    """Otherwise the only way to tell an invisible trajectory from a missing
    one is to read the browser's console -- which is asking the user to be an
    instrument."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    for report in ('"received "', '"drawing"', '"setPositions did not draw"',
                   '"no setPositions on the page"'):
        assert report in watcher, f"the page never reports {report}"
    assert 'send("gfnplay", text)' in watcher
    assert 'verb+":"+play.serial' in watcher, "the line has to name its verb"

    handler = source.split("def on_submit_cmd")[1].split("\n    def ")[0]
    assert "verb == 'gfnplay'" in handler
    assert "state['gfn_play_note'] = str(payload)" in handler
    assert "Trajectory: {payload}" in handler


def test_the_fullscreen_copy_is_not_seen_next_to_the_original(editor):
    """Both lines carry the same text, so both visible prints it twice."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    assert "mol_status_fs.layout.display = 'none'" in source
    assert ".submit-fs-overlay .submit-fs-member-status {" in source
    assert "display: block !important;" in source


def test_fullscreen_is_not_told_to_enter_coordinates(editor):
    """The copy reports work: a spinner, a trajectory, a result.

    In fullscreen there is a structure on screen, so a prompt to load one is a
    permanent line saying nothing.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    setter = source.split("def _set_mol_status")[1].split("\n    def ")[0]
    assert "'enter XYZ' in str(line)" in setter
    assert "mol_status_fs.value = '' if prompt else rendered_html" in setter


def test_the_player_arrives_with_the_startup_scripts(tmp_path):
    """run_js clears its output before displaying, so a script sent at click
    time can be replaced before the page has run it.  That is how the player
    came to be missing while everything it depends on was working -- and why
    the page never reported anything at all.  add_init_js is the channel the
    explorer's own JS arrives on."""
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_submit
    from delfin.dashboard.context import DashboardContext

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    through_run_js: list[str] = []
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = through_run_js.append
    tab_submit.create_tab(ctx)

    startup = "\n".join(ctx.init_js_parts)
    assert "__delfinGfnPlay" in startup, "the player is not on the page"
    assert not any("__delfinGfnPlay" in s for s in through_run_js), (
        "it must not go through the channel that clears itself"
    )
    assert '"gfnplay"' in startup, "it has to be able to report back"


def test_the_dashboard_runs_without_a_clock_because_it_has_a_switch(editor):
    """A clock that stops a converging optimisation at an arbitrary second is
    worse than one that never stops it.  Optimise is a switch, so the person
    watching decides when a run has gone on long enough."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "timeout=None" in handler
    assert "should_stop=_stopped" in handler, "without it nothing could stop it"


@_needs_xtb
def test_without_a_timeout_the_hand_is_what_stops_it():
    import threading

    stretched = "3\nstretched\nO 0 0 0\nH 1.40 0 0\nH -0.40 1.35 0\n"
    stop = [False]
    threading.Timer(0.05, lambda: stop.__setitem__(0, True)).start()

    result = gfn.optimize_with_gfn(
        stretched, "gfnff", should_stop=lambda: stop[0], timeout=None)

    assert "stopped" in result["status"] or result["ok"]


def test_the_path_is_handed_over_while_it_is_still_being_walked(editor):
    """The trajectory only played once the run was switched off, and frame 1
    stood there while it worked -- because xtbopt.log was read at the end.
    xtb writes it as it goes, so it can be read as it goes."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "def _push_frames" in handler
    assert "on_frames=_push_frames" in handler

    runner = open(
        __import__("delfin.dashboard.gfn_optimize", fromlist=["x"]).__file__,
        encoding="utf-8").read()
    loop = runner.split("while running.poll() is None:")[1][:4200]
    assert "read_trajectory(folder)" in loop, "the log is not read while running"
    assert "on_frames(walking)" in loop
    # and the switch is checked before any of that, every pass
    assert loop.index("should_stop()") < loop.index("read_trajectory(folder)")


@_needs_xtb
def test_frames_arrive_during_a_long_run_not_after_it():
    import threading
    import time

    zn = pathlib.Path(__file__).parent / "fixtures" / "does-not-exist"
    del zn
    stretched = "6\nbent\nC 0 0 0\nC 1.9 0 0\nH -0.6 0.9 0\nH -0.6 -0.9 0\n" \
                "H 2.5 0.9 0\nH 2.5 -0.9 0\n"
    seen: list[int] = []
    stop = [False]
    threading.Timer(3.0, lambda: stop.__setitem__(0, True)).start()

    result = gfn.optimize_with_gfn(
        stretched, "gfnff", timeout=None, should_stop=lambda: stop[0],
        on_frames=lambda frames: seen.append(len(frames)))

    # a molecule this small may well finish inside one poll; what must hold is
    # that nothing is reported twice and the final result agrees
    assert seen == sorted(seen)
    if seen:
        assert seen[-1] <= len(result["frames"])


@_needs_xtb
def test_the_switch_still_works_on_a_run_that_is_reporting_frames():
    """Reading the log is the expensive part and it grows, so parsing it on
    every pass crowded the stop check out: the trajectory appeared and the
    switch stopped working.  Stopping comes first and every time round."""
    import threading
    import time

    zn_like = "3\nstretched\nO 0 0 0\nH 1.45 0 0\nH -0.45 1.35 0\n"
    stop = [False]
    threading.Timer(0.3, lambda: stop.__setitem__(0, True)).start()

    started = time.perf_counter()
    result = gfn.optimize_with_gfn(
        zn_like, "gfnff", timeout=None, should_stop=lambda: stop[0],
        on_frames=lambda frames: None)
    elapsed = time.perf_counter() - started

    assert elapsed < 5, "the stop check has to run whatever else is going on"
    assert result["ok"] or "stopped" in result["status"]


def test_the_playback_speeds_up_when_it_falls_behind(editor):
    """xtb computes faster than a fixed frame rate can show.

    75 frames arrive in 0.4 s and would take 4 s at 55 ms each, so the picture
    trails the calculation and keeps trailing it further.  A backlog is played
    faster; the whole path is still shown.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "function stepMs()" in watcher
    assert "play.queue.length" in watcher
    assert "stepMs()" in watcher.split("var t=(now-play.started)")[1][:20]


def test_stopping_keeps_the_frame_that_was_on_screen(editor):
    """Stop means what the user was looking at.

    xtb runs ahead of the picture -- it produces frames faster than any frame
    rate shows them -- so keeping its newest geometry hands back something
    nobody saw and did not choose.  The page says which frame it was showing;
    that one is kept and the rest are discarded.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "play.shown=(play.shown||0)+1" in watcher, "nothing counts what was shown"
    assert '"stopped at frame "' in watcher

    cmd = source.split("def on_submit_cmd")[1].split("\n    def ")[0]
    assert "state['gfn_shown_frame']" in cmd

    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "gfn_shown_frame" in handler
    assert "stopped at the frame on" in handler


def test_a_backlog_is_skipped_rather_than_played_out(editor):
    """A queue allowed to grow puts the picture permanently behind the run."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "play.queue.slice(-20)" in watcher


def test_the_page_stops_the_picture_without_asking_the_kernel(editor):
    """Waiting to be told the switch went off costs a round trip, and the
    playback ran on for the length of it.  ipywidgets marks a pressed toggle
    with mod-active; reading that is instant."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    assert "submit_optimize_btn.add_class('submit-optimize-switch')" in source

    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "function switchIsOn()" in watcher
    assert 'classList.contains("mod-active")' in watcher
    # and it is checked before anything is drawn, every frame
    body = watcher.split("function frame(now)")[1]
    assert body.index("switchIsOn()") < body.index("read(now);")


# ---------------------------------------------------------------------------
# the controls, split and cleaned up
# ---------------------------------------------------------------------------
def test_controls_that_cannot_work_under_gfn_are_taken_away(editor, monkeypatch):
    """Greying them out invites the question of why they are dead.

    Strength is how many steps the browser's field takes per animation frame,
    and that field does not run under GFN.  Relax and Settle both keep a
    meaning -- relax the structure and follow a dragged atom, and tidy up when
    it is let go -- and both are carried out by the method that is chosen,
    never by anything else.

    With an xtb on the machine, that is: whether Relax is there at all is a
    different question, and the machine running the tests answers it either
    way, so it is said here rather than inherited.
    """
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "xtb_available", lambda: True)
    editor["submit_ff_dd"].value = "uff"

    editor["submit_ff_dd"].value = "gfnff"
    assert editor["submit_strength_slider"].layout.display == "none"
    assert editor["submit_relax_btn"].layout.display in (None, "")
    assert editor["submit_settle_btn"].layout.display in (None, "")

    editor["submit_ff_dd"].value = "uff"
    assert editor["submit_strength_slider"].layout.display == ""


def test_optimise_and_optimise_all_are_two_switches(editor):
    import ipywidgets as w

    one = editor["submit_optimize_btn"]
    every = editor["submit_optimize_all_btn"]
    assert isinstance(one, w.ToggleButton) and isinstance(every, w.ToggleButton)
    assert one.description == "Optimize" and every.description == "all"


def test_only_one_of_them_runs_at_a_time(editor):
    """A login node is shared; two sets of xtb processes is how it is noticed."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "other.value = False" in handler
    # and the frames of a set are walked in one loop, not started side by side
    assert "for position, item in enumerate(targets):" in handler
    assert "threading.Thread" in handler and handler.count("threading.Thread") == 1


def test_only_all_takes_the_whole_set(editor):
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "if every_frame else []" in handler


@_needs_xtb
def test_the_energy_is_reported_like_the_force_field_reports_one(editor):
    import time

    editor["coords_widget"].value = "3\nw\nO 0 0 0\nH 1.05 0 0\nH -0.28 1.02 0\n"
    editor["submit_ff_dd"].value = "gfnff"

    editor["submit_optimize_btn"].value = True
    deadline = time.time() + 30
    while time.time() < deadline and editor["submit_optimize_btn"].value:
        time.sleep(0.05)

    energy = editor["editor_state"].get("gfn_energy")
    assert energy is not None and energy < 0, "no energy came back"

    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "E = {energy:.6f} Eh" in handler
    assert "kcal/mol" in handler, "the tab speaks kcal/mol elsewhere"


# ---------------------------------------------------------------------------
# undo answers for what the user did, not for what the optimiser did
# ---------------------------------------------------------------------------
def test_a_whole_optimisation_is_one_step_of_undo(editor):
    """The playback writes coordinates dozens of times a second.

    If each of those were a step, Undo would walk back through an optimisation
    frame by frame and the operation before it would be unreachable.  Undo
    answers for dashboard events: a run is one.
    """
    from delfin.dashboard import tab_submit
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    source = open(tab_submit.__file__, encoding="utf-8").read()
    undo = source.split("def on_submit_manip_undo")[1].split("\n    def ")[0]
    assert "state.pop('pre_optimize_frames'" in undo, (
        "the geometries from before the run are what one undo restores"
    )

    # and the frames the page is given never enter the browser's own stack
    editor_js = submit_manip_bootstrap_js()
    body = editor_js[editor_js.index("function setPositions("):][:900]
    assert "snapshotForUndo" not in body, (
        "a played frame is not something the user did"
    )


# ---------------------------------------------------------------------------
# a pipe nobody reads
# ---------------------------------------------------------------------------
def test_xtb_talks_to_a_file_and_not_into_a_pipe():
    """A pipe holds 64 KiB and then blocks whoever is writing to it.

    The watching loop reads nothing until the process has ended, so xtb waited
    for the loop and the loop waited for xtb: an optimisation that takes half a
    second never finished at all.  It is the switch's own path -- the plain
    call drains its pipes -- so it was every stoppable run, which is every run
    the dashboard starts.
    """
    source = open(gfn.__file__, encoding="utf-8").read()
    runner = source.split("def optimize_with_gfn")[1].split("\ndef ")[0]
    assert "stdout=sink" in runner and "stderr=subprocess.STDOUT" in runner
    assert "subprocess.PIPE" not in runner, "a pipe is what deadlocked"
    assert "record.read_text" in runner, "the output is read back off disk"


@_needs_xtb
def test_a_run_that_says_more_than_a_pipe_holds_still_ends(tmp_path):
    """The regression itself, at the size that showed it.

    A GFN2 optimisation of a decane says 77 276 bytes, which is past what a
    pipe will hold.  Through the stoppable path -- the one the switch uses --
    it used to hang for as long as anyone was willing to wait.
    """
    import time as _time

    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles("CCCCCCCCCC"))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    conf = mol.GetConformer()
    rows = []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        rows.append(f"{atom.GetSymbol()} {p.x:.6f} {p.y:.6f} {p.z:.6f}")
    decane = f"{len(rows)}\ndecane\n" + "\n".join(rows) + "\n"

    started = _time.perf_counter()
    # should_stop given, so it takes the path that watches rather than waits
    result = gfn.optimize_with_gfn(
        decane, "gfn2", charge=0, uhf=0, timeout=120,
        should_stop=lambda: False, on_frames=lambda _frames: None,
    )
    spent = _time.perf_counter() - started

    assert result["ok"] is True, result["status"]
    assert spent < 60, f"it took {spent:.0f} s; it used not to end at all"
    assert len(result["frames"]) > 20, "the whole path has to come back"


# ---------------------------------------------------------------------------
# a hand on the structure while xtb is minimising it
# ---------------------------------------------------------------------------
def test_the_playback_lets_go_of_the_picture_while_an_atom_is_dragged(editor):
    """Otherwise the drag cannot happen at all.

    The player writes every atom's position once per animation frame, so an
    atom being pulled is put back where xtb had it sixty times a second.  The
    page decides this itself: asking the kernel costs a round trip, and the
    picture fights the hand for the length of it.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "function grabbed()" in watcher
    assert "_submitManipStateByScope" in watcher, "it reads the drag off the page"
    assert 'drag.kind==="translate"' in watcher
    assert 'send(held?"gfngrab":"gfnfree","")' in watcher
    assert "if(play.held&&!followIsOn()){" in watcher, (
        "with nothing following, the drag has the picture to itself"
    )
    assert "play.queue=[]; play.last=null;" in watcher, (
        "the queued frames belong to a structure that has just been changed"
    )


def test_the_grab_ends_the_run_and_the_release_starts_the_next_one(editor):
    """The kernel is told at the grab, not at the release.

    A GFN2 run is seconds long; every one of them would otherwise be spent
    minimising a structure the user is in the middle of changing.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def on_submit_cmd")[1].split("\n    def ")[0]
    assert "verb == 'gfngrab'" in handler and "_interrupt_gfn()" in handler
    assert "verb == 'gfnfree'" in handler and "_arm_gfn_restart()" in handler

    # Set, Hold and a bond edit arrive as coordinates, and count the same
    sync = source.split("def on_submit_manip_sync")[1].split("\n    def ")[0]
    assert "if drag_ended:" in sync
    assert "_interrupt_gfn()" in sync and "_arm_gfn_restart()" in sync

    optimise = source.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0]
    assert "state.get('optimize_interrupted') is token" in optimise, (
        "an interrupted run must not write the geometry it reached"
    )
    assert "state.get('gfn_run') != run_id" in optimise, (
        "nor draw the path it had walked over the structure now on screen"
    )


@_needs_xtb
def test_a_moved_atom_is_what_the_next_run_starts_from(editor, monkeypatch):
    """The whole of point five, driven the way the browser drives it.

    xtb really runs; what is recorded is the geometry each run was handed.  The
    first is the structure as loaded, the second is the one with an atom moved
    -- a run that started again from the old coordinates would be optimising
    something the user had already changed.
    """
    import time as _time

    from delfin.dashboard import tab_submit

    refs = editor
    state = refs["editor_state"]
    xyz = (
        "9\npropane\n"
        "C -1.26 0.00 0.00\nC 0.00 0.86 0.00\nC 1.26 0.00 0.00\n"
        "H -2.15 0.63 0.00\nH -1.30 -0.64 0.88\nH 0.00 1.50 0.89\n"
        "H 0.00 1.50 -0.89\nH 2.15 0.63 0.00\nH 1.30 -0.64 0.88\n"
    )
    refs["coords_widget"].value = xyz
    state["current_xyz_for_copy"] = {"content": xyz}
    refs["submit_ff_dd"].value = "gfnff"

    handed: list[str] = []
    real = tab_submit._gfn.optimize_with_gfn

    def recording(text, method, **kwargs):
        handed.append(text)
        return real(text, method, **kwargs)

    monkeypatch.setattr(tab_submit._gfn, "optimize_with_gfn", recording)

    refs["submit_optimize_btn"].value = True
    deadline = _time.time() + 30
    while _time.time() < deadline and not handed:
        _time.sleep(0.02)

    # an atom is picked up, and the page says so before the run can finish
    refs["submit_cmd_sync"].value = "gfngrab:1:"
    assert state.get("optimize_run") is None, "the run was not ended"
    assert state.get("optimize_interrupted") is not None

    body = [line for line in xyz.splitlines()[2:] if line.strip()]
    parts = body[0].split()
    parts[1] = f"{float(parts[1]) - 0.7:.6f}"
    body[0] = " ".join(parts)
    moved = f"{len(body)}\nDELFIN drag-end\n" + "\n".join(body) + "\n"
    refs["submit_manip_sync"].value = moved
    refs["submit_cmd_sync"].value = "gfnfree:2:"

    deadline = _time.time() + 30
    while _time.time() < deadline and len(handed) < 2:
        _time.sleep(0.02)
    assert len(handed) >= 2, "the run never started again"

    assert gfn.atom_lines(handed[0])[0].split()[1].startswith("-1.26")
    assert gfn.atom_lines(handed[1])[0].split()[1].startswith("-1.96"), (
        "the second run was handed the structure from before the drag"
    )

    deadline = _time.time() + 60
    while _time.time() < deadline and refs["submit_optimize_btn"].value:
        _time.sleep(0.05)
    assert refs["submit_optimize_btn"].value is False, "the switch stayed down"
    assert state.get("optimize_run") is None


# ---------------------------------------------------------------------------
# what Hold and Set mean to xtb
# ---------------------------------------------------------------------------
def test_a_held_value_becomes_a_constraint_block_counted_from_one():
    """xtb numbers atoms from one; the editor numbers them from zero."""
    built = gfn.constraint_input([
        {"kind": "distance", "atoms": [0, 1], "value": 1.3, "mode": "pull"},
        {"kind": "angle", "atoms": [1, 0, 2], "value": 109.5, "mode": "pull"},
        {"kind": "dihedral", "atoms": [0, 1, 2, 3], "value": 180.0,
         "mode": "pull"},
    ], atoms=4)

    lines = built["text"].splitlines()
    assert lines[0] == "$constrain"
    assert lines[-1] == "$end"
    assert "  distance: 1, 2, 1.300000" in lines
    assert "  angle: 2, 1, 3, 109.500000" in lines
    assert "  dihedral: 1, 2, 3, 4, 180.000000" in lines
    assert built["held"] == 3


def test_a_pull_and_a_fix_are_two_different_force_constants():
    """Measured on a propane asked for C-C 1.700 A: 0.5 gives 1.6704, which is
    a compromise, and 20 gives 1.6992, which is the value."""
    pull = gfn.constraint_input(
        [{"kind": "distance", "atoms": [0, 1], "value": 1.7, "mode": "pull"}])
    fix = gfn.constraint_input(
        [{"kind": "distance", "atoms": [0, 1], "value": 1.7, "mode": "fix"}])

    assert pull["force"] == gfn.PULL_FORCE_CONSTANT < gfn.FIX_FORCE_CONSTANT
    assert fix["force"] == gfn.FIX_FORCE_CONSTANT
    assert f"force constant={gfn.PULL_FORCE_CONSTANT}" in pull["text"]
    assert f"force constant={gfn.FIX_FORCE_CONSTANT}" in fix["text"]


def test_one_force_constant_holds_the_whole_set_and_it_is_said():
    """A second `force constant=` line is read and ignored -- measured: a block
    asking 0.25 for a distance and 10 for an angle held both at 0.25, digit for
    digit.  So a set with anything exact in it is exact throughout, and a pull
    that has stopped negotiating has to say why."""
    both = gfn.constraint_input([
        {"kind": "distance", "atoms": [0, 1], "value": 1.7, "mode": "pull"},
        {"kind": "angle", "atoms": [0, 1, 2], "value": 100.0, "mode": "fix"},
    ], atoms=3)

    assert both["text"].count("force constant=") == 1
    assert both["force"] == gfn.FIX_FORCE_CONSTANT
    assert both["mixed"] is True
    assert "one force constant for the whole set" in gfn.held_note(both)

    alone = gfn.constraint_input(
        [{"kind": "distance", "atoms": [0, 1], "value": 1.7, "mode": "pull"}])
    assert alone["mixed"] is False


def test_a_held_value_naming_an_atom_that_is_not_there_is_dropped_and_said():
    """xtb would stop with a parse error about a line the user never typed."""
    built = gfn.constraint_input([
        {"kind": "distance", "atoms": [0, 40], "value": 1.7, "mode": "pull"},
        {"kind": "distance", "atoms": [0, 1], "value": 1.7, "mode": "pull"},
        {"kind": "angle", "atoms": [0, 1], "value": 100.0, "mode": "pull"},
    ], atoms=9)

    assert built["held"] == 1, "only the one that names atoms it has"
    assert len(built["dropped"]) == 2
    assert "dropped" in gfn.held_note(built)


def test_nothing_held_writes_no_input_file_at_all():
    assert gfn.constraint_input([])["text"] == ""
    assert gfn.constraint_input([])["force"] is None
    assert gfn.held_note(gfn.constraint_input([])) == ""


def test_holding_an_atom_where_it_is_is_not_asked_of_xtb():
    """Measured on xtb 6.7.1, and the reason the editor does not try.

    `$constrain atoms:` naming one atom changes nothing -- the geometry comes
    back identical to the free one -- and `$fix atoms:` is worse than nothing:
    three fixed carbons of a propane came back with their C-C at 1.4555 A
    instead of 1.5255 under GFN-FF, and at 0.4623 A under GFN2.
    """
    source = open(gfn.__file__, encoding="utf-8").read()
    assert "$fix" not in gfn.constraint_input(
        [{"kind": "distance", "atoms": [0, 1], "value": 1.7, "mode": "fix"}]
    )["text"]
    assert "0.4623" in source, "the measurement that decided it has to be here"


@_needs_xtb
def test_a_pull_negotiates_and_a_fix_is_met():
    """The whole of point six, against the program itself."""
    propane = (
        "9\npropane\n"
        "C -1.26 0.00 0.00\nC 0.00 0.86 0.00\nC 1.26 0.00 0.00\n"
        "H -2.15 0.63 0.00\nH -1.30 -0.64 0.88\nH 0.00 1.50 0.89\n"
        "H 0.00 1.50 -0.89\nH 2.15 0.63 0.00\nH 1.30 -0.64 0.88\n"
    )

    def carbon_carbon(xyz):
        rows = [line.split() for line in gfn.atom_lines(xyz)]
        first = [float(v) for v in rows[0][1:4]]
        second = [float(v) for v in rows[1][1:4]]
        return sum((a - b) ** 2 for a, b in zip(first, second)) ** 0.5

    free = gfn.optimize_with_gfn(propane, "gfnff")
    assert free["ok"], free["status"]
    assert abs(carbon_carbon(free["xyz"]) - 1.49) < 0.05, "free C-C is ~1.49 A"

    asked = [{"kind": "distance", "atoms": [0, 1], "value": 1.70}]
    pull = gfn.optimize_with_gfn(
        propane, "gfnff", constraints=[dict(asked[0], mode="pull")])
    fix = gfn.optimize_with_gfn(
        propane, "gfnff", constraints=[dict(asked[0], mode="fix")])

    assert pull["ok"] and fix["ok"]
    pulled, fixed = carbon_carbon(pull["xyz"]), carbon_carbon(fix["xyz"])
    assert 1.60 < pulled < 1.69, f"a pull settles short of the value: {pulled}"
    assert abs(fixed - 1.70) < 0.005, f"a fix meets it: {fixed}"
    assert "force constant" in fix["status"]


@_needs_xtb
def test_what_the_editor_holds_is_what_the_optimisation_holds(editor):
    """Held on screen and given up the moment GFN is chosen is the bug this
    closes: the angle is asked for and the angle is what comes back."""
    import time as _time

    refs = editor
    state = refs["editor_state"]
    propane = (
        "9\npropane\n"
        "C -1.26 0.00 0.00\nC 0.00 0.86 0.00\nC 1.26 0.00 0.00\n"
        "H -2.15 0.63 0.00\nH -1.30 -0.64 0.88\nH 0.00 1.50 0.89\n"
        "H 0.00 1.50 -0.89\nH 2.15 0.63 0.00\nH 1.30 -0.64 0.88\n"
    )
    refs["coords_widget"].value = propane
    state["current_xyz_for_copy"] = {"content": propane}
    refs["submit_ff_dd"].value = "gfnff"
    state["constraints"] = [
        {"kind": "angle", "atoms": [0, 1, 2], "value": 100.0, "mode": "fix"},
    ]

    refs["submit_optimize_btn"].value = True
    deadline = _time.time() + 60
    while _time.time() < deadline and refs["submit_optimize_btn"].value:
        _time.sleep(0.05)

    rows = [line.split() for line in gfn.atom_lines(refs["coords_widget"].value)]
    points = [[float(v) for v in row[1:4]] for row in rows]
    first = [points[0][n] - points[1][n] for n in range(3)]
    third = [points[2][n] - points[1][n] for n in range(3)]
    import math

    dot = sum(a * b for a, b in zip(first, third))
    norms = math.dist(points[0], points[1]) * math.dist(points[2], points[1])
    angle = math.degrees(math.acos(max(-1.0, min(1.0, dot / norms))))

    assert abs(angle - 100.0) < 1.0, f"asked for 100 deg, got {angle:.2f}"
    assert "held value" in refs["mol_status"].value, (
        "the status has to say what was held while it ran"
    )


# ---------------------------------------------------------------------------
# the molecule follows the atom being dragged
# ---------------------------------------------------------------------------
def test_the_page_hands_the_geometry_over_while_the_mouse_is_down(player_js):
    """A drag that only reports at the release cannot be followed."""
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    assert "function followIsOn()" in player_js
    assert "submit-gfn-follow" in player_js, "the switch is read off the page"
    assert 'pushXyz(scope,"drag-follow")' in player_js
    assert "var since=now-(play.pushed||0);" in player_js, (
        "the pace is the machine's, and it has to be measured to be followed"
    )
    # and the manip API has to offer that push at all
    assert "pushXyz: pushXyzToPython" in submit_manip_bootstrap_js()


def test_the_dragged_atom_belongs_to_the_cursor_not_to_the_answer(player_js):
    """The coordinates come back describing where the atom was when they were
    sent, and the cursor has moved on since.  Written back, the atom would be
    pulled to where it was a fifth of a second ago, sixty times a second."""
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    assert "function heldSerials()" in player_js
    assert "setPositions(scope,out,heldSerials())" in player_js

    editor_js = submit_manip_bootstrap_js()
    body = editor_js[editor_js.index("function setPositions("):][:1200]
    assert "keepSerials" in body
    assert "ffIndicesOf(viewer, keepSerials)" in body
    assert "ffReadPositions(viewer)" in body, (
        "the kept atoms have to come from where they are now"
    )


def test_a_truncated_trail_says_where_in_the_run_it_starts(editor):
    """The tail is sent rather than the whole path -- every write is a message.

    Counting from the front of the message would replay what has been shown and
    then stop showing anything new: the player counts frames of the run, not
    frames of the message.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "var from=(data&&data.from)||0;" in watcher
    assert "play.seen=from+frames.length;" in watcher

    handler = source.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0]
    assert "'from': first" in handler
    follow = source.split("def _gfn_follow_step")[1].split("\n    def ")[0]
    assert "'from': len(frames) - len(trail)" in follow


def test_the_follow_runs_one_process_at_a_time_and_takes_the_newest(editor):
    """A hand moves faster than xtb answers.  A queue of answers about where
    the atom used to be is worse than no answer at all."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    follow = source.split("def _gfn_follow_step")[1].split("\n    def ")[0]
    assert "state.get('gfn_follow_busy')" in follow
    assert "state.pop('gfn_follow_xyz', None)" in follow, "the newest wins"
    assert "relax_steps(" in follow and "_GFN_FOLLOW_CYCLES" in follow
    assert "constraints=constraints" in follow, (
        "a value held on screen is held while the molecule follows too"
    )


@_needs_xtb
def test_the_rest_of_the_molecule_follows_the_atom_that_is_dragged(editor):
    """Point seven, played through the tab the way the browser plays it.

    The carbon is walked out 0.36 A in six pushes, each one carrying the
    relaxed rest from the answer before it -- which is what the page does when
    it writes every atom but the one under the cursor.
    """
    import json as _json
    import math
    import time as _time

    refs = editor
    state = refs["editor_state"]
    propane = (
        "9\npropane\n"
        "C -1.26 0.00 0.00\nC 0.00 0.86 0.00\nC 1.26 0.00 0.00\n"
        "H -2.15 0.63 0.00\nH -1.30 -0.64 0.88\nH 0.00 1.50 0.89\n"
        "H 0.00 1.50 -0.89\nH 2.15 0.63 0.00\nH 1.30 -0.64 0.88\n"
    )

    def points(xyz):
        return [[float(v) for v in line.split()[1:4]]
                for line in gfn.atom_lines(xyz)]

    refs["coords_widget"].value = propane
    state["current_xyz_for_copy"] = {"content": propane}
    refs["submit_ff_dd"].value = "gfnff"
    refs["submit_relax_btn"].value = True
    assert refs["submit_relax_btn"].value is True, "xtb was not found"

    refs["submit_cmd_sync"].value = "gfngrab:1:"
    assert state.get("gfn_follow") is True

    current, seen = propane, 0
    for _step in range(6):
        rows = [line.split() for line in gfn.atom_lines(current)]
        rows[0][1] = f"{float(rows[0][1]) - 0.06:.6f}"     # the cursor moves
        current = (f"{len(rows)}\nDELFIN drag-follow\n"
                   + "\n".join(" ".join(r) for r in rows) + "\n")
        refs["submit_manip_sync"].value = current
        deadline = _time.time() + 20
        while _time.time() < deadline:
            data = _json.loads(refs["submit_gfn_frame"].value or "{}")
            if len(data.get("frames") or []) > seen:
                seen = len(data["frames"])
                break
            _time.sleep(0.01)
        else:
            raise AssertionError("no answer came back for a push")
        # the page keeps the dragged atom and writes the rest
        flat = _json.loads(refs["submit_gfn_frame"].value)["frames"][-1]
        rows = [line.split() for line in gfn.atom_lines(current)]
        for i in range(1, len(rows)):
            rows[i][1:4] = [f"{flat[3 * i + n]:.6f}" for n in range(3)]
        current = (f"{len(rows)}\nDELFIN drag-follow\n"
                   + "\n".join(" ".join(r) for r in rows) + "\n")

    refs["submit_cmd_sync"].value = "gfnfree:2:"
    assert state.get("gfn_follow") is False

    began, ended = points(propane), points(current)
    assert abs(ended[0][0] - (began[0][0] - 0.36)) < 1e-4, (
        "the dragged atom did not stay where the cursor put it"
    )
    assert math.dist(began[1], ended[1]) > 0.05, "the neighbour did not follow"
    assert math.dist(began[2], ended[2]) > 0.05, "nor did the far carbon"
    stretched = math.dist(ended[0], ended[1])
    assert stretched < math.dist(began[0], began[1]) + 0.36, (
        "the rest of the molecule did not close any of the gap"
    )


# ---------------------------------------------------------------------------
# xtb that is not there yet
# ---------------------------------------------------------------------------
def test_the_installer_is_delfins_own_and_is_asked_for_xtb_alone():
    """Naming the tool keeps it to that one: with no arguments the script
    fetches crest, dftb+ and the stda bundle behind it as well."""
    script = gfn.install_script()
    assert script is not None and script.name == "install_qm_tools.sh"

    command = gfn.install_command()
    assert command[0] == "bash" and command[-1] == "xtb"
    assert str(script) in command


def test_the_offer_appears_only_when_it_is_needed(editor, monkeypatch):
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_xtb", lambda: None)
    monkeypatch.setattr(tab_submit._gfn, "xtb_available", lambda: False)

    assert editor["submit_xtb_install_btn"].layout.display == "none"

    editor["submit_ff_dd"].value = "gfnff"
    assert editor["submit_xtb_install_btn"].layout.display == ""
    assert "Install xtb" in editor["mol_status"].value

    editor["submit_ff_dd"].value = "uff"
    assert editor["submit_xtb_install_btn"].layout.display == "none", (
        "nothing to install for a field that runs in the browser"
    )


def test_the_offer_stays_away_when_xtb_is_already_there(editor, monkeypatch):
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "xtb_available", lambda: True)

    editor["submit_ff_dd"].value = "gfn2"

    assert editor["submit_xtb_install_btn"].layout.display == "none"


def test_the_first_press_says_what_the_second_one_would_do(editor, monkeypatch):
    """A few hundred megabytes through conda is not a thing to start on one
    click, and on a cluster the right answer is often the module instead."""
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_xtb", lambda: None)
    monkeypatch.setattr(tab_submit._gfn, "xtb_available", lambda: False)
    ran: list = []
    monkeypatch.setattr(tab_submit._gfn, "install_xtb",
                        lambda **kwargs: ran.append(kwargs) or {"ok": True})

    editor["submit_ff_dd"].value = "gfnff"
    editor["submit_xtb_install_btn"].click()

    assert ran == [], "the first press must not install anything"
    assert editor["submit_xtb_install_btn"].layout.display == "none"
    assert editor["submit_xtb_confirm_btn"].layout.display == ""
    assert editor["submit_xtb_cancel_btn"].layout.display == ""
    said = editor["mol_status"].value
    assert "install_qm_tools.sh xtb" in said, "the command has to be on screen"
    assert "megabytes" in said and "module load xtb" in said

    editor["submit_xtb_cancel_btn"].click()
    assert ran == []
    assert editor["submit_xtb_install_btn"].layout.display == ""
    assert editor["submit_xtb_confirm_btn"].layout.display == "none"


def test_the_second_press_runs_it_and_says_where_xtb_ended_up(editor, monkeypatch):
    import time as _time

    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_xtb", lambda: None)
    monkeypatch.setattr(tab_submit._gfn, "xtb_available", lambda: False)

    def fake_install(on_line=None, **_kwargs):
        if on_line:
            on_line("[qm_tools] link xtb")
        return {"ok": True, "binary": "/somewhere/bin/xtb", "lines": [],
                "status": "xtb installed at /somewhere/bin/xtb in 0.4 min."}

    monkeypatch.setattr(tab_submit._gfn, "install_xtb", fake_install)

    editor["submit_ff_dd"].value = "gfnff"
    editor["submit_xtb_install_btn"].click()
    editor["submit_xtb_confirm_btn"].click()

    deadline = _time.time() + 10
    while (_time.time() < deadline
           and "/somewhere/bin/xtb" not in editor["mol_status"].value):
        _time.sleep(0.02)

    assert editor["editor_state"].get("xtb_installing") is False
    assert "/somewhere/bin/xtb" in editor["mol_status"].value
    assert editor["mol_status"].value.count("Installing xtb...") == 0, (
        "the spinner has to go when the install is over"
    )


def test_an_install_that_produced_no_xtb_is_a_failure_not_a_shrug(monkeypatch):
    """Reporting success without a binary sends the user back to a button that
    already said it worked."""
    class _Fake:
        returncode = 0
        stdout = iter(["[qm_tools] ERROR: micromamba/conda not found\n"])

        def wait(self, timeout=None):
            return 0

    monkeypatch.setattr(gfn.subprocess, "Popen", lambda *a, **k: _Fake())
    monkeypatch.setattr(gfn, "find_xtb", lambda: None)

    outcome = gfn.install_xtb()

    assert outcome["ok"] is False
    assert "micromamba/conda not found" in outcome["status"]


@_needs_xtb
def test_the_installer_runs_and_ends_with_an_xtb_the_dashboard_can_find():
    """Run for real.  With an xtb already on the machine the script links that
    one instead of downloading -- which is the path a user whose cluster
    provides xtb takes, and the one that has to work without the network.
    """
    outcome = gfn.install_xtb(timeout=900)

    assert outcome["ok"] is True, outcome["status"]
    assert outcome["binary"] and pathlib.Path(outcome["binary"]).is_file()
    assert gfn.find_xtb() == outcome["binary"]


# ---------------------------------------------------------------------------
# the method that is chosen is the method that acts
# ---------------------------------------------------------------------------
def test_nothing_of_the_browsers_field_may_run_under_gfn(editor):
    """A dozen handlers install UFF parameters -- Hold, a polyhedron, a
    hybridisation, a bond edit.  Any one of them put a UFF relaxation under a
    molecule whose method box said GFN: it settled on release, and the geometry
    that reached the coordinate box was UFF's."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    enable = source.split("def _enable_live_forcefield")[1].split("\n    def ")[0]
    head = enable.split("\"\"\"", 2)[2]
    assert "if _gfn.is_gfn_method(submit_ff_dd.value):" in head
    assert head.index("_stop_browser_field()") < head.index("xyz =")

    toggle = source.split("def on_submit_relax_toggle")[1].split("\n    def ")[0]
    gfn_branch = toggle.split("if _gfn.is_gfn_method(submit_ff_dd.value):")[1]
    assert "_stop_browser_field()" in gfn_branch, (
        "a field left running from a UFF session goes on relaxing"
    )

    changed = source.split("def on_submit_ff_changed")[1].split("\n    def ")[0]
    assert "_stop_browser_field()" in changed


def test_the_follow_uses_the_method_that_is_on_screen(editor):
    """GFN-FF whatever the box said is a picture of a calculation nobody asked
    for -- and, from the outside, indistinguishable from the right one."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    began = source.split("def _begin_gfn_follow")[1].split("\n    def ")[0]
    assert "state['gfn_follow_method'] = str(submit_ff_dd.value)" in began

    follow = source.split("def _gfn_follow_step")[1].split("\n    def ")[0]
    assert "method=method" in follow, "the chosen method has to reach xtb"
    assert "'gfnff'" not in follow, "no method may be baked in here"
    assert "{label} is following the drag" in follow, (
        "and the status has to name the one that ran"
    )


def test_relax_steps_runs_the_method_it_is_given():
    """It was GFN-FF only, which is what let the wrong one run quietly."""
    import inspect

    signature = inspect.signature(gfn.relax_steps)
    assert "method" in signature.parameters
    assert signature.parameters["method"].default == "gfnff"

    source = inspect.getsource(gfn.relax_steps)
    assert "optimize_with_gfn(\n        xyz_text, method," in source


def test_dynamik_opt_goes_away_when_there_is_no_xtb_behind_it(editor, monkeypatch):
    """A switch that cannot do anything is worse than one that is not there:
    it is pressed, nothing happens, and nothing says why."""
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "xtb_available", lambda: False)

    editor["submit_ff_dd"].value = "gfnff"
    assert editor["submit_relax_btn"].layout.display == "none"
    assert editor["submit_relax_btn"].value is False

    monkeypatch.setattr(tab_submit._gfn, "xtb_available", lambda: True)
    editor["submit_ff_dd"].value = "gfn2"
    assert editor["submit_relax_btn"].layout.display == ""


def test_switching_methods_back_and_forth_re_arms_the_right_engine(editor):
    """Switching is something a user does constantly. It must not cost a press
    each time, nor leave the previous engine running under the new choice."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    changed = source.split("def on_submit_ff_changed")[1].split("\n    def ")[0]
    assert "submit_relax_btn.value = False" not in changed.split(
        "usable = not gfn")[0], "the switch keeps its position"
    assert "elif submit_relax_btn.value:" in changed
    assert "startAutoOptimize" in changed, "going back to UFF has to re-arm it"


def test_settle_under_gfn_is_the_chosen_method_tidying_up(editor):
    """Settle means: what reaches the coordinate box has relaxed around where
    the atom was put, not wherever the cursor stopped."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    settle = source.split("def _gfn_settle_now")[1].split("\n    def ")[0]
    assert "method = str(submit_ff_dd.value)" in settle
    assert "_GFN_SETTLE_CYCLES" in settle and "on_frames=_push" in settle
    assert "coords_widget.value = (" in settle, "the result has to land"

    handler = source.split("def on_submit_cmd")[1].split("\n    def ")[0]
    assert "_arm_gfn_settle()" in handler, "a release is what triggers it"

    toggle = source.split("def on_submit_settle_toggle")[1].split("\n    def ")[0]
    assert "if _gfn.is_gfn_method(submit_ff_dd.value):" in toggle
    assert toggle.index("return") < toggle.index("_ensure_manip_bootstrap()"), (
        "the browser must not be told to settle with a field it does not have"
    )


@_needs_xtb
def test_letting_go_settles_with_the_method_that_is_chosen(editor):
    """Driven the way the page drives it: a release with Settle on."""
    import time as _time

    refs = editor
    state = refs["editor_state"]
    strained = (
        "9\npropane, one carbon pulled out\n"
        "C -1.66 0.00 0.00\nC 0.00 0.86 0.00\nC 1.26 0.00 0.00\n"
        "H -2.15 0.63 0.00\nH -1.30 -0.64 0.88\nH 0.00 1.50 0.89\n"
        "H 0.00 1.50 -0.89\nH 2.15 0.63 0.00\nH 1.30 -0.64 0.88\n"
    )
    refs["coords_widget"].value = strained
    state["current_xyz_for_copy"] = {"content": strained}
    refs["submit_ff_dd"].value = "gfnff"
    refs["submit_settle_btn"].value = True

    refs["submit_cmd_sync"].value = "gfnfree:1:"

    deadline = _time.time() + 60
    while _time.time() < deadline and "Settled with" not in refs["coords_widget"].value:
        _time.sleep(0.05)

    assert "Settled with GFN-FF" in refs["coords_widget"].value, (
        refs["mol_status"].value
    )
    rows = [line.split() for line in gfn.atom_lines(refs["coords_widget"].value)]
    first = [float(v) for v in rows[0][1:4]]
    second = [float(v) for v in rows[1][1:4]]
    length = sum((a - b) ** 2 for a, b in zip(first, second)) ** 0.5
    assert length < 1.75, f"the stretched C-C was not tidied up: {length:.3f}"


# ---------------------------------------------------------------------------
# the atom under the cursor stays under the cursor
# ---------------------------------------------------------------------------
def test_the_held_atoms_are_put_back_where_they_were_sent():
    """xtb pulls a dragged atom most of the way home in five cycles -- 244 mA
    of a 250 mA drag -- and the answer outlives the drag.  Applied after the
    release, when nothing is being held any more, it took the atom with it."""
    sent = ("3\nsent\nC 5.000000 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n")
    answered = ("3\nanswered\nC 1.100000 0.2 0.1\nH 1.4 0.1 0.0\n"
                "H 0.1 1.3 0.0\n")

    kept = gfn.hold_atoms_at(answered, sent, [0])

    rows = [line.split() for line in gfn.atom_lines(kept)]
    assert rows[0][1:4] == ["5.000000", "0.0", "0.0"], "the held atom moved"
    assert rows[1][1:4] == ["1.4", "0.1", "0.0"], "the rest must keep the answer"

    both = gfn.hold_atoms_at(answered, sent, [0, 2])
    rows = [line.split() for line in gfn.atom_lines(both)]
    assert rows[0][1] == "5.000000" and rows[2][1:4] == ["0.0", "1.0", "0.0"], (
        "a selection of several atoms is dragged as one, and held as one"
    )

    assert gfn.hold_atoms_at(answered, sent, []) == answered
    # a structure that changed size under it is left alone rather than mangled
    assert gfn.hold_atoms_at(answered, "2\nx\nH 0 0 0\nH 1 0 0\n", [0]) == answered


def test_the_page_names_the_atoms_the_hand_is_on(player_js):
    """Python cannot know which atoms are being dragged; the page can."""
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    push = submit_manip_bootstrap_js().split("function pushXyzToPython(")[1][:1400]
    assert "held=" in push and "ffIndicesOf(viewer, targets)" in push

    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    sync = source.split("def on_submit_manip_sync")[1].split("\n    def ")[0]
    assert "word.startswith('held=')" in sync
    assert "_gfn_follow_step(new_xyz, holding)" in sync

    follow = source.split("def _gfn_follow_step")[1].split("\n    def ")[0]
    assert "_gfn.hold_atoms_at(" in follow


def test_the_follow_is_paced_by_the_machine_not_by_a_clock(player_js):
    """A fixed fifth of a second threw nine tenths of GFN-FF away: it answers
    a small molecule in under twenty milliseconds.  Measured in a real engine,
    one second of dragging went from 5 answers to 16."""
    assert "answered=play.seen>(play.pushedAt||0)" in player_js
    assert "if(since>500||(answered&&since>50)){" in player_js, (
        "as soon as the last answer landed, with a floor and a ceiling"
    )
    # and drawn over exactly as long as the next one takes to arrive
    assert "if(play.follow&&play.gap) return play.gap;" in player_js
    assert "play.gap=Math.min(600," in player_js


@_needs_xtb
def test_a_dragged_atom_comes_back_exactly_where_it_was_put(editor):
    """The whole of the spring-back, driven the way the page drives it."""
    import json as _json
    import time as _time

    refs = editor
    state = refs["editor_state"]
    propane = (
        "9\npropane\n"
        "C -1.26 0.00 0.00\nC 0.00 0.86 0.00\nC 1.26 0.00 0.00\n"
        "H -2.15 0.63 0.00\nH -1.30 -0.64 0.88\nH 0.00 1.50 0.89\n"
        "H 0.00 1.50 -0.89\nH 2.15 0.63 0.00\nH 1.30 -0.64 0.88\n"
    )
    refs["coords_widget"].value = propane
    state["current_xyz_for_copy"] = {"content": propane}
    refs["submit_ff_dd"].value = "gfnff"
    refs["submit_relax_btn"].value = True
    refs["submit_cmd_sync"].value = "gfngrab:1:"

    rows = [line.split() for line in gfn.atom_lines(propane)]
    rows[0][1] = f"{float(rows[0][1]) - 0.30:.6f}"
    pushed = ("9\nDELFIN drag-follow held=0\n"
              + "\n".join(" ".join(r) for r in rows) + "\n")
    refs["submit_manip_sync"].value = pushed

    deadline = _time.time() + 30
    while _time.time() < deadline and not (
            _json.loads(refs["submit_gfn_frame"].value or "{}").get("frames")):
        _time.sleep(0.01)
    frames = _json.loads(refs["submit_gfn_frame"].value)["frames"]
    assert frames, refs["mol_status"].value

    answered = frames[-1]
    assert abs(answered[0] - float(rows[0][1])) < 1e-6, (
        f"the held atom came back at {answered[0]:.4f}, not "
        f"{float(rows[0][1]):.4f} -- that is the spring back"
    )
    moved = max(
        abs(answered[3 * i + n] - float(rows[i][1 + n]))
        for i in range(1, 9) for n in range(3)
    )
    assert moved > 0.005, "and everything else has to have followed"


@_needs_xtb
def test_holding_a_value_moves_the_structure_to_it_there_and_then(editor):
    """Point of the whole thing, driven the way the buttons drive it: Hold is
    pressed, Optimise is not, and the angle is what it was asked to be."""
    import math
    import time as _time

    refs = editor
    state = refs["editor_state"]
    propane = (
        "9\npropane\n"
        "C -1.26 0.00 0.00\nC 0.00 0.86 0.00\nC 1.26 0.00 0.00\n"
        "H -2.15 0.63 0.00\nH -1.30 -0.64 0.88\nH 0.00 1.50 0.89\n"
        "H 0.00 1.50 -0.89\nH 2.15 0.63 0.00\nH 1.30 -0.64 0.88\n"
    )

    def carbon_angle(xyz):
        points = [[float(v) for v in line.split()[1:4]]
                  for line in gfn.atom_lines(xyz)]
        first = [points[0][n] - points[1][n] for n in range(3)]
        third = [points[2][n] - points[1][n] for n in range(3)]
        dot = sum(a * b for a, b in zip(first, third))
        norms = (math.dist(points[0], points[1])
                 * math.dist(points[2], points[1]))
        return math.degrees(math.acos(max(-1.0, min(1.0, dot / norms))))

    refs["coords_widget"].value = propane
    state["current_xyz_for_copy"] = {"content": propane}
    refs["submit_ff_dd"].value = "gfnff"
    refs["submit_relax_btn"].value = True
    deadline = _time.time() + 60
    while _time.time() < deadline and "converged after" not in refs["mol_status"].value:
        _time.sleep(0.05)

    state["picked"] = [0, 1, 2]
    refs["submit_internal_value"].value = 95.0
    refs["submit_hold_mode"].value = "fix"
    refs["submit_hold_btn"].click()

    deadline = _time.time() + 60
    while _time.time() < deadline and "converged after" not in refs["mol_status"].value:
        _time.sleep(0.05)

    held = carbon_angle(refs["coords_widget"].value)
    assert abs(held - 95.0) < 1.5, (
        f"asked for 95 deg and got {held:.2f} without pressing Optimise"
    )
    refs["submit_relax_btn"].value = False


# ---------------------------------------------------------------------------
# a relaxation that will not end, and a spin that cannot exist
# ---------------------------------------------------------------------------
def test_a_multiplicity_the_electrons_cannot_make_is_refused():
    """xtb does not refuse it.  Water asked for a doublet came back with the
    singlet's energy to the last digit -- -5.070543980552 either way -- and a
    confidently mislabelled answer is the one failure this module exists to
    prevent."""
    refused = gfn.optimize_with_gfn(_WATER, "gfn2", charge=0, uhf=1)

    assert refused["ok"] is False
    assert "even number of electrons" in refused["status"]
    assert "cannot make M = 2" in refused["status"]
    assert refused["xyz"] == _WATER

    # and the ones it can make are not refused
    for uhf in (0, 2):
        assert "cannot make" not in gfn.optimize_with_gfn(
            _WATER, "gfn2", charge=0, uhf=uhf)["status"]
    # an odd electron count is the other way round
    assert gfn.optimize_with_gfn(_WATER, "gfn2", charge=1, uhf=0)["ok"] is False


def test_how_far_the_structure_moved_is_asked_of_the_coordinates():
    """Whether a relaxation is still getting anywhere is a different question
    from whether xtb calls it converged."""
    before = "2\nx\nH 0.0 0.0 0.0\nH 1.0 0.0 0.0\n"
    after = "2\nx\nH 0.0 0.0 0.0\nH 1.1 0.0 0.0\n"

    assert abs(gfn.largest_shift(before, after) - 0.1) < 1e-9
    assert gfn.largest_shift(before, before) == 0.0
    assert gfn.largest_shift(before, "3\nx\nH 0 0 0\nH 1 0 0\nH 2 0 0\n") == 0.0


def test_the_relaxation_ends_three_ways_and_says_which(editor):
    """Converged, standing still, or out of rounds.  A structure that has
    stopped improving and one that is finished look identical, and only one of
    them is worth pressing Optimise on."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    settle = source.split("def _gfn_settle_now")[1].split("\n    def ")[0]
    assert "rounds < _GFN_SETTLE_ROUNDS" in settle
    assert "moved > _GFN_SETTLE_STILL" in settle
    assert "stopped without converging" in settle
    assert "Held values that cannot all be met at once" in settle


def test_the_whole_relaxation_is_one_run_not_one_per_round(editor):
    """A new run number resets the player: it drops what it had not drawn and
    applies the next frame outright instead of moving to it.  With a round
    every few tenths of a second that is a twitch per round."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    settle = source.split("def _gfn_settle_now")[1].split("\n    def ")[0]
    assert "if rounds == 1:" in settle, "only the first round takes a number"
    assert "state['gfn_settle_offset'] = 0" in settle
    assert "offset + len(walked) - len(trail)" in settle, (
        "and the frames of later rounds carry on from where the last stopped"
    )


def test_the_live_relaxation_asks_the_same_spin_optimise_did(editor):
    """Optimise scanning multiplicities while the live relaxation ran the box's
    fixed one is two answers about two different molecules -- and pressing
    Optimise afterwards then moves the structure again, for no visible reason.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    picked = source.split("def _gfn_uhf_now")[1].split("\n    def ")[0]
    assert "submit_gfn_autospin.value" in picked
    assert "state.get('gfn_scanned_uhf')" in picked

    handler = source.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0]
    assert "state['gfn_scanned_uhf'] = int(outcome['uhf'])" in handler

    for name in ("_gfn_follow_step", "_gfn_settle_now"):
        body = source.split(f"def {name}")[1].split("\n    def ")[0]
        assert "uhf = _gfn_uhf_now()" in body, name


def test_the_status_counts_the_atoms_the_hand_is_on(editor):
    """Grabbing an atom that is part of a selection drags the whole selection,
    so one left over from earlier makes every drag move everything -- which
    reads as the molecule fighting itself and is otherwise invisible."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    follow = source.split("def _gfn_follow_step")[1].split("\n    def ")[0]
    assert "holding {len(holding)} atoms" in follow


def test_letting_go_gives_the_relaxation_back_what_the_drag_took(editor):
    """While Dynamik Opt is down, carrying on after a release has nothing to do
    with Settle: the switch means the structure is being kept relaxed.

    It used to stop dead at the release unless Settle happened to be on as
    well, which is why pressing Optimise afterwards still had work to do.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def on_submit_cmd")[1].split("\n    def ")[0]
    free = handler.split("verb == 'gfnfree'")[1].split("return")[0]
    assert "_arm_gfn_takeup()" in free, "the relaxation has to carry on"
    assert "_arm_gfn_settle()" in free, "and Settle alone still tidies one"


@_needs_xtb
def test_a_release_relaxes_to_convergence_with_settle_switched_off(editor):
    """The complaint this closes: it reacted briefly, stopped where the drag
    left it, and pressing Optimise then went on optimising."""
    import time as _time

    refs = editor
    state = refs["editor_state"]
    propane = (
        "9\npropane\n"
        "C -1.26 0.00 0.00\nC 0.00 0.86 0.00\nC 1.26 0.00 0.00\n"
        "H -2.15 0.63 0.00\nH -1.30 -0.64 0.88\nH 0.00 1.50 0.89\n"
        "H 0.00 1.50 -0.89\nH 2.15 0.63 0.00\nH 1.30 -0.64 0.88\n"
    )
    refs["coords_widget"].value = propane
    state["current_xyz_for_copy"] = {"content": propane}
    refs["submit_ff_dd"].value = "gfnff"
    refs["submit_settle_btn"].value = False          # deliberately off
    refs["submit_relax_btn"].value = True
    deadline = _time.time() + 60
    while _time.time() < deadline and "converged after" not in refs["mol_status"].value:
        _time.sleep(0.05)

    refs["submit_cmd_sync"].value = "gfngrab:1:"
    now = state["current_xyz_for_copy"]["content"]
    rows = [line.split() for line in gfn.atom_lines(now)]
    rows[0][2] = f"{float(rows[0][2]) + 0.45:.6f}"
    body = "\n".join(" ".join(r) for r in rows)
    refs["submit_manip_sync"].value = f"9\nDELFIN drag-follow held=0\n{body}\n"
    _time.sleep(0.3)
    refs["submit_manip_sync"].value = f"9\nDELFIN drag-end\n{body}\n"
    refs["submit_cmd_sync"].value = "gfnfree:1:"

    deadline = _time.time() + 90
    while _time.time() < deadline and "converged after" not in refs["mol_status"].value:
        _time.sleep(0.05)
    assert "converged after" in refs["mol_status"].value, refs["mol_status"].value

    # and there is nothing left for Optimise to find
    settled = state["current_xyz_for_copy"]["content"]
    again = gfn.optimize_with_gfn(settled, "gfnff")
    assert gfn.largest_shift(settled, again["xyz"]) < 0.01, (
        "Optimise still had work to do after the relaxation said it was done"
    )
    refs["submit_relax_btn"].value = False
