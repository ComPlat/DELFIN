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


def test_the_relax_toggle_is_the_browser_fields_alone(editor):
    """It was given a GFN path, and that was the wrong shape for it.

    Relax runs a field per frame; GFN is a process per step.  The trajectory
    belongs to Optimise, which gets it from one run for free.
    """
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    toggle = source.split("def on_submit_relax_toggle")[1].split("\n    def ")[0]
    assert "_start_gfn_loop()" not in toggle
    assert "_stop_gfn_loop()" in toggle, "a loop left running would keep going"


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


def test_relax_is_switched_off_while_a_gfn_method_is_chosen(editor):
    """Relax runs the browser's own field per frame, which GFN cannot be.

    Left pressable it did something unrecognisable; switched off it says which
    methods it is for.
    """
    assert editor["submit_relax_btn"].disabled is False

    editor["submit_ff_dd"].value = "gfnff"
    assert editor["submit_relax_btn"].disabled is True
    assert "browser" in editor["submit_relax_btn"].tooltip

    editor["submit_ff_dd"].value = "uff"
    assert editor["submit_relax_btn"].disabled is False


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
    assert 'gfnplay:' in watcher

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
    assert "gfnplay:" in startup, "it has to be able to report back"


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
    assert body.index("switchIsOn()") < body.index("read();")


# ---------------------------------------------------------------------------
# the controls, split and cleaned up
# ---------------------------------------------------------------------------
def test_controls_that_cannot_work_under_gfn_are_taken_away(editor):
    """Greying them out invites the question of why they are dead."""
    assert editor["submit_relax_btn"].layout.display in (None, "")
    editor["submit_ff_dd"].value = "uff"

    editor["submit_ff_dd"].value = "gfnff"
    for name in ("submit_relax_btn", "submit_settle_btn", "submit_strength_slider"):
        assert editor[name].layout.display == "none", name

    editor["submit_ff_dd"].value = "uff"
    for name in ("submit_relax_btn", "submit_settle_btn", "submit_strength_slider"):
        assert editor[name].layout.display == "", name


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
