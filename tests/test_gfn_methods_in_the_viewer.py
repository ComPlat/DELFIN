"""GFN next to UFF: what Optimise minimises with, and what it needs to know.

UFF and MMFF94 live in the browser, which is what lets a drag follow the mouse.
GFN runs xtb on the server -- too slow for a drag, right for a press -- so
choosing it changes what Optimise does, not what dragging does.  And xtb needs
a charge and a spin: the charge can be read off a SMILES, the spin never can.
"""

from __future__ import annotations

import os
import pathlib
import shutil

import pytest

from delfin.dashboard import gfn_optimize as gfn
from delfin.dashboard.context import DashboardContext
from editor_source import (
    EDITOR_SOURCE, FULLSCREEN_CSS, FULLSCREEN_JS, SUBMIT_SOURCE,
)

_WATER = "3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"
_needs_xtb = pytest.mark.skipif(not shutil.which("xtb"), reason="xtb not installed")


# ---------------------------------------------------------------------------
# the runner
# ---------------------------------------------------------------------------


def test_the_methods_offered_are_the_ones_xtb_knows():
    assert set(gfn.GFN_METHODS) == {"gfnff", "gfn2", "gfn1", "gxtb"}
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
    assert values[:5] == ["uff", "mmff94", "gfnff", "gfn2", "gxtb"]
    # And the MOPAC ones behind them: measured against literature bond
    # lengths, PM6 is closer than GFN2 on small organics (5.0 against
    # 11.3 mA), and PM6-D3H4 keeps that while binding the water dimer that
    # plain PM6 lets come apart.
    assert values[5:] == ["pm6d3h4", "pm6", "pm7"]


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

    source = SUBMIT_SOURCE
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
    # The places are gathered in one list now, because the first of them is no
    # longer taken on trust: one that cannot optimise is stepped over.
    places = source.split("def _xtb_candidates")[1].split("\ndef ")[0]
    assert "find_tool_executable" in places, (
        "DELFIN's own resolver knows about qm_tools and XTBHOME; ask it first"
    )
    assert "shutil.which" in places
    assert "sys.prefix" in places
    finder = source.split("def find_xtb")[1].split("\ndef ")[0]
    assert "_xtb_candidates()" in finder and "judge_xtb" in finder


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

    source = SUBMIT_SOURCE
    body = source.split("def _set_ff_notes")[1].split("\n    def ")[0]
    assert "_server_method()" in body
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

    source = SUBMIT_SOURCE
    follow = source.split("def _gfn_follow_step")[1].split("\n    def ")[0]
    assert "while state.get('gfn_follow')" in follow, (
        "it has to stop when the drag does"
    )

    for name in ("on_submit_relax_toggle", "on_submit_ff_changed",
                 "on_submit_cmd"):
        body = source.split(f"def {name}")[1].split("\n    def ")[0]
        assert "_end_gfn_follow()" in body, name


def test_a_held_value_is_taken_up_without_pressing_optimise(editor):
    """Setting a value, watching nothing happen and having to press Optimise
    for it is the switch claiming to be live and not being it."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    for name in ("on_submit_hold", "on_submit_hold_mode", "_apply_internal_now"):
        body = source.split(f"def {name}")[1].split("\n    def ")[0]
        assert "_arm_gfn_takeup(" in body, name

    takeup = source.split("def _arm_gfn_takeup")[1].split("\n    def ")[0]
    assert "_gfn_live_is_on()" in takeup, "only while the switch is down"
    assert "_arm_gfn_settle(forced=True)" in takeup, (
        "asked for by a hand, so it is not a tidy-up that may be skipped"
    )


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
    body = editor_js[editor_js.index("function setPositions("):][:1800]
    assert "ffWritePositions" in body
    assert "redrawHighlights" in body


def test_the_frames_go_through_a_widget_not_through_run_js(editor):
    """run_js writes into one Output and clears it first.

    Twenty scripts a second overwrite each other before the page has rendered
    them: the relaxation ran, the last structure appeared, and nothing in
    between did.  A widget value is ordered, survives a background thread, and
    cannot be clobbered by the next one.
    """
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
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

    A run of them, not a long tail of what the page already holds: the reader
    runs at 60 Hz and the writer at 20, so eight cover a 400 ms gap with room
    to spare. Four hundred of them was 5.1 MB of JSON per push at 400 atoms,
    serialised in 196 ms on the kernel's own thread; eight is 0.12 MB and
    5.3 ms.
    """
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    follow = source.split("def _gfn_follow_step")[1].split("\n    def ")[0]
    assert "frames.append(" in follow
    assert "frames[-40:]" in follow, "one frame per write can be missed"

    handler = source.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0]
    # The same concern, kept, and the fixed eight gone with it.  Eight was a
    # guess at how many frames fit between two reads, and xtb beats it easily
    # -- a benzene runs 23 cycles in a fraction of a second -- so everything
    # past the last eight was never sent and the viewer showed a sample of the
    # path.  The window starts where the previous window started instead:
    # every frame goes out twice, which is the same insurance, and it is still
    # bounded at two reads' worth rather than growing with the path.
    assert "walked[-8:]" not in handler
    assert "walked[start:]" in handler
    assert "state['gfn_push_start']" in handler
    assert "walked[-400:]" not in handler, 'a tail of what the page holds'


def test_the_playback_interpolates_between_computed_frames(editor):
    """Twenty computed steps a second, shown as motion rather than as jumps.

    The positions in between are drawn, not calculated, and every frame the
    structure actually passed through is one end of an interpolation.
    """
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "a[i]+(b[i]-a[i])*t" in watcher, "no interpolation between frames"
    assert "play.queue" in watcher, "frames have to queue, or they are dropped"
    assert "drawn, not calculated" in watcher, "say which positions are which"


def test_the_bootstrap_is_on_the_page_before_anything_pushes(editor):
    """Without it there is no __delfinSubmitManip, and every push is a no-op --
    which is exactly what "Relaxing..." and nothing moving looked like."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    toggle = source.split("def on_submit_relax_toggle")[1].split("\n    def ")[0]
    gfn_branch = toggle.split("if _server_method():")[1]
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

    source = SUBMIT_SOURCE
    handler = source.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0]
    apply_step = handler.split("def _apply()")[1]
    assert "_write_coords(" in apply_step
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


def test_fullscreen_gets_the_status_line_because_it_gets_the_picture():
    """It travels inside the stack, so there is nothing to lend and nothing
    to copy.

    Fullscreen relocates its members by hand and ipywidgets knows nothing
    about such a move, so a status line borrowed for the big view did not come
    back to the small one. The answer used to be a second widget carrying the
    same text -- one that stays, one that goes. The line lies on the picture
    now and the picture is what travels: it leaves and returns with the thing
    it is about, in both tabs, and neither view is ever without it.
    """
    from delfin.dashboard import tab_orca_builder, tab_submit

    setter = SUBMIT_SOURCE.split(
        "def _set_mol_status")[1].split("\n    def ")[0]
    assert "mol_status.value = rendered_html" in setter
    assert "mol_status_fs.value = '' if prompt else rendered_html" in setter, (
        "a tab that still places the twin has to get the same text"
    )

    for module in (tab_submit, tab_orca_builder):
        source = open(module.__file__, encoding='utf-8').read()
        stack = source.split('delfin-structure-viewer-stack')[0]
        stack = stack.rsplit('widgets.Box(', 1)[1]
        assert 'mol_status' in stack, (
            f'{module.__name__}: the line has to be in the stack, or the big '
            f'view has none')
        assert 'mol_output' in stack, f'{module.__name__}: and so does the view'
        assert "add_class('delfin-structure-fs-member')" in source

    # And the picture itself is not a member on its own: lifted straight out
    # of the stack, it would leave the line behind on the page.
    builder = open(tab_orca_builder.__file__, encoding='utf-8').read()
    marks = builder.split('orca_mol_output.add_class')
    assert not any(m.startswith("('delfin-structure-fs-member')")
                   for m in marks[1:]), 'the stack is the member, not the view'


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

    source = SUBMIT_SOURCE
    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "outcome.get('frames')" in handler
    assert "submit_gfn_frame" in handler
    assert "_install_gfn_frame_watcher" in handler
    # both ways of showing a result rebuild the viewer, and both would tear
    # the playback down
    # The set goes back the way it came -- to the tab that keeps blocks as
    # blocks, or to the isomer stepper -- and showing one re-renders, which is
    # why a running playback keeps the results without being shown.
    assert "show=not played[0]" in handler, "the isomer path re-renders too"
    assert "_offer_isomers(results, show=not played[0])" in handler, (
        "the set goes back to the tab whatever the picture is doing -- held "
        "back, every optimised geometry was lost in a tab that keeps blocks")


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
    # And it says what xtb said, not that a run ended. "Error termination.
    # Backtrace:" is an announcement; the reason is a different line.
    said = result["status"].lower()
    assert "error termination" not in said, result["status"]
    assert "backtrace" not in said, result["status"]
    assert len(result["status"]) > 40, "there is more to say than this"
    assert result.get("output"), "the output is kept, or nobody can look"


def test_each_run_is_told_apart_so_a_short_one_still_plays(editor):
    """The player counted the frames it had seen, and the count carried over.

    A run with fewer frames than the one before it therefore played nothing --
    which is why the playback looked like it worked only sometimes.
    """
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "run!==play.run" in watcher, "a new run has to reset the count"
    assert "play.seen=0" in watcher

    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "'run': r" in handler, "every payload has to name its run"
    assert "state['gfn_run']" in handler


def test_leaving_fullscreen_puts_every_member_back():
    """In every tab, not only the one that was burned by it.

    Fullscreen moves the widgets into an overlay and the overlay is removed on
    exit -- so a member whose recorded parent was replaced in the meantime has
    nowhere to go back to and is carried out of the page with it.

    The Submit tab grew a rescue for this; the other three had a fullscreen of
    their own and did not. Driven in chromium with the box under the members
    replaced while the overlay was open, before and after the two were made
    one::

        tab                       before   after
        Submit                    back     back
        ORCA Builder              lost     back
        Calculations browser      lost     back
        Remote archive            lost     back

    "lost" is every member of the molecule panel at once -- toolbar, status
    line and viewer -- gone from the tab with nothing to say they existed.
    """
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    exit_body = FULLSCREEN_JS.split('function exitFullscreen')[1].split(
        '\n    function ')[0]
    assert 'insertBefore' in exit_body and 'appendChild' in exit_body
    assert 'isConnected' in exit_body, 'an orphaned member is a lost control'
    assert 'home.appendChild(member)' in exit_body

    # One rescue for every tab: there is one exitFullscreen in the dashboard.
    assert FULLSCREEN_JS.count('function exitFullscreen') == 1
    assert 'exitFullscreen' not in submit_manip_bootstrap_js(), (
        'the editor carried a second fullscreen of its own, and a fix to one '
        'was never a fix to the other')


def test_the_finished_geometry_does_not_tear_down_the_playback(editor):
    """Writing the coordinates the ordinary way rebuilds the viewer.

    Done while the trajectory was playing, that destroyed it milliseconds
    after it started -- which is why only the end of the optimisation was ever
    seen.  The playback's last frame is that geometry already.
    """
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "played = [False]" in handler
    assert "played[0] = True" in handler
    apply_body = handler.split("def _apply")[1]
    assert "if played[0]:" in apply_body
    # The flag rides with the write now, in the one function that puts a
    # geometry in the box -- raised only when the value really changes, or it
    # stays up and swallows the next redraw instead.
    assert "drawn=played[0]" in apply_body
    # the flag has to be set before the write that would re-render
    # Two writes live in _apply now -- one cuts a run a hand interrupted, the
    # other ends a finished one -- and this is about the second.  Measured
    # from where the interrupt block returns, so the anchor cannot drift onto
    # the wrong one again.
    finished = apply_body.split("'stopped where you took hold'")[-1]
    assert (finished.index("if played[0]:")
            < finished.index("_write_coords("))


def test_the_playback_finds_its_field_in_fullscreen_but_only_its_own(player_js):
    """In every root carrying the scope, and nowhere else.

    Fullscreen moves the toolbar into an overlay carrying the same scope class,
    so one root is not enough -- that is what the loop over roots is for, and
    the playback showed nothing in the big view before it was there.

    Behind the loop stood the whole document, which was written down as safe
    while there was one editor per dashboard.  There are two.  A page-wide
    lookup answers with the first field in document order, which belongs to
    whichever editor was written out first -- so an editor that could not find
    its own would play the other one's trajectory into its viewer.

    Two players driven in chromium on one page, each field naming its own run,
    read back as the run each player believes it is playing::

                              before        after
        both fields present   A-run B-run   A-run B-run
        A's toolbar in an     A-run B-run   A-run B-run
          overlay
        B's field removed     A-run A-run   A-run None

    The third row is the point, and the second is what must not break to get
    it: B is told there is none rather than handed A's.
    """
    reader = player_js.split("function read(")[1].split("function ")[0]

    assert "var field=inScope(" in reader, (
        "the reader walks every root carrying the scope, which is the helper")
    assert "document.querySelector(" not in reader, (
        "no page-wide fallback: with two editors that is the other one's field")

    found = player_js.split("function inScope(selector){")[1].split("function ")[0]
    assert 'document.querySelectorAll("."+scope)' in found
    assert "roots[i].querySelector(selector)" in found
    assert "document.querySelector(selector)" not in found


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

    source = SUBMIT_SOURCE
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

    source = SUBMIT_SOURCE
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    for report in ('"received "', '"drawing"', '"setPositions did not draw"',
                   '"no setPositions on the page"'):
        assert report in watcher, f"the page never reports {report}"
    assert 'send("gfnplay", text)' in watcher
    assert 'verb+":"+play.serial' in watcher, "the line has to name its verb"

    handler = source.split("def on_submit_cmd")[1].split("\n    def ")[0]
    assert "verb == 'gfnplay'" in handler
    assert "state['gfn_play_note'] = str(payload)" in handler
    assert "_gfn_status_lines()" in handler, (
        "and it says it in the shape the follow step says its half in")


@pytest.fixture
def bare_editor(tmp_path):
    """One structure editor, over a coordinate box of its own.

    Built here rather than reached through a tab: what the status line does is
    the part's, and the tab hands out only the widgets it places.
    """
    pytest.importorskip("ipywidgets")
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = lambda _script: None
    state = {}
    part = structure_editor.build(
        ctx,
        state=state,
        coords_widget=widgets.Textarea(value=_WATER),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None,
    )
    return part, state


def test_the_two_ends_of_a_drag_write_one_message_not_two(bare_editor):
    """A drag is reported from both ends, and they used to fight over the line.

    The kernel counts the follow steps it has run; the page counts the frames
    it has drawn. Each wrote the status line for itself and wiped the other:
    one row with a spinner, then two rows without, then one again, several
    times a second for as long as the drag lasted. The status line stands
    above the viewer in the column, so the structure the user was aiming at
    stepped up and down with it.

    Measured in chromium at 1440x900 against the tab's own stylesheet, one
    drag, the two ends reporting in turn::

                          before                     after
        follow step 14    18 px  1 row,  spinner     18 px  1 row, spinner
        page drew 14      35 px  2 rows, no spinner  18 px  1 row, spinner
        follow step 15    18 px  1 row,  spinner     18 px  1 row, spinner
        page drew 15      35 px  2 rows, no spinner  18 px  1 row, spinner

        movement per report:  17 px            ->    0 px

    The page's count is gone from the line entirely: it was there to tell an
    invisible trajectory from a missing one, and "received 41 frames" is the
    case where nothing is wrong. What the playback could not do still gets
    said, on the end of the line that is already there.
    """
    part, state = bare_editor
    part.submit_ff_dd.value = "gxtb"
    # The switch is detached from its own observer: that one refuses to stay
    # on where xtb is not installed, which is a different question.
    part.submit_relax_btn.unobserve_all()
    part.submit_relax_btn.value = True
    assert part._begin_gfn_follow(), "the follow did not start"
    assert not state["gfn_play_note"], (
        "whatever the page said about the last drag is not about this one")

    said = "g-xTB is following the drag: 15 step(s), 973 ms each."
    state["gfn_last_status"] = said
    part._set_mol_status(*part._gfn_status_lines(said), spinner=True)
    from_kernel = part.mol_status.value

    part.submit_cmd_sync.value = "gfnplay:7:received 15 frames"
    from_page = part.mol_status.value

    for name, html_value in (("kernel", from_kernel), ("page", from_page)):
        assert html_value.count("<br>") == 0, f"{name}: more than one row"
        assert "delfin-busy" in html_value, f"{name}: no spinner"
        assert "is following the drag" in html_value, name
    # The page's count is not a row. It said the playback was working, several
    # times a second, under a viewer that moved every time it appeared.
    assert "received 15 frames" not in from_page
    assert from_page == from_kernel, (
        "the page's report must not change what the line looks like at all")

    # A fault does get said -- that is the case the row was added for -- and it
    # goes on the end of the line rather than under it.
    part.submit_cmd_sync.value = "gfnplay:8:setPositions did not draw"
    faulted = part.mol_status.value
    assert "setPositions did not draw" in faulted
    assert faulted.count("<br>") == 0, "still one row"


def test_the_status_line_lies_on_the_picture_and_costs_it_no_room(bare_editor):
    """Not a row above the viewer, where its height was the picture's place.

    Messages are not the same length -- one row while the structure follows the
    hand, two when the optimisation reports an energy and what it held, more
    when it refuses and explains itself. In a row above the picture, every one
    of those moved the structure the user was aiming at. Giving that row a
    fixed height stopped the movement by spending two rows of height, empty
    most of the time.

    On the picture it costs nothing: it is positioned against the stack, grows
    upwards into the view when there is more to say, and disappears when there
    is nothing. Measured in chromium, a 660 px panel and a 300 px picture::

        message           viewer top   viewer height   line height
        nothing yet       92           304             0
        follow step       92           304             26
        run finished      92           304             44
        a long refusal    92           304             61

    The first two columns are the point: the picture does not move and does not
    change size, whatever is being said over it. The third is what used to be
    taken out of the layout.
    """
    part, _state = bare_editor
    for status in (part.mol_status, part.mol_status_fs):
        assert 'delfin-structure-status-over' in status._dom_classes
        assert not status.layout.height, (
            'a height here is layout the picture pays for')

    rules = FULLSCREEN_CSS.split('.delfin-structure-viewer-stack {')[1]
    assert 'position: relative' in rules.split('}')[0], (
        'the stack is what the line is positioned against')
    over = FULLSCREEN_CSS.split(
        '.delfin-structure-viewer-stack > .delfin-structure-status-over {')[1]
    over = over.split('}')[0]
    assert 'position: absolute' in over
    assert 'bottom:' in over, 'anchored to the bottom, so it grows upwards'
    assert 'pointer-events: none' in over, (
        'the structure underneath is being dragged')


def test_letting_go_of_an_atom_does_one_thing_or_the_other_and_says_which(
        bare_editor):
    """Auto: down to a minimum on release, or stop where the hand left it.

    It used to be neither, reliably. A drag that interrupted a running Optimise
    brought that run back and the structure went down to a minimum; the same
    drag with no run behind it got Settle's short tidy and nothing else. One
    gesture, two outcomes, decided by a switch on the other side of the
    toolbar -- which from the user's seat is it finishing the structure only
    sometimes.

    Off is the other mode, asked for by name: it stops, you move what else you
    want to move, and Optimize takes it down when the structure is what you
    meant.
    """
    part, state = bare_editor
    assert part.submit_auto_btn.value is True, 'the mode that finishes is first'
    assert part.submit_auto_btn in part.submit_manip_toolbar.children

    # A run was going when the atom was picked up. On, it comes back.
    part.submit_optimize_btn.unobserve_all()
    part.submit_optimize_btn.value = True
    state['optimize_interrupted'] = 'token'
    part._arm_gfn_restart()
    assert state.get('gfn_restart_armed') is True
    assert part.submit_optimize_btn.value is True, 'it is still running'


def test_with_auto_off_nothing_carries_on_and_the_switch_says_so(bare_editor):
    """Off, the run does not come back -- and the switch goes up with it.

    Left lit over a structure nobody is minimising, the switch says a
    calculation is running that is not.
    """
    part, state = bare_editor
    part.submit_auto_btn.value = False
    part.submit_optimize_btn.unobserve_all()
    part.submit_optimize_btn.value = True
    state['optimize_interrupted'] = 'token'

    part._arm_gfn_restart()

    assert state.get('gfn_restart_armed') is None, 'nothing was armed'
    assert state.get('optimize_interrupted') is None, 'and nothing is pending'
    assert part.submit_optimize_btn.value is False, 'the switch must not lie'
    assert 'press Optimize' in part.mol_status.value, (
        'the way back has to be said, or the structure just stops')


def test_a_drag_with_no_run_behind_it_still_reaches_a_minimum(bare_editor):
    """The half that was missing: nothing to resume is not nothing to do.

    _minimise_now is what the wait lands on. It presses the switch rather than
    calling the handler, because the switch has to be seen to be down for as
    long as the run lasts -- it is what the user presses to stop it again.
    """
    part, state = bare_editor
    part.submit_optimize_btn.unobserve_all()
    assert part.submit_optimize_btn.value is False

    part._minimise_now()
    assert part.submit_optimize_btn.value is True, 'Auto is on; one should start'


def test_with_auto_off_a_drag_starts_nothing_at_all(bare_editor):
    part, _state = bare_editor
    part.submit_auto_btn.value = False
    part.submit_optimize_btn.unobserve_all()

    part._minimise_now()

    assert part.submit_optimize_btn.value is False, 'Auto is off; none should'


def test_the_restart_is_gated_where_every_path_goes_through_it():
    """Three things ask for the restart, not one.

    A release asks for it; so does the arrival of the coordinates a drag
    produced, which happens first; so does adding or removing an atom. Gating
    only the release would have left the switch working for a drag that moved
    nothing and not for one that did.
    """
    source = EDITOR_SOURCE
    assert source.count('_arm_gfn_restart()') >= 3
    guard = source.split('def _arm_gfn_restart')[1].split('\n    def ')[0]
    assert 'submit_auto_btn.value' in guard
    assert '_stand_down_after_interrupt()' in guard


def test_one_press_keeps_running_until_there_is_nothing_left_to_do(bare_editor):
    """A press was one xtb run, and one run is often not a minimum.

    xtb's optimiser has a cycle limit of its own. At the limit it hands back
    the geometry it reached and reports that it did not converge -- and that
    was taken as the end: the switch went up over a structure better than it
    had been and not at a minimum, and the user pressed again. Pressing again
    worked because it is a NEW process from that geometry, with a fresh cycle
    budget and a fresh optimiser history. Nothing was stuck; the run had run
    out of room.

    Measured on a pulled-about propane under GFN-FF at three cycles a run,
    each row one press::

        press  converged  C-C      energy / Eh   largest shift
        start  -          1.735
        1      False      1.509    -1.514088     0.253 A
        2      False      1.523    -1.515532     0.066 A
        3      True       1.523    -1.515536     0.002 A

    The pressing is done for the user now, and it ends the three ways a settle
    has always ended: converged, no longer moving, or out of rounds. The last
    two matter -- held values that cannot all be met at once never converge,
    and a run that will not end is a process per round for as long as the
    switch is down.
    """
    part, _state = bare_editor
    part.submit_optimize_btn.unobserve_all()
    part.submit_optimize_btn.value = True
    asks = dict(converged=False, moved=0.25, rounds=1, failed=False,
                every_frame=False, still=0)

    assert part._optimise_carries_on(**asks) is True, (
        'a run that stopped for want of cycles is exactly the case to carry on')
    assert part._optimise_carries_on(**{**asks, 'converged': True}) is False
    assert part._optimise_carries_on(**{**asks, 'failed': True}) is False, (
        'a run that produced no geometry ends the press')
    # One quiet round is not proof.  A run that reached its cycle limit in a
    # flat stretch moves almost nothing, and its successor -- a new process
    # with a fresh cycle budget and a fresh optimiser history -- can still
    # take a real step.  Stopping on the first made "no longer moving" mean
    # "did not move much just then", and the user pressed again and watched it
    # move: the very thing these rounds were written to stop.
    assert part._optimise_carries_on(**{**asks, 'moved': 0.0, 'still': 1}) is True, (
        'one quiet round earns the fresh optimiser one more go')
    assert part._optimise_carries_on(**{**asks, 'moved': 0.0, 'still': 2}) is False, (
        'two in a row is a structure that has stopped improving')
    assert part._optimise_carries_on(**{**asks, 'rounds': 99}) is False, (
        'held values that cannot all be met never converge; it has to end')

    part.submit_optimize_btn.value = False
    assert part._optimise_carries_on(**asks) is False, (
        'the switch is up: the user asked it to stop')

    # And an engine that cannot say whether it converged counts as finished,
    # or it would be run again for ever.
    assert "outcome.get('converged', True)" in EDITOR_SOURCE


def test_a_run_that_stopped_short_is_not_counted_among_the_failures(
        bare_editor):
    """It is a geometry, and it is the reason the rounds exist.

    Written into the same list as a run that produced nothing -- which is
    where it went, because both end up on the status line -- it made the
    decision above answer False for the only situation it was written for.
    The rounds were there, the tests were green, and the user went on pressing
    Optimize because nothing ever carried on.
    """
    body = EDITOR_SOURCE.split(
        'def on_submit_optimize(change=None, every_frame=False)')[1]
    body = body.split('\n    def on_submit')[0]

    assert 'results, failures, unfinished = [], [], []' in body
    # Just that branch: the one beside it is a run that produced nothing, and
    # that one belongs among the failures.
    note = body.split("if 'before converging' in note:")[1]
    note = note.split('\n                else:')[0]
    assert 'unfinished.append(' in note
    assert 'failures.append(' not in note, (
        'counted as a failure, it stops the rounds it is supposed to start')
    # Both are still said, so a press that gives up still explains itself.
    assert '(failures + unfinished)[:2]' in body


def test_carrying_a_press_on_does_not_take_a_second_undo_point(bare_editor):
    """Undo after it should reach the geometry the user pressed on.

    Each round is another call into the same handler, so without this the
    undo history filled with one entry per round and Undo walked back through
    the optimisation a step at a time instead of leaving it.
    """
    part, state = bare_editor
    part.coords_widget.value = _WATER
    state['history'] = []

    part.submit_optimize_btn.unobserve_all()
    state['optimize_carrying_on'] = True
    part.submit_optimize_btn.value = True
    part.on_submit_optimize(None)
    carried = len(state.get('history') or [])

    state['optimize_carrying_on'] = False
    part.on_submit_optimize(None)
    pressed = len(state.get('history') or [])

    assert carried == 0, 'a round in the middle of a press is not a new start'
    assert pressed > carried, 'and a press of its own still is'


def test_the_fullscreen_copy_is_not_seen_next_to_the_original(editor):
    """Both lines carry the same text, so both visible prints it twice."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    assert "mol_status_fs.layout.display = 'none'" in source
    rule = FULLSCREEN_CSS.split(
        '.delfin-structure-fs-overlay .delfin-structure-fs-status {')[1]
    assert 'display: block !important;' in rule.split('}')[0]


def test_fullscreen_is_not_told_to_enter_coordinates(editor):
    """The copy reports work: a spinner, a trajectory, a result.

    In fullscreen there is a structure on screen, so a prompt to load one is a
    permanent line saying nothing.
    """
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
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

    source = SUBMIT_SOURCE
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

    source = SUBMIT_SOURCE
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

    source = SUBMIT_SOURCE
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "function stepMs()" in watcher
    assert "play.queue.length" in watcher
    assert "var ms=stepMs();" in watcher
    assert "var t=(now-play.started)/ms;" in watcher


def test_stopping_keeps_the_frame_that_was_on_screen(editor):
    """Stop means what the user was looking at.

    xtb runs ahead of the picture -- it produces frames faster than any frame
    rate shows them -- so keeping its newest geometry hands back something
    nobody saw and did not choose.  The page says which frame it was showing;
    that one is kept and the rest are discarded.
    """
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "play.shown=(play.shown||0)+1" in watcher, "nothing counts what was shown"
    assert '"stopped at frame "' in watcher

    cmd = source.split("def on_submit_cmd")[1].split("\n    def ")[0]
    assert "state['gfn_shown_frame']" in cmd

    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "gfn_shown_frame" in handler
    assert "stopped at the frame on" in handler


def test_the_backlog_is_the_path_and_is_kept(editor):
    """A queue allowed to grow puts the picture behind the run -- and being
    behind is what the pace control asks for.

    Twenty was the bound, and it was the second of two rules quietly sampling
    the trajectory.  Driven in a browser against a sixty-frame path arriving
    in bursts of twenty, it drew thirty-five and never showed the other
    twenty-five: the oldest of every burst, dropped to keep the bound.  See
    tests/test_the_player_draws_the_whole_path.py, which measures it rather
    than reading it.

    Where the picture has got to is where the user is: taking hold of an atom
    keeps the frame on screen and drops what was computed past it, so the
    frames in front of it are what is being chosen between.  A ceiling stays,
    far above any real path, because a queue with no ceiling is a leak.
    """
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "play.queue.slice(-20)" not in watcher
    assert "play.queue.slice(-100000)" in watcher


def test_the_page_stops_the_picture_without_asking_the_kernel(editor):
    """Waiting to be told the switch went off costs a round trip, and the
    playback ran on for the length of it.  ipywidgets marks a pressed toggle
    with mod-active; reading that is instant."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
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

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    editor["submit_ff_dd"].value = "uff"

    editor["submit_ff_dd"].value = "gfnff"
    assert editor["submit_strength_slider"].layout.display == "none"
    assert editor["submit_relax_btn"].layout.display in (None, "")
    # Settle went with the strength slider.  On the server it had become the
    # ordinary optimisation run on release, which is what Auto does, and two
    # switches for one behaviour is one too many -- the spare one was the one
    # still relaxing structures with everything visible switched off.
    assert editor["submit_settle_btn"].layout.display == "none"

    editor["submit_ff_dd"].value = "uff"
    assert editor["submit_strength_slider"].layout.display == ""
    assert editor["submit_settle_btn"].layout.display in (None, "")


def test_optimise_and_optimise_all_are_two_switches(editor):
    import ipywidgets as w

    one = editor["submit_optimize_btn"]
    every = editor["submit_optimize_all_btn"]
    assert isinstance(one, w.ToggleButton) and isinstance(every, w.ToggleButton)
    assert one.description == "Optimize" and every.description == "all"


def test_only_one_of_them_runs_at_a_time(editor):
    """A login node is shared; two sets of xtb processes is how it is noticed."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    handler = source.split("def on_submit_optimize(change=None, every_frame=False)")[1].split("\n    def ")[0]
    assert "other.value = False" in handler
    # and the frames of a set are walked in one loop, not started side by side
    assert "for position, item in enumerate(targets):" in handler
    assert "threading.Thread" in handler and handler.count("threading.Thread") == 1


def test_only_all_takes_the_whole_set(editor):
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
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

    source = SUBMIT_SOURCE
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

    source = SUBMIT_SOURCE
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

    source = SUBMIT_SOURCE
    watcher = source.split("def _install_gfn_frame_watcher")[1].split("\n    def ")[0]
    assert "function grabbed()" in watcher
    assert "_submitManipStateByScope" in watcher, "it reads the drag off the page"
    assert 'drag.kind==="translate"' in watcher
    assert 'send(held?"gfngrab":"gfnfree",' in watcher
    assert 'held?String(play.shown||0):""' in watcher, (
        "and it says which frame the picture stood on")
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

    source = SUBMIT_SOURCE
    handler = source.split("def on_submit_cmd")[1].split("\n    def ")[0]
    assert "verb == 'gfngrab'" in handler and "_interrupt_gfn()" in handler
    # A release goes through _after_release, which is where the one question
    # about what a release means is asked.
    assert "verb == 'gfnfree'" in handler and "_after_release()" in handler
    release = EDITOR_SOURCE.split("def _after_release")[1].split("\n    def ")[0]
    assert "_arm_gfn_restart()" in release

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
        "11\npropane\n"
        "C 1.16 0.48 -0.22\nC 0.13 -0.61 0.01\nC -1.26 -0.02 0.21\n"
        "H 2.15 0.03 -0.37\nH 0.92 1.07 -1.11\nH 1.22 1.16 0.63\n"
        "H 0.41 -1.20 0.89\nH 0.11 -1.29 -0.85\nH -1.99 -0.82 0.38\n"
        "H -1.28 0.64 1.08\nH -1.58 0.55 -0.66\n"
    )
    refs["coords_widget"].value = xyz
    state["current_xyz_for_copy"] = {"content": xyz}
    refs["submit_ff_dd"].value = "gfnff"

    handed: list[str] = []
    real = tab_submit._gfn.optimize_with_gfn

    def recording(text, method, **kwargs):
        # The single cycle that reads GFN-FF's bonding is not an optimisation
        # of anything; counting it here would make the first real run look
        # like the second.
        if kwargs.get("max_steps") != 1:
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
    parts[1] = f"{float(parts[1]) + 0.7:.6f}"
    body[0] = " ".join(parts)
    moved = f"{len(body)}\nDELFIN drag-end\n" + "\n".join(body) + "\n"
    refs["submit_manip_sync"].value = moved
    refs["submit_cmd_sync"].value = "gfnfree:2:"

    deadline = _time.time() + 30
    while _time.time() < deadline and len(handed) < 2:
        _time.sleep(0.02)
    assert len(handed) >= 2, "the run never started again"

    first = float(gfn.atom_lines(handed[0])[0].split()[1])
    second = float(gfn.atom_lines(handed[1])[0].split()[1])
    assert abs(second - (first + 0.7)) < 1e-4, (
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
        "11\npropane\n"
        "C 1.16 0.48 -0.22\nC 0.13 -0.61 0.01\nC -1.26 -0.02 0.21\n"
        "H 2.15 0.03 -0.37\nH 0.92 1.07 -1.11\nH 1.22 1.16 0.63\n"
        "H 0.41 -1.20 0.89\nH 0.11 -1.29 -0.85\nH -1.99 -0.82 0.38\n"
        "H -1.28 0.64 1.08\nH -1.58 0.55 -0.66\n"
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
        "11\npropane\n"
        "C 1.16 0.48 -0.22\nC 0.13 -0.61 0.01\nC -1.26 -0.02 0.21\n"
        "H 2.15 0.03 -0.37\nH 0.92 1.07 -1.11\nH 1.22 1.16 0.63\n"
        "H 0.41 -1.20 0.89\nH 0.11 -1.29 -0.85\nH -1.99 -0.82 0.38\n"
        "H -1.28 0.64 1.08\nH -1.58 0.55 -0.66\n"
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

    source = SUBMIT_SOURCE
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

    source = SUBMIT_SOURCE
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
        "11\npropane\n"
        "C 1.16 0.48 -0.22\nC 0.13 -0.61 0.01\nC -1.26 -0.02 0.21\n"
        "H 2.15 0.03 -0.37\nH 0.92 1.07 -1.11\nH 1.22 1.16 0.63\n"
        "H 0.41 -1.20 0.89\nH 0.11 -1.29 -0.85\nH -1.99 -0.82 0.38\n"
        "H -1.28 0.64 1.08\nH -1.58 0.55 -0.66\n"
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
    # It names the directory the resolver reads on the way, or it would
    # install beside its own file -- inside the package -- while the Settings
    # tab filled the user's copy, and only one of the two was searched.
    assert command[-1] == "xtb"
    assert "bash" in command and command[-2].endswith("install_qm_tools.sh")
    if command[0] == "env":
        assert any(part.startswith("DELFIN_QM_TOOLS_ROOT=") for part in command)
    assert str(script) in command


def test_the_offer_appears_only_when_it_is_needed(editor, monkeypatch):
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: None)
    monkeypatch.setattr(tab_submit._gfn, "find_xtb", lambda: None)
    monkeypatch.setattr(tab_submit._gfn, "xtb_available", lambda: False)

    assert editor["submit_xtb_install_btn"].layout.display == "none"

    editor["submit_ff_dd"].value = "gfnff"
    assert editor["submit_xtb_install_btn"].layout.display == ""
    assert "was not found" in editor["mol_status"].value

    editor["submit_ff_dd"].value = "uff"
    assert editor["submit_xtb_install_btn"].layout.display == "none", (
        "nothing to install for a field that runs in the browser"
    )


def test_the_offer_stays_away_when_xtb_is_already_there(editor, monkeypatch):
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")

    editor["submit_ff_dd"].value = "gfn2"

    assert editor["submit_xtb_install_btn"].layout.display == "none"


def test_the_first_press_says_what_the_second_one_would_do(editor, monkeypatch):
    """A few hundred megabytes through conda is not a thing to start on one
    click, and on a cluster the right answer is often the module instead."""
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: None)
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

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: None)
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

    source = SUBMIT_SOURCE
    enable = source.split("def _enable_live_forcefield")[1].split("\n    def ")[0]
    head = enable.split("\"\"\"", 2)[2]
    assert "if _server_method():" in head
    assert head.index("_stop_browser_field()") < head.index("xyz =")

    toggle = source.split("def on_submit_relax_toggle")[1].split("\n    def ")[0]
    gfn_branch = toggle.split("if _server_method():")[1]
    assert "_stop_browser_field()" in gfn_branch, (
        "a field left running from a UFF session goes on relaxing"
    )

    changed = source.split("def on_submit_ff_changed")[1].split("\n    def ")[0]
    assert "_stop_browser_field()" in changed


def test_the_follow_uses_the_method_that_is_on_screen(editor):
    """GFN-FF whatever the box said is a picture of a calculation nobody asked
    for -- and, from the outside, indistinguishable from the right one."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
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

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: None)

    editor["submit_ff_dd"].value = "gfnff"
    assert editor["submit_relax_btn"].layout.display == "none"
    assert editor["submit_relax_btn"].value is False

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    editor["submit_ff_dd"].value = "gfn2"
    assert editor["submit_relax_btn"].layout.display == ""


def test_switching_methods_back_and_forth_re_arms_the_right_engine(editor):
    """Switching is something a user does constantly. It must not cost a press
    each time, nor leave the previous engine running under the new choice."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    changed = source.split("def on_submit_ff_changed")[1].split("\n    def ")[0]
    assert "submit_relax_btn.value = False" not in changed.split(
        "usable = not gfn")[0], "the switch keeps its position"
    assert "elif submit_relax_btn.value:" in changed
    assert "startAutoOptimize" in changed, "going back to UFF has to re-arm it"


def test_settle_under_gfn_is_the_chosen_method_tidying_up(editor):
    """Settle means: what reaches the coordinate box has relaxed around where
    the atom was put, not wherever the cursor stopped."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    settle = source.split("def _gfn_settle_now")[1].split("\n    def ")[0]
    assert "method = str(submit_ff_dd.value)" in settle
    assert "max_steps=None" in settle, (
        "it is the ordinary optimisation, run on the frame that is on screen"
    )
    assert "on_frames=_push" in settle
    assert "_write_coords(" in settle, "the result has to land"

    # Not a release any more -- under a server method Settle is gone, and
    # going to a minimum when an atom is let go is Auto's, with one switch on
    # screen saying so.  What still reaches this machinery is a value the user
    # has just set or held: the structure has to arrange itself around it, and
    # that is asked for by the hand rather than by a switch.
    takeup = EDITOR_SOURCE.split("def _arm_gfn_takeup")[1].split("\n    def ")[0]
    assert "_arm_gfn_settle(forced=True)" in takeup

    # And the browser is told, before anything returns, that it is not the one
    # settling.  Leaving early instead -- on the grounds that a server method
    # has no browser field to settle with -- left the page settling on release
    # with a field installed under whatever method had been chosen before, so
    # a structure went on relaxing with every switch on screen switched off.
    toggle = source.split("def on_submit_settle_toggle")[1].split("\n    def ")[0]
    assert (toggle.index("setSettleOnRelease")
            < toggle.index("if _server_method():")), (
        "the browser hears about it before the server branch bows out"
    )
    assert "settling = active and not _server_method()" in toggle


@_needs_xtb


def test_letting_go_goes_to_a_minimum_with_the_method_that_is_chosen(editor):
    """Driven the way the page drives it: a release with Auto and Dynamik Opt
    on, which is the one thing on a server method that acts on a release."""
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
    # Auto and Dynamik Opt, not Settle: under a server method Settle is gone,
    # because what it ran there was the whole minimisation rather than a
    # tidy-up, and that is what these two already say they do.
    refs["submit_relax_btn"].value = True
    refs["submit_auto_btn"].value = True

    refs["submit_cmd_sync"].value = "gfnfree:1:"

    deadline = _time.time() + 60
    while (_time.time() < deadline
           and "Optimised in DELFIN viewer" not in refs["coords_widget"].value):
        _time.sleep(0.05)

    assert "Optimised in DELFIN viewer" in refs["coords_widget"].value, (
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

    source = SUBMIT_SOURCE
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
    assert "if(since>500||(answered&&since>floor)){" in player_js, (
        "as soon as the last answer landed, with a floor and a ceiling"
    )
    assert "var floor=Math.max(16,Math.min(120,(play.gap||60)/2));" in player_js, (
        "and the floor is the machine's own answer time, not a constant"
    )
    # and drawn over exactly as long as the next one takes to arrive -- but
    # only while they arrive one at a time; a burst keeps the backlog rules
    assert "if(play.follow&&play.gap&&n<=3) return play.gap;" in player_js
    assert "play.gap=play.gap?(play.gap*0.6+measured*0.4):measured;" in player_js, (
        "averaged, or the drawing speed jumps about as much as the arrivals do"
    )


@_needs_xtb


def test_a_dragged_atom_comes_back_exactly_where_it_was_put(editor):
    """The whole of the spring-back, driven the way the page drives it."""
    import json as _json
    import time as _time

    refs = editor
    state = refs["editor_state"]
    propane = (
        "11\npropane\n"
        "C 1.16 0.48 -0.22\nC 0.13 -0.61 0.01\nC -1.26 -0.02 0.21\n"
        "H 2.15 0.03 -0.37\nH 0.92 1.07 -1.11\nH 1.22 1.16 0.63\n"
        "H 0.41 -1.20 0.89\nH 0.11 -1.29 -0.85\nH -1.99 -0.82 0.38\n"
        "H -1.28 0.64 1.08\nH -1.58 0.55 -0.66\n"
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
        "11\npropane\n"
        "C 1.16 0.48 -0.22\nC 0.13 -0.61 0.01\nC -1.26 -0.02 0.21\n"
        "H 2.15 0.03 -0.37\nH 0.92 1.07 -1.11\nH 1.22 1.16 0.63\n"
        "H 0.41 -1.20 0.89\nH 0.11 -1.29 -0.85\nH -1.99 -0.82 0.38\n"
        "H -1.28 0.64 1.08\nH -1.58 0.55 -0.66\n"
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

    source = SUBMIT_SOURCE
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

    source = SUBMIT_SOURCE
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

    source = SUBMIT_SOURCE
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

    source = SUBMIT_SOURCE
    follow = source.split("def _gfn_follow_step")[1].split("\n    def ")[0]
    assert "holding {len(holding)} atoms" in follow


def test_a_run_that_ends_lands_the_picture_on_its_last_frame(player_js):
    """Otherwise the viewer keeps whatever it had drawn while the kernel keeps
    the geometry it computed, and the two drift apart.

    Driven in a real JS engine: a relaxation walking an atom from 0 to 29 was
    at 11.9 when the next drag started a new run, and the viewer stayed there.
    The next drag then hands over a structure that is behind, and the
    relaxation nobody saw is walked again -- which is every earlier drag being
    made a second time.
    """
    assert "show(play.last,play.queue[play.queue.length-1],1);" in player_js
    landing = player_js.split("if(run!==play.run){")[1].split("play.run=run;")[0]
    assert "if(play.queue.length){" in landing, (
        "the frames of the ending run have to be landed before they are dropped"
    )


def test_a_burst_is_not_played_at_the_pace_of_a_followed_hand(player_js):
    """A follow arrives one answer at a time and is drawn over the gap between
    them; a relaxation arrives in bursts, and those keep the backlog rules."""
    step = player_js.split("function stepMs(){")[1].split("\n  }")[0]
    assert "if(n>60) return 8;" in step
    assert step.index("if(n>10) return 35;") < step.index("play.gap"), (
        "the backlog rules come first"
    )
    assert "play.follow&&play.gap&&n<=3" in step


@_needs_xtb


def test_a_run_shorter_than_the_reading_interval_still_hands_its_path_over():
    """The watching loop reads the log five times a second at most.

    A settle of a small molecule is twenty milliseconds, so it handed over
    nothing at all: the picture never saw the relaxation it is a picture of,
    kept the geometry the drag had left, and gave that one to the next drag --
    which walked the whole path again, in front of the user, as though every
    earlier drag were being redone.
    """
    import time as _time

    strained = (
        "9\npropane, strained\n"
        "C -1.26 0.00 0.30\nC 0.00 1.16 0.00\nC 1.26 0.00 0.00\n"
        "H -2.15 0.63 0.00\nH -1.30 -0.64 0.88\nH 0.00 1.50 0.89\n"
        "H 0.00 1.50 -0.89\nH 2.15 0.63 0.00\nH 1.30 -0.64 0.88\n"
    )
    handed: list = []
    began = _time.perf_counter()
    result = gfn.optimize_with_gfn(
        strained, "gfnff", max_steps=40, timeout=60,
        should_stop=lambda: False, on_frames=handed.append,
    )
    spent = _time.perf_counter() - began

    assert result["ok"] is True, result["status"]
    assert spent < 0.2, f"the run has to be shorter than the interval: {spent:.3f} s"
    assert handed, "a run that fast handed over nothing at all"
    assert len(handed[-1]) == len(result["frames"]), (
        "and what it handed over has to be the whole path"
    )


def test_the_path_is_handed_over_once_at_the_end_whatever_the_clock_did():
    source = open(gfn.__file__, encoding="utf-8").read()
    runner = source.split("def optimize_with_gfn")[1].split("\ndef ")[0]
    assert "if on_frames is not None and len(frames) > sent:" in runner
    assert runner.index("sent = 0") < runner.index("while running.poll()")


def test_letting_go_does_not_walk_the_atom_back_through_the_drag(player_js):
    """The queued follow answers describe the drag: each carries the dragged
    atom where the hand had it when the frame was computed, and while the hand
    was down they were drawn around it.  Drawn after the release, with nothing
    held any more, they walk that atom back through the whole drag.

    Driven in a real JS engine: three answers queued at x = 1, 2, 3, the hand
    lets go, and the only thing drawn is 3 -- where the hand left it.
    """
    release = player_js.split("play.queue=[]; play.last=null;\n      }")[1]
    release = release.split("send(held")[0]
    assert "show(play.last,play.queue[play.queue.length-1],1);" in release, (
        "it has to land on the newest before dropping the rest"
    )
    assert "play.queue=[]; play.last=null;" in release


def test_a_live_run_shows_the_frame_that_is_current_not_the_way_there(player_js):
    """The frames before the newest describe where the structure was on the way
    here -- the drag that has just happened, or a relaxation of it.  Playing
    those is showing the user their own past.

    Driven in a real JS engine: a live run carrying a path of twelve frames
    draws frame twelve and nothing else.  An Optimise run still carries its
    path, because showing what the optimiser walked through is what that button
    is for.
    """
    assert "if(play.follow&&play.seen===0&&frames.length>1){" in player_js
    assert "frames=[frames[frames.length-1]];" in player_js
    assert "from=from+frames.length-1;" in player_js, (
        "the count has to stay honest, or the next payload replays"
    )


def test_a_past_drag_has_no_hold_on_the_present(editor):
    """Switching the toggle off and on again finished the job that letting go
    should have -- which is what a leftover from an earlier drag looks like.

    Every timer and every worker belongs to the structure it was started for.
    A drag, a Hold, a Set: each makes a new one, and whatever does not belong
    to the current generation does nothing at all.
    """
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    fresh = source.split("def _gfn_new_generation")[1].split("\n    def ")[0]
    assert "state['gfn_generation'] = int(state.get('gfn_generation', 0)) + 1" in fresh
    for flag in ("gfn_settle_forced", "gfn_settle_rounds", "gfn_settle_again"):
        assert f"state['{flag}']" in fresh, flag

    for name in ("_begin_gfn_follow", "_arm_gfn_takeup"):
        body = source.split(f"def {name}")[1].split("\n    def ")[0]
        assert "_gfn_new_generation()" in body, name

    armed = source.split("def _arm_gfn_settle")[1].split("\n    def ")[0]
    assert "if int(state.get('gfn_generation', 0)) != generation:" in armed, (
        "a timer for a structure that is history must not fire"
    )
    settle = source.split("def _gfn_settle_now")[1].split("\n    def ")[0]
    assert "int(state.get('gfn_generation', 0)) != generation" in settle, (
        "nor may a finished run write its geometry over a newer one"
    )


@_needs_xtb
def test_gfnff_loses_the_bond_it_was_never_told_to_keep(tmp_path):
    """The cliff, measured, and the reason the perception has to be kept.

    GFN-FF works its topology out of whatever geometry it is handed.  Pulling
    a carbon of a propane out relaxes back cleanly to 1.49 A while the bond is
    still seen -- and the moment it is not, the same relaxation pushes the two
    apart instead.  On a drag, where every step is a fresh perception, that is
    a molecule that can fall apart between one answer and the next.
    """
    propane = (
        "11\npropane\n"
        "C 1.16 0.48 -0.22\nC 0.13 -0.61 0.01\nC -1.26 -0.02 0.21\n"
        "H 2.15 0.03 -0.37\nH 0.92 1.07 -1.11\nH 1.22 1.16 0.63\n"
        "H 0.41 -1.20 0.89\nH 0.11 -1.29 -0.85\nH -1.99 -0.82 0.38\n"
        "H -1.28 0.64 1.08\nH -1.58 0.55 -0.66\n"
    )

    def pulled(by):
        rows = [line.split() for line in gfn.atom_lines(propane)]
        rows[0][1] = f"{float(rows[0][1]) + by:.6f}"
        return (f"{len(rows)}\npulled\n"
                + "\n".join(" ".join(r) for r in rows) + "\n")

    def carbon_carbon(xyz):
        rows = [line.split() for line in gfn.atom_lines(xyz)]
        first = [float(v) for v in rows[0][1:4]]
        second = [float(v) for v in rows[1][1:4]]
        return sum((a - b) ** 2 for a, b in zip(first, second)) ** 0.5

    far = pulled(0.9)
    assert carbon_carbon(far) > 2.2

    # perceived afresh: the bond is not there any more, and it is pushed apart
    fresh = gfn.relax_steps(far, cycles=5)
    assert carbon_carbon(fresh["xyz"]) > carbon_carbon(far), (
        "this is the failure the kept topology exists for"
    )

    # perceived once from the intact molecule, then kept: it pulls back
    keep = tmp_path / "topo"
    gfn.relax_steps(propane, cycles=5, topology=keep)
    assert (keep / "gfnff_topo").is_file(), "nothing was kept"
    held = gfn.relax_steps(far, cycles=5, topology=keep)
    assert carbon_carbon(held["xyz"]) < 1.8, (
        f"kept topology should pull back: {carbon_carbon(held['xyz']):.3f} A"
    )


def test_the_bonding_is_kept_for_the_molecule_it_was_perceived_from(editor):
    """It belongs to one molecule: an atom added or taken away makes it another
    one, and the count is what says so."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    keeper = source.split("def _gfn_topology_dir")[1].split("\n    def ")[0]
    assert "kept.get('who') == who" in keeper
    dropper = source.split("def _drop_gfn_topology")[1].split("\n    def ")[0]
    assert "shutil.rmtree" in dropper

    edit = source.split("def _apply_structure")[1].split("\n    def ")[0]
    assert "_drop_gfn_topology()" in edit, "a structural edit is another molecule"

    follow = source.split("def _gfn_follow_step")[1].split("\n    def ")[0]
    assert "topology=_gfn_topology_dir(" in follow
    settle = source.split("def _gfn_settle_now")[1].split("\n    def ")[0]
    assert "_gfn_topology_dir(xyz)" in settle
    assert "topology=perceived" in settle
    # And asking for one makes a directory, so the engine that has no topology
    # does not ask: MOPAC works the bonding out for itself every run.
    assert "if _gfn.is_gfn_method(method) else None" in settle


def test_only_gfnff_has_a_topology_to_keep():
    """GFN2 works the bonding out from the wavefunction every time."""
    source = open(gfn.__file__, encoding="utf-8").read()
    runner = source.split("def optimize_with_gfn")[1].split("\ndef ")[0]
    assert "keeping = topology is not None and key == 'gfnff'" in runner


def test_optimise_keeps_the_bonding_the_editor_has_been_working_with(editor):
    """Pressing Optimise on a structure that has been pulled about used to
    re-perceive the topology and shove the molecule apart -- the same cliff the
    drag had, in the one place a user reaches for when the drag went wrong."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    handler = source.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0]
    assert "topology=perceived" in handler
    assert "None if every_frame" in handler, (
        "a set of isomers is a set of different molecules, and one perception "
        "cannot describe them all"
    )


def test_optimise_supersedes_the_live_relaxation(editor):
    """Two xtb runs writing the same coordinate box is how the first press came
    to look like it had only worked out an energy: the live one landed after it
    and put its own, older geometry back, so a second press was needed."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    settle = source.split("def _gfn_settle_now")[1].split("\n    def ")[0]
    assert "or state.get('optimize_run') is not None" in settle, (
        "the live relaxation has to stop when Optimise starts"
    )
    assert settle.count("state.get('optimize_run') is not None") >= 3, (
        "it must also not start one, nor write its geometry over Optimise's"
    )


def test_dragging_moves_along_the_surface_and_optimise_goes_down_it(editor):
    """The division the editor settled on, and the reason for it.

    Dragging an atom moves the structure along the potential surface and
    leaves it where the hand left it -- an atom pulled off the place it was
    just put is the opposite of placing it.  Optimise is what goes downhill.
    Settle asks for that on every release, and is a choice rather than the
    default under GFN.
    """
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    release = EDITOR_SOURCE.split("def _after_release")[1].split("\n    def ")[0]
    assert "_arm_gfn_takeup()" not in release, (
        "letting go must not walk the structure downhill unasked"
    )
    assert "_arm_gfn_settle()" in release, "Settle is how that is asked for"
    # And with Auto on it is asked for: going to a minimum on release is a
    # mode with a switch on it, not something that happens unannounced.
    assert "submit_auto_btn.value" in release

    toggle = source.split("def on_submit_relax_toggle")[1].split("\n    def ")[0]
    gfn_branch = toggle.split("if _server_method():")[1]
    assert "_arm_gfn_takeup()" not in gfn_branch, (
        "switching it on arms the follow; it does not optimise"
    )

    controls = source.split("def _refresh_method_controls")[1].split("\n    def ")[0]
    assert "submit_settle_btn.value = False" in controls


def test_choosing_gfn_leaves_the_structure_where_the_hand_puts_it(editor):
    editor["submit_settle_btn"].value = True
    editor["submit_ff_dd"].value = "gfnff"
    assert editor["submit_settle_btn"].value is False, (
        "tidying on every release is a thing to ask for, not the default"
    )
    assert editor["submit_settle_btn"].layout.display == "none", (
        "and it goes with it: on the server Auto is the switch for what "
        "happens on release, and Settle sitting beside it did the same thing "
        "under a name that did not say so"
    )


def test_a_held_value_is_never_dropped_because_something_else_was_running(editor):
    """A value held while a relaxation was in flight did nothing at all until
    the toggle was switched off and on again."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    settle = source.split("def _gfn_settle_now")[1].split("\n    def ")[0]
    assert "state['gfn_settle_pending'] = True" in settle
    assert "state.pop('gfn_settle_pending', False)" in settle

    armed = source.split("def _arm_gfn_settle")[1].split("\n    def ")[0]
    assert "if not forced and state.get('optimize_interrupted')" in armed, (
        "a value the user just held is not a tidy-up to be skipped"
    )
    takeup = source.split("def _arm_gfn_takeup")[1].split("\n    def ")[0]
    assert "_arm_gfn_settle(forced=True)" in takeup


def test_the_bonding_is_read_before_a_hand_is_laid_on_the_molecule(editor):
    """Left to the first calculation that needs it, the perception is made
    from a geometry that has already been pulled about -- and past the point
    where a bond is still recognised, the topology kept for the whole session
    is one with the bond missing.  Measured: a propane whose C-C had been
    pulled to 2.1 A perceived that way came back at 3.57 A."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    began = source.split("def _begin_gfn_follow")[1].split("\n    def ")[0]
    assert "state['gfn_topology_source']" in began

    keeper = source.split("def _gfn_topology_dir")[1].split("\n    def ")[0]
    assert "gfn_topology_source" in keeper
    assert "relax_steps(seed, cycles=1" in keeper, (
        "the perception is made from it there and then"
    )


@_needs_xtb


def test_the_whole_cycle_end_to_end(editor):
    """Drag, let go, drag again, let go, press Optimise once.

    This is the mode Auto is off for: the structure stays where each drag left
    it, the bond survives being pulled well past where GFN-FF would otherwise
    stop seeing it, and one press takes it to the minimum. Placing an atom and
    having it pulled off the place you placed it is what that mode exists to
    avoid; with Auto on it is the other way round, and that is the point of
    there being a switch.
    """
    import time as _time

    refs = editor
    state = refs["editor_state"]
    propane = (
        "11\npropane\n"
        "C 1.16 0.48 -0.22\nC 0.13 -0.61 0.01\nC -1.26 -0.02 0.21\n"
        "H 2.15 0.03 -0.37\nH 0.92 1.07 -1.11\nH 1.22 1.16 0.63\n"
        "H 0.41 -1.20 0.89\nH 0.11 -1.29 -0.85\nH -1.99 -0.82 0.38\n"
        "H -1.28 0.64 1.08\nH -1.58 0.55 -0.66\n"
    )

    def carbon_carbon(xyz):
        rows = [line.split() for line in gfn.atom_lines(xyz)]
        a = [float(v) for v in rows[0][1:4]]
        b = [float(v) for v in rows[1][1:4]]
        return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5

    refs["coords_widget"].value = propane
    state["current_xyz_for_copy"] = {"content": propane}
    refs["submit_ff_dd"].value = "gfnff"
    refs["submit_auto_btn"].value = False        # the mode under test
    refs["submit_relax_btn"].value = True
    assert refs["submit_relax_btn"].value is True, "xtb was not found"

    for step, (atom, pull) in enumerate(((0, 0.7), (2, 0.6))):
        refs["submit_cmd_sync"].value = f"gfngrab:{step}:"
        rows = [line.split()
                for line in gfn.atom_lines(state["current_xyz_for_copy"]["content"])]
        rows[atom][1] = f"{float(rows[atom][1]) + pull:.6f}"
        body = "\n".join(" ".join(r) for r in rows)
        refs["submit_manip_sync"].value = f"{len(rows)}\nDELFIN drag-follow held={atom}\n{body}\n"
        _time.sleep(0.5)
        refs["submit_manip_sync"].value = f"{len(rows)}\nDELFIN drag-end\n{body}\n"
        refs["submit_cmd_sync"].value = f"gfnfree:{step}:"
        _time.sleep(1.0)

    pulled_about = state["current_xyz_for_copy"]["content"]
    assert carbon_carbon(pulled_about) > 1.9, (
        "the drag should have left it stretched, not tidied it up"
    )

    refs["submit_optimize_btn"].value = True              # exactly one press
    deadline = _time.time() + 120
    while _time.time() < deadline and refs["submit_optimize_btn"].value:
        _time.sleep(0.05)

    minimised = state["current_xyz_for_copy"]["content"]
    assert carbon_carbon(minimised) < 1.6, (
        f"one press has to reach the minimum: {carbon_carbon(minimised):.3f} A"
    )
    again = gfn.optimize_with_gfn(minimised, "gfnff")
    assert gfn.largest_shift(minimised, again["xyz"]) < 0.02, (
        "and there must be nothing left for a second press"
    )


def test_with_auto_on_letting_go_reaches_the_minimum_without_being_asked(editor):
    """The other mode, driven the same way and end to end.

    The same drag, and no press: letting go is the ask. It used to depend on
    whether an optimisation happened to be running when the atom was picked up
    -- one had something to come back to and the other did not, so the same
    gesture finished the structure or left it strained.
    """
    import time as _time

    refs = editor
    state = refs["editor_state"]
    propane = (
        "11\npropane\n"
        "C 1.16 0.48 -0.22\nC 0.13 -0.61 0.01\nC -1.26 -0.02 0.21\n"
        "H 2.15 0.03 -0.37\nH 0.92 1.07 -1.11\nH 1.22 1.16 0.63\n"
        "H 0.41 -1.20 0.89\nH 0.11 -1.29 -0.85\nH -1.99 -0.82 0.38\n"
        "H -1.28 0.64 1.08\nH -1.58 0.55 -0.66\n"
    )

    def carbon_carbon(xyz):
        rows = [line.split() for line in gfn.atom_lines(xyz)]
        a = [float(v) for v in rows[0][1:4]]
        b = [float(v) for v in rows[1][1:4]]
        return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5

    refs["coords_widget"].value = propane
    state["current_xyz_for_copy"] = {"content": propane}
    refs["submit_ff_dd"].value = "gfnff"
    refs["submit_relax_btn"].value = True
    assert refs["submit_relax_btn"].value is True, "xtb was not found"
    refs["submit_auto_btn"].value = True       # nothing was running: the mode
    assert state.get("optimize_interrupted") is None

    refs["submit_cmd_sync"].value = "gfngrab:0:"
    rows = [line.split()
            for line in gfn.atom_lines(state["current_xyz_for_copy"]["content"])]
    rows[0][1] = f"{float(rows[0][1]) + 0.7:.6f}"
    body = "\n".join(" ".join(r) for r in rows)
    refs["submit_manip_sync"].value = (
        f"{len(rows)}\nDELFIN drag-follow held=0\n{body}\n")
    _time.sleep(0.5)
    refs["submit_manip_sync"].value = f"{len(rows)}\nDELFIN drag-end\n{body}\n"
    refs["submit_cmd_sync"].value = "gfnfree:0:"

    # Letting go is the whole instruction. Nothing is pressed below this line.
    # The wait is generous because two xtb runs stand between the release and
    # the answer, and this suite is run on machines that are busy with other
    # things -- a short deadline here fails for want of a core, which reads as
    # the release having started nothing.
    deadline = _time.time() + 600
    while _time.time() < deadline and not refs["submit_optimize_btn"].value:
        _time.sleep(0.05)
    assert refs["submit_optimize_btn"].value, "letting go started nothing"
    while _time.time() < deadline and refs["submit_optimize_btn"].value:
        _time.sleep(0.05)

    reached = carbon_carbon(state["current_xyz_for_copy"]["content"])
    assert reached < 1.6, f"it did not reach a minimum on its own: {reached:.3f} A"
    refs["submit_relax_btn"].value = False


# ---------------------------------------------------------------------------
# a structure optimised in water is not one optimised in vacuum
# ---------------------------------------------------------------------------


def test_the_solvents_are_the_ones_this_xtb_is_parametrised_for():
    """Asked of the binary, not taken from a manual: every name here came back
    parametrised from xtb 6.7.1, for GFN2 and GFN-FF alike.

    Which model wraps them is a separate choice -- ALPB, GBSA or ddCOSMO --
    and it is no longer hard-coded here; see test_solvation_models.py.
    """
    assert gfn.SOLVENTS[""] == "none (gas phase)"
    for name in ("water", "dmso", "thf", "toluene", "ethanol", "ch2cl2"):
        assert name in gfn.SOLVENTS
    # aliases xtb also takes, deliberately not offered twice
    assert "h2o" not in gfn.SOLVENTS and "n-hexane" not in gfn.SOLVENTS

    source = open(gfn.__file__, encoding="utf-8").read()
    runner = source.split("def optimize_with_gfn")[1].split("\ndef ")[0]
    # The flags come from the shared table, so that MOPAC can be told about
    # the same liquids; hard-coding --alpb here is what made the model
    # unchooseable in the first place.
    assert "_solvents.xtb_flags(model, wet)" in runner
    assert "'--alpb'" not in runner


def test_a_solvent_it_does_not_know_is_refused_with_the_list():
    refused = gfn.optimize_with_gfn(_WATER, "gfn2", solvent="schnaps")

    assert refused["ok"] is False
    assert "is not a solvent that is known here" in refused["status"]
    assert "acetone" in refused["status"], "the ones it does know have to be named"


def test_which_solvent_a_result_is_about_is_part_of_the_result():
    """A geometry optimised in the gas phase and one optimised in water are two
    answers to two different questions, and one that does not say which invites
    them to be compared."""
    assert gfn.solvent_note("water") == " In water (ALPB)."
    assert gfn.solvent_note("ch2cl2") == " In dichloromethane (ALPB)."
    assert gfn.solvent_note("") == ""
    assert gfn.solvent_note(None) == ""


def test_the_solvent_reaches_every_way_of_running_it(editor):
    """Optimise, the follow and the relaxation after it are the same
    calculation asked for in three ways; a solvent that reached only one of
    them would make them disagree."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    for name in ("on_submit_optimize(change=None, every_frame=False)",
                 "_gfn_follow_step", "_gfn_settle_now"):
        body = source.split(f"def {name}")[1].split("\n    def ")[0]
        assert "submit_gfn_solvent.value" in body, name
        assert "solvent=wet" in body, name


def test_the_solvent_box_appears_with_the_method(editor):
    assert editor["submit_gfn_solvent"].layout.display == "none"

    editor["submit_ff_dd"].value = "gfn2"
    assert editor["submit_gfn_solvent"].layout.display == ""
    assert editor["submit_gfn_solvent"].value == "", "gas phase unless asked"

    editor["submit_gfn_solvent"].value = "water"
    assert "water" in editor["mol_status"].value
    assert "worth optimising again" in editor["mol_status"].value, (
        "the structure on screen was not optimised in it"
    )

    editor["submit_ff_dd"].value = "uff"
    assert editor["submit_gfn_solvent"].layout.display == "none"


@_needs_xtb


def test_a_solvent_changes_the_answer_and_the_answer_says_so(editor):
    import time as _time

    refs = editor
    state = refs["editor_state"]
    refs["coords_widget"].value = _WATER
    state["current_xyz_for_copy"] = {"content": _WATER}
    refs["submit_ff_dd"].value = "gfn2"

    energies = {}
    for solvent in ("", "water"):
        refs["submit_gfn_solvent"].value = solvent
        refs["submit_optimize_btn"].value = True
        deadline = _time.time() + 120
        while _time.time() < deadline and refs["submit_optimize_btn"].value:
            _time.sleep(0.05)
        energies[solvent] = state["gfn_energy"]

    assert energies[""] is not None and energies["water"] is not None
    assert energies["water"] < energies[""], (
        "water in water has to be stabilised by the solvent"
    )
    kcal = (energies[""] - energies["water"]) * 627.5094740631
    assert 3 < kcal < 30, f"solvation energy of {kcal:.1f} kcal/mol is not credible"
    assert "In water (ALPB)" in refs["mol_status"].value


@_needs_xtb


def test_gfnff_takes_a_solvent_too():
    """It is the method the drag and the release use, so if it could not be
    solvated the live half of the editor would be answering a different
    question from the pressed half."""
    dry = gfn.optimize_with_gfn(_WATER, "gfnff")
    wet = gfn.optimize_with_gfn(_WATER, "gfnff", solvent="water")

    assert dry["ok"] and wet["ok"]
    assert wet["energy"] != dry["energy"]
    assert wet["solvent"] == "water" and dry["solvent"] == ""


def test_the_lines_can_be_asked_to_follow_the_distances(editor):
    """3Dmol works its bond list out once, when the model is built.  Moving
    atoms does not touch it, so a bond pulled apart goes on being drawn as a
    bond -- the picture keeps the connectivity the structure had rather than
    the one it has.

    Off by default: it costs a rebuild per frame, and in a crowded
    coordination sphere the perception is at its limit and the lines flicker.
    """
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    editor_js = submit_manip_bootstrap_js()
    assert "setDynamicBonds: setDynamicBonds" in editor_js

    # every drag frame and every set of coordinates comes through here, which
    # is what makes the lines follow during a manipulation and not only after.
    # The asking and the drawing are two functions now -- a drag asks about
    # twice per mouse event and a mouse reports faster than a screen refreshes,
    # so the drawing happens once a frame however often it is asked for.
    redraw = editor_js[editor_js.index("function drawHighlightsNow("):][:600]
    assert "if (state.dynamicBonds) perceiveBonds(viewer);" in redraw

    setter = editor_js[editor_js.index("function setDynamicBonds("):][:900]
    assert "redrawHighlights(scopeKey)" in setter, (
        "switching it on has to show the truth at once, not at the next frame"
    )

    assert editor["submit_dyn_bonds_btn"].value is False
    editor["submit_dyn_bonds_btn"].value = True
    assert "follow the distances" in editor["mol_status"].value
    assert "Bond and Unbond" in editor["mol_status"].value, (
        "the picture and what the calculation holds together are two questions"
    )


def test_the_bond_perception_is_the_covalent_radii_and_says_where_it_stops():
    """3Dmol has no re-perception to ask for -- rebuildBonds does not exist in
    this build -- so the lines are worked out here, from the radii.

    Driven in a real JS engine: a C-C is drawn out to 1.92 A and gone at 1.95;
    a double bond that survives keeps its order, because the distance says
    whether there is a bond and not what kind; two atoms on top of each other
    are a mistake rather than a bond; and Pt-P at 2.30 A is drawn, which a rule
    tuned on carbon would have missed.
    """
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    editor_js = submit_manip_bootstrap_js()
    assert "var COVALENT_RADII" in editor_js
    assert "var BOND_TOLERANCE = 0.40;" in editor_js
    body = editor_js[editor_js.index("function perceiveBonds("):][:2200]
    assert "was[i + '-' + j] || 1" in body, "a double bond has to stay double"
    assert "d2 < 0.16" in body, "overlapping atoms are not a bond"
    for element in ("Pt", "Fe", "Re", "Pd", "Ru"):
        assert f"{element}:" in editor_js, f"no covalent radius for {element}"


def test_the_reach_is_added_not_multiplied_so_metals_keep_one_sphere():
    """A factor grows the reach with the radii, so a metal gets a
    proportionally huge one and starts drawing lines to the second
    coordination sphere.

    Measured on a manganese, radii 1.50 and 0.76: a factor of 1.30 reaches
    2.94 A and picks up a carbon at 2.90, which is second sphere.  Adding 0.40
    reaches 2.66 and does not, while every first-sphere contact still is --
    Mn-N 2.25, Mn-O 2.15, Mn-Cl 2.40, Pt-P 2.30.

    Checked across 35 metals from scandium to uranium against N, O, P, S, Cl
    and C donors: no first-sphere bond falls outside the reach, and the
    tightest margin to the atom behind the donor is 0.56 A (Sc-O).  It holds
    for any element, because the second sphere is one covalent bond further
    off -- an N-H is 1.02 A, the shortest thing that can be behind a donor --
    and the tolerance is 0.40.
    """
    from delfin.dashboard.molecule_viewer import submit_manip_bootstrap_js

    editor_js = submit_manip_bootstrap_js()
    body = editor_js[editor_js.index("function perceiveBonds("):][:2400]
    assert "var reach = radii[i] + radii[j] + BOND_TOLERANCE;" in body, (
        "multiplying is what reached into the second sphere"
    )
    assert "* BOND_TOLERANCE" not in body

    radii = editor_js[editor_js.index("var COVALENT_RADII"):]
    radii = radii[:radii.index("};")]
    import re

    table = {m.group(1): float(m.group(2))
             for m in re.finditer(r"([A-Z][a-z]?): ([0-9.]+)", radii)}
    tolerance = 0.40
    donors = {"N": 1.02, "O": 0.96, "P": 1.84, "S": 1.34, "C": 1.09}
    metals = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
              "Zr", "Mo", "Ru", "Rh", "Pd", "Ag", "W", "Re", "Os", "Ir",
              "Pt", "Au", "Hg", "La", "Ce", "Nd", "Eu", "Gd", "Yb", "U"]
    for metal in metals:
        assert metal in table, f"no covalent radius for {metal}"
        for donor, behind in donors.items():
            reach = table[metal] + table[donor] + tolerance
            bond = table[metal] + table[donor]
            assert reach >= bond, f"{metal}-{donor} bond would not be drawn"
            assert reach < bond + behind, (
                f"{metal}-{donor} reaches past the atom behind the donor"
            )


# ---------------------------------------------------------------------------
# g-xTB, which is not the xtb beside it
# ---------------------------------------------------------------------------
def test_gxtb_is_looked_for_under_its_own_name():
    """It ships as a statically linked xtb 6.7.1 of its own -- the method is
    not in a released tblite, and tblite 0.4.0 offers gfn1, gfn2 and ipea1.

    An ordinary xtb **accepts** --gxtb and silently runs GFN2: measured here,
    the energy came back at -5.070383095219 either way, and the Hamiltonian
    line said GFN2-xTB.  So the binary is never assumed to be the one beside
    it.
    """
    assert gfn.GFN_METHODS["gxtb"]["flags"] == ["--gxtb"]
    assert gfn.GFN_METHODS["gxtb"]["binary"] == "gxtb"
    assert gfn.GFN_METHODS["gxtb"]["reports"] is None, (
        "g-xTB prints no Hamiltonian line of its own; that is the check"
    )

    source = open(gfn.__file__, encoding="utf-8").read()
    finder = source.split("def find_gxtb")[1].split("\ndef ")[0]
    assert "DELFIN_GXTB_BINARY" in finder
    assert "xtb-gxtb" in finder


def test_an_xtb_that_ignored_the_flag_is_caught(monkeypatch):
    """The one way this can go wrong silently, so the one thing checked."""
    monkeypatch.setattr(gfn, "find_gxtb", gfn.find_xtb)   # the wrong binary

    if gfn.find_xtb() is None:
        pytest.skip("no xtb to be wrong with")
    refused = gfn.optimize_with_gfn(_WATER, "gxtb")

    assert refused["ok"] is False
    assert "ignored the flag" in refused["status"]
    assert "GFN2-xTB" in refused["status"], "it has to name what ran instead"


def test_gxtb_is_offered_no_solvent_because_it_takes_none():
    """--alpb, --gbsa and --cpcmx all stop this build dead, and --cosmo runs
    but only writes a file: it moved water in water *up* by 1.6 kcal/mol,
    which is the wrong direction for a solvation energy."""
    assert gfn.GFN_METHODS["gxtb"].get("solvation") is False

    refused = gfn.optimize_with_gfn(_WATER, "gxtb", solvent="water")
    assert refused["ok"] is False
    assert "no implicit solvation" in refused["status"]
    assert "GFN2-xTB or GFN-FF" in refused["status"], (
        "and it has to say which method does have one"
    )


def test_the_installer_can_be_asked_for_gxtb():
    assert gfn.install_command("gxtb")[-1] == "gxtb"
    assert gfn.install_command()[-1] == "xtb", "xtb stays the default"

    script = open(gfn.install_script(), encoding="utf-8").read()
    assert "install_gxtb()" in script
    # The dispatch moved into install_one, where each tool is attempted on its
    # own so that one failing does not take the rest of the list with it.
    assert "gxtb|g-xtb)     install_gxtb ;;" in script
    assert "sha256sum" in script, (
        "a binary fetched from the network and run on a user's structures is "
        "worth the one line it costs to check it is the one that was published"
    )
    assert "link_into_bin \"${target}/bin/xtb\" xtb-gxtb" in script


def test_the_offer_names_the_program_the_method_needs(editor, monkeypatch):
    """g-xTB is not the xtb beside it, so offering the wrong one would install
    something that changes nothing."""
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: None)

    editor["submit_ff_dd"].value = "gfnff"
    assert editor["submit_xtb_install_btn"].description == "Install xtb"

    editor["submit_ff_dd"].value = "gxtb"
    assert editor["submit_xtb_install_btn"].description == "Install g-xTB"
    assert editor["editor_state"]["xtb_install_tool"] == "gxtb"


def test_the_solvent_box_is_gone_under_a_method_that_has_no_solvation(editor):
    editor["submit_ff_dd"].value = "gfn2"
    editor["submit_gfn_solvent"].value = "water"
    assert editor["submit_gfn_solvent"].layout.display == ""

    editor["submit_ff_dd"].value = "gxtb"
    assert editor["submit_gfn_solvent"].layout.display == "none"
    assert editor["submit_gfn_solvent"].value == "", (
        "a control that can only produce a refusal is worse than no control"
    )


_needs_gxtb = pytest.mark.skipif(
    __import__("delfin.dashboard.gfn_optimize", fromlist=["x"]).find_gxtb() is None,
    reason="the g-xTB build is not installed",
)


@_needs_gxtb
def test_gxtb_optimises_and_gets_the_geometry_right():
    """It approximates wB97M-V/def2-TZVPPD, and on a propane it does.

    Experiment: C-C 1.526 A, C-C-C 112.4 deg, C-H 1.09 A.  Measured here --
    GFN-FF 1.5213 / 112.51, GFN2 1.5238 / 111.67, g-xTB 1.5162 / 112.49 -- so
    it is at least as good as the others and takes 0.21 s where GFN2 takes
    0.03; a full optimisation is 0.29 s at 32 atoms and 2.19 s at 92.
    """
    import math

    propane = (
        "11\npropane\n"
        "C 1.16 0.48 -0.22\nC 0.13 -0.61 0.01\nC -1.26 -0.02 0.21\n"
        "H 2.15 0.03 -0.37\nH 0.92 1.07 -1.11\nH 1.22 1.16 0.63\n"
        "H 0.41 -1.20 0.89\nH 0.11 -1.29 -0.85\nH -1.99 -0.82 0.38\n"
        "H -1.28 0.64 1.08\nH -1.58 0.55 -0.66\n"
    )
    result = gfn.optimize_with_gfn(propane, "gxtb")

    assert result["ok"] is True, result["status"]
    assert result["converged"] is True
    assert result["frames"], "the path has to come back, as for any other method"
    points = [[float(v) for v in line.split()[1:4]]
              for line in gfn.atom_lines(result["xyz"])]
    assert abs(math.dist(points[0], points[1]) - 1.526) < 0.05, "C-C"
    first = [points[0][n] - points[1][n] for n in range(3)]
    third = [points[2][n] - points[1][n] for n in range(3)]
    dot = sum(a * b for a, b in zip(first, third))
    norms = math.dist(points[0], points[1]) * math.dist(points[2], points[1])
    angle = math.degrees(math.acos(max(-1.0, min(1.0, dot / norms))))
    assert abs(angle - 112.4) < 3.0, f"C-C-C came out at {angle:.2f} deg"

    # and its energy is a different quantity from a GFN one, not a variant
    assert result["energy"] < -50, "g-xTB totals include the core electrons"


@_needs_gxtb
def test_gxtb_holds_what_the_editor_holds():
    """Charge, spin and held values reach it the same way they reach GFN2 --
    it is the same command line, which is the whole point of it being an xtb.
    """
    held = [{"kind": "distance", "atoms": [0, 1], "value": 1.05, "mode": "fix"}]
    result = gfn.optimize_with_gfn(_WATER, "gxtb", constraints=held)

    assert result["ok"] is True, result["status"]
    rows = [line.split() for line in gfn.atom_lines(result["xyz"])]
    first = [float(v) for v in rows[0][1:4]]
    second = [float(v) for v in rows[1][1:4]]
    length = sum((a - b) ** 2 for a, b in zip(first, second)) ** 0.5
    assert abs(length - 1.05) < 0.01, f"the held O-H came out at {length:.3f} A"


@_needs_xtb
def test_xtb_uses_its_own_parameters_and_not_whatever_is_lying_around(tmp_path,
                                                                     monkeypatch):
    """A stray parameter file in the home directory killed every GFN2 run.

    xtb reads ``param_gfn2-xtb.txt`` from XTBPATH or from the home directory
    *instead of* the parameters compiled into it. A truncated one there gives

        no basis found for atom 1 Z= 8
        ERROR STOP / Error termination. Backtrace:

    on every GFN2 run, while GFN-FF -- which does not read them -- goes on
    working. That pair is the fingerprint, and it is the report that came in
    from a user: an optimisation that failed with a backtrace and nothing to
    act on. Pointing XTBPATH at the share directory beside the binary restores
    the right answer with the bad file still in place.
    """
    from delfin.dashboard import gfn_optimize as gfn

    real = gfn.parameter_home(gfn.find_xtb())
    if not real:
        pytest.skip('this xtb keeps its parameters compiled in')

    home = tmp_path / 'home'
    home.mkdir()
    source = pathlib.Path(real) / 'param_gfn2-xtb.txt'
    (home / 'param_gfn2-xtb.txt').write_bytes(source.read_bytes()[:400])
    monkeypatch.setenv('HOME', str(home))
    monkeypatch.delenv('XTBPATH', raising=False)

    water = '3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n'

    # As it was: xtb reads the broken file and stops.
    monkeypatch.setattr(gfn, 'parameter_home', lambda _binary: None)
    broken = gfn.optimize_with_gfn(water, method='gfn2', max_steps=1)
    assert not broken['ok'], 'the poisoned parameter file did not bite'
    # And the message names the cause rather than the announcement: the
    # complaint is printed to stdout and lands *after* the terminators in the
    # merged output, which is why reading upwards from them found the
    # citation list.
    assert 'no basis found' in str(broken['status']).lower(), broken['status']

    # And GFN-FF is untouched, which is what makes the pair a fingerprint.
    monkeypatch.undo()
    monkeypatch.setenv('HOME', str(home))
    monkeypatch.delenv('XTBPATH', raising=False)
    monkeypatch.setattr(gfn, 'parameter_home', lambda _binary: None)
    assert gfn.optimize_with_gfn(water, method='gfnff', max_steps=1)['ok']

    # As it is now: the parameters that belong to the binary are used.
    monkeypatch.undo()
    monkeypatch.setenv('HOME', str(home))
    monkeypatch.delenv('XTBPATH', raising=False)
    fixed = gfn.optimize_with_gfn(water, method='gfn2', max_steps=1)
    assert fixed['ok'], fixed['status']
    assert fixed['energy'] == pytest.approx(-5.07, abs=0.05)


@_needs_xtb
def test_the_parameter_directory_is_the_one_beside_the_binary():
    from delfin.dashboard import gfn_optimize as gfn

    found = gfn.parameter_home(gfn.find_xtb())
    if found:
        assert (pathlib.Path(found) / 'param_gfn2-xtb.txt').is_file()
    # Nothing sensible to point at is not an error: g-xTB has them built in.
    assert gfn.parameter_home('/nowhere/at/all/xtb') is None
    assert gfn.parameter_home(None) is None


def test_an_interactive_run_may_use_more_than_one_core():
    """One core was the rule, on the grounds that a login node is shared.

    It is -- but a press of Optimise is one person waiting, not a job. Measured
    through the dashboard on a 74-atom cholesterol, five optimisation steps:
    719 ms on one core, 569 on four, and flat beyond that. GFN-FF gains
    nothing, 69 ms either way, because at that size it is mostly the program
    starting up -- 12 ms of it.
    """
    assert gfn.INTERACTIVE_CORES == 4
    assert gfn.interactive_cores() >= 1
    assert gfn.interactive_cores() <= (os.cpu_count() or 1)


def test_the_core_count_can_be_set_and_never_exceeds_the_machine(monkeypatch):
    monkeypatch.setenv('DELFIN_GFN_CORES', '2')
    assert gfn.interactive_cores() == 2
    monkeypatch.setenv('DELFIN_GFN_CORES', '100000')
    assert gfn.interactive_cores() == (os.cpu_count() or 1)
    monkeypatch.setenv('DELFIN_GFN_CORES', 'nonsense')
    assert gfn.interactive_cores() == min(gfn.INTERACTIVE_CORES, os.cpu_count() or 1)


def test_the_trajectory_is_shown_as_fast_as_it_is_computed():
    """The limit was the watching, not the calculating.

    The log was read five times a second, so a GFN-FF optimisation put five
    frames on screen per second while it was computing a step every fourteen
    milliseconds. Asking whether the file grew costs 0.002 ms and reading it
    1.73 ms for a 211 KiB log of 37 frames -- so the question is asked every
    pass and only the reading is limited.

    Measured end to end: at 74 atoms one delivery became three and the run
    itself finished in 0.14 s instead of 0.30; at 405 atoms, 27 deliveries
    became 64 of the 66 frames that exist, which is every frame as it is
    computed.
    """
    assert gfn.WATCH_INTERVAL <= 0.01
    assert gfn.FRAME_READ_INTERVAL <= 0.05

    import inspect

    body = inspect.getsource(gfn.optimize_with_gfn)
    # The size is checked outside the rate limit, the reading inside it.
    watch = body.split('xtb writes the log as it optimises')[1]
    assert watch.index('log.stat().st_size') < watch.index('FRAME_READ_INTERVAL')
    assert 'time.sleep(WATCH_INTERVAL)' in body
    assert '0.2' not in watch.split('read_trajectory')[0]


def test_an_empty_comment_line_does_not_cost_an_atom():
    """XYZ's second line is a comment, and it is allowed to be empty.

    The coordinate lines were picked out by dropping every blank line and then
    skipping two -- so with an empty comment the first atom was skipped as
    though it were the comment. A water went to xtb as two hydrogens and came
    back as two hydrogens, at -0.982686 Eh instead of -5.070544, and every
    other sign said the run had succeeded.

    A named block in the ORCA Builder writes exactly that: name, count, empty
    comment, atoms. So the editor there optimised a molecule short of its first
    atom, silently, and the same would happen to anyone pasting such a block
    into the Submit tab.
    """
    from delfin.dashboard.gfn_optimize import atom_lines

    shapes = {
        'blank comment': '3\n\nO 0 0 0\nH 0.9 0 0\nH -0.3 0.85 0\n',
        'named comment': '3\nwater\nO 0 0 0\nH 0.9 0 0\nH -0.3 0.85 0\n',
        'no header at all': 'O 0 0 0\nH 0.9 0 0\nH -0.3 0.85 0\n',
        'blank line first': '\n3\nwater\nO 0 0 0\nH 0.9 0 0\nH -0.3 0.85 0\n',
        'blanks in the middle': '3\nwater\nO 0 0 0\n\nH 0.9 0 0\nH -0.3 0.85 0\n',
        'trailing blanks': '3\nwater\nO 0 0 0\nH 0.9 0 0\nH -0.3 0.85 0\n\n\n',
        # The count is not trusted either: xtb reads that many and ignores the
        # rest, so a header that disagrees with the block must not decide.
        'a count that lies': '89\nwater\nO 0 0 0\nH 0.9 0 0\nH -0.3 0.85 0\n',
    }
    for shape, xyz in shapes.items():
        assert len(atom_lines(xyz)) == 3, shape
        assert atom_lines(xyz)[0].split()[0] == 'O', shape


# ---------------------------------------------------------------------------
# one layout, and a header that is true about the coordinates under it
# ---------------------------------------------------------------------------


def test_every_coordinate_the_editor_writes_is_in_one_layout():
    """Which layout a box was in used to say where it had been written from.

    Fourteen decimals meant xtb, eight the frame a stopped run was showing,
    and six the browser -- whose model is serialised in JavaScript.  Two
    geometries that differed only in that were two different histories, and
    counting decimals is not how a user should have to find that out.
    """
    from delfin.dashboard import structure_editor as se

    line = se.xyz_line("C", 4.64514437834514, 1.05045715164728,
                       -0.41772688977231)
    assert line == (
        "C            4.64514437834514        1.05045715164728"
        "       -0.41772688977231"
    )
    assert len(line) == 77
    # Two letters and a wider number keep the same columns.
    assert len(se.xyz_line("Mn", -12.5, 0.0, 123.25)) == 77
    # The count comes from the lines, not from a header that may disagree.
    doc = se.xyz_document(["O 0 0 0", "H 0.96 0 0"], "a note")
    assert doc.splitlines()[0] == "2"
    assert doc.splitlines()[1] == "a note"
    assert doc.splitlines()[2].startswith("O    ")


def test_a_line_it_cannot_read_is_left_alone_rather_than_dropped():
    """The box is the user's: a column someone is relying on is not this
    function's to throw away."""
    from delfin.dashboard import structure_editor as se

    kept = se.xyz_body(["O 0 0 0", "this is not an atom", "", "H 1 x 0"])
    assert kept[0].startswith("O    ")
    assert kept[1] == "this is not an atom"
    assert kept[2] == "H 1 x 0"
    assert len(kept) == 3, "blank lines go, unreadable ones stay"


def test_the_browser_serialises_in_the_same_layout(player_js):
    """Six decimals from the page and fourteen from xtb is the same box
    changing shape according to which side last wrote it."""
    from delfin.dashboard import tab_submit

    js = tab_submit.submit_manip_bootstrap_js()
    assert "toFixed(14)" in js
    assert "toFixed(6)" not in js.split("function serializeXyz")[1][:800]
    assert "function xyzColumn" in js


def test_a_geometry_from_the_browser_does_not_wear_the_optimised_label():
    """These coordinates are what the hand has just done to the structure.

    Carrying the old comment over put "Optimised in DELFIN viewer" above a
    geometry that had been dragged out of shape since -- so the word said
    where the box was last written from rather than what was in it, and a
    benzene with a hydrogen 2.66 A off its carbon read as the result of an
    optimisation.
    """
    from delfin.dashboard import structure_editor as se

    assert se._is_editor_comment("Optimised in DELFIN viewer")
    assert se._is_editor_comment("Settled with GFN2-xTB")
    assert se._is_editor_comment("stopped at the frame on screen")
    assert se._is_editor_comment("DELFIN drag-end held=3,4")
    # A name off a file or a note the user typed is theirs, and stays.
    assert not se._is_editor_comment("cholesterol, from the CSD")
    assert not se._is_editor_comment("Converted from SMILES (isomer: E)")

    sync = EDITOR_SOURCE.split("def on_submit_manip_sync")[1].split("\n    def ")[0]
    assert "_is_editor_comment(kept)" in sync
    assert "'Edited in DELFIN viewer'" in sync


def test_a_run_that_produced_no_geometry_does_not_relabel_the_box():
    """It hands the input straight back, and writing that under "Optimised"
    labelled the geometry the user made as the answer to a question that was
    never answered."""
    apply_body = SUBMIT_SOURCE.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0].split("def _apply()")[1]
    guard = apply_body.split("elif failures:")[1].split("else:")[0]
    assert "coords_widget.value" not in guard
    assert "pass" in guard
    # And it does not say it optimised something either.
    assert "f'{label} could not optimise it.' if not done else" in apply_body


# ---------------------------------------------------------------------------
# what the chosen method can do, and what becomes of the held values
# ---------------------------------------------------------------------------


def test_the_toolbar_shows_what_the_method_has_and_nothing_else(editor, monkeypatch):
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    dd = editor["submit_ff_dd"]

    dd.value = "uff"
    assert editor["submit_gfn_charge"].layout.display == "none"
    assert editor["submit_gfn_mult"].layout.display == "none"
    assert editor["submit_gfn_autospin"].layout.display == "none"
    assert editor["submit_strength_slider"].layout.display == ""
    assert editor["submit_settle_btn"].layout.display in (None, "")
    # Auto is the other way round: going down to a minimum on release is
    # refused outright for a browser method, so the switch has nothing to do.
    # Its value is left as it stands -- nothing reads it here, and switching
    # it off would cost the user their setting for stepping through UFF.
    assert editor["submit_auto_btn"].layout.display == "none"

    dd.value = "gfn2"
    assert editor["submit_gfn_charge"].layout.display == ""
    assert editor["submit_gfn_mult"].layout.display == ""
    assert editor["submit_gfn_autospin"].layout.display == ""
    assert editor["submit_strength_slider"].layout.display == "none"
    assert editor["submit_settle_btn"].layout.display == "none"
    assert editor["submit_auto_btn"].layout.display == ""

    dd.value = "pm7"
    assert editor["submit_gfn_charge"].layout.display == ""
    # Scanning the multiplicity is xtb's; MOPAC is given the one on screen.
    assert editor["submit_gfn_autospin"].layout.display == "none"
    assert editor["submit_settle_btn"].layout.display == "none"
    assert editor["submit_strength_slider"].layout.display == "none"
    assert editor["submit_auto_btn"].layout.display == ""


def test_a_multiplicity_of_zero_or_below_is_refused_by_the_box(editor):
    """M = 2S+1, so the smallest there is is 1.

    The box took 0 and -3 as readily as 2, and neither reached xtb as itself:
    the conversion to unpaired electrons floors at zero, so both were quietly
    run as a singlet -- the confident wrong answer arrived at from a number
    the user could see was not what they typed.
    """
    box = editor["submit_gfn_mult"]
    assert box.min == 1
    box.value = 0
    assert box.value == 1
    box.value = -3
    assert box.value == 1
    box.value = 3
    assert box.value == 3


def test_held_values_survive_a_change_of_method_and_are_spoken_for(editor, monkeypatch):
    """The list describes the molecule, not the program.

    What each engine will do with it differs completely, and an engine that is
    handed none of them has to say so: a held bond that stays in the list and
    stops being held is the list describing a thing that is not happening.
    """
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    state = editor["editor_state"]
    held = [{"kind": "distance", "atoms": [0, 1], "value": 1.6, "mode": "fix"}]
    state["constraints"] = list(held)

    # What each of them will do with the one held value, said at the moment
    # the method changes rather than left for the run to reveal.
    expected = {
        "gfn2": "will hold all 1",
        "gfnff": "will hold all 1",
        "pm7": "holds them its own way",
        "uff": "pulls towards a held value rather than fixing it",
    }
    for method, phrase in expected.items():
        editor["submit_ff_dd"].value = method
        assert state["constraints"] == held, (
            f"{method} must not eat the held values"
        )
        assert phrase in editor["mol_status"].value, method


def test_mopac_holds_a_value_by_fixing_the_atoms_that_name_it(editor):
    """MOPAC's Cartesian input carries an optimisation flag per coordinate,
    and that is the only constraint it has.  Measured on a propane with PM7:
    free, C1 and C3 relax from 3.000 A to 2.523 and each moves 0.302 A; with
    their flags at zero each moves 0.0000 A and the distance stays 3.000."""
    state = editor["editor_state"]
    state["constraints"] = [
        {"kind": "distance", "atoms": [0, 1], "value": 1.6, "mode": "fix"}
    ]
    editor["submit_ff_dd"].value = "pm7"
    said = editor["mol_status"].value
    assert "holds them its own way" in said
    # And says what that is not: it fixes the atoms, not the value between
    # them, so those atoms also stop turning and moving.
    assert "fixes atoms, not the value between them" in said


def test_a_pull_is_refused_by_mopac_rather_than_frozen(editor):
    """One flag, on or off, so there is nothing to negotiate with.  Freezing a
    pull would hold as exact a value the user asked to have argued with."""
    from delfin.dashboard import mopac_optimize as mopac

    pull = [{"kind": "distance", "atoms": [0, 1], "value": 1.6, "mode": "pull"}]
    fix = [{"kind": "distance", "atoms": [0, 1], "value": 1.6, "mode": "fix"}]

    assert mopac.freeze_flags(pull)["frozen"] == set()
    assert mopac.freeze_flags(pull)["pulls"] == 1
    assert mopac.freeze_flags(pull)["held"] == 0
    assert "not honoured" in mopac.freeze_note(mopac.freeze_flags(pull))

    assert mopac.freeze_flags(fix)["frozen"] == {0, 1}
    assert mopac.freeze_flags(fix)["held"] == 1

    # An atom this structure does not have is dropped, not passed on.
    assert mopac.freeze_flags(fix, atoms=1)["dropped"]
    assert mopac.freeze_flags(fix, atoms=1)["frozen"] == set()

    state = editor["editor_state"]
    state["constraints"] = pull
    editor["submit_ff_dd"].value = "pm7"
    assert "not honoured" in editor["mol_status"].value


def test_with_every_switch_off_letting_go_of_an_atom_runs_nothing(editor, monkeypatch):
    """GFN2 chosen, Dynamik Opt off, Optimise off -- and it optimised anyway.

    The gate on the release path read Settle's widget, and Settle defaulted to
    on.  Under a server method that switch is not on the toolbar at all, so it
    could be left on from an earlier method or restored with a structure and
    never be seen again -- and what it then ran was not a tidy-up but the whole
    minimisation, this path having stopped capping its cycles when it was made
    to match Optimise.  A structure went down to a minimum with nothing on
    screen claiming to be running.

    Driven the way the page drives it: a release arrives as a command, and
    nothing may be armed by it.
    """
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    state = editor["editor_state"]
    editor["submit_ff_dd"].value = "gfn2"
    editor["submit_relax_btn"].value = False
    editor["submit_auto_btn"].value = False
    # The toolbar has taken Settle away and switched it off with it.
    assert editor["submit_settle_btn"].value is False
    assert editor["submit_settle_btn"].layout.display == "none"

    # Forced back on behind the toolbar's back -- the state the bug was
    # reached from, and the gate must not care what it says under GFN.
    editor["submit_settle_btn"].value = True
    state["gfn_settle_forced"] = False
    state.pop("gfn_settle_armed", None)
    state.pop("gfn_minimise_armed", None)

    # The page sends "verb:serial:payload"; the serial only makes the same
    # command twice read as two changes.
    editor["submit_cmd_sync"].value = "gfnfree:1:"

    assert not state.get("gfn_settle_armed"), (
        "a release with everything off must not arm a relaxation"
    )
    assert not state.get("gfn_minimise_armed"), (
        "and must not arm a minimisation either -- that is Auto's"
    )
    assert state.get("optimize_run") is None

    # And the gate says it in one line: on a server method only a hand -- a
    # value set, a value held -- arms this at all.  It used to read Settle's
    # widget, which under GFN is not on the toolbar to be turned off.
    arm = EDITOR_SOURCE.split("def _arm_gfn_settle")[1].split("\n    def ")[0]
    assert "if not state.get('gfn_settle_forced'):" in arm
    assert "submit_settle_btn.value" not in arm


def test_with_auto_and_dynamik_opt_on_a_release_does_go_to_a_minimum(editor, monkeypatch):
    """The other half: switched on, it is asked for, and it happens."""
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    state = editor["editor_state"]
    editor["submit_ff_dd"].value = "gfn2"
    editor["submit_relax_btn"].value = True
    editor["submit_auto_btn"].value = True
    state.pop("gfn_minimise_armed", None)

    editor["submit_cmd_sync"].value = "gfnfree:2:"

    assert state.get("gfn_minimise_armed") is True, (
        "with both switches on, letting go is what starts the minimisation"
    )


def test_a_restored_structure_cannot_hand_back_a_switch_the_method_lacks(editor, monkeypatch):
    """A memory made under one method is applied under another.

    The restore writes the switches one by one, and the method box may already
    hold the value it is being given -- so its handler never fires, and a
    switch that method does not have would sit there on.
    """
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    editor["submit_ff_dd"].value = "gfn2"
    editor["submit_settle_btn"].value = True      # as a restore would leave it

    apply_body = EDITOR_SOURCE.split("def _apply_controls")[1].split("\n    def ")[0]
    assert "_refresh_method_controls()" in apply_body


def test_retuning_a_held_value_reaches_the_engine_without_pressing_hold():
    """Three ways to alter a held value, and one of them went nowhere.

    Hold arms the take-up, and so does changing pull to fix; typing a new
    number into the box did not.  Under the browser's field that never showed,
    because its restraints ride along with the parameters that are handed over
    again -- but under GFN nothing runs between drags, so the change is what
    has to start it.  The list said one thing, the structure went on standing
    at another, and pressing Hold again was what made it happen.
    """
    for name in ("on_submit_hold", "on_submit_hold_mode",
                 "on_submit_constraint_retune"):
        body = EDITOR_SOURCE.split(f"def {name}")[1].split("\n    def ")[0]
        assert "_arm_gfn_takeup(" in body, name
        assert "_enable_live_forcefield()" in body, name


def test_xtb_holds_a_value_with_one_force_constant_for_the_whole_set():
    """pull and fix are two force constants, not two mechanisms.

    xtb takes one ``force constant=`` for the whole ``$constrain`` block -- a
    second is read and ignored -- so a set with anything exact in it is held
    exact throughout, and the caller is told so rather than left to wonder why
    a pull stopped negotiating.
    """
    pull = [{"kind": "distance", "atoms": [0, 1], "value": 1.6, "mode": "pull"}]
    fix = [{"kind": "distance", "atoms": [0, 1], "value": 1.6, "mode": "fix"}]

    assert gfn.constraint_input(pull)["force"] == gfn.PULL_FORCE_CONSTANT
    assert gfn.constraint_input(fix)["force"] == gfn.FIX_FORCE_CONSTANT
    assert not gfn.constraint_input(pull)["mixed"]
    assert not gfn.constraint_input(fix)["mixed"]

    both = gfn.constraint_input(pull + [
        {"kind": "angle", "atoms": [0, 1, 2], "value": 109.5, "mode": "fix"}])
    assert both["held"] == 2
    assert both["force"] == gfn.FIX_FORCE_CONSTANT, (
        "one constant for the block, and anything exact in it decides"
    )
    assert both["mixed"] is True
    assert "held as firmly as the exact values" in gfn.held_note(both)

    # Only internal coordinates.  xtb's own $fix atoms: is not asked for --
    # three fixed carbons of a propane came back at 0.4623 A under GFN2.
    assert set(gfn.CONSTRAINT_ATOMS) == {"distance", "angle", "dihedral"}
    assert "$fix" not in gfn.constraint_input(fix)["text"]
    assert gfn.constraint_input(
        [{"kind": "atom", "atoms": [0], "value": 0.0, "mode": "fix"}]
    )["held"] == 0


def test_the_optimisation_owns_the_coordinate_box_while_it_runs(editor, monkeypatch):
    """The browser's own field reports where it has got to twice a second and
    once more when it is switched off.

    Those arrived with no reason on them, took neither the edit path nor any
    other, and wrote the coordinate box regardless -- so the picture's idea of
    the geometry landed on top of the calculation's, and which of the two the
    user was left with came down to which wrote last.  That is the shape the
    "Optimised in DELFIN viewer" header was found over a structure no run had
    produced.
    """
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    state = editor["editor_state"]
    editor["submit_ff_dd"].value = "gfn2"
    before = editor["coords_widget"].value

    token = object()
    state["optimize_run"] = token
    editor["submit_manip_sync"].value = (
        "3\nEdited in DELFIN viewer\n"
        "O 9.000000 0.000000 0.000000\n"
        "H 0.960000 0.000000 0.000000\n"
        "H -0.240000 0.930000 0.000000\n"
    )
    assert editor["coords_widget"].value == before, (
        "a report from the field must not land on a running optimisation"
    )
    assert state["optimize_run"] is token, "and must not stop it either"

    # A hand still takes it back: an edit interrupts the run and lands.
    editor["submit_manip_sync"].value = (
        "3\nDELFIN drag-end\n"
        "O 8.000000 0.000000 0.000000\n"
        "H 0.960000 0.000000 0.000000\n"
        "H -0.240000 0.930000 0.000000\n"
    )
    assert "8.00000000000000" in editor["coords_widget"].value
    assert state["optimize_run"] is None, "the run is about a structure that is gone"


def test_undoing_a_drag_in_the_browser_counts_as_an_edit():
    """It went back with no reason on it, so it was neither an edit nor
    anything else: the optimisation ran on over a geometry that had just been
    taken back."""
    from delfin.dashboard import tab_submit

    js = tab_submit.submit_manip_bootstrap_js()
    assert "pushXyzToPython(scopeKey, 'undo')" in js
    # And the field's own reports say which they are, so they can be told apart.
    assert "pushXyzToPython(scopeKey, 'field')" in js
    assert "pushXyzToPython(scopeKey);" not in js, "every push carries a reason"

    sync = EDITOR_SOURCE.split("def on_submit_manip_sync")[1].split("\n    def ")[0]
    assert "note.startswith('DELFIN undo')" in sync
    assert "_interrupt_gfn()" in sync


def test_the_follow_does_not_move_a_structure_an_optimisation_owns(editor, monkeypatch):
    """Two xtb processes arranging the same molecule around each other.

    Picking an atom up interrupts a running optimisation before the follow
    begins, so the ordinary way round never reaches this.  Pressing Optimise
    while an atom is already held does: the follow went on asking for a run per
    push while the press was running.  It is the same collision the coordinate
    box was given an owner for, one step earlier in the same path, and the
    press is the later of the two things the user asked for.

    Driven on gfn_follow_busy, which _gfn_follow_step sets before it starts its
    thread -- the run itself happens off this one, so waiting for it to appear
    is a race and its absence would prove nothing.
    """
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    state = editor["editor_state"]
    editor["submit_ff_dd"].value = "gfn2"
    editor["submit_relax_btn"].value = True

    def push():
        state["gfn_follow"] = True
        state["gfn_follow_method"] = "gfn2"
        state.pop("gfn_follow_busy", None)
        editor["submit_manip_sync"].value = (
            f"3\nDELFIN drag-follow held=0 n={state.get('optimize_run') is None}\n"
            "O 0.100000 0.000000 0.000000\n"
            "H 0.960000 0.000000 0.000000\n"
            "H -0.240000 0.930000 0.000000\n"
        )
        return bool(state.get("gfn_follow_busy"))

    state["optimize_run"] = object()          # a press is running
    assert push() is False, 'the follow must not run against a press'

    # The other way round, so the check above is worth something: with nothing
    # optimising, the same push does start a follow.
    state["optimize_run"] = None
    assert push() is True, 'and it still follows when nothing owns the structure'

    follow = EDITOR_SOURCE.split("def _gfn_follow_step")[1].split("\n    def ")[0]
    assert "if state.get('optimize_run') is not None:" in follow


def test_one_quiet_round_is_not_a_finished_structure(bare_editor):
    """"No longer moving" used to mean "did not move much just then".

    A run is one xtb --opt and its optimiser has a cycle limit of its own; one
    that reaches that limit in a flat stretch hands back a geometry that moved
    almost nothing.  Its successor is a new process with a fresh cycle budget
    and a fresh optimiser history, which is the whole reason these rounds
    exist -- so stopping on the first quiet one ended the press exactly where
    another round would still have got somewhere, and the user pressed again
    and watched it move.
    """
    part, _state = bare_editor
    part.submit_optimize_btn.unobserve_all()
    part.submit_optimize_btn.value = True
    quiet = dict(converged=False, moved=0.0, rounds=1, failed=False,
                 every_frame=False)

    assert part._optimise_carries_on(**{**quiet, 'still': 1}) is True
    assert part._optimise_carries_on(**{**quiet, 'still': 2}) is False
    assert part._OPTIMISE_STILL_ROUNDS == 2

    # It never outranks the endings that are certain.
    assert part._optimise_carries_on(
        **{**quiet, 'still': 1, 'converged': True}) is False
    assert part._optimise_carries_on(
        **{**quiet, 'still': 1, 'failed': True}) is False
    assert part._optimise_carries_on(
        **{**quiet, 'still': 1, 'rounds': 99}) is False

    # And a round that moved is not quiet at all, whatever came before it.
    assert part._optimise_carries_on(
        **{**quiet, 'moved': 0.25, 'still': 0}) is True


def test_the_quiet_count_starts_over_with_each_press(bare_editor):
    """A press that ended on two quiet rounds must not hand its count to the
    next one, or the next press stops after a single round."""
    apply_body = SUBMIT_SOURCE.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0]
    # Reset beside the round counter, on a press that is not carrying on.
    assert "state['optimize_still'] = 0" in apply_body
    assert "state['optimize_still'] = still if carry_on else 0" in apply_body
    # Counted from the movement of this round, not from what xtb called it.
    assert "if moved_now <= _OPTIMISE_STILL else 0" in apply_body


def test_the_viewer_is_sent_every_frame_of_the_path(editor, monkeypatch):
    """The whole path, not a sample of it.

    The frame field is one slot, not a queue: a write that lands before the
    page has read the one before it replaces it.  The eight-frame tail was the
    insurance against that -- it re-sent recent frames so a missed read still
    caught them -- but eight is a *fixed* number and xtb makes frames far
    faster than the page is asked to look.  A benzene runs 23 cycles in a
    fraction of a second and a 149-atom chain 260, so everything between two
    reads beyond the last eight was never sent at all.

    The window starts where the previous window started now, so every frame
    goes out in two consecutive writes -- the same insurance, and nothing
    skipped however fast they arrive.
    """
    import json
    import time as _time

    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    editor["submit_ff_dd"].value = "gfn2"

    total, atoms = 60, 3
    path = [[float(i) + 0.001 * k for k in range(atoms * 3)]
            for i in range(total)]

    def fake(xyz, method, **kw):
        on = kw.get("on_frames")
        # Bursts of twenty: far more than the eight the tail used to carry,
        # and exactly the shape a real run has between two reads.
        for cut in (20, 40, 60):
            if on:
                on(path[:cut])
        return {"ok": True, "xyz": xyz, "energy": -1.0, "converged": True,
                "frames": path, "status": "converged in 1.0 s"}

    monkeypatch.setattr(tab_submit._gfn, "optimize_with_gfn", fake)

    seen: list[str] = []
    editor["submit_gfn_frame"].observe(
        lambda c: seen.append(c["new"]), names="value")

    editor["submit_optimize_btn"].value = True
    deadline = _time.time() + 30
    while _time.time() < deadline and editor["submit_optimize_btn"].value:
        _time.sleep(0.02)

    carried: dict[int, int] = {}
    for text in seen:
        try:
            payload = json.loads(text)
        except ValueError:
            continue
        if not payload.get("frames"):
            continue
        first = int(payload.get("from") or 0)
        for k in range(len(payload["frames"])):
            carried[first + k] = carried.get(first + k, 0) + 1

    missing = sorted(set(range(total)) - set(carried))
    assert not missing, f"frames never sent to the viewer: {missing}"
    # And each one carried twice, which is the insurance the tail used to give.
    thin = sorted(i for i, n in carried.items() if n < 2)
    assert not thin, f"frames sent only once, a missed read loses them: {thin}"


def test_a_frame_is_not_sent_twice_over_when_nothing_is_new(editor, monkeypatch):
    """A write that carries nothing new is a message for no reason, and the
    field is read sixty times a second."""
    import json
    import time as _time

    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    editor["submit_ff_dd"].value = "gfn2"
    path = [[float(i)] * 9 for i in range(12)]

    def fake(xyz, method, **kw):
        on = kw.get("on_frames")
        if on:
            on(path)          # the same trajectory, three times over
            on(path)
            on(path)
        return {"ok": True, "xyz": xyz, "energy": -1.0, "converged": True,
                "frames": path, "status": "converged in 1.0 s"}

    monkeypatch.setattr(tab_submit._gfn, "optimize_with_gfn", fake)
    seen: list[str] = []
    editor["submit_gfn_frame"].observe(
        lambda c: seen.append(c["new"]), names="value")

    editor["submit_optimize_btn"].value = True
    deadline = _time.time() + 30
    while _time.time() < deadline and editor["submit_optimize_btn"].value:
        _time.sleep(0.02)

    with_frames = [t for t in seen
                   if json.loads(t).get("frames")] if seen else []
    assert len(with_frames) <= 2, (
        f"nothing new arrived, yet {len(with_frames)} writes carried frames")


def test_the_playback_pace_is_the_users_and_not_the_backlogs(player_js):
    """The pacing sped the playback up whenever the queue grew.

    That was right while a trailing picture counted as a fault.  It is the
    behaviour a slow setting exists to ask for -- the whole path is in the page
    now, so where the picture has got to is where the user is, and grabbing an
    atom there keeps that frame and drops what was computed past it.  A rule
    that quietly sped it back up would take that away, so the setting wins.
    """
    step = player_js.split("function stepMs()")[1].split("function ")[0]
    assert "if(play.pace) return play.pace;" in step
    # And it is asked first, before any of the backlog rules.
    assert step.index("play.pace") < step.index("n>60")


def test_the_speed_slider_reaches_the_page_in_the_unit_it_counts_in(tmp_path):
    """Frames a second on the slider, milliseconds a frame on the page."""
    pytest.importorskip("ipywidgets")
    from delfin.dashboard.context import DashboardContext
    from delfin.dashboard import tab_submit

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    sent: list[str] = []
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = sent.append
    _widget, refs = tab_submit.create_tab(ctx)
    refs["coords_widget"].value = _WATER

    slider = refs["submit_play_speed"]
    sent.clear()
    slider.value = 25                       # 25 frames a second
    said = "".join(sent)
    assert ".pace=40;" in said, said[-400:]  # 1000/25

    sent.clear()
    slider.value = 4                        # slow: a quarter of a second each
    assert ".pace=250;" in "".join(sent)

    sent.clear()
    slider.value = 1                        # the floor: a second a frame
    assert ".pace=1000;" in "".join(sent)


def test_the_pace_is_offered_only_where_a_path_is_walked(editor, monkeypatch):
    """The browser's own field draws its frames as it computes them; there is
    no queue to pace.  A server engine walks a path worth watching."""
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    dd = editor["submit_ff_dd"]

    dd.value = "uff"
    assert editor["submit_play_speed"].layout.display == "none"
    assert editor["submit_strength_slider"].layout.display == ""

    dd.value = "gfn2"
    assert editor["submit_play_speed"].layout.display == ""
    assert editor["submit_strength_slider"].layout.display == "none"

    dd.value = "pm7"
    assert editor["submit_play_speed"].layout.display == ""


def test_the_pace_control_is_on_the_toolbar(editor):
    toolbar = editor["submit_manip_toolbar"]
    assert editor["submit_play_speed"] in set(toolbar.children)


def test_the_grab_says_which_frame_the_picture_stood_on():
    """The message went out with no frame on it.

    So the kernel knew a hand had arrived and not where the picture stood, and
    the geometry the user had taken hold of lived only in the browser until
    some later drag pushed it back.  Taking hold and letting go without moving
    left the box and the picture apart -- the box holding what was there
    before the run, the picture holding the frame that was grabbed.
    """
    watcher = EDITOR_SOURCE.split("def _install_gfn_frame_watcher")[1].split(
        "\n    def ")[0]
    assert 'held?String(play.shown||0):""' in watcher, (
        "the grab has to name the frame it happened on")

    cmd = EDITOR_SOURCE.split("def on_submit_cmd")[1].split("\n    def ")[0]
    grab = cmd.split("if verb == 'gfngrab':")[1].split("if verb ==")[0]
    assert "state['gfn_shown_frame'] = int(str(payload).strip())" in grab
    # and it is read before the run is cut, or the cut has nothing to go on
    assert grab.index("gfn_shown_frame") < grab.index("_interrupt_gfn()")


def test_an_interrupted_run_leaves_the_frame_that_was_on_screen():
    """Not where xtb had got to, which nobody saw.

    The whole point of cutting at the grab is that the frames past the picture
    were computed for a structure the hand is changing.  Writing none of them
    left the coordinates behind the picture instead of level with it.
    """
    apply_body = SUBMIT_SOURCE.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0].split("def _apply()")[1]
    cut = apply_body.split("if state.get('optimize_interrupted') is token:")[1]
    cut = cut.split("return")[0]
    assert "shown = state.get('gfn_shown_frame')" in cut
    assert "walked = trail[0]" in cut
    assert "'stopped where you took hold'" in cut
    # The path is kept somewhere a cut-short run can still reach it: reading
    # the loop's own name there is a crash in the one place written for a run
    # that ended early.
    whole = SUBMIT_SOURCE.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0]
    assert "trail = [None]" in whole
    assert "trail[0] = outcome['frames']" in whole


def test_a_write_that_changes_nothing_lowers_the_render_flag(bare_editor):
    """The picture stopped following the box, and cutting the coordinates out
    and pasting them back was the cure.

    manip_inflight tells the host that the playback has drawn this geometry
    already, so it must not redraw.  The host clears it in update_view, and
    traitlets only calls that when the value actually changes.  Raised before a
    write that changes nothing -- a structure already at its minimum, written
    back by a second press of Optimise, which is the ordinary case -- it stayed
    raised, and the *next* genuine change was the one swallowed.  Cut and paste
    is two real changes, which is why that cleared it.
    """
    part, state = bare_editor

    same = part.coords_widget.value
    state["manip_inflight"] = False
    assert part._write_coords(same, drawn=True) is False, (
        "nothing changed, so nothing was written")
    assert state["manip_inflight"] is False, (
        "a write that changes nothing must not leave the flag raised")

    other = (same or "") + "\nH 9.0 9.0 9.0"
    assert part._write_coords(other, drawn=True) is True
    assert state["manip_inflight"] is True, "the playback has drawn this one"


def test_every_geometry_the_editor_writes_goes_through_the_one_door():
    """Three places wrote the box and raised the flag by hand, and each was
    one no-op write away from stopping the picture for the rest of the
    session."""
    for call in ("_write_coords(xyz_document(lines, f'Settled with {label}')",
                 "xyz_document(lines, 'Optimised in DELFIN viewer')",
                 "'stopped where you took hold'), drawn=True)"):
        assert call in EDITOR_SOURCE, call
    # No geometry is written beside a hand-raised flag any more.
    assert "state['manip_inflight'] = True\n                    coords_widget" \
        not in EDITOR_SOURCE
    assert EDITOR_SOURCE.count("_write_coords(") >= 4


def test_the_path_is_played_at_a_pace_that_can_be_watched(editor):
    """Twelve a second: a path being walked, not a jump.

    It sat at the top of the range for a while, on the theory that a viewer
    taking longer to arrive was the complaint.  It was not -- that was a flag
    left raised over a write that changed nothing, fixed where it was.  At
    sixty a second a twenty-three-frame optimisation is over in 0.4 s, which
    is the single jerk the playback exists to replace; at twelve it is 1.9 s,
    and a 260-frame one is 21.7 s.  Both ends of that are a setting away.
    """
    slider = editor["submit_play_speed"]
    assert slider.value == 12
    assert slider.description == "Speed"
    # One a second at the bottom -- already a second of looking at each
    # geometry -- and sixty at the top, which is not a speed so much as
    # "live": the frames only exist as fast as xtb makes them, so on a large
    # structure the picture drains the queue and waits, and that is the
    # calculation watched as it happens.  Whole steps: a tenth of a frame a
    # second is a distinction nobody makes.
    assert (slider.min, slider.max) == (1, 60)
    assert slider.step == 1
    assert slider.value < slider.max, "the default cannot be the fastest there is"


def test_the_frame_a_hand_arrived_on_belongs_to_one_run(editor, monkeypatch):
    """A number left over from an earlier grab is a plausible index into the
    next run's path.

    The page reports it only when a hand arrives or the switch goes up, so it
    survives whole runs otherwise -- and an edit interrupting a later run would
    cut that run at a frame nobody ever saw, taken from a trajectory that no
    longer exists.  The index is only ever about the run it came from.
    """
    import time as _time

    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    state = editor["editor_state"]
    editor["submit_ff_dd"].value = "gfn2"

    state["gfn_shown_frame"] = 13          # from a grab during an earlier run
    monkeypatch.setattr(
        tab_submit._gfn, "optimize_with_gfn",
        lambda xyz, method, **kw: {"ok": True, "xyz": xyz, "energy": -1.0,
                                   "converged": True, "frames": [],
                                   "status": "converged in 1.0 s"})
    editor["submit_optimize_btn"].value = True
    deadline = _time.time() + 20
    while _time.time() < deadline and editor["submit_optimize_btn"].value:
        _time.sleep(0.02)

    assert state.get("gfn_shown_frame") is None, (
        "a new run starts with no frame remembered from the last one")

    start = SUBMIT_SOURCE.split(
        "def on_submit_optimize(change=None, every_frame=False)"
    )[1].split("\n    def ")[0]
    assert "state.pop('gfn_shown_frame', None)" in start


# ---------------------------------------------------------------------------
# what this structure can do at this temperature
# ---------------------------------------------------------------------------


def test_the_ceiling_is_eyring_turned_around():
    """"Possible at this temperature" has two halves, and the second is
    usually left out.

    A rate is k = (kB T / h) exp(-dG/RT), so a barrier crossed within tau may
    be at most R T ln(kB T tau / h).  At 298 K one second buys 17.5 kcal/mol
    and a year 27.7 -- the same temperature, four times the chemistry.  A
    ceiling quoted without a timescale is not a number.
    """
    from delfin.dashboard.structure_editor import thermal_ceiling as ceiling

    table = {
        (298.15, 1.0): 17.5, (298.15, 3600.0): 22.3,
        (298.15, 3.15576e7): 27.7,
        (198.15, 3600.0): 14.7, (373.15, 3600.0): 28.1,
        (773.15, 3600.0): 59.3,
    }
    for (T, tau), expected in table.items():
        got = ceiling(T, tau)
        assert abs(got - expected) < 0.1, f"{T} K / {tau} s: {got:.1f}"

    # Hotter and longer both buy more, and neither is optional.
    assert ceiling(373.15, 3600.0) > ceiling(298.15, 3600.0)
    assert ceiling(298.15, 86400.0) > ceiling(298.15, 1.0)
    # Never negative, however short the time.  A picosecond is the floor --
    # below a tenth of that a molecule has not finished one vibration, so
    # there is no crossing to speak of -- and everything under it lands there.
    assert ceiling(298.15, 1e-20) == ceiling(298.15, 1e-12)
    assert 0.0 <= ceiling(298.15, 1e-20) < 2.0


def test_the_budget_forbids_tearing_a_ring_and_allows_a_real_distortion():
    """Measured with GFN2 on a benzene, the ring bond stretched with
    everything else relaxed -- which is the quantity the follow reports.

        1.45 A   +0.5 kcal/mol      1.75 A  +33.3
        1.55 A   +8.4               1.90 A  +54.5
        1.65 A  +20.0

    At 298 K within an hour the structure has 22.3, so the wall falls at about
    1.66 A.  The ring cannot be pulled open, and that falls out of the energy
    rather than out of a rule about aromatics -- while 1.55 A, a real
    distortion, stays available.
    """
    from delfin.dashboard.structure_editor import thermal_ceiling as ceiling

    room = ceiling(298.15, 3600.0)
    measured = {1.45: 0.5, 1.55: 8.4, 1.65: 20.0, 1.75: 33.3, 1.90: 54.5}
    allowed = {r for r, cost in measured.items() if cost <= room}
    assert allowed == {1.45, 1.55, 1.65}, allowed

    # And at 773 K the same bond reaches 1.9 A, which is the temperature at
    # which benzene really does come apart.
    assert measured[1.90] <= ceiling(773.15, 3600.0)


def test_a_single_point_leaves_the_geometry_alone():
    """The anchor is the energy of the structure as it stands.

    It shares every other decision with the optimisation -- which binary,
    which parameters, how many cores, the solvent, the charge, the spin --
    which is why it is a flag on that function rather than a second one that
    would drift from it.
    """
    source = open(_gfn_source(), encoding="utf-8").read()
    body = source.split("def optimize_with_gfn(")[1].split("\ndef ")[0]
    assert "optimise: bool = True" in body
    assert "*(['--opt'] if optimise else [])" in body
    # No --opt writes no xtbopt.xyz and no path: reading them would hand back
    # None and read as a failure.
    assert "if optimise else xyz_text" in body
    assert "read_trajectory(folder) if optimise else []" in body
    # And a single point has no geometry to converge, so it is not counted
    # among the runs that stopped short.
    assert "(not optimise) or 'GEOMETRY OPTIMIZATION CONVERGED' in output" in body


def _gfn_source():
    from delfin.dashboard import gfn_optimize
    return gfn_optimize.__file__


def test_the_thermal_controls_belong_to_xtb_alone(editor, monkeypatch):
    """MOPAC's follow reports a heat of formation, which is not the quantity
    an xtb anchor is in and cannot be differenced against it."""
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, "find_binary", lambda _m=None: "/x/xtb")
    dd = editor["submit_ff_dd"]

    dd.value = "gfn2"
    assert editor["submit_thermal_btn"].layout.display == ""
    dd.value = "pm7"
    assert editor["submit_thermal_btn"].layout.display == "none"
    dd.value = "uff"
    assert editor["submit_thermal_btn"].layout.display == "none"


def test_the_anchor_belongs_to_one_structure(bare_editor):
    """An energy measured on one molecule says nothing about another.

    Carried over quietly, the budget would report the difference between two
    unrelated energies as though it were a distortion of the second.  The
    anchor names its molecule the way the perception and the GFN-FF topology
    do -- by the element column -- so it stops answering by itself rather than
    by a callback firing.  structure_changed is called from one place, and a
    structure can arrive without it.
    """
    part, state = bare_editor
    part.submit_ff_dd.value = "gfn2"
    part.submit_temperature.value = 298.15
    part.submit_timescale.value = 3600.0

    here = part.coords_widget.value
    state["thermal_e0"] = -15.877561
    state["thermal_for"] = part._structure_fingerprint(here)

    anchor, ceiling = part._thermal_budget()
    assert anchor == -15.877561
    assert abs(ceiling - 22.3) < 0.1
    assert part._thermal_note(-15.870000), "it reports while it is about this one"

    # Another molecule, and the same stored anchor answers for nothing.
    part.coords_widget.value = (
        "2\ncarbon monoxide\nC 0.000 0.000 0.000\nO 1.128 0.000 0.000\n")
    anchor, ceiling = part._thermal_budget()
    assert anchor is None, "the budget cannot outlive the molecule"
    assert abs(ceiling - 22.3) < 0.1, "the ceiling is about T, not about this"
    assert part._thermal_note(-15.870000) == ""

    # T and the timescale are not per-structure: they are how the user is
    # working, like the method and the solvent.
    controls = EDITOR_SOURCE.split("def _structure_controls")[1].split(
        "\n    def ")[0]
    assert "submit_temperature" not in controls
    assert "submit_thermal_btn" not in controls


def test_the_budget_line_goes_on_the_row_that_is_already_there():
    """The status row stands above the viewer, and a second row moves the
    picture while the user is aiming an atom -- measured once already."""
    follow = EDITOR_SOURCE.split("def _gfn_follow_step")[1].split(
        "\n    def ")[0]
    assert "state['thermal_now'] = outcome.get('energy')" in follow, (
        "the follow already computes the energy; it has to be read")
    assert "said = f'{said} {spent}'" in follow
    # One line handed to the status, never two.
    assert "_gfn_status_lines(said)" in follow
