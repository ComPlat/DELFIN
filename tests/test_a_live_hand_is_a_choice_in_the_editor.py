"""The live hand is wired in as a third choice, opt-in, beside pull and move.

The steered engine (:func:`climb.steer`, tested for its physics in
test_a_live_hand_steers_one_path) is put in the editor without touching the two
hands that were there: a new value in the hand box, a branch in the follow
worker that runs it, and the budget taught that it is a force hand like the
pull.  Pull and move fall through to exactly the block they always did, which
is why this cannot regress them.

These read the editor's source rather than build a dashboard, the way the
other hand-wiring tests do -- the wiring is the claim, and it is in the text.
"""
from __future__ import annotations

import pathlib

_EDITOR = (pathlib.Path(__file__).resolve().parents[1]
           / 'delfin' / 'dashboard' / 'structure_editor.py').read_text()


def test_the_box_offers_a_live_hand():
    """A third value, and it names the engine rather than the gesture."""
    assert "('live dynamics', 'live')" in _EDITOR


def test_the_follow_worker_runs_the_steered_engine_for_it():
    """The live hand branches out of the follow loop before the pull's own
    perception, runs its own answer, and continues -- so nothing below it
    reads a hand it was not written for.  The branch and the call are kept
    free of any climb name (a helper decides whether the engine can run), so
    the follow step itself names no climb machinery -- see
    test_the_climb_can_be_helped."""
    follow = _EDITOR.split('def _gfn_follow_step(')[1].split(
        "_start_background(_work, 'The relaxation under the hand')")[0]
    branch = follow.split('if _live_hand_runs(method):')[1][:2000]
    assert '_live_answer(current, holding, began' in branch, branch
    assert 'driven=driven' in branch, branch
    # The drive hand works out its coordinate before the call; the live hand
    # passes driven None.
    assert 'if _hand_drives():' in branch, branch
    assert '_coordinate_wish(current, holding)' in branch, branch
    assert 'climb' not in follow, 'the follow step must name no climb'
    # The live answer is a self-contained function at editor scope, not a
    # change threaded through the pull's path and not buried in the follow
    # step (so the source slices that guard the follow step do not read it).
    assert 'def _live_answer(current, holding, began, method' in _EDITOR
    assert 'def _live_answer(' not in follow
    # The coordinate the drive hand forces is worked out from the selection,
    # off editor scope and not in the follow step.
    assert 'def _coordinate_wish(' in _EDITOR
    assert 'def _coordinate_wish(' not in follow
    # Its engine only runs where climb can drive it: GFN2/GFN1/GFN-FF.
    runs = _EDITOR.split('def _live_hand_runs(')[1].split('def ')[0]
    assert 'in _climb.CLIMB_METHODS' in runs


def test_the_budget_treats_the_live_hand_as_a_force():
    """A force hand is priced, and both engines are forces -- so the budget is
    live for the live hand and its temperature box is on screen for it."""
    live = _EDITOR.split('def _thermal_live(')[1].split('def ')[0]
    assert 'return bool(submit_thermal_btn.value) and _hand_is_a_force()' in live
    force = _EDITOR.split('def _hand_is_a_force(')[1].split('def ')[0]
    assert 'return _hand_pulls() or _hand_is_live()' in force
    # The temperature box follows the same rule.
    assert 'shown = bool(xtb and _hand_is_a_force())' in _EDITOR


def test_the_pull_and_placement_engines_are_untouched():
    """The two old hands still key on _hand_pulls, which stayed strictly the
    pull -- the live hand did not widen it and so cannot have changed them."""
    pulls = _EDITOR.split('def _hand_pulls(')[1].split('def ')[0]
    assert "return str(submit_hand_dd.value) == 'pull'" in pulls
    # The pull's own force is still a pull-only thing.
    force = _EDITOR.split('def _pull_force(')[1].split('def ')[0]
    assert 'if not _hand_pulls():' in force


def test_the_steered_hands_are_offered_only_where_an_engine_can_drive_them():
    """The live and drive hands are the server's steered dynamics, so they need
    xtb gradients: offered under a steerable method, and not under MOPAC (no
    gradients) or a browser method (no server engine).  A choice made there
    survives a detour through either and comes back."""
    block = _EDITOR.split("driving = ('drive coordinate', 'drive')")[1]
    block = block.split("finally:")[0]
    # Only the steerable branch lists the steered hands: one options list names
    # both, and the two others (MOPAC, browser) do not.
    assert 'submit_hand_dd.options = [pulling, moving, living, driving]' in block
    assert block.count('living, driving]') == 1
    assert 'options = [moving]' in block            # MOPAC: placement only
    assert block.count('options = [pulling, moving]') == 1   # browser: neither
    # MOPAC and the browser both remember a steered choice rather than dropping
    # it, so switching method and back does not cost it.
    assert "state['hand_was'] = had" in block
    assert "had in ('pull', 'live', 'drive')" in block
    assert "had in ('live', 'drive')" in block


def test_the_drive_hand_forces_the_selected_coordinate():
    """The drive hand is the one that can force a stiff torsion.  It reads the
    coordinate from the selection -- 2/3/4 picked atoms a bond/angle/torsion --
    and drives it with climb.steer_coordinate through the same answer the live
    hand uses, only with a coordinate named."""
    assert "('drive coordinate', 'drive')" in _EDITOR
    drives = _EDITOR.split('def _hand_drives(')[1].split('def ')[0]
    assert "return str(submit_hand_dd.value) == 'drive'" in drives
    # It is a force hand, so the budget prices it.
    force = _EDITOR.split('def _hand_is_a_force(')[1].split('def ')[0]
    assert '_hand_drives()' in force
    # The coordinate/kind comes from the selection, the same 2/3/4 mapping the
    # constraints use.
    wish = _EDITOR.split('def _coordinate_wish(')[1].split('def ')[0]
    assert "_CONSTRAINT_KINDS.get(len(atoms))" in wish
    assert "state.get('picked')" in wish
    # And the engine driven is the coordinate one, not the Cartesian spring.
    answer = _EDITOR.split('def _live_answer(')[1][:6000]
    assert 'steer_coordinate(' in answer
    assert 'if driven is not None:' in answer


def test_the_smearing_the_pull_uses_is_reused_by_the_live_hand():
    """A closed frontier gap fails a gradient the same way it fails a
    relaxation, and the same electronic temperature rescues it -- the live
    answer reaches for _smearing_for and retries warmed, exactly as the pull's
    relax_steps path does."""
    answer = _EDITOR.split('def _live_answer(')[1].split('def _work():')[0]
    assert '_smearing_for(method)' in answer
    assert '_smearing_for(method, failed=True)' in answer
    assert 'etemp=' in answer


def test_the_live_hand_lags_like_the_pull():
    """The page's own force field lags a grabbed atom behind the cursor only
    when the share it is told is above zero; at zero it places the atom
    rigidly.  The live hand sends the share like the pull, so it lags -- which
    is the whole point of a live drag -- rather than sticking to the cursor.
    The drive hand stays at zero (it forces a coordinate, not the atom)."""
    share = _EDITOR.split('def _hand_share(')[1].split('def ')[0]
    assert 'if not (_hand_pulls() or _hand_is_live()):' in share
