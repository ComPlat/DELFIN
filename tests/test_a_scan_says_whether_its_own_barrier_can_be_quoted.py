"""A driven scan walks itself back, and says whether the two legs are one path.

A relaxed scan holds one coordinate and relaxes everything else, and nothing
makes the rest follow continuously.  Where it does not, the profile depends on
which way the walk went and neither direction is the path.  Jonsson, Mills and
Jacobsen said it in 1998: "the path generated may be discontinuous and the
procedure may depend on the direction of the drag (hysteresis effects).  In
particular, some atomic coordinates may 'slip' near the saddle point region
and the saddle point configuration will then be missed."  Bofill and Quapp
(Mol. Phys. 2019) give the condition it holds under -- no turning point and no
valley-ridge inflection -- and there is no way to test that from one leg.

So the editor walks the second leg.  Measured here, butadiene and ethylene
under GFN2 with one forming C-C driven from 3.40 A to 1.60 and back, 0.1 A at
a time, everything relaxed at every point:

                                 forward   backward   ORCA's converged saddle
    apparent barrier               +7.25     +11.72         +6.74 kcal/mol
    where the maximum is            2.20 A     2.90 A         2.315 A
    the *undriven* forming C-C
       at that maximum              2.92 A     1.76 A         2.315 A
    the mode that goes the
       wrong way                 -128 cm-1  -267 cm-1      -394 cm-1

    largest gap at the same coordinate:  23.77 kcal/mol, at 2.20 A

The saddle is symmetric -- both forming bonds at 2.315 A -- and neither leg of
the scan is anywhere near it: one arrives with the undriven bond a whole
Angstrom too long and the other with it half an Angstrom too short, and each
misses the crossing on its own side.  Both maxima carry exactly one imaginary
frequency, at gradient norms far from converged, so that test says nothing
here.  A user reading the barrier off either leg alone is out by more than
half.

The null result is why this is on by default rather than a curiosity.  Eleven
scans were walked out and back:

    ethane C-C stretch                     0.000 kcal/mol apart
    propane C-C-C angle                    0.001
    butanol, butane, glycol torsions       0.002, 0.004, 0.006
    SN2 Cl-/CH3Cl                          0.032
    water dimer O-O                        0.803
    --------------------------------------------- RT ln 10 = 1.364
    N-methylacetamide torsion             14.2
    Diels-Alder                           23.8
    formate/water H transfer              60.8
    methylcyclopropane ring opening      129.1

Seven of them are scans worth quoting and the editor had no way to say so.
And the two groups are not two kinds of coordinate: an amide torsion jumps and
an SN2 -- one bond made as another breaks -- does not.  There is no rule to
apply instead of walking the second leg.

The threshold has a meaning rather than a round shape: RT ln 10, which at
298 K is 1.36 kcal/mol, is the difference in barrier that is a factor of ten
in rate.  Two legs nearer than that are the same answer to any use a barrier
is put to; two legs further apart than it are two different answers.  It is a
temperature, so it is worked out at the temperature the walk was priced at.

And where the two legs disagree, the second difference of the energy names the
step it happened at and an undriven internal coordinate names the culprit --
which is what a user needs, because the answer to a scan that jumped is to arm
the coordinate that slipped as well and walk both together.
"""

from __future__ import annotations

import math
import shutil

import pytest

from delfin.dashboard import gfn_optimize as gfn
from editor_source import EDITOR_SOURCE

_needs_xtb = pytest.mark.skipif(not shutil.which('xtb'),
                                reason='xtb not installed')

#: The Diels-Alder, measured.  (driven C-C in Angstrom, kcal/mol against the
#: first point of the forward leg), in the order each leg was walked.
_FWD = [(3.4, 0.0), (3.3, -0.012), (3.2, 0.012), (3.1, 0.086), (3.0, 0.217),
        (2.9, 0.436), (2.8, 0.792), (2.7, 1.113), (2.6, 1.803), (2.5, 2.763),
        (2.4, 4.025), (2.3, 5.544), (2.2, 7.254), (2.1, -24.522),
        (2.0, -33.445), (1.9, -42.303), (1.8, -50.897), (1.7, -58.214),
        (1.6, -62.955)]
_BACK = [(1.6, -62.955), (1.7, -60.088), (1.8, -52.375), (1.9, -42.538),
         (2.0, -32.218), (2.1, -22.805), (2.2, -16.511), (2.3, -9.386),
         (2.4, -3.336), (2.5, 1.561), (2.6, 5.364), (2.7, 8.212),
         (2.8, 10.282), (2.9, 11.725), (3.0, 0.242), (3.1, 0.101),
         (3.2, 0.029), (3.3, 0.001), (3.4, 0.002)]

#: An ethylene glycol OCCO torsion turned right round and back, measured the
#: same way.  A scan that is a path, and the editor should be able to say so.
_GLYCOL_GRID = [-180.0 + 20.0 * n for n in range(19)]
_GLYCOL_FWD = [0.0, 0.3713, 1.4653, 2.4043, 2.2417, 1.0501, 0.0517, 0.0278,
               0.7163, 1.2479, 1.0614, 0.4744, 0.2265, 0.5998, 1.3291, 1.837,
               1.5226, 0.5859, -0.002]
_GLYCOL_BACK = [0.0, 0.588, 1.5242, 1.839, 1.331, 0.6015, 0.2277, 0.4763,
                1.0629, 1.2499, 0.7182, 0.0298, 0.0564, 1.0576, 2.2452,
                2.4068, 1.4673, 0.3735, 0.0003]


# ---------------------------------------------------------------------------
# what a factor of ten in rate is worth
# ---------------------------------------------------------------------------


def test_the_threshold_is_a_rate_and_not_a_round_number():
    assert math.isclose(gfn.a_rate_apart(298.15), 1.3642, abs_tol=5e-4)
    # It is a temperature: a scan run hot is judged at the temperature it was
    # run, because the same kcal/mol buys less there.
    assert gfn.a_rate_apart(500.0) > gfn.a_rate_apart(298.15)
    assert gfn.a_rate_apart(77.0) < gfn.a_rate_apart(298.15)
    # And nonsense does not raise: an unset box is 298 K.
    assert gfn.a_rate_apart(None) == gfn.a_rate_apart(298.15)


# ---------------------------------------------------------------------------
# the two legs compared
# ---------------------------------------------------------------------------


def test_the_diels_alder_legs_disagree_by_more_than_the_barrier_itself():
    found = gfn.paths_disagree(_FWD, _BACK)
    assert found is not None
    assert math.isclose(found['gap'], 23.766, abs_tol=0.01)
    assert found['at'] == 2.2
    assert found['points'] == 19
    # And it is worse than the barrier either leg reports.
    assert found['gap'] > max(one[1] for one in _FWD)
    assert found['gap'] > gfn.a_rate_apart(298.15)


def test_a_torsion_that_is_a_path_says_so():
    there = list(zip(_GLYCOL_GRID, _GLYCOL_FWD))
    back = list(zip(reversed(_GLYCOL_GRID), _GLYCOL_BACK))
    found = gfn.paths_disagree(there, back)
    assert found is not None
    assert found['gap'] < 0.01, found
    assert found['gap'] < gfn.a_rate_apart(298.15)
    assert found['points'] == 19


def test_it_compares_at_the_coordinate_and_not_at_the_two_maxima():
    """The two maxima are not two views of one place.

    Forward puts its highest point at 2.20 A and backward at 2.90.  Comparing
    the heights would be comparing two different geometries and would report
    4.5 kcal/mol; comparing where the driven coordinate agrees reports 23.8,
    which is the disagreement that is really there.
    """
    heights = abs(max(one[1] for one in _FWD) - max(one[1] for one in _BACK))
    assert heights < 5.0
    assert gfn.paths_disagree(_FWD, _BACK)['gap'] > 20.0


def test_a_return_leg_that_was_cut_short_is_compared_over_what_it_walked():
    found = gfn.paths_disagree(_FWD, _BACK[:6])
    assert found is not None
    assert found['points'] == 6


def test_nothing_to_compare_is_not_an_answer():
    assert gfn.paths_disagree([], _BACK) is None
    assert gfn.paths_disagree(_FWD, []) is None
    assert gfn.paths_disagree([(1.0, 0.0)], [(1.0, 0.0)]) is None
    # Two legs that share no coordinate value have nothing to say about each
    # other, and saying it at whichever point happened to be nearest would be
    # inventing a comparison.
    away = [(v + 50.0, e) for v, e in _BACK]
    assert gfn.paths_disagree(_FWD, away) is None


# ---------------------------------------------------------------------------
# where it jumped, and what moved
# ---------------------------------------------------------------------------


def test_the_fall_is_found_on_the_second_difference():
    found = gfn.where_a_walk_jumped([one[1] for one in _FWD])
    assert found is not None
    # Between 2.20 and 2.10 A, which is where the pair falls into the product.
    assert found['step'] == 13
    assert math.isclose(found['fell'], -31.776, abs_tol=0.01)
    assert found['times'] > 100
    assert math.isclose(found['scale'], 0.256, abs_tol=0.01)


def test_the_backward_leg_jumped_too_and_in_a_different_place():
    found = gfn.where_a_walk_jumped([one[1] for one in _BACK])
    assert found is not None
    assert found['step'] == 14           # between 2.90 and 3.00 A
    assert found['times'] > gfn.WALK_JUMPED_TIMES


def test_a_torsion_is_never_called_a_jump():
    assert gfn.where_a_walk_jumped(_GLYCOL_FWD) is None
    assert gfn.where_a_walk_jumped(_GLYCOL_BACK) is None


def test_a_curve_that_is_out_of_line_but_small_is_not_a_jump():
    """The floor, and why it is not decoration.

    Measured on a water dimer whose O-O was stretched from 2.70 to 4.50 A and
    back: the largest second difference is 21 times the median on the way out
    and 38 times on the way back -- past the ratio on its own -- and it is
    0.62 kcal/mol.  Two waters turning over inside a hydrogen bond is not a
    scan that jumped, and the two legs agree to 0.803 kcal/mol, which is
    inside a factor of ten in rate.
    """
    flat = [0.0, 0.001, 0.002, 0.003, 0.004, 0.62, 0.005, 0.006, 0.007]
    assert max(flat) < gfn.WALK_JUMPED_AT_LEAST
    assert gfn.where_a_walk_jumped(flat) is None


def test_a_path_too_short_to_have_a_middle_is_not_judged():
    assert gfn.where_a_walk_jumped([0.0, 5.0, 0.0]) is None
    assert gfn.where_a_walk_jumped([]) is None
    assert gfn.where_a_walk_jumped([0.0, 'x', 1.0, 2.0]) is None


def test_it_names_the_step_the_fall_is_in_and_not_the_one_beside_it():
    """A step shows in the second difference at both of its ends.

    Which of the two is the fall is decided by the plain difference, so the
    step named is the one the energy actually moved in.
    """
    made = [0.0, 0.1, 0.2, 0.3, -30.0, -30.1, -30.2, -30.3]
    found = gfn.where_a_walk_jumped(made)
    assert found is not None and found['step'] == 4


def test_the_culprit_is_an_internal_coordinate_nobody_was_driving():
    """Named by its two atoms, so it can be armed."""
    before = ('4\nx\nC 0 0 0\nC 2.2 0 0\nC 0 2.9 0\nC 2.2 2.9 0\n')
    after = ('4\nx\nC 0 0 0\nC 2.1 0 0\nC 0 1.55 0\nC 2.1 1.55 0\n')
    found = gfn.what_else_moved(before, after, [[0, 1]])
    assert found is not None
    # 0-2 is the pair that moved 1.35 A; 0-1 is the driven one and is not
    # eligible however much it moved.
    assert found['pair'] == (0, 2)
    assert math.isclose(found['moved'], 1.35, abs_tol=1e-6)
    assert gfn.pair_named(found['pair'], ['C', 'C', 'C', 'C']) == 'C1-C3'


def test_every_atom_of_a_driven_leg_is_left_out_of_the_comparison():
    """An angle drives three atoms and a torsion four, and all of them move."""
    assert set(gfn._pairs_of([[1, 2, 3]])) == {(1, 2), (1, 3), (2, 3)}
    assert (0, 3) in gfn._pairs_of([[0, 1, 2, 3]])
    # A single leg given flat rather than as a list of legs reads the same,
    # because leaving the driven pair in would make it the culprit on every
    # step of every scan.
    assert gfn._pairs_of([0, 1]) == [(0, 1)]


def test_two_geometries_that_cannot_be_compared_say_nothing():
    assert gfn.what_else_moved('', '', []) is None
    assert gfn.what_else_moved('2\nx\nC 0 0 0\nC 1 0 0\n',
                               '3\nx\nC 0 0 0\nC 1 0 0\nC 2 0 0\n', []) is None


# ---------------------------------------------------------------------------
# and what the editor does with all that
# ---------------------------------------------------------------------------


def test_the_return_leg_retraces_the_points_the_walk_really_held():
    """Not the range that was armed.

    A scan stops at the next minimum by default, so the values it walked are a
    part of what was set up -- and it is that part the second leg has to
    retrace, or the two are not comparable at any coordinate.
    """
    assert 'drove.append(' in EDITOR_SOURCE
    assert 'for back_n, values in enumerate(' in EDITOR_SOURCE
    assert 'reversed(drove[:-1])' in EDITOR_SOURCE


def test_it_starts_from_where_the_walk_ended_and_not_from_the_box():
    """After coming back to a minimum the box holds the minimum, not the end."""
    assert 'here = standing if standing is not None else walked' in EDITOR_SOURCE


def test_it_does_not_run_after_a_stop_a_collapse_or_a_push():
    assert ("if (not pushing and bool(submit_scan_back.value)"
            in EDITOR_SOURCE)
    assert "and len(drove) > 2 and not state.get('scan_stop')" in EDITOR_SOURCE
    assert "and state.get('scan_crowded') is None" in EDITOR_SOURCE


def test_the_second_leg_is_watched_like_the_first():
    """Down the frame channel, not into the box: a write to the box rebuilds
    the viewer from nothing, and the picture is what makes a jump obvious."""
    assert EDITOR_SOURCE.count("'follow': 1,") >= 2
    assert 'is walking it back: step ' in EDITOR_SOURCE


def test_the_jump_is_looked_for_on_a_walk_and_not_on_a_push():
    """A push means to fall through its crossing and prices it afterwards.

    The same test on a push would fire on the thing a push is for, and on a
    path whose points are not evenly spaced besides -- ``_across`` inserts its
    own between two of them, and a second difference is about even steps.
    """
    where = EDITOR_SOURCE.index('_gfn.where_a_walk_jumped(')
    assert 'if not pushing:' in EDITOR_SOURCE[where - 300:where]


def test_the_verdict_says_which_of_the_three_things_happened():
    assert 'Both legs agree to {gap["gap"]:.2f} kcal/mol' in EDITOR_SOURCE
    assert 'legs disagree by {gap["gap"]:.1f} kcal/mol' in EDITOR_SOURCE
    assert 'Nothing walked it back' in EDITOR_SOURCE
    assert 'arm that as well and ' in EDITOR_SOURCE
    # And it is said beside the barrier it is about rather than after the
    # temperature, because a caveat read after the number is a caveat nobody
    # applied.
    assert 'bits += _phrases(_scan_can_be_quoted(T))' in EDITOR_SOURCE


def _an_editor(text):
    """One structure editor over a coordinate box of its own.

    Driven the way a user drives it, because whether the second leg runs is a
    fact about the widgets and the state they answer to rather than about what
    the source says it means to do.
    """
    import pathlib
    import tempfile

    pytest.importorskip('ipywidgets')
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room = pathlib.Path(tempfile.mkdtemp())
    for name in ('calc', 'archive', 'office'):
        (room / name).mkdir()
    ctx = DashboardContext(calc_dir=room / 'calc', archive_dir=room / 'archive',
                           office_dir=room / 'office')
    ctx.run_js = lambda _script: None
    return structure_editor.build(
        ctx, state={}, coords_widget=widgets.Textarea(value=text),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)


_BUTANE = """14
butane
C -1.9026 -0.2415  0.0000
C -0.6178  0.5776  0.0000
C  0.6178 -0.5776  0.0000
C  1.9026  0.2415  0.0000
H -1.9741 -0.8776  0.8834
H -1.9741 -0.8776 -0.8834
H -2.7573  0.4382  0.0000
H -0.5834  1.2225  0.8813
H -0.5834  1.2225 -0.8813
H  0.5834 -1.2225  0.8813
H  0.5834 -1.2225 -0.8813
H  1.9741  0.8776  0.8834
H  1.9741  0.8776 -0.8834
H  2.7573 -0.4382  0.0000
"""


@_needs_xtb
def test_the_editor_really_walks_the_second_leg_and_says_they_agree():
    """The whole thing, on the real editor.

    A butane C-C-C-C torsion under GFN-FF, walked out and back.  Nothing slips
    in an alkane torsion, so what the editor has to be able to do here is say
    that the profile is trustworthy -- the check is not worth having if it can
    only ever be bad news.
    """
    import time

    part = _an_editor(_BUTANE)
    part.submit_ff_dd.value = 'gfnff'
    part.submit_scan_how.value = 'hold'
    part.submit_scan_back.value = True
    part.submit_scan_whole.value = True     # the whole curve, not the first well
    part.state['scan_legs'] = [
        {'kind': 'dihedral', 'atoms': [0, 1, 2, 3], 'from': 180.0, 'to': 100.0,
         'steps': 6, 'structure': None},
    ]
    part.on_submit_scan_run()
    began = time.time()
    while part.state.get('scan_run') and time.time() - began < 600:
        time.sleep(0.05)
    assert not part.state.get('scan_run'), 'the scan never finished'

    both = part._scan_two_legs()
    assert len(both['there']) >= 5
    # The second leg retraces the first, so it has the points the walk really
    # held minus the one it starts standing on -- plus that one as its anchor.
    assert len(both['back']) == len(both['there'])
    assert both['back'][0][0] == pytest.approx(both['there'][-1][0])
    assert both['back'][-1][0] == pytest.approx(both['there'][0][0])
    # And they agree, which is the answer for a torsion of an alkane.
    assert both['disagree'] is not None
    assert both['disagree']['gap'] < gfn.a_rate_apart(298.15), both['disagree']
    assert both['jumped'] is None
    said = ' '.join(part.state['mol_status_lines'])
    assert 'legs agree to' in said, said
    assert 'the profile is the path' in said


@_needs_xtb
def test_the_second_leg_is_not_walked_when_the_toggle_is_up():
    """Off is one press, for a large system where the second leg is felt."""
    import time

    part = _an_editor(_BUTANE)
    part.submit_ff_dd.value = 'gfnff'
    part.submit_scan_how.value = 'hold'
    part.submit_scan_back.value = False
    part.submit_scan_whole.value = True
    part.state['scan_legs'] = [
        {'kind': 'dihedral', 'atoms': [0, 1, 2, 3], 'from': 180.0, 'to': 140.0,
         'steps': 4, 'structure': None},
    ]
    part.on_submit_scan_run()
    began = time.time()
    while part.state.get('scan_run') and time.time() - began < 600:
        time.sleep(0.05)
    assert not part.state.get('scan_run')
    both = part._scan_two_legs()
    assert both['there'] and not both['back']
    assert both['disagree'] is None
    said = ' '.join(part.state['mol_status_lines'])
    assert 'Nothing walked it back' in said, said


@_needs_xtb
def test_a_torsion_walked_out_and_back_really_does_agree(tmp_path):
    """The null result, run live.

    A butane C-C-C-C torsion under GFN-FF, out five points and back over the
    same five.  Nothing slips in an alkane torsion, so the two legs are the
    same curve -- and this is the case the editor has to be able to call
    trustworthy, or the check is only ever bad news.
    """
    butane = ('14\nbutane\n'
              'C -1.9026 -0.2415  0.0000\nC -0.6178  0.5776  0.0000\n'
              'C  0.6178 -0.5776  0.0000\nC  1.9026  0.2415  0.0000\n'
              'H -1.9741 -0.8776  0.8834\nH -1.9741 -0.8776 -0.8834\n'
              'H -2.7573  0.4382  0.0000\nH -0.5834  1.2225  0.8813\n'
              'H -0.5834  1.2225 -0.8813\nH  0.5834 -1.2225  0.8813\n'
              'H  0.5834 -1.2225 -0.8813\nH  1.9741  0.8776  0.8834\n'
              'H  1.9741  0.8776 -0.8834\nH  2.7573 -0.4382  0.0000\n')
    grid = [180.0, 160.0, 140.0, 120.0, 100.0]
    here, there = butane, []
    for want in grid:
        got = gfn.optimize_with_gfn(
            here, 'gfnff', max_steps=60, timeout=None,
            topology=tmp_path / 'topo',
            constraints=[{'kind': 'dihedral', 'atoms': [0, 1, 2, 3],
                          'mode': 'fix', 'value': want}])
        if not got['ok']:
            pytest.skip(got['status'])
        here = got['xyz']
        there.append((want, float(got['energy'])))
    zero = there[0][1]
    there = [(v, (e - zero) * 627.5094740631) for v, e in there]
    back = [there[-1]]
    for want in reversed(grid[:-1]):
        got = gfn.optimize_with_gfn(
            here, 'gfnff', max_steps=60, timeout=None,
            topology=tmp_path / 'topo',
            constraints=[{'kind': 'dihedral', 'atoms': [0, 1, 2, 3],
                          'mode': 'fix', 'value': want}])
        assert got['ok'], got['status']
        here = got['xyz']
        back.append((want, (float(got['energy']) - zero) * 627.5094740631))
    found = gfn.paths_disagree(there, back)
    assert found is not None and found['points'] == len(grid)
    assert found['gap'] < gfn.a_rate_apart(298.15), found
    assert gfn.where_a_walk_jumped([e for _, e in there]) is None


def test_the_toggle_belongs_to_walking_and_is_on_by_default():
    """It retraces the same points, so it belongs to walking a value.

    Driven rather than read.  This asserted the exact line that hid it, which
    is a check on how the rule is spelled and not on what it does -- and the
    rule grew a second case: a pull is a ramp of forces with the coordinate
    left out entirely, so it has no grid to retrace either, and the switch was
    being offered for a second leg that cannot exist.
    """
    assert 'submit_scan_back = widgets.ToggleButton(' in EDITOR_SOURCE
    assert "value=True, description='Walk it back'" in EDITOR_SOURCE

    part = _an_editor(_BUTANE)
    part.submit_ff_dd.value = 'gfn2'
    part.state['picked'] = [0, 1]
    part.submit_scan_way.value = 'to'
    part.submit_scan_to.value = 1.2
    part.on_submit_scan(None)
    part.submit_scan_gear.value = True

    part.submit_scan_how.value = 'hold'
    part._refresh_scan()
    assert part.submit_scan_back.layout.display == '', 'a walk can retrace'
    assert part.submit_scan_back.value is True, 'and it does so by default'

    for way in ('push', 'load'):
        part.submit_scan_how.value = way
        part._refresh_scan()
        assert part.submit_scan_back.layout.display == 'none', way
