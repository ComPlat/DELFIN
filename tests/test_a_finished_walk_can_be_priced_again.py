"""A finished profile can be priced again, and the walk says what it walked into.

Two checks on an answer, and both exist because a barrier that comes out of
this editor is a number somebody is going to use.

**The barrier is probably too small, and by how much is known.**  GFN2 has no
Fock exchange at all, which puts it at the far end of the self-interaction
axis; a transition state is stretched bonds and near-degenerate weakly bound
orbitals, which Bursch, Mewes, Hansen and Grimme (Angew. Chem. Int. Ed. 2022,
https://doi.org/10.1002/anie.202205735) say "often leads to a systematic
underestimation of their electronic energy and, in turn, barrier heights".
The published mean absolute deviations say the error is structured rather than
random -- 10.22 kcal/mol on pericyclic barriers, 8.12 on BHDIV10, 13.05 on BH9
aggregate, against 1.17 on BHROT27, where the coordinate is a rotation and no
bond is broken.  A torsion profile is nearly right; a bond-breaking one is
not.

Measured here through the editor's own controls -- a butadiene and an ethylene
walked together under GFN2 over nineteen points, then the new press::

    GFN2-xTB   highest +6.3 kcal/mol at 2.29 A, ending -64.0
    g-xTB      highest +21.7 kcal/mol at 2.21 A, ending -42.7

a shift of +15.4 against a published pericyclic deviation of 10.22, and the
top moved as well as rose.  One reaction is not a benchmark; the direction and
the size are what the literature says to expect.  It cost 0.40 s a point on a
quiet machine, against minutes to walk the scan again under another method --
which is why it is a press on a finished walk rather than a box set before
one, and why the geometries the walk produced are now kept instead of dropped
once their frame had been drawn.

What comes back is not a barrier, and the three reasons ship in the answer
rather than in a docstring: the geometries are still GFN2's and that method is
out by about 0.2 A on a half-broken bond against 0.03 A on an ordinary one;
g-xTB is a preprint method whose installed build calls itself a development
version differing from the paper; and the points are wherever the walk left
them.

**And whether one arrangement of the electrons was ever the structure.**
``xtb --fod`` runs one more SCF at 5000 K and counts how far the occupations
are from two-and-zero -- Grimme and Hansen, Angew. Chem. Int. Ed. 54 (2015)
12308, https://doi.org/10.1002/anie.201501887.  Measured here under GFN2::

    ethane C-C at 1.54 A      N_FOD 0.000   gap 15.08 eV
    ethane C-C at 3.50 A      N_FOD 1.726   gap  0.24 eV
    ethylene twisted 45 deg   N_FOD 0.075   gap  3.42 eV
    ethylene twisted 60 deg   N_FOD 0.251   gap  2.33 eV
    ethylene twisted 90 deg   N_FOD 2.000   gap  0.00 eV

which is the frontier-gap rule this editor already has, agreeing: 60 degrees
is where that rule fires and where N_FOD has just left zero.  So it is not a
second opinion but the same one with a scale on it, and it is reported as a
*change* across the walk and never as a threshold -- because against ORCA's
own FOD at TPSS/def2-TZVP the same measurement is right to three percent on
ozone (0.426 against 0.439) and four times too large on Ni(CO)4 (0.335 against
0.085).  A fixed cutoff would tell everyone working on metal complexes that
their structures are diradicals.

Driven whole on the user's 57-atom manganese complex at charge +1, closed
shell, the Mn-Br bond walked out over six points: it starts at 2.65 electrons
in half-filled orbitals before anything has moved -- which an absolute rule
would call a triradical -- and ends at 3.42, most of it on the metal.  The
+0.77 is the finding and the 2.65 is what a threshold would have reported
instead.  That run also settles the cost question: the walk took 1097 s on a
loaded machine and pricing the same six geometries again took 35.8 s, three
percent of it, because a constrained optimisation grows with the molecule far
faster than a single point does.

It is asked of GFN2 whatever walked the scan.  Two of the four methods on the
toolbar cannot be asked: GFN-FF has no electrons, and g-xTB takes ``--fod``,
converges, terminates normally and prints nothing at all.  And under g-xTB the
gap rule is deaf, because its thresholds were calibrated on GFN2 and the two
are not the same scale -- measured on an ethylene turned about its double
bond, where 90 degrees is a diradical no closed shell describes::

    twist       0     30     60     90
    GFN2     5.47   4.37   2.32   0.00 eV
    g-xTB   12.17  10.55   7.94   5.12 eV

so a g-xTB profile never trips the rule, under the one method a barrier is
most likely to be quoted from.  Handing the two geometries to GFN2 gives every
walk the same check for two single points -- measured at 0.14 s each on the
sixteen-atom complex whose whole scan is 17 s, and 1.7 s each on a 57-atom
manganese complex whose optimisation steps are 11 to 33 s.
"""

from __future__ import annotations

import pathlib
import shutil
import tempfile
import time

import pytest

from delfin.dashboard import gfn_optimize as gfn


_needs_xtb = pytest.mark.skipif(not shutil.which('xtb'),
                                reason='xtb not installed')
_needs_gxtb = pytest.mark.skipif(gfn.find_gxtb() is None,
                                 reason='the g-xTB build is not installed')

_ETHANE = """8
ethane
C            0.00000000000000        0.00000000000000        0.00000000000000
C            1.53000000000000        0.00000000000000        0.00000000000000
H           -0.36000000000000        1.02000000000000        0.00000000000000
H           -0.36000000000000       -0.51000000000000        0.88400000000000
H           -0.36000000000000       -0.51000000000000       -0.88400000000000
H            1.89000000000000        0.51000000000000        0.88400000000000
H            1.89000000000000        0.51000000000000       -0.88400000000000
H            1.89000000000000       -1.02000000000000        0.00000000000000
"""

#: Another molecule, for the press that must stop being offered.  Told apart
#: by the element column, the way the scan's two ends and the held values are.
_WATER = """3
water
O  0.0000  0.0000  0.0000
H  0.7570  0.5860  0.0000
H -0.7570  0.5860  0.0000
"""


def _an_editor(text=_ETHANE):
    """One structure editor over a coordinate box of its own.

    The real part, driven the way a user drives it: which control is on screen
    is a fact about the widgets and not about what the source says it intends.
    """
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
    state = {}
    part = structure_editor.build(
        ctx, state=state, coords_widget=widgets.Textarea(value=text),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)
    part.submit_manip_toolbar.layout.display = 'flex'
    return part, state


def _shown(widget):
    return (widget.layout.display or '') != 'none'


def _walk_the_ethane(part, state, method='gfn2', steps=6, to=3.50):
    """Ethane's C-C, walked out until the bond is not one any more."""
    part.submit_ff_dd.value = method
    part.submit_scan_how.value = 'hold'
    part.submit_scan_energy.value = 'E'
    part.submit_scan_whole.value = True
    part.submit_scan_steps.value = steps
    part.submit_pick_sync.value = '0,1'
    part.submit_scan_way.value = 'to'
    part.submit_scan_to.value = to
    part.submit_scan_add_btn.click()
    part.submit_scan_run_btn.click()
    began = time.time()
    while state.get('scan_run') and time.time() - began < 600:
        time.sleep(0.05)
    time.sleep(0.2)


# -- what the toolbar says can be done ---------------------------------------

def test_nothing_is_offered_before_there_is_a_walk_to_price():
    """The absence is the statement: there is no profile, so there is no press.

    The visible set of controls is this editor's answer to "what can I do
    now", so a press offering to price a walk that has not happened would be
    the toolbar claiming something it cannot do.
    """
    part, _state = _an_editor()
    assert not _shown(part.submit_scan_price_btn)
    assert part.submit_scan_price_btn.disabled


@_needs_xtb
def test_the_press_arrives_with_the_finished_profile():
    """A walk that has finished is the moment, and the row says so by changing.

    The same shape as the two ends f1be8954 put on the saddle press: what a
    walk made possible arrives on screen when the walk ends, rather than
    waiting to be looked for.
    """
    part, state = _an_editor()
    _walk_the_ethane(part, state)
    assert state.get('scan_walk'), 'the walk kept no geometries'
    assert len(state['scan_walk']['points']) >= 5
    assert _shown(part.submit_scan_price_btn)
    assert not part.submit_scan_price_btn.disabled
    assert 'g-xTB' in part.submit_scan_price_btn.description


@_needs_xtb
def test_a_profile_of_another_molecule_is_not_this_molecules_profile():
    """The walk outlives the structure it was made on, and must not be offered.

    A set of geometries is only a profile *of* something, and geometries
    belonging to another molecule re-priced under this one's name would be a
    second opinion about a structure nobody is looking at.  Told apart by the
    element column, which is what makes two geometries the same molecule --
    the rule the scan's two ends and the thermal anchor already answer to.
    """
    part, state = _an_editor()
    _walk_the_ethane(part, state)
    assert _shown(part.submit_scan_price_btn)

    part.coords_widget.value = _WATER
    state['current_xyz_for_copy'] = {'content': _WATER}
    part._refresh_scan()
    assert not _shown(part.submit_scan_price_btn)
    # And it is kept rather than thrown away, so going back offers it again.
    assert state.get('scan_walk')


@_needs_xtb
@_needs_gxtb
def test_the_best_method_in_the_list_is_offered_no_second_opinion():
    """Walked with g-xTB there is nothing better here, so there is no press.

    A box or a button offering to price a g-xTB profile with g-xTB is a
    question with one wrong answer.
    """
    part, state = _an_editor()
    _walk_the_ethane(part, state, method='gxtb', steps=4)
    assert state.get('scan_walk'), 'the walk kept no geometries'
    assert state['scan_walk']['method'] == 'gxtb'
    assert not _shown(part.submit_scan_price_btn)


@_needs_xtb
def test_the_press_goes_with_the_rest_of_the_scan_under_a_browser_method():
    """A scan is xtb's, and so is pricing one again."""
    part, state = _an_editor()
    _walk_the_ethane(part, state)
    assert _shown(part.submit_scan_price_btn)
    part.submit_ff_dd.value = 'uff'
    assert not _shown(part.submit_scan_price_btn)


# -- the second opinion itself ------------------------------------------------

@_needs_xtb
@_needs_gxtb
def test_the_second_profile_prices_the_walk_and_does_not_move_it():
    """Every point is re-priced where it stood, and the energies change.

    The whole claim of the press is "the same geometries, a better method", so
    the coordinate of every point has to come back identical and the energies
    have to not.  Measured on ethane's C-C walked to 3.50 A.
    """
    part, state = _an_editor()
    _walk_the_ethane(part, state)
    walked = [one[0] for one in state['scan_walk']['points']]

    part.submit_scan_price_btn.click()
    began = time.time()
    while state.get('reprice_run') and time.time() - began < 600:
        time.sleep(0.05)
    time.sleep(0.2)

    again = state.get('scan_repriced')
    assert again, 'nothing was priced'
    assert again['label'] == 'g-xTB'
    assert again['walked_with'] == 'GFN2-xTB'
    assert [one[0] for one in again['points']] == walked
    assert again['points'][0][1] == pytest.approx(0.0)
    # A homolysis is exactly the case the deviations above are about, so the
    # two profiles must not agree.
    was = state['scan_walk']['points'][-1][1]
    assert abs(again['points'][-1][1] - was) > 1.0
    # And the button is a press again rather than a Stop.
    assert part.submit_scan_price_btn.description.startswith('Price')


@_needs_xtb
@_needs_gxtb
def test_the_verdict_says_what_to_do_and_the_tooltip_says_what_it_is():
    """Two different readers, two different moments, two different texts.

    What a second opinion *is* -- a screen and not a barrier, taken on
    structures the first method found -- is true of every press and changes
    nothing about any of them. It belongs on the button, which is what
    somebody deciding whether to press it is looking at.

    What comes back is a decision, and it is different every time: either
    the two methods agree and the walk's own number stands, or they do not
    and the top has to be optimised again. That is the sentence under the
    picture, and putting the unchanging half there as well turned a verdict
    into five paragraphs that were read past.
    """
    part, state = _an_editor()
    _walk_the_ethane(part, state)
    part.submit_scan_price_btn.click()
    began = time.time()
    while state.get('reprice_run') and time.time() - began < 600:
        time.sleep(0.05)
    time.sleep(0.2)

    said = ' '.join(state.get('mol_status_lines') or ())
    # The numbers, and then what they mean for the walk.
    assert 'prices the same' in said
    assert 'kcal/mol' in said
    assert ('optimise the top again' in said
            or "walk's own number stands" in said), said
    # And none of the standing caveats, which changed no decision.
    for essay in ('preprint', 'development version', 'never the geometries',
                  'Three things this is not'):
        assert essay not in said, essay
    # They are on the button instead, read before the press rather than
    # after it.
    tip = part.submit_scan_price_btn.tooltip
    assert 'screen, not a barrier to quote' in tip
    assert 'still the ones GFN2 found' in tip


# -- how much of it was never a closed shell ---------------------------------

def test_a_method_with_no_wavefunction_is_refused_rather_than_left_silent():
    """Both ways of not answering are silent, so both are refused in words.

    GFN-FF prints a runtime warning, no NFOD and no cube, and then says normal
    termination.  g-xTB is worse: it takes the flag, converges, terminates
    normally and prints nothing, so a caller reading only the exit code would
    report a molecule with no static correlation rather than a question that
    was never put.
    """
    for method in ('gfnff', 'gxtb'):
        got = gfn.optimize_with_gfn(_ETHANE, method, optimise=False, fod=True)
        assert not got['ok']
        assert 'closed shell' in got['status']
    assert gfn.can_measure_fod('gfn2')
    assert gfn.can_measure_fod('gfn1')
    assert not gfn.can_measure_fod('gfnff')
    assert not gfn.can_measure_fod('gxtb')


def test_it_is_a_single_point_by_construction():
    """The extra SCF runs at 5000 K, so optimising under it is relaxing on a
    surface nobody wants."""
    got = gfn.optimize_with_gfn(_ETHANE, 'gfn2', optimise=True, fod=True)
    assert not got['ok']
    assert '5000 K' in got['status']


def test_nothing_printed_is_not_the_same_as_nothing_there():
    """A clean molecule prints zero; a method that cannot be asked prints
    nothing, and the two must not come back the same."""
    from delfin.dashboard.gfn_optimize import _read_fod

    assert _read_fod('') is None
    assert _read_fod('nothing to see here') is None
    got = _read_fod('NFOD :     0.0000\n')
    assert got is not None and got['total'] == 0.0


def test_the_per_atom_table_is_read_in_this_editors_atom_numbering():
    """xtb counts from one and the editor counts from zero, and the two are
    one line apart wherever this is reported."""
    from delfin.dashboard.gfn_optimize import _read_fod

    got = _read_fod(
        'NFOD :     2.6031\n'
        '\n'
        ' Loewdin FODpop      n(s)   n(p)   n(d)\n'
        '     1Br    0.0848   0.000  0.083  0.002\n'
        '     8Mn    1.5533   0.001  0.002  1.550\n'
        '\n'
        'Wiberg/Mayer (AO) data.\n'
    )
    assert got['total'] == pytest.approx(2.6031)
    assert got['on'] == [('Br', 0, 0.0848), ('Mn', 7, 1.5533)]


def test_it_is_said_as_a_change_and_only_when_it_moved():
    """No absolute threshold, ever -- the reason is Ni(CO)4.

    Measured against ORCA's FOD at TPSS/def2-TZVP: ozone 0.426 against 0.439,
    three percent on the textbook diradicaloid, and Ni(CO)4 0.335 against
    0.085, four times too large on a closed-shell metal carbonyl.  So there is
    no value that means "bad", and what is reported is the difference between
    two points of one walk on one molecule.
    """
    assert gfn.fod_moved(0.03, 0.05) == ''
    assert gfn.fod_moved(0.335, 0.335) == ''
    said = gfn.fod_moved(0.00, 1.73)
    assert 'from 0.00 to 1.73' in said
    assert 'Read the change, not the number' in said
    assert 'four times high' in said
    # An absolute number never speaks on its own.
    assert gfn.fod_moved(2.60, 2.60) == ''


def test_an_atom_is_named_only_where_that_is_true():
    """A bond broken symmetrically has half on each end, and naming one of
    them says something false about which end is the trouble.

    Measured: the ethane homolysis puts 0.842 of 1.726 on each carbon, and the
    manganese complex puts 1.553 of 2.603 on the metal -- 1.550 of it d.
    """
    even = gfn.fod_moved(0.00, 1.726,
                         [('C', 0, 0.8421), ('C', 1, 0.8421)])
    assert 'most of it on' not in even
    lopsided = gfn.fod_moved(0.03, 2.603,
                             [('Mn', 7, 1.5533), ('Br', 0, 0.0848)])
    assert 'most of it on Mn7' in lopsided


@_needs_xtb
def test_the_gap_and_the_count_agree_about_where_the_trouble_is():
    """The rule this stands beside fires at 60 degrees; N_FOD leaves zero there.

    That agreement is the whole reason the count is worth having: it is not a
    second opinion, it is the same one with a scale on it -- so a walk into a
    region the gap warns about must not come back with nothing to count.
    """
    part, state = _an_editor()
    _walk_the_ethane(part, state)
    said = state.get('scan_depth') or ''
    assert 'frontier gap' in said, 'the gap said nothing about a broken bond'
    assert 'half-filled orbitals' in said, 'nothing was counted'
    # And both are in what the user is shown, not only in the state.
    assert 'half-filled orbitals' in ' '.join(state['mol_status_lines'])


@_needs_xtb
@_needs_gxtb
def test_a_walk_under_a_method_that_cannot_be_asked_is_checked_anyway():
    """It is a probe, so it is GFN2's job whatever walked the scan.

    Under g-xTB this is not a corner case but the main one: the gap rule is
    calibrated on GFN2 and cannot fire there at all -- 5.12 eV at a perfect
    diradical against GFN2's 0.00 -- so without this a g-xTB barrier goes out
    with no check on it of any kind.  And the sentence names who was asked,
    because a number from a method other than the one that walked is a
    different claim.
    """
    part, state = _an_editor()
    _walk_the_ethane(part, state, method='gxtb', steps=4)
    said = state.get('scan_depth') or ''
    assert 'half-filled orbitals' in said
    assert 'GFN2-xTB was asked, on the geometries g-xTB walked' in said
