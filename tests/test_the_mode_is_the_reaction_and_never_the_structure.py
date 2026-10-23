"""What there is to do with a transition state once the editor has found one.

Four ways up were built and nothing at the top.  A search converged, said "one
imaginary mode at -394 cm-1", and that number was the whole of what the editor
could tell anybody about the structure it had just made.  The two questions a
chemist actually has about a saddle both went unanswered.

**What reaction is this.**  The imaginary mode *is* the reaction coordinate:
the atoms it moves are the atoms the reaction moves.  Drawing it is therefore
the single most informative thing there is to show about a saddle, and it costs
one xtb Hessian -- measured, 0.41 s at sixteen atoms.  It is drawn down the
frame channel the trajectory player already owns, because a mode is frames like
any other path, and it is drawn at an amplitude that is measured rather than
chosen: on the sixteen-atom Diels-Alder saddle, whose two forming bonds stand
at 2.315 A, 0.35 A of largest atomic displacement swings both of them between
1.63 A and 3.01 A, and the tightest contact anywhere in that swing is 0.81 of a
bond where two thirds is the editor's floor for "squeezed".  That same
structure first tears at 0.70, so the default is half of where this molecule
breaks -- and it is cut further on a structure whose mode drives two atoms
together sooner.

**Does it connect what I meant.**  One imaginary mode makes a structure *a*
transition state; whether it is *the* one is a different question, and the
standard way to ask it is to push the structure a little way down the mode each
way and relax.  Measured on this machine, on the Diels-Alder saddle:

    displace and relax    1.0 s   ring closed at 1.54 A, -70.6 kcal/mol
                                  two pieces 3.3 A apart, -6.7 kcal/mol
    ! XTB2 IRC          207.0 s   2.83 A and -5.3, 1.54 A and -70.0

Both find the same two ends.  The IRC converges in the valley rather than at
the bottom of it, so it needs those same two relaxations afterwards to name
what it reached -- and 207 s is longer than :func:`saddle.seconds_for` allows
the saddle search itself, on the smallest case anybody runs, on a machine that
on a cluster is the login node.  So the press is the cheap one.  ORCA's IRC
does work over ``ExtOpt`` as well, which is how g-xTB is driven: run on the
same saddle it terminated normally, converged both ways onto the same two
ends, and took 960 s.

And the last part, which is the one that would make all of it worthless: **an
animated frame must never reach the coordinate box**.  A geometry displaced
along an eigenvector is a picture and not a structure anybody chose.  The
animation begins and ends on the structure itself, and the run carries the
frame to put the picture back to if a hand lands in the middle of it -- drawn
on the page, before the grab is reported, because the drag pushes the viewer's
own coordinates to the kernel ten times a second and a round trip is longer
than that.
"""

from __future__ import annotations

import json
import math
import re
import tempfile
from pathlib import Path

import numpy as np
import pytest

from delfin.dashboard import climb, gfn_optimize as gfn
from delfin.dashboard.context import DashboardContext

_needs_xtb = pytest.mark.skipif(
    gfn.find_binary('gfn2') is None, reason='xtb not installed')

#: The converged Diels-Alder transition state, as ORCA's OptTS on XTB2 left it
#: from the path finder's estimate: the two forming bonds at 2.315 A and one
#: imaginary mode at -393.5 cm-1.  Atoms 0-9 are the butadiene, 10-15 the
#: ethene, so the bonds the reaction makes are C0-C10 and C3-C11.
_SADDLE = """16
ORCA OptTS on XTB2 from the estimate, E = -17.812259376910
  C          -1.43478368987847     -0.11555756879252      0.39801753593899
  C          -0.70800000464735      0.96029290831639      0.01351101273233
  C           0.70792543036907      0.96066325147865      0.01382699497627
  C           1.43510098205032     -0.11480687774322      0.39865617430983
  H          -2.50480092399673     -0.04135148267055      0.50582952622170
  H          -1.04673865581999     -1.11397983254932      0.31729931422258
  H          -1.20115753813432      1.91631252668957     -0.10370598026690
  H           1.20063543422368      1.91694042137411     -0.10317121370013
  H           2.50503091371526     -0.04004089983163      0.50694578141932
  H           1.04761410021811     -1.11343199175848      0.31776831289286
  C          -0.67562925889034     -0.26587715375425      2.58016852501963
  C           0.67503524104879     -0.26553006909586      2.58047503090081
  H          -1.22616409779150      0.63483074833459      2.78770928087516
  H          -1.23028498550611     -1.18636553719529      2.63888406495727
  H           1.22501381825161      0.63546069529481      2.78826237108677
  H           1.23013805516997     -1.18573310998389      2.63943734003805
"""

#: Ammonia forced flat, at the N-H distance where D3h is stationary under GFN2.
#: It is the inversion transition state, and it is the case where both ends of
#: the mode are the same molecule -- which the sentences have to say rather
#: than dress up as a reaction.
_FLAT_AMMONIA = """4
planar ammonia, D3h stationary under GFN2
N     0.000000000000     0.000000000000     0.000000000000
H     0.992730000000     0.000000000000     0.000000000000
H    -0.496365000000     0.859729166847     0.000000000000
H    -0.496365000000    -0.859729166847     0.000000000000
"""

#: Methane forced square planar: a saddle of higher order, and the cheapest
#: honest one there is.  Two modes go the wrong way under GFN2, at -3262 and
#: -768 cm-1, and following either of them down gives the same tetrahedral
#: methane both ways -- which is the case where a verdict that named a
#: reaction would be inventing one.
_FLAT_METHANE = """5
methane forced square planar
C   0.000000   0.000000   0.000000
H   1.090000   0.000000   0.000000
H  -1.090000   0.000000   0.000000
H   0.000000   1.090000   0.000000
H   0.000000  -1.090000   0.000000
"""

_ETHANE = """8
ethane
C  0.000  0.000  0.000
C  1.520  0.000  0.000
H -0.380  1.020  0.000
H -0.380 -0.510  0.880
H -0.380 -0.510 -0.880
H  1.900  1.020  0.000
H  1.900 -0.510  0.880
H  1.900 -0.510 -0.880
"""

#: Words that name a role or a kind of chemistry.  None of them may appear in
#: anything the descent says: every sort of system is computed in this editor,
#: and a sentence built for one of them is wrong about the rest.  Which end of
#: a mode is the reactant is not a question about the saddle at all.
_ASSUMES_CHEMISTRY = (
    'reactant', 'product', 'substrate', 'nucleophile', 'electrophile',
    'leaving group', 'adduct', 'diene', 'ligand', 'catalyst', 'educt',
)


def _rows(text):
    return np.asarray(climb._elements(text)['angstrom'], dtype=float)


def _pieces(text):
    return climb.pieces_in(text)


# ---------------------------------------------------------------------------
# the picture
# ---------------------------------------------------------------------------


@_needs_xtb
def test_the_mode_comes_back_in_the_metric_a_picture_is_drawn_in():
    """Cartesian and scaled, not mass-weighted.

    xtb's Hessian is read out of the ``hessian`` file rather than off the
    printout because the file carries the mode *shapes*, and those come out in
    the mass-weighted metric the eigenvalues live in.  A picture is not drawn
    in that metric: left as they are, a hydrogen moves a third of what it
    really does and a proton transfer looks like the heavy atoms doing the
    work.  Dividing by the root mass puts it back, and normalising so that the
    furthest-moving atom moves exactly one makes an amplitude in Angstrom mean
    Angstrom.

    Measured here on the Diels-Alder saddle: the softest mode is -393.5 cm-1,
    and in it the two ethene carbons move 1.00 while the four hydrogens on the
    ethene move 0.04 to 0.07 -- which is the reaction and not a wag.
    """
    got = climb.modes_of(_SADDLE, 'gfn2')
    assert got is not None, 'no xtb to take a Hessian with'
    assert float(got['cm'][0]) == pytest.approx(-393.5, abs=2.0)
    assert climb.imaginary_among(got['cm']) == [0]

    way = got['ways'][:, 0].reshape(-1, 3)
    moved = np.linalg.norm(way, axis=1)
    assert float(moved.max()) == pytest.approx(1.0)
    # The two carbons that form the bonds are what the mode moves.
    assert moved[10] > 0.9 and moved[11] > 0.9
    assert moved[0] > 0.9 and moved[3] > 0.9
    # And the ethene's own hydrogens barely move at all.
    assert max(float(moved[i]) for i in (12, 13, 14, 15)) < 0.2


def test_the_animation_begins_and_ends_on_the_structure_exactly():
    """The one property that makes it safe to draw at all.

    Every frame between the ends is a geometry nobody chose and nothing
    computed.  So the picture has to come back: the last frame is a whole
    number of swings, and the sine of that is nought, so the animation ends on
    the structure to the last decimal it is written with.  Whatever else goes
    wrong, a run that plays out leaves the viewer showing what the coordinate
    box holds.

    No xtb here -- this is arithmetic on a geometry and a direction.
    """
    here = _rows(_ETHANE)
    way = np.zeros_like(here)
    way[0] = (1.0, 0.0, 0.0)
    way[1] = (-1.0, 0.0, 0.0)

    frames = climb.mode_pictures(here, way, amplitude=0.3)
    assert len(frames) == climb.MODE_PICTURES * climb.MODE_SWINGS + 1

    home = [round(float(one), 4) for one in here.reshape(-1)]
    assert frames[0] == home
    assert frames[-1] == home
    # And it really does move in between, or it is not an animation.
    quarter = np.asarray(frames[climb.MODE_PICTURES // 4]).reshape(-1, 3)
    assert math.dist(quarter[0], quarter[1]) == pytest.approx(
        math.dist(here[0], here[1]) - 0.6, abs=1e-6)


def test_the_amplitude_is_cut_where_the_picture_would_tear():
    """Big enough to read, never big enough to break the molecule.

    "Tear" is the floor a drag is already held to: two atoms inside two thirds
    of the bond they would make are not on any path at any temperature.  A
    changed bond graph is deliberately *not* the test -- the whole point of the
    picture is that bonds appear and disappear in it, and on the Diels-Alder
    saddle the graph changes at an amplitude of 0.2 A while nothing is tight
    until 0.70.

    Driven on two hydrogens 1.0 A apart, with a mode along the line between
    them: at the default amplitude one end of the swing has them 0.30 A apart,
    which is half the floor, so the amplitude is cut until it does not.

    Both ends are asked, and there is no such thing as a mode that only opens
    things up: a mode is symmetric, so the direction that pulls two atoms
    apart on the way out pushes them together on the way back.  Which is why
    a structure keeps the whole amplitude only when *neither* end of the swing
    is tight -- the same two hydrogens 1.4 A apart do.
    """
    tight = "2\ntwo hydrogens\nH 0.0 0.0 0.0\nH 1.0 0.0 0.0\n"
    here = _rows(tight)
    along = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])

    cut = climb.amplitude_that_fits(['H', 'H'], here, along)
    assert cut < climb.MODE_AMPLITUDE, 'a torn picture was allowed'
    for sign in (1.0, -1.0):
        moved = climb.displaced_along(here, along, sign * cut)
        text = climb.xyz_document(['H', 'H'], moved, 'x')
        assert gfn.closest_contact(text)[0] >= gfn._TOO_CLOSE

    roomy = _rows("2\ntwo hydrogens\nH 0.0 0.0 0.0\nH 1.4 0.0 0.0\n")
    assert climb.amplitude_that_fits(['H', 'H'], roomy, along) == \
        climb.MODE_AMPLITUDE


@_needs_xtb
def test_the_diels_alder_mode_is_drawn_with_the_bond_made_and_gone():
    """The measurement the default amplitude was chosen from.

    Not "it moves": at 0.35 A of largest atomic displacement the two forming
    bonds of this saddle swing between 1.63 A -- a formed C-C bond -- and
    3.01 A, which is no bond at all, and nothing in the molecule comes closer
    than 0.81 of a bond anywhere in the swing.  That is a picture nobody has to
    have explained to them, and it is inside the floor by a wide margin.
    """
    got = climb.modes_of(_SADDLE, 'gfn2')
    assert got is not None
    here = np.asarray(got['angstrom'])
    way = got['ways'][:, 0].reshape(-1, 3)
    assert climb.amplitude_that_fits(got['symbols'], here, way) == \
        climb.MODE_AMPLITUDE

    frames = climb.mode_pictures(here, way, amplitude=climb.MODE_AMPLITUDE)
    reach = [math.dist(np.asarray(one).reshape(-1, 3)[0],
                       np.asarray(one).reshape(-1, 3)[10]) for one in frames]
    assert min(reach) == pytest.approx(1.63, abs=0.05)
    assert max(reach) == pytest.approx(3.01, abs=0.05)
    tightest = min(
        gfn.closest_contact(climb.xyz_document(
            got['symbols'], np.asarray(one).reshape(-1, 3), 'x'))[0]
        for one in frames)
    assert tightest > gfn._TOO_CLOSE
    assert tightest == pytest.approx(0.81, abs=0.03)


# ---------------------------------------------------------------------------
# where the mode leads
# ---------------------------------------------------------------------------


@_needs_xtb
def test_a_diels_alder_saddle_comes_down_to_the_two_ends_it_should():
    """The known case, which is what makes the rest of it believable.

    A Diels-Alder transition state has to come down to the diene and the
    dienophile on one side and to the ring on the other.  Measured here, from
    the converged saddle, at 0.3 A down the imaginary mode each way and
    relaxed under GFN2: one direction gives **one piece** with six C-C bonds
    and the two forming bonds closed at 1.54 A, 70.6 kcal/mol below the
    saddle; the other gives **two separate pieces** with four C-C bonds and
    the fragments 3.3 A apart, 6.7 below.  The whole thing costs about a
    second at this size.

    The step size does not decide the answer, which was checked rather than
    assumed: 0.1, 0.2, 0.3 and 0.5 A all give the same two ends to two
    decimals.
    """
    got = climb.follow_the_mode_down(_SADDLE, 'gfn2')

    assert got['ok'], got['status']
    assert got['order'] == 1
    ends = [end for end in got['ends'] if end['ok']]
    assert len(ends) == 2, 'both directions have to relax to something'

    closed = [end for end in ends if end['pieces'] == 1]
    apart = [end for end in ends if end['pieces'] == 2]
    assert len(closed) == 1 and len(apart) == 1, \
        'one end is the ring and the other is the two fragments'

    ring = _rows(closed[0]['xyz'])
    assert math.dist(ring[0], ring[10]) == pytest.approx(1.54, abs=0.03)
    assert math.dist(ring[3], ring[11]) == pytest.approx(1.54, abs=0.03)
    assert closed[0]['kcal'] == pytest.approx(-70.6, abs=1.5)
    # The two bonds it has that the saddle did not are exactly those two.
    assert {tuple(one) for one in closed[0]['made']} == {(0, 10), (3, 11)}

    loose = _rows(apart[0]['xyz'])
    assert math.dist(loose[0], loose[10]) > 3.0
    assert apart[0]['kcal'] == pytest.approx(-6.7, abs=1.5)
    assert not apart[0]['made'] and not apart[0]['broke']

    said = ' '.join(got['lines'])
    assert 'do not have the same bonds' in said
    assert 'C0-C10' in said and 'C3-C11' in said


@_needs_xtb
def test_both_ways_to_the_same_molecule_is_said_and_not_dressed_up():
    """Planar ammonia, where the mode is an inversion and not a reaction.

    Both ends of the umbrella mode are the same molecule: the same bonds, the
    same energy, and 0.32 A RMSD apart once the best rotation that is not a
    reflection has been taken out.  A verdict that reported "a reaction was
    found" here would be inventing one, and one that reported the two ends as
    different molecules would be wrong the other way.  What it has to do is say
    what it measured.
    """
    got = climb.follow_the_mode_down(_FLAT_AMMONIA, 'gfn2')

    assert got['ok'], got['status']
    assert got['order'] == 1
    said = ' '.join(got['lines'])
    assert 'same bonds and the same energy' in said
    assert 'two arrangements of one structure' in said
    ends = [end for end in got['ends'] if end['ok']]
    assert ends[0]['kcal'] == pytest.approx(ends[1]['kcal'], abs=0.05)
    assert climb.turned_onto(ends[0]['there'], ends[1]['there']) > 0.1


@_needs_xtb
def test_a_saddle_of_higher_order_says_so_before_it_says_anything_else():
    """Square-planar methane: two modes going the wrong way, not one.

    The two ends of the deepest of them are worth reporting -- they are where
    that mode leads -- but a structure that is a maximum in two directions is
    not a point any single step passes through, and a report that started with
    "one way it relaxed to..." would read as a reaction.  So the count is said
    first and the ends are named as what they are: the two ways along *one* of
    the modes.

    Measured under GFN2: -3262 and -768 cm-1, and both directions of the
    deepest come back to the same tetrahedral methane, 124.0 kcal/mol below
    and 1.12 A RMSD apart from each other.
    """
    got = climb.follow_the_mode_down(_FLAT_METHANE, 'gfn2')

    assert got['ok'], got['status']
    assert got['order'] == 2
    assert got['lines'][0].startswith('This structure has two modes going the '
                                      'wrong way')
    assert 'not a transition state' in got['lines'][0]
    assert 'two ways along one of them' in got['lines'][0]
    # And it does not call a pair of equivalent minima a reaction.
    assert 'two arrangements of one structure' in ' '.join(got['lines'])
    ends = [end for end in got['ends'] if end['ok']]
    assert ends[0]['kcal'] == pytest.approx(-124.0, abs=2.0)
    assert ends[0]['pieces'] == 1 and ends[1]['pieces'] == 1


@_needs_xtb
def test_a_minimum_has_no_reaction_coordinate_and_says_so():
    """Pressed on something that is not a saddle at all.

    There is no mode going the wrong way, so both directions would come
    straight back and there is nothing to report.  Said rather than run: two
    relaxations that both return the structure they were given look exactly
    like a reaction that connects a thing to itself.

    On a real minimum, which is not the same thing as a relaxed structure and
    was worth finding out: the hand-typed ethane in this file is eclipsed, and
    relaxing it under GFN2 keeps the symmetry it was typed with and lands on
    the *torsional* transition state -- one mode going the wrong way at
    -292 cm-1.  Water is a minimum whichever way it is written down.
    """
    water = "3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"
    settled = gfn.optimize_with_gfn(water, 'gfn2')
    assert settled['ok'], settled['status']
    got = climb.follow_the_mode_down(settled['xyz'], 'gfn2')

    assert not got['ok']
    assert got['order'] == 0
    assert 'no mode' in got['status'].lower()
    assert 'minimum' in got['status']
    assert not got['ends']


@_needs_xtb
def test_the_verdict_names_no_kind_of_chemistry():
    """Universality, checked on the sentences themselves rather than intended.

    This editor is used on organics, on metal complexes, on radicals, on
    clusters -- and on structures that are none of those.  So nothing the
    descent says may name a role: there is no reactant and no product, because
    which end of a mode is which is not a question the saddle answers, and
    there is no reaction type because the saddle does not know one.  What is
    said instead is what was measured: how many separate pieces, which bonds,
    how many kcal/mol.
    """
    for structure in (_SADDLE, _FLAT_AMMONIA, _FLAT_METHANE):
        got = climb.follow_the_mode_down(structure, 'gfn2')
        said = (' '.join(got['lines']) + ' ' + str(got['status'])).lower()
        for word in _ASSUMES_CHEMISTRY:
            assert word not in said, f'{word!r} in {said!r}'


def test_the_sentences_themselves_name_no_kind_of_chemistry():
    """The same rule, read off the source, so a branch that never ran here is
    covered too.

    The measured test above only sees the sentences the two cases it runs
    produce.  A branch for a relaxation that failed, or for two ends that are
    the same place, is not exercised by either -- and that is exactly where a
    hurried sentence about a reactant would survive.
    """
    source = Path(climb.__file__).read_text(encoding='utf-8')
    block = source.split('def _end_said(')[1]
    block += source.split('def follow_the_mode_down(')[1]
    quoted = ' '.join(re.findall(r"'([^']*)'", block)).lower()
    for word in _ASSUMES_CHEMISTRY:
        assert word not in quoted, word


# ---------------------------------------------------------------------------
# the editor
# ---------------------------------------------------------------------------


@pytest.fixture
def editor(tmp_path):
    """One structure editor over a coordinate box of its own."""
    pytest.importorskip('ipywidgets')
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor

    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / 'calc',
                           archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path / 'office')
    ctx.run_js = lambda _script: None
    state = {}

    def build(text=_SADDLE):
        box = widgets.Textarea(value=text)
        part = structure_editor.build(
            ctx, state=state, coords_widget=box, viewer_height=560,
            schedule_ui_update=lambda func, *a, **k: func(*a, **k),
            update_view=lambda *a, **k: None,
            get_smiles_charge=lambda *a, **k: None)
        state['current_xyz_for_copy'] = {'content': text}
        part._set_manip_toolbar_enabled(True)
        return part

    return build


def _visible(widget):
    return str(getattr(widget.layout, 'display', '') or '') != 'none'


def _said(part):
    return ' '.join(re.sub('<[^>]+>', ' ', part.mol_status.value or '').split())


def _one_imaginary(deepest=-393.5):
    """What a search hands over when it has found a transition state."""
    return {'count': 1, 'modes': [deepest], 'real': [119.3, 173.4]}


def _fake_modes(text, deepest=-393.5, others=(119.3, 173.4)):
    """Modes for a structure without paying xtb for them.

    What the two presses need is a frequency list and one direction per mode;
    where those came from is xtb's business and is measured above.  Made up
    here so that the editor's own rules -- what appears, what is drawn, what
    never reaches the box -- can be driven without a chemistry package.
    """
    here = _rows(text)
    ways = np.zeros((here.size, 1 + len(others)))
    ways[0, 0] = 1.0
    ways[3, 0] = -1.0
    for n in range(len(others)):
        ways[1, n + 1] = 1.0
    return {'cm': np.array([deepest, *others]), 'ways': ways,
            'symbols': [line.split()[0] for line in gfn.atom_lines(text)],
            'angstrom': here}


def test_the_two_presses_are_absent_until_there_is_a_saddle(editor):
    """An absence is a statement, so it has to be the true one.

    The visible controls are the editor's answer to "what can I do now".  A
    press that draws the imaginary mode of a structure that has none would be
    the row saying something untrue, so neither of these is on screen until a
    search has reported a mode going the wrong way -- and both go the moment
    the box stops holding the structure it was reported on.
    """
    part = editor()
    part.submit_ff_dd.value = 'gfn2'
    assert not _visible(part.submit_mode_btn)
    assert not _visible(part.submit_ends_btn)

    part._note_the_saddle(part.coords_widget.value, _one_imaginary())
    part._refresh_saddle_controls()
    assert _visible(part.submit_mode_btn)
    assert _visible(part.submit_ends_btn)
    # One mode is not a choice, so the box that picks between them is not on
    # screen -- the same rule the two boxes beside the press already follow.
    assert not _visible(part.submit_mode_dd)

    # And a different structure takes them away without anything having to
    # remember to.  The saddle is kept against the geometry, not against a
    # flag, so every way of changing what is being worked on clears it at
    # once.  Both places the editor reads the structure from are moved, which
    # is what a host does when it loads one.
    part.coords_widget.value = _ETHANE
    part.state['current_xyz_for_copy'] = {'content': _ETHANE}
    part._refresh_saddle_controls()
    assert not _visible(part.submit_mode_btn)
    assert not _visible(part.submit_ends_btn)


def test_the_mode_box_appears_only_where_there_is_more_than_one(editor):
    """Which is a saddle of higher order, and the one place it earns its keep.

    The verdict on a second-order saddle already tells the user to move the
    structure along the *second* of its modes and climb again.  Until this
    there was no way to look at the second of anything.
    """
    part = editor()
    part.submit_ff_dd.value = 'gfn2'
    part._note_the_saddle(part.coords_widget.value,
                          {'count': 2, 'modes': [-52.0, -26.0], 'real': []})
    part._refresh_saddle_controls()

    assert _visible(part.submit_mode_dd)
    assert [value for _label, value in part.submit_mode_dd.options] == [0, 1]
    assert part.submit_mode_dd.options[1][0] == '-26 cm-1'
    part.submit_mode_dd.value = 1
    assert part._which_mode() == 1


def test_the_animation_is_drawn_and_never_kept(editor):
    """The coordinate box across a whole animation, and what the run carries.

    The editor's rule is that what is in the viewer is what is being worked on.
    A geometry displaced along an eigenvector breaks that: nobody chose it and
    no method computed it.  So the whole animation goes down the frame channel
    -- which is drawn, not kept -- and the box is not written at all.

    Three things are checked because three different mistakes would each be
    invisible on its own: the box is byte-for-byte what it was; the frames
    really do move, so this is not a picture of nothing; and the run carries
    the geometry to return to, which is the box's own and is also the last
    frame, so the animation running out and the animation being cut short land
    in the same place.
    """
    part = editor()
    part.submit_ff_dd.value = 'gfn2'
    before = part.coords_widget.value
    part._note_the_saddle(before, _one_imaginary())
    part._refresh_saddle_controls()

    part._draw_the_mode(_fake_modes(before), 0)

    assert part.coords_widget.value == before, 'the animation reached the box'
    payload = json.loads(part.submit_gfn_frame.value)
    assert payload['final'] == 1, 'the player would throw the queue away'
    assert payload['pace'] == climb.MODE_PACE_MS
    assert payload['from'] == 0
    frames = payload['frames']
    assert len(frames) == climb.MODE_PICTURES * climb.MODE_SWINGS + 1

    home = [round(float(one), 4) for one in _rows(before).reshape(-1)]
    assert payload['home'] == home
    assert frames[0] == home and frames[-1] == home
    # And it moved: the extremes of the swing are a long way from the box.
    furthest = max(max(abs(a - b) for a, b in zip(one, home)) for one in frames)
    assert furthest > 0.1, 'nothing was drawn away from the structure'

    # Nothing was recorded as a change either, or Undo would offer to take a
    # picture back.
    assert not any('mode' in str(step.get('what', '')).lower()
                   for step in (part.state.get('history') or ()))


def test_a_hand_during_the_animation_hands_the_kernel_nothing(editor):
    """The grab, from the kernel's side.

    The page reports which frame the picture stands on when a hand lands, and
    for a real run that number is how the kernel cuts the path -- the frame on
    screen becomes the structure.  An animation has no path to cut and must
    not gain one: it never puts a path down, so the message names a frame of a
    walk the kernel is not holding and nothing is written.

    The other half of the same rule is on the page and is measured in a
    browser below: the picture is put back to the run's ``home`` frame before
    the grab is reported at all, so the coordinates the drag pushes back are
    the ones the box already holds.
    """
    part = editor()
    part.submit_ff_dd.value = 'gfn2'
    before = part.coords_widget.value
    part._note_the_saddle(before, _one_imaginary())
    part._draw_the_mode(_fake_modes(before), 0)
    run = json.loads(part.submit_gfn_frame.value)['run']

    # Halfway through the swing, which is the worst moment there is: the
    # structure on screen is as far from the box as it ever gets.
    part.submit_cmd_sync.value = f'gfngrab:1:{climb.MODE_PICTURES // 4},{run}'

    assert part.coords_widget.value == before, \
        'a displaced frame was handed to the kernel'
    assert part.state.get('gfn_stopped_path') is None, \
        'the animation put a path down for a Stop to cut'

    part.submit_cmd_sync.value = 'gfnfree:2:'
    assert part.coords_widget.value == before


def test_every_search_that_reaches_a_saddle_offers_what_to_do_with_it(editor):
    """All four routes up, and the same two presses at the top of each.

    There are four ways to a transition state in this editor and they end in
    three different places in the source.  A capability that arrived on one of
    them and not the others would be the worst kind of inconsistency: the same
    structure, reached two ways, offering different things.

    Read off the source rather than driven, because driving all four means
    ORCA four times; what each of them *reaches* is measured next door.
    """
    from delfin.dashboard import structure_editor

    source = Path(structure_editor.__file__).read_text(encoding='utf-8')
    for where in ('def _saddle_from_here():',
                  'def _path_then_orca(ends):',
                  'def _climb_now():'):
        block = source.split(where)[1].split('\n    def _')[0]
        assert '_note_the_saddle(' in block, where
        assert '_the_mode_is_offered(' in block, where

    # And it is written down after the geometry, because what it records is
    # checked against the structure the box holds: said first, it would name a
    # saddle nobody is standing on yet and the presses would stay hidden.
    for where in ('def _saddle_from_here():', 'def _climb_now():'):
        block = source.split(where)[1].split('\n    def _')[0]
        assert block.index('_write_coords(kept') < \
            block.index('_note_the_saddle(')


def test_the_verdict_says_the_presses_have_arrived(editor):
    """A control appearing is a statement, and it is one somebody has to
    notice.

    This row has already had the other failure once: a capability arrived, the
    controls it needed were folded into an existing box, and nothing on screen
    said so -- "wo sind die hin".  So the sentence at the end of a saddle
    search names the two presses, once, at the moment they become possible.
    """
    part = editor()
    part.submit_ff_dd.value = 'gfn2'

    said = ' '.join(part._the_mode_is_offered(_one_imaginary()))
    assert 'Show the mode' in said and 'Follow it down' in said
    assert 'reaction coordinate' in said

    # Two modes: the same two presses, and the box that says which.
    said = ' '.join(part._the_mode_is_offered({'count': 2, 'modes': [-52, -26]}))
    assert 'Show the mode' in said and 'the box beside them' in said

    # A minimum offers nothing, and the sentence says nothing.
    assert part._the_mode_is_offered({'count': 0, 'modes': []}) == []

    # And under a method whose modes xtb cannot be asked for, the presses are
    # not on screen -- so the sentence names the methods that would put them
    # there rather than two buttons that are not there.
    part.submit_ff_dd.value = 'gxtb'
    said = ' '.join(part._the_mode_is_offered(_one_imaginary()))
    assert 'Show the mode' not in said
    assert 'GFN2' in said
    part._note_the_saddle(part.coords_widget.value, _one_imaginary())
    part._refresh_saddle_controls()
    assert not _visible(part.submit_mode_btn)


def test_a_method_that_cannot_be_asked_for_modes_says_so_rather_than_running(
        editor):
    """Pressed under g-xTB, which has no ``--hess`` of its own.

    ORCA can drive g-xTB to a saddle through ExtOpt, so a g-xTB transition
    state is a thing this editor can produce; xtb cannot be asked what its
    modes look like, because g-xTB is a build of its own and an ordinary xtb
    accepts ``--gxtb`` and silently runs GFN2.  Refused with a sentence that
    says what to do, rather than answered with GFN2's modes under g-xTB's name.
    """
    part = editor()
    part.submit_ff_dd.value = 'gxtb'
    part._note_the_saddle(part.coords_widget.value, _one_imaginary())
    part.on_submit_show_mode()

    assert 'GFN2' in _said(part)
    assert part.submit_gfn_frame.value == '', 'it drew something anyway'


# ---------------------------------------------------------------------------
# and the page
# ---------------------------------------------------------------------------


_PAGE = """<!doctype html><html><body>
<div class="__SCOPE__">
  <div class="submit-gfn-frame"><textarea></textarea></div>
  <div class="submit-cmd-sync"><input></div>
  <button class="submit-optimize-switch">Optimize</button>
</div>
<script>
window.__drawn = [];
window.__pushed = [];
window._submitManipStateByScope = {};
window.__delfinSubmitManip = {
  setPositions: function(scope, out, held) {
    window.__drawn.push(out.slice());
    return true;
  },
  pushXyz: function(scope, why) {
    window.__pushed.push(window.__drawn[window.__drawn.length - 1]);
    return true;
  }
};
</script>
</body></html>"""


def _player_script():
    """The whole self-executing player block, as the page receives it."""
    from delfin.dashboard import tab_submit

    folder = Path(tempfile.mkdtemp())
    for name in ('calc', 'archive', 'office'):
        (folder / name).mkdir()
    ctx = DashboardContext(calc_dir=folder / 'calc',
                           archive_dir=folder / 'archive',
                           office_dir=folder / 'office')
    ctx.run_js = lambda _script: None
    tab_submit.create_tab(ctx)
    startup = '\n'.join(ctx.init_js_parts)
    at = startup.index('window.__delfinGfnPlay')
    return startup[startup.rindex('(function(){', 0, at):]


@pytest.fixture(scope='module')
def browser():
    pytest.importorskip('ipywidgets')
    playwright = pytest.importorskip('playwright.sync_api')
    with playwright.sync_playwright() as p:
        try:
            engine = p.chromium.launch()
        except Exception as exc:                        # no browser installed
            pytest.skip(f'chromium not available: {exc}')
        yield engine
        engine.close()


def test_a_hand_on_an_animation_puts_the_picture_back_before_it_pushes(browser):
    """The half of "never kept" that lives on the page, driven in Chromium.

    A drag pushes the *viewer's* coordinates back to the kernel about ten
    times a second as the geometry the next round starts from.  So it is not
    enough for the kernel to refuse a displaced frame: the browser must not be
    standing on one when the hand lands, and asking the kernel to put it back
    is a round trip longer than the interval between two pushes.

    The run carries the frame to return to, and the player draws it inside the
    same animation frame that notices the grab -- before ``gfngrab`` is sent
    and before any push.  Measured here: the structure is 5.0 away from home
    at the moment the hand arrives, the last thing drawn is home exactly, and
    what the drag pushes is home.
    """
    script = _player_script()
    scope = script.split('var scope=')[1].split(';')[0].strip().strip('"')
    home = [0.0] * 9
    path = [[float(one)] * 9 for one in range(1, 11)] + [home]

    page = browser.new_page()
    try:
        page.set_content(_PAGE.replace('__SCOPE__', scope))
        page.evaluate(script)
        page.evaluate(
            """([sel, text]) => { document.querySelector(sel).value = text; }""",
            ['.submit-gfn-frame textarea',
             json.dumps({'run': 1, 'from': 0, 'frames': path, 'final': 1,
                         'home': home, 'pace': 120})])
        page.wait_for_timeout(700)
        away = page.evaluate('window.__drawn[window.__drawn.length-1]')
        assert away[0] > 1.0, 'nothing was drawn away from home'

        # A hand lands, the way the manipulator says one has.
        page.evaluate(
            """(s) => { window._submitManipStateByScope[s] =
                 {drag: {kind: 'translate', movedEnough: true, targets: []}}; }""",
            scope)
        page.wait_for_timeout(400)

        drawn = page.evaluate('window.__drawn')
        assert drawn[-1] == home, \
            f'the picture was left displaced at {drawn[-1][0]}'
        assert page.evaluate('window.__drawn.length') == len(drawn), \
            'it went on drawing under the hand'
        sent = page.evaluate(
            """(s) => document.querySelector('.submit-cmd-sync input').value""",
            scope)
        assert sent.startswith('gfngrab:'), sent
    finally:
        page.close()


def test_a_climb_that_has_been_superseded_says_nothing():
    """Reported from a real session: "Was soll NEB-TS hier jetzt machen?"

    The journal answers it. A hand-climb runs in rounds and takes a fresh run
    number for each one, so between two rounds its own switch reads as free
    and another press may start something else -- which is what happened::

         969.2s  run 47 claimed by climb
        1001.7s  run 48 claimed by climb
        1005.0s  press To the saddle
        1005.0s  run 49 claimed by band
                 "relaxing a band of 8 images between the two ends..."
        1010.8s  "The climb stopped at the frame you were looking at,
                  0 step(s) in (41.6 s). Press Climb to TS again"

    The band was started and was doing exactly what it said. Then the finished
    round of the climb wrote its goodbye over the band's line, and from the
    user's seat NEB-TS had announced itself and immediately reported a climb
    stopping. It is not the band that was wrong; it is the row.

    The rule is the one the frame channel has always followed: what has been
    superseded says nothing.
    """
    from editor_source import EDITOR_SOURCE

    stopped = EDITOR_SOURCE.split(
        'The climb stopped at the frame you were looking ', 1)[1].split(
        '\n                    return', 1)[0]
    assert '_frame_run_is_current(run)' in stopped, stopped
    assert '_set_mol_status(*walked_said, said)' in stopped
    # But being superseded is not enough on its own, and getting that wrong
    # cost a red CI: a Stop moves the run number too -- deliberately, so the
    # page drops what it was playing -- and a climb takes a fresh number for
    # every round, so it supersedes itself constantly. Silencing on that
    # silenced every Stop the user pressed. It is the walker that has to be a
    # different one.
    for allowed in ('press', 'abandoned', 'climb'):
        assert repr(allowed) in stopped, allowed


# ---------------------------------------------------------------------------
# The two structures the press reached, on screen
# ---------------------------------------------------------------------------
#
# Follow it down finds two minima, describes them in sentences, and used to
# throw the geometries away: two structures came out, one box holds one, and
# the one being worked on is the saddle.  What the sentences cannot say is
# what the two ends look like, and that is the reason anybody follows a mode
# down -- "bei Follow it down will ich natuerlich auch beide enden im viewer
# sehen koennen", said twice.
#
# The saddle is the first entry in the box rather than something to find
# again with Undo, because it is what the user was working on and the two
# ends are a fact about it.


def _put(part, text):
    """Put a structure where the editor reads it from.

    Both places, the way the tab does it: the box, and the key the Submit tab
    keeps the current structure in and every force field here reads
    (``_current_xyz``).  A test that wrote only one of them would be driving
    an editor no host produces.
    """
    part.coords_widget.value = text
    part.state['current_xyz_for_copy'] = {'content': text}


def _two_ends(part, saddle=None):
    """What a Follow it down leaves behind when it has reached both ends."""
    saddle = saddle or part.coords_widget.value
    rows = [line for line in saddle.splitlines()[2:] if line.strip()]

    def moved(by):
        out = []
        for line in rows:
            bits = line.split()
            out.append(f'{bits[0]} {float(bits[1]) + by:.6f} '
                       f'{bits[2]} {bits[3]}')
        return '{}\nan end\n'.format(len(out)) + '\n'.join(out)

    part.state['down_ends'] = [
        {'which': 'One way', 'xyz': moved(-0.4), 'kcal': -70.6},
        {'which': 'The other way', 'xyz': moved(0.4), 'kcal': -6.7},
    ]
    part.state['down_saddle'] = saddle
    part._refresh_the_down_ends()
    return part.state['down_ends']


def test_the_ends_box_is_absent_until_a_press_has_reached_them(editor):
    part = editor(_SADDLE)
    assert not _visible(part.submit_down_dd), 'nothing has been followed down'
    _two_ends(part)
    assert _visible(part.submit_down_dd)
    assert [value for _label, value in part.submit_down_dd.options] == [
        'top', '0', '1']


def test_either_end_goes_on_screen_and_the_saddle_comes_back(editor):
    """Which is the whole of the ask. The numbers say what the two ends are;
    only the picture says what they look like."""
    part = editor(_SADDLE)
    ends = _two_ends(part)
    saddle = part.coords_widget.value

    part.submit_down_dd.value = '0'
    assert (part._geometry_key(part.coords_widget.value)
            == part._geometry_key(ends[0]['xyz']))
    assert 'One way' in _said(part) and '-70.6' in _said(part)

    part.submit_down_dd.value = '1'
    assert (part._geometry_key(part.coords_widget.value)
            == part._geometry_key(ends[1]['xyz']))

    part.submit_down_dd.value = 'top'
    assert part._geometry_key(part.coords_widget.value) == part._geometry_key(
        saddle), 'the saddle is one press away, not a count of Undos'


def test_the_box_reads_whichever_of_the_three_is_on_screen(editor):
    """It follows the structure rather than leading it, so the row cannot say
    "one way" over a picture of the saddle."""
    part = editor(_SADDLE)
    ends = _two_ends(part)

    _put(part, ends[1]['xyz'])
    part._refresh_the_down_ends()
    assert part.submit_down_dd.value == '1'

    _put(part, part.state['down_saddle'])
    part._refresh_the_down_ends()
    assert part.submit_down_dd.value == 'top'


def test_looking_at_both_ends_is_one_press_of_undo(editor):
    """Nothing is lost by going there. A run of looking is one step in the
    history, the way a sweep of an arrow key is one rather than two hundred."""
    part = editor(_SADDLE)
    _two_ends(part)
    saddle = part.coords_widget.value
    depth = len(part.state.get('history') or ())

    part.submit_down_dd.value = '0'
    part.submit_down_dd.value = '1'
    part.submit_down_dd.value = '0'
    assert len(part.state.get('history') or ()) == depth + 1

    part.on_submit_manip_undo()
    assert part._geometry_key(part.coords_widget.value) == part._geometry_key(
        saddle)


def test_the_ends_belong_to_the_saddle_and_not_to_the_molecule(editor):
    """The walk and the scan's two ends answer to the element column: they
    describe a molecule and outlive any one geometry of it. These describe one
    stationary point. Optimise an end or drag it, and there is no longer a
    saddle on screen for two ends to be the ends of."""
    part = editor(_SADDLE)
    _two_ends(part)
    assert _visible(part.submit_down_dd)

    rows = [line for line in _SADDLE.splitlines()[2:] if line.strip()]
    bits = rows[0].split()
    rows[0] = f'{bits[0]} {float(bits[1]) + 0.9:.6f} {bits[2]} {bits[3]}'
    _put(part, '{}\nEdited in DELFIN viewer\n'.format(
        len(rows)) + '\n'.join(rows))
    part._refresh_the_down_ends()
    assert not _visible(part.submit_down_dd)
