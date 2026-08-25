"""A nudged elastic band, for when the fast answer is not believed.

ORCA's ``! NEB-TS``: a chain of images relaxed onto the way between two ends
with energy-weighted springs, the highest turned into a climbing image, and
that one handed to the transition-state optimiser.  Asgeirsson and Jonsson,
JCTC 17, 4929 (2021).

It is the arbiter and not the default, and the reason was measured rather than
inherited.  On the sixteen-atom Diels-Alder, from the two ends a scan leaves:
the band converged in 18 band iterations and 17 climbing steps -- **203
gradients** -- and reached an imaginary mode of **-393.60 cm-1**, where
``path_to_saddle`` on the same two ends within the same hour reached
**-393.53**.  The same saddle to 0.07 cm-1.  So the band buys a second opinion
that shares no machinery with the first, and buys nothing else; the chain goes
first because it is cheaper, and this sits beside it because it is different.

**The wall time was not what this repository had recorded**, and the
correction is the point of running it again.  416 s was the serial number.
Measured here, the same band, same box, same hour: 272 s on one process and
**39.4 s on eight**, a factor of 6.9 -- the images are independent gradients
and ORCA computes them at once, which is why ``NProcs`` above ``NImages`` is
the one setting that cannot help.  Later rounds under a load average of 909
and 912 on 384 cores gave 106 s and 156 s for the band against 208 s and 279 s
for the chain, so on a many-core box the band is not the slow route in wall
time at all; it is the parallel one.  Seconds here are about the machine and
are quoted with the load they were taken at.  Gradients are not, which is why
the gradient count leads.

**Two hazards, handled rather than discovered.**  Both fail in the way that
costs the whole timeout and returns nothing, and both are a second's
arithmetic to refuse:

* Two ends whose interpolation strands a fragment.  Measured on a proton four
  Angstrom from where it is going: ORCA announced "Two fragments detected in
  initial structure. Preparing initial structure... Done", failed three
  interior images with "the Numerical calculation ISN'T COMPLETE", and stopped
  on ``ERROR: GBW-File s_im1.gbw not found``.  Both endpoints computed
  perfectly, which is what makes it hard to see coming, and ORCA's own
  fragment handling did not rescue it.  ``path_comes_apart`` walks the
  straight line between the two ends and counts the pieces at each point.
* Two ends whose atoms are in different orders.  A band joins atom 1 to atom
  1.  The editor's usual way of making a second end preserves the order by
  construction -- morphing a structure moves atoms and does not renumber them
  -- but ``Adjust H`` fills and trims hydrogens, and a hydrogen taken out
  moves every atom after it up one.  ORCA's own NEB tutorial names
  hand-prepared endpoints in the wrong order as the commonest failure.  The
  check is on the element column, so it catches that and does **not** catch
  two atoms of the same element exchanged -- measured, and pinned below,
  because a check whose limit is not written down gets mistaken for cover.

**It reaches g-xTB.**  ``! ExtOpt NEB-TS`` over the interface
:mod:`delfin.dashboard.gxtb_engrad` answers converged on the same Diels-Alder
in **716 s** and 225 gradients, over a barrier of 22.9 kcal/mol -- so the most
accurate method the editor has is not shut out of the band, and it wants half
an hour rather than ten minutes because every gradient there is a separate
process.
"""

import pytest

from delfin.dashboard import saddle
from editor_source import EDITOR_SOURCE

_needs_orca = pytest.mark.skipif(saddle.find_orca() is None,
                                 reason='ORCA not installed')


def _has_xtb():
    from delfin.dashboard import gfn_optimize as gfn
    return gfn.find_binary('gfn2') is not None


_needs_xtb = pytest.mark.skipif(not _has_xtb(), reason='xtb not installed')

_REACTANT = """16
butadiene and ethylene, relaxed
C           -1.51175969445403       -0.02647380860866       -0.06557436490527
C           -0.72674009436823        1.04324898539530       -0.13156348350727
C            0.72671613615743        1.04323688164824       -0.13200313736446
C            1.51174712621305       -0.02646817123412       -0.06592265675999
H           -2.58463144839309        0.06203170569435       -0.06842720691043
H           -1.11830994961107       -1.02635341307553       -0.00678112167498
H           -1.18275196117662        2.02412528656961       -0.19063949165969
H            1.18273106412604        2.02410195941484       -0.19113584108773
H            2.58463111371453        0.06202058682634       -0.06829703412017
H            1.11831903097964       -1.02629431053540       -0.00619815580858
C           -0.65812656382745       -0.39453753897557        3.15439731529169
C            0.65814270849405       -0.39445452845154        3.15457854734613
H           -1.22868654033087        0.51904761763282        3.14158888465065
H           -1.22983917433790       -1.30795394609327        3.16689639740950
H            1.22859345191931        0.51920014666317        3.14189155138364
H            1.22996179489533       -1.30780145287077        3.16718979771703
"""

_PRODUCT = """16
cyclohexene, relaxed
C           -1.36290194815897       -0.06313150970649        0.75107481477849
C           -0.66218260952155        1.09869132587048        0.12233880799005
C            0.66173055731454        1.09913795139967        0.12248552655879
C            1.36310557821274       -0.06226438766873        0.75126945964802
H           -2.43534601907205        0.11530592358156        0.83579350839840
H           -1.21914017994592       -0.94894767327948        0.12363419262958
H           -1.24659833279792        1.90123497817873       -0.30035668134291
H            1.24569990280881        1.90207650485640       -0.30007507370404
H            2.43541532174521        0.11689111711109        0.83621469395190
H            1.22007694226306       -0.94810437965944        0.12369377517092
C           -0.76948705831100       -0.34196841669303        2.14125370859433
C            0.76954534349931       -0.34170006236609        2.14127653272432
H           -1.13126031423455        0.42887199800881        2.82393080033168
H           -1.14957836569781       -1.29937167674247        2.50173435160042
H            1.13098517132138        0.42903787478979        2.82424151569521
H            1.14993301057964       -1.29908356768101        2.50149006697273
"""

#: A proton four Angstrom from where it is going, on the far side of nothing.
#: Both ends are a hydrogen fluoride beside a fluoride and both compute; the
#: straight line between them leaves the proton bonded to neither, which is the
#: shape of the failure the pre-screen exists for.
_STRANDED_FROM = """3
HF beside F
F 0.000 0.000 0.000
F 0.000 4.000 0.000
H 0.930 0.000 0.000
"""
_STRANDED_TO = """3
F beside HF
F 0.000 0.000 0.000
F 0.000 4.000 0.000
H 0.930 4.000 0.000
"""


def test_the_pieces_of_a_geometry_are_counted_by_the_graph_it_already_has():
    """One union-find over the bonds the viewer draws, and nothing else.

    Whether a structure has come apart has to be answerable at every point of
    an interpolation before anything is submitted, so it is arithmetic on the
    coordinates rather than a calculation.
    """
    assert saddle.pieces(frozenset(), 3) == 3
    assert saddle.pieces(frozenset({(0, 1), (1, 2)}), 3) == 1
    assert saddle.pieces(frozenset({(0, 1)}), 3) == 2
    # Pairs naming atoms that are not there do not join anything.
    assert saddle.pieces(frozenset({(0, 9)}), 2) == 2


def test_a_way_that_strands_a_fragment_is_refused_before_orca_is_started():
    """The failure costs the whole timeout, so it is refused in a second.

    Endpoints whose interpolation pulls a piece away compute perfectly at both
    ends and kill every image between them, which is exactly why it has to be
    looked for rather than waited out.  Named in words a chemist can act on:
    how far across it happens, how many pieces, and which bond went.
    """
    said = saddle.path_comes_apart(_STRANDED_FROM, _STRANDED_TO)
    assert said
    assert 'comes apart in the middle' in said
    assert 'pieces' in said
    # And the pair that really is two ends of one reaction is not refused.
    assert saddle.path_comes_apart(_REACTANT, _PRODUCT) == ''


def test_two_ends_whose_atoms_are_in_different_orders_are_refused():
    """A band joins atom 1 to atom 1, and Adjust H is what breaks that.

    ORCA's own NEB tutorial names hand-prepared endpoints in the wrong order
    as the commonest way of getting nothing back, and the editor has a control
    that can cause it: filling or trimming a hydrogen moves every atom after
    it.  So the message names that control rather than leaving the user to
    work out what they did.
    """
    # A hydrogen where a carbon was, which is what a trimmed hydrogen does to
    # every atom after it.
    swapped = _REACTANT.splitlines()
    swapped[2], swapped[6] = swapped[6], swapped[2]
    swapped = '\n'.join(swapped) + '\n'
    said = saddle.same_atoms(_REACTANT, swapped)
    assert 'same order' in said
    assert 'Adjust H' in said

    short = _REACTANT.splitlines()
    short = '15\none H fewer\n' + '\n'.join(short[2:-1]) + '\n'
    said = saddle.same_atoms(_REACTANT, short)
    assert '16 and 15 atoms' in said
    assert 'Adjust H' in said

    assert saddle.same_atoms(_REACTANT, _PRODUCT) == ''


def test_the_order_check_does_not_catch_two_atoms_of_the_same_element():
    """The limit of the cheap check, pinned so it is not mistaken for cover.

    Measured: exchanging two carbons between the two ends leaves the element
    column identical and this passes.  The band then runs and converges onto
    whatever that mapping describes, which is a real saddle of a different
    reaction.  There is no cheap test for it -- telling two orderings of one
    molecule apart is the atom-mapping problem -- so what the message says is
    the defence: build the second end from the first rather than separately.
    """
    swapped = _REACTANT.splitlines()
    swapped[2], swapped[12] = swapped[12], swapped[2]     # two carbons
    assert saddle.same_atoms(_REACTANT, '\n'.join(swapped) + '\n') == ''
    source = open(saddle.__file__, encoding='utf-8').read()
    assert 'two atoms of the same element exchanged' in source


def test_the_barrier_is_read_from_the_last_table_and_from_the_right_column():
    """ORCA prints its path summary twice and the two are not the same shape.

    The one written when the band converges carries a distance column that the
    one written after the climb does not, so a fixed column index read a
    distance as an energy -- caught on the g-xTB run, where the first table has
    six columns and the second five.  The forces are the last two either way,
    and the rise is the one in front of them.

    And the rise that matters is the row marked ``TS``, not the highest image:
    measured on the Diels-Alder, the highest image says 6.98 kcal/mol and the
    sharpened saddle says 6.74.
    """
    with_distance = (
        'PATH SUMMARY\n'
        'Image Dist.(Ang.)    E(Eh)   dE(kcal/mol)  max(|Fp|)  RMS(Fp)\n'
        '  0     0.000    -234.53571      0.00       0.01096   0.00302\n'
        '  1     0.695    -234.53309      1.64       0.00426   0.00151\n')
    assert saddle._band_profile(with_distance)['barrier'] == 1.64

    both = with_distance + (
        'PATH SUMMARY FOR NEB-TS\n'
        'Image     E(Eh)   dE(kcal/mol)  max(|Fp|)  RMS(Fp)\n'
        '  0    -17.82300     0.00       0.00016   0.00006\n'
        '  5    -17.81188     6.98       0.00105   0.00050 <= CI\n'
        ' TS    -17.81226     6.74       0.00003   0.00001 <= TS\n'
        '  9    -17.92484   -63.90       0.00020   0.00008\n'
        'Number of SCF / gradient calculations:\n'
        '     NEB                               ... 186  91.6%\n'
        '     TS optimization                   ... 17  8.4%\n')
    read = saddle._band_profile(both)
    assert read['barrier'] == 6.74
    assert read['reaction'] == -63.90
    assert read['gradients'] == 203


def test_a_band_asks_for_no_more_processes_than_it_has_images():
    """``NProcs <= NImages``, because past that there is nothing to compute.

    The images are independent gradients and that is the whole of the
    parallelism: measured on the sixteen-atom Diels-Alder, the same band took
    272 s on one process and 39.4 s on eight.  A ninth would have nothing to
    do and would take a core off somebody else on the login node.
    """
    source = open(saddle.__file__, encoding='utf-8').read()
    assert 'ranks = 1 if own_program is not None else min(_share(cores), band)' \
        in source
    assert 'NImages {band}' in source


@_needs_orca
@_needs_xtb
@pytest.mark.slow
def test_the_band_reaches_the_same_saddle_the_fast_chain_does():
    """The measurement the whole design rests on, run rather than quoted.

    Measured on this box: 203 gradients, a barrier of 6.74 kcal/mol, and an
    imaginary mode of -393.60 cm-1 against the chain's -393.53 on the same two
    ends within the same hour.  The same saddle to 0.07 cm-1, for about twice
    the work -- which is the argument for keeping it as the second entry in
    the box and not the first.

    Marked slow because it is minutes: 39.4 s on an eight-process run on an
    idle-ish box, and 106 to 156 s at a load average of 909 on 384 cores.
    """
    seen = []
    got = saddle.neb_to_saddle(_REACTANT, _PRODUCT, 'gfn2', cores=8,
                               on_frame=lambda w, e: seen.append(len(w)))
    assert got.get('ok'), got.get('status')
    modes = (got.get('imaginary') or {}).get('modes') or []
    assert modes and -420.0 < modes[0] < -370.0, modes
    assert got['verdict']['first_order']
    assert 4.0 < got['barrier'] < 10.0, got['barrier']
    # It streamed while it ran, which is the whole reason it is allowed to
    # take minutes.
    assert seen and seen[-1] > 20


def test_the_band_is_another_how_and_not_another_button():
    """One press, two boxes, and a band is an answer to the second one.

    "Where the start comes from" and "how it gets there" are the two questions
    the transition-state control asks, and a band is plainly another *how*
    from a two-ended start.  A fourth press beside the others would be the
    third time this row learnt that lesson.

    Offered after ORCA rather than before, because the order of the list is
    the recommendation, and only from two ends -- there is no band between one
    structure and itself.
    """
    source = EDITOR_SOURCE
    assert "out.append(('through NEB-TS', 'neb'))" in source
    # Inside the `start != 'here'` branch, with the by-hand and walk-only
    # entries, and after the ORCA entry.
    assert source.index("out.append(('through ORCA', 'orca'))") < \
        source.index("out.append(('through NEB-TS', 'neb'))") < \
        source.index("out.append(('by hand', 'hand'))")
    assert "if how == 'neb':" in source
    assert '_band_between(ends)' in source


def test_a_band_streams_and_the_same_press_stops_it():
    """It runs for minutes, so it has to be watchable and interruptible.

    ORCA writes every band it accepts to ``<base>_MEP_ALL_trj`` and then every
    step of the climb to ``<base>_trj``, both with energies on the comment
    line, so one channel plays the whole run in the order it happened.  The
    press says Stop while it runs, and stopping keeps what was reached --
    which matters most here, because a band is the longest thing this press
    starts.
    """
    source = open(saddle.__file__, encoding='utf-8').read()
    assert "trails = (folder / 'in_MEP_ALL_trj.xyz', folder / 'in_trj.xyz')" \
        in source
    assert 'if should_stop is not None and should_stop():' in source

    editor = EDITOR_SOURCE
    assert "if state.get('band_run'):" in editor
    assert "state['band_stop'] = True" in editor
    assert "should_stop=lambda: bool(state.get('band_stop'))" in editor
    assert "state['band_frame_run']" in editor


def test_a_band_on_g_xtb_is_given_the_time_a_process_per_gradient_needs():
    """It reaches g-xTB, and wants half an hour rather than ten minutes.

    ``! ExtOpt NEB-TS`` over :mod:`delfin.dashboard.gxtb_engrad` converged on
    the same Diels-Alder in 716 s and 225 gradients -- past the ten minutes the
    methods ORCA drives itself are given, which would have stopped it a few
    iterations from the answer and bought nothing.  Every gradient there is a
    separate process with an ORCA start-up in front of it.
    """
    assert saddle.neb_seconds_for('gxtb') > 716.0
    assert saddle.neb_seconds_for('gfn2') == saddle.NEB_SECONDS
    assert 'gxtb' in saddle.SADDLE_METHODS
    source = open(saddle.__file__, encoding='utf-8').read()
    # One process through ExtOpt, for the reason optimise_to_saddle gives.
    assert 'ranks = 1 if own_program is not None' in source
