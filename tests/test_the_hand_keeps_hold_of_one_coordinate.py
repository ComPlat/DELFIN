"""What the hand is driving is decided once, not argued about every answer.

The perception worked out afresh, on every answer, which internal coordinates
a drag is changing.  Between candidates that say the same thing it was decided
by noise, and a different answer is not a smaller answer -- it is a different
experiment, so the structure jumps to suit it.

Two things made it noisy, and both are here.

* **Everything the hand holds travels; not everything it holds has anything to
  say.**  A page names every atom under the cursor, so a hand on a selected
  methyl names its hydrogens too, and each of them claimed a coordinate of its
  own.  Those coordinates are H...H contacts across the gap the drag has
  opened, every one of them changed by about the same amount, and which is
  largest depends on where the methyl happens to be pointing.

* **A hand pulls rather than places**, so the structure never arrives where
  the cursor is.  Every coordinate spanning that standing gap therefore reads
  as strongly changed, and the ranking between them carries no information at
  all.  The one already being held has to keep its place unless something
  beats it clearly.

Measured on an ethane, the far methyl drawn 1.2 A out and the hand then held
still to a thousandth of an angstrom for twenty answers:

    the set changed          3 of 20 answers  ->  0 of 20
    largest atom jump        1.287 A          ->  0.002 A
    C-C over the whole hold  0.116 A of span  ->  0.002 A

The jump is the whole of what the user sees.  It is not the chemistry and it
is not the optimiser: it is the question changing.
"""
from __future__ import annotations

import pytest

from delfin.dashboard import climb as _climb
from delfin.dashboard import gfn_optimize as gfn

# ``have_fast_gradients`` is climb's, not gfn_optimize's -- see the same guard
# in test_keep_bonds_keeps_the_structure_the_user_takes_away, and the CI run it
# took down when it named the wrong module.
_needs_xtb = pytest.mark.skipif(
    gfn.find_xtb() is None and not _climb.have_fast_gradients(),
    reason='no xtb to relax with')

_ETHANE = """8
ethane
C  0.000000  0.000000  0.762900
C  0.000000  0.000000 -0.762900
H -0.505000  0.874000  1.162900
H -0.505000 -0.874000  1.162900
H  1.010000  0.000000  1.162900
H  0.505000  0.874000 -1.162900
H  0.505000 -0.874000 -1.162900
H -1.010000  0.000000 -1.162900
"""

#: One carbon and two partners, one of them a little more the way it is going.
#: The scores are 0.480 and 0.432 Angstrom, which is eleven per cent apart --
#: inside the margin, so whichever is held keeps its place.
_A_NEAR_TIE = (
    "3\ntwo partners, one a little more the way the hand is going\n"
    "C   0.00 0.00 0.00\n"
    "F   0.80 0.00 3.00\n"
    "F  -1.60 0.00 3.00\n"
)

#: The same with the second partner further round: 0.480 against 0.363, which
#: is a third better and takes the place whatever was in it.
_A_CLEAR_WIN = _A_NEAR_TIE.replace("-1.60", "-2.60")


def _pushed(text, atom=0, dz=0.50):
    """*text* with one atom moved along z, as a hand would move it."""
    rows = gfn.atom_lines(text)
    said = rows[atom].split()
    rows[atom] = (f"{said[0]}  {float(said[1]):.4f} {float(said[2]):.4f} "
                  f"{float(said[3]) + dz:.4f}")
    return f"{len(rows)}\nmoved\n" + "\n".join(rows)


def _holding(atoms):
    """One held entry, in the shape the editor keeps its last answer in."""
    return [{'kind': 'distance', 'atoms': list(atoms), 'value': 0.0,
             'mode': 'drag'}]


# ---------------------------------------------------------------------------
# Who speaks for the hand
# ---------------------------------------------------------------------------

def test_a_grabbed_methyl_speaks_through_its_carbon():
    """The hydrogens ride on it, whether or not the page named them."""
    where = [(float(a), float(b), float(c)) for _s, a, b, c in
             (line.split() for line in gfn.atom_lines(_ETHANE))]
    from delfin.atom_mapping import cov_radius
    radius = [cov_radius(line.split()[0]) for line in gfn.atom_lines(_ETHANE)]
    assert gfn.speaking_for_the_drag(where, radius, [1, 5, 6, 7]) == [1]


def test_a_hydrogen_the_hand_has_on_its_own_still_speaks():
    """Its carbon is not in the grab, so the contact is the hand's own."""
    where = [(float(a), float(b), float(c)) for _s, a, b, c in
             (line.split() for line in gfn.atom_lines(_ETHANE))]
    from delfin.atom_mapping import cov_radius
    radius = [cov_radius(line.split()[0]) for line in gfn.atom_lines(_ETHANE)]
    assert gfn.speaking_for_the_drag(where, radius, [5]) == [5]
    # Two of them on a carbon nobody is holding: two hands, two statements.
    assert gfn.speaking_for_the_drag(where, radius, [2, 3]) == [2, 3]


def test_atoms_holding_nothing_but_each_other_all_speak():
    """Silencing the whole grab would leave the drag with no coordinate.

    Two carbons approaching two others is a cycloaddition and each half is
    part of the coordinate; a hydrogen molecule brought up to a metal is the
    same shape.  Neither is riding on anything.
    """
    pair = "2\nhydrogen\nH 0.00 0.00 0.00\nH 0.00 0.00 0.74\n"
    where = [(float(a), float(b), float(c)) for _s, a, b, c in
             (line.split() for line in gfn.atom_lines(pair))]
    from delfin.atom_mapping import cov_radius
    radius = [cov_radius(line.split()[0]) for line in gfn.atom_lines(pair)]
    assert gfn.speaking_for_the_drag(where, radius, [0, 1]) == [0, 1]


def test_a_dragged_methyl_is_held_by_one_coordinate_and_it_is_the_bond():
    """Not by three, two of which are hydrogens facing hydrogens.

    This is the scored branch, which is where a drag lives after its first
    answer -- the opening one had this rule already, and the drag then spent
    the rest of the gesture without it.
    """
    was = _ETHANE
    now = _pushed(_pushed(_pushed(_pushed(was, 1, 0.6), 5, 0.6), 6, 0.6),
                  7, 0.6)
    held = gfn.contacts_holding(now, [1, 5, 6, 7], most=3, was=was)
    assert [set(one['atoms']) for one in held] == [{0, 1}], held


# ---------------------------------------------------------------------------
# And it keeps hold of it
# ---------------------------------------------------------------------------

def test_a_held_coordinate_keeps_its_place_against_a_near_tie():
    """Eleven per cent is not a better description of the same gesture."""
    now = _pushed(_A_NEAR_TIE)
    free = gfn.contacts_holding(now, [0], most=3, was=_A_NEAR_TIE)
    assert [set(one['atoms']) for one in free] == [{0, 1}], free
    kept = gfn.contacts_holding(now, [0], most=3, was=_A_NEAR_TIE,
                                holding=_holding([0, 2]))
    assert [set(one['atoms']) for one in kept] == [{0, 2}], kept


def test_a_clearly_better_coordinate_still_takes_the_place():
    """Otherwise the first answer of a drag would decide the whole of it."""
    now = _pushed(_A_CLEAR_WIN)
    kept = gfn.contacts_holding(now, [0], most=3, was=_A_CLEAR_WIN,
                                holding=_holding([0, 2]))
    assert [set(one['atoms']) for one in kept] == [{0, 1}], kept


def test_the_margin_is_the_share_the_perception_already_uses():
    """One number for "these two say the same thing", not two.

    :data:`gfn_optimize._MOVED_SHARE` decides which coordinates are kept
    alongside the leading one because they are part of the same drag; a
    challenger inside that band is by the perception's own reckoning not
    describing a different gesture, so it has no claim on the place.
    """
    assert gfn._HELD_KEEPS_ITS_PLACE == gfn._MOVED_SHARE


# ---------------------------------------------------------------------------
# The whole of it, with a hand that stands still
# ---------------------------------------------------------------------------

@_needs_xtb
def test_a_hand_that_stands_still_is_answered_the_same_way_every_time():
    """The measurement the two changes above exist for.

    A hand that keeps advancing hides this completely: once the leash clamps,
    repeated wishes are deduplicated and the trace goes flat.  A hand that has
    arrived somewhere and is merely being held -- which is what a hand does
    most of the time -- asks the same question over and over, and the answer
    was allowed to be different every time.

    Judged on the largest distance any atom moved between two answers, after
    the overall shift and turn are taken out, so a molecule drifting as a
    whole does not read as one shaking.  Watching a bond length instead sees
    nothing: what jumped was the methyl turning.
    """
    import time

    import numpy as np

    helper = pytest.importorskip('test_the_budget_prices_a_relaxed_path')

    # That module's ethane, not this one's: its C-C lies along x and _shifted
    # moves along x, so the hand is drawing the methyl straight off the bond.
    # Laid the other way the same numbers describe a sideways swing, which is
    # a different gesture and does not show this at all.
    start = gfn.optimize_with_gfn(helper._ETHANE, 'gfn2', max_steps=400,
                                  timeout=300)
    assert start.get('ok'), start.get('status')
    begin = start['xyz']
    methyl = {1, 5, 6, 7}

    def rows(xyz):
        return np.asarray(gfn.coordinates_of(xyz), float).reshape(-1, 3)

    def moved(a, b):
        """The largest an atom moved, once the whole molecule is laid back."""
        a, b = a - a.mean(0), b - b.mean(0)
        left, _s, right = np.linalg.svd(a.T @ b)
        turn = left @ np.diag(
            [1, 1, np.sign(np.linalg.det(left @ right))]) @ right
        return float(np.linalg.norm(a @ turn - b, axis=1).max())

    part = helper._a_part(begin)
    part.submit_ff_dd.value = 'gfn2'
    part.submit_relax_btn.value = True
    part.submit_hand_dd.value = 'pull'
    part.submit_temperature.value = 298.15
    part.submit_thermal_btn.value = True
    began = time.time()
    while part.state.get('thermal_e0') is None and time.time() - began < 300:
        time.sleep(0.05)
    assert part.state.get('thermal_e0') is not None, 'no budget to price with'
    part._begin_gfn_follow()
    part._arm_thermal_leash()

    # Out to where a hand would stop, and then held there.
    for n in range(1, 7):
        part.submit_manip_sync.value = helper._drag_message(
            helper._shifted(begin, methyl, 1.2 * n / 6),
            'DELFIN drag-follow held=1,5,6,7')
        helper._quiet(part.state)

    rng = np.random.default_rng(3)
    jumps, named = [], []
    last = rows(part.coords_widget.value)
    for _n in range(20):
        # A real hand standing still still trembles a thousandth.
        part.submit_manip_sync.value = helper._drag_message(
            helper._shifted(begin, methyl, 1.2 + float(rng.normal(0, 0.001))),
            'DELFIN drag-follow held=1,5,6,7')
        helper._quiet(part.state)
        here = rows(part.coords_widget.value)
        jumps.append(moved(last, here))
        named.append(tuple(sorted(
            tuple(sorted(one['atoms']))
            for one in (part.state.get('thermal_holding') or ()))))
        last = here

    assert len(set(named)) == 1, f'the question changed: {named}'
    assert max(jumps) < 0.1, (
        f'an atom moved {max(jumps):.3f} A under a hand that was not moving; '
        f'jumps {["%.3f" % one for one in jumps]}')
