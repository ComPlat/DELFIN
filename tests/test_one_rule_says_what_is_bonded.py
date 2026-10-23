"""What the graph holds the user to and what the hand perceives are one rule.

They were two.  :func:`gfn_optimize.bond_graph` judged at
:data:`gfn_optimize.BOND_STARTS_AT`, 1.3 times the two covalent radii, and
``_is_a_bond`` -- which every perception in that file leans on -- carried its
own default of 1.25 with nothing written down beside it.  Between them lay a
band where a step that broke the bond was still refused while the perception
had already stopped seeing one, and for a C-C that band, 1.900 to 1.976 A, is
where a hand pulling a bond apart spends its time.

The cost is not the disagreement itself but a coordinate: a torsion is built
through bonds, so inside the band the gesture could not be named as the turn
it was.  Priced live under GFN2 on butane, its C2-C3 out at 1.947 A and the end
carbon swung: +0.32 kcal/mol and the bonding held with one rule, +25.63 and the
step refused for breaking C4-H14 with two.  The measurement is written up
beside :data:`gfn_optimize.BOND_STARTS_AT`; what is checked here is the rule
itself, which needs no binary.
"""
from __future__ import annotations

import math

import pytest

from delfin.atom_mapping import cov_radius
from delfin.dashboard import gfn_optimize as gfn


def _chain(last: float) -> str:
    """Four carbons with hydrogens, the last C-C stretched to *last*."""
    rows = [('C', 0.0, 0.0, 0.0), ('C', 1.53, 0.0, 0.0), ('C', 2.04, 1.44, 0.0)]
    rows.append(('C', 2.04 + last * 0.34, 1.44 + last * 0.94, 0.0))
    for _, x, y, _z in list(rows):
        rows.append(('H', x, y - 0.30, 1.05))
        rows.append(('H', x, y + 0.30, -1.05))
    return f'{len(rows)}\nchain\n' + ''.join(
        f'{s} {x:.6f} {y:.6f} {z:.6f}\n' for s, x, y, z in rows)


def _swung(xyz_text: str, index: int, by: float) -> str:
    """One atom moved across its bond rather than along it."""
    head, name, *body = xyz_text.splitlines()
    rows = [one.split() for one in body if one.strip()]
    rows[index][3] = f'{float(rows[index][3]) + by:.6f}'
    return f'{head}\n{name}\n' + ''.join(' '.join(one) + '\n' for one in rows)


def _named(said):
    return [(one['kind'], tuple(one['atoms'])) for one in said]


def _turns_about(said, axis):
    """Whether the one coordinate named is a turn about *axis*.

    Which atom names the far end of a torsion is not the question: every one
    on the same axis is the same turn and takes the grabbed atom to the same
    place, so a chain carbon and a hydrogen on it are interchangeable here and
    the scores are all but tied between them.
    """
    return (len(said) == 1 and said[0]['kind'] == 'dihedral'
            and tuple(said[0]['atoms'][:3]) == tuple(axis))


def test_the_perception_judges_a_bond_by_the_number_the_graph_holds():
    """One rule, so there is no band for a coordinate to fall into."""
    assert (gfn._is_a_bond.__defaults__ or ())[-1] == gfn.BOND_STARTS_AT


@pytest.mark.parametrize('one,two', [
    ('C', 'C'), ('C', 'H'), ('O', 'H'), ('C', 'Cl'), ('Pd', 'Br'), ('Cs', 'I'),
])
def test_the_two_never_disagree_about_a_pair(one, two):
    """Across the whole approach, and on the pairs whose radii differ most.

    A multiplied threshold puts the two rules further apart the larger the
    radii are -- 0.076 A for a C-C, 0.192 for a Cs-I -- so a heavy pair is
    where a second number would show first.
    """
    reach = cov_radius(one) + cov_radius(two)
    for share in (0.6, 0.9, 1.0, 1.2, 1.24, 1.25, 1.28, 1.3, 1.31, 1.5, 2.0):
        span = share * reach
        text = f'2\npair\n{one} 0 0 0\n{two} 0 0 {span:.6f}\n'
        # Both asked of the same parsed numbers, so a pair sitting exactly
        # on the threshold cannot be told apart by rounding rather than rule.
        flat = gfn.coordinates_of(text)
        where = [tuple(flat[:3]), tuple(flat[3:6])]
        graph = (0, 1) in gfn.bond_graph(text)
        hand = gfn._is_a_bond(where, [cov_radius(one), cov_radius(two)], 0, 1)
        assert graph == hand, (one, two, share, span, graph, hand)


def test_a_swing_about_a_stretched_bond_is_still_named_as_the_turn():
    """The band's cost, at the distance it used to be paid at.

    1.939 A is inside the old band: the graph holds the chain to that C-C, so
    a step which breaks it is refused, and the perception used to see nothing
    there to turn about.  What it named instead was the distance from the
    swung carbon to a hydrogen on a different one -- an unrelated contact
    pinned while the turn that was happening stayed free.
    """
    text = _chain(1.939)
    at = [tuple(map(float, one.split()[1:4])) for one in gfn.atom_lines(text)]
    assert math.dist(at[2], at[3]) == pytest.approx(1.939, abs=5e-3)
    # The premise: the graph is holding the user to this bond.
    assert (2, 3) in gfn.bond_graph(text)
    said = gfn.contacts_holding(_swung(text, 3, 0.25), [3], was=text)
    assert _turns_about(said, (3, 2, 1)), _named(said)


def test_and_either_side_of_that_band_it_always_was():
    """So the test above is about the band and not about the gesture."""
    for short in (1.849, 1.75, 1.60):
        text = _chain(short)
        said = gfn.contacts_holding(_swung(text, 3, 0.25), [3], was=text)
        assert _turns_about(said, (3, 2, 1)), (short, _named(said))
    # And past where either rule calls it a bond there is no turn to name,
    # which is right: a fragment on its way off is a distance.
    text = _chain(2.049)
    assert (2, 3) not in gfn.bond_graph(text)
    said = _named(gfn.contacts_holding(_swung(text, 3, 0.25), [3], was=text))
    assert said and all(kind == 'distance' for kind, _ in said), said
