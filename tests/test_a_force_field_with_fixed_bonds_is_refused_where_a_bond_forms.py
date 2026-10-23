"""GFN-FF is offered for the scans it can walk and refused for the ones it cannot.

GFN-FF perceives its bonding once, from the geometry it is first handed, and
then holds the molecule together with it.  xtb's own documentation draws the
consequence: "GFN-FF can only break bonds, dissociation reactions will
therefore usually work fine, while association reactions are likely to fail."
ORCA's GOAT bars it from uphill steps for the same reason.

So the boundary is a *direction*, not a method, and it was established by
measurement rather than by refusing the method wholesale -- GFN-FF is a
hundred times faster than GFN2 and is the right engine for a great many scans.
Butadiene and ethylene, one forming C-C driven in 0.1 A steps with everything
relaxed at every point, under GFN2 and under GFN-FF:

    a bond FORMED, 3.40 -> 1.60 A

        GFN2      crosses at +7.3 kcal/mol at 2.20 A and lands at -63.0 in
                  the product; the *other* forming C-C closes to 1.53 A
                  without being asked
        GFN-FF    climbs to +94.1 kcal/mol and crosses nothing; the other
                  forming C-C ends at 3.39 A.  No reaction, and 87 kcal/mol
                  of error in the one number a scan exists to produce

    the same bond BROKEN, 1.53 -> 3.40 A

        GFN2      +79.0 kcal/mol at 2.85 A
        GFN-FF    +84.3 kcal/mol -- 7 % apart on a bond dissociation, which
                  is a force field being a force field

The editor's topology cache makes the first of those worse rather than better.
The cache exists so that a drag cannot fall apart between one frame and the
next, and it works; but with the topology pinned the false profile is smooth
and monotonic, which is what a wall of repulsion looks like when it is drawn
carefully.  Unpinned, the same scan gives +108.6 with its maximum in a
different place -- visibly wrong, and therefore less dangerous.

Refused before the run, in the way an unparametrised solvent is, because xtb
does not refuse it: it converges, it reports a number, and the number is
repulsion wearing the name of a barrier.
"""

from __future__ import annotations

import math
import shutil

import pytest

from delfin.dashboard import gfn_optimize as gfn
from editor_source import EDITOR_SOURCE

_needs_xtb = pytest.mark.skipif(not shutil.which('xtb'),
                                reason='xtb not installed')

#: Butadiene and ethylene, the forming C-C pair 3.40 A apart.  Atoms 1 and 2
#: are one forming bond and atoms 0 and 5 the other; the scan drives the first
#: and nobody drives the second, which is the whole point of the measurement.
_PAIR = """16
butadiene and ethylene, one forming C-C at 3.40 A
C    -1.664289   -1.281060   -0.752904
C    -2.185424   -0.117961   -0.424283
C     0.355880    1.565209    1.081937
C     1.148259    1.329226    0.042004
C     1.725647    0.046361   -0.324306
C     1.557030   -1.103104    0.320139
H    -1.088669   -1.414591   -1.653679
H    -1.785440   -2.156569   -0.136069
H    -2.061675    0.757232   -1.040193
H    -2.759310    0.017138    0.477951
H    -0.031590    2.549441    1.282094
H     0.062789    0.791023    1.769230
H     1.407068    2.150058   -0.615584
H     2.346407    0.062614   -1.211750
H     0.951475   -1.182611    1.205918
H     2.021841   -2.012403   -0.020504
"""


def _armed(atoms, to, kind='distance'):
    return [{'kind': kind, 'atoms': list(atoms), 'to': float(to)}]


def test_driving_two_atoms_together_across_a_bond_is_refused():
    said = gfn.gfnff_refusal(_PAIR, _armed((1, 2), 1.60))
    assert said
    assert 'GFN-FF' in said and 'association' in said
    # Named by element and number, and told which way it was being driven.
    assert 'C2-C3' in said
    assert 'GFN2' in said and 'g-xTB' in said


def test_driving_the_same_two_atoms_apart_is_not():
    """Two atoms 3.40 A apart, driven further apart, cross no bond at all."""
    assert gfn.gfnff_refusal(_PAIR, _armed((1, 2), 4.50)) == ''
    # And stopping short of the line is not crossing it either.
    assert gfn.gfnff_refusal(_PAIR, _armed((1, 2), 2.40)) == ''


def test_a_bond_being_stretched_is_never_refused():
    """The ordinary dissociation scan, which is what the method is good at."""
    ethane = ('8\nethane\nC 0 0 0\nC 1.53 0 0\n'
              'H -0.36 1.02 0\nH -0.36 -0.51 0.88\nH -0.36 -0.51 -0.88\n'
              'H 1.89 -1.02 0\nH 1.89 0.51 0.88\nH 1.89 0.51 -0.88\n')
    assert gfn.gfnff_refusal(ethane, _armed((0, 1), 3.50)) == ''
    # And squeezing a bond that is already a bond is not making one either.
    assert gfn.gfnff_refusal(ethane, _armed((0, 1), 1.30)) == ''


def test_an_angle_or_a_torsion_makes_no_bond_and_is_left_alone():
    assert gfn.gfnff_refusal(_PAIR, _armed((1, 2, 3), 90.0, 'angle')) == ''
    assert gfn.gfnff_refusal(
        _PAIR, _armed((1, 2, 3, 4), 180.0, 'dihedral')) == ''


def test_the_line_is_the_one_the_rest_of_the_editor_draws_bonds_at():
    """One threshold for what a bond is, not a second one invented here."""
    found = gfn.gfnff_would_form(_PAIR, _armed((1, 2), 1.60))
    assert found is not None
    # Two carbons: 1.3 times the sum of the covalent radii.
    from delfin.atom_mapping import cov_radius
    want = gfn.BOND_STARTS_AT * 2.0 * cov_radius('C')
    assert math.isclose(found['bonds_at'], want, rel_tol=1e-9)


def test_nothing_is_refused_for_a_structure_that_is_not_there():
    assert gfn.gfnff_refusal('', _armed((1, 2), 1.60)) == ''
    assert gfn.gfnff_refusal(_PAIR, []) == ''
    assert gfn.gfnff_refusal(_PAIR, _armed((1, 99), 1.60)) == ''


def test_the_scan_asks_before_it_starts_and_only_of_gfnff():
    """Said before the run, because xtb will not say it after one."""
    assert "if str(method).strip().lower() == 'gfnff':" in EDITOR_SOURCE
    assert '_gfn.gfnff_refusal(xyz, legs)' in EDITOR_SOURCE


# ---------------------------------------------------------------------------
# and the other route where the question is decidable
# ---------------------------------------------------------------------------


def test_a_pair_of_ends_that_makes_a_bond_is_refused_too():
    """Two ends say what the reaction is; one geometry does not.

    Measured with xtb's own path finder given the separated pair and
    cyclohexene: GFN2 reports the product 68.0 kcal/mol below the start and
    GFN-FF reports it 34.3 *above*, because it never sees the product as a
    molecule -- it prices cyclohexene as a strained contact between the two
    things it still believes are there.  The sign of a reaction energy is not
    a detail.
    """
    ethene = ('6\nethene\nC 0 0 0\nC 1.33 0 0\n'
              'H -0.57 0.94 0\nH -0.57 -0.94 0\n'
              'H 1.90 0.94 0\nH 1.90 -0.94 0\n')
    ethane = ethene.replace('C 1.33 0 0', 'C 1.53 0 0')
    apart = ethene.replace('C 1.33 0 0', 'C 3.60 0 0')
    # A bond that appears between the first end and the second.
    said = gfn.gfnff_pair_refusal(apart, ethene)
    assert said and 'C1-C2' in said and 'cannot make one' in said
    # The same pair the other way round breaks one, and is left alone.
    assert gfn.gfnff_pair_refusal(ethene, apart) == ''
    # And two ends of the same molecule ask nothing of the topology.
    assert gfn.gfnff_pair_refusal(ethene, ethane) == ''


def test_the_two_ends_press_asks_before_it_walks():
    """All three ways from a pair go through the path finder, so all three
    are refused together."""
    assert '_gfn.gfnff_pair_refusal(ends[0], ends[1])' in EDITOR_SOURCE
    where = EDITOR_SOURCE.index('_gfn.gfnff_pair_refusal(')
    # Before the branch that chooses between them.
    assert EDITOR_SOURCE.index("if how == 'orca':\n            _path_then_orca") \
        > where


@_needs_xtb
def test_it_really_does_answer_the_wrong_question(tmp_path):
    """The measurement the refusal is made of, run again.

    Three points rather than the nineteen the table above was taken over --
    enough to show that GFN-FF climbs where GFN2 falls, which is the whole
    disagreement.  The two engines are asked exactly the same thing.
    """
    here = _PAIR
    there = _PAIR
    rose = []
    fell = []
    for want in (2.60, 2.20, 1.80):
        held = [{'kind': 'distance', 'atoms': [1, 2], 'mode': 'fix',
                 'value': want}]
        one = gfn.optimize_with_gfn(here, 'gfn2', max_steps=60, timeout=None,
                                    constraints=held)
        two = gfn.optimize_with_gfn(there, 'gfnff', max_steps=60, timeout=None,
                                    constraints=held,
                                    topology=tmp_path / 'topo')
        if not (one['ok'] and two['ok']):
            pytest.skip(one.get('status') or two.get('status'))
        here, there = one['xyz'], two['xyz']
        fell.append(float(one['energy']))
        rose.append(float(two['energy']))
    # GFN2 goes over and down into the product; GFN-FF only ever goes up.
    assert fell[-1] < fell[0], fell
    assert rose[-1] > rose[0], rose
    # And the other forming bond follows under GFN2 and does not under GFN-FF.
    where = gfn.coordinates_of(here)
    other = gfn.coordinates_of(there)

    def span(c, i, j):
        return math.dist(c[3 * i:3 * i + 3], c[3 * j:3 * j + 3])

    assert span(where, 0, 5) < span(other, 0, 5)
