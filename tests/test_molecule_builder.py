"""Structural edits have to produce chemistry, not just coordinates.

Every expectation here is a molecule anyone can check by name: a placed carbon
is methane, growing twice gives propane, a bond raised to double turns ethane
into ethene at 1.329 A.
"""
from __future__ import annotations

import math
from collections import Counter

import pytest

rdkit = pytest.importorskip('rdkit')

from delfin.dashboard import molecule_builder as B  # noqa: E402
from delfin.dashboard import molecule_forcefield as mff  # noqa: E402


def _empty():
    return B.Structure([], [], {})


def _formula(structure):
    return dict(Counter(structure.symbols))


def _carbons(structure):
    return [i for i, s in enumerate(structure.symbols) if s == 'C']


def _terms(structure):
    xyz = B.to_xyz(structure)
    perceived = mff.perceive_molecule(xyz)
    mff.apply_bond_edits(perceived, dict(structure.bonds))
    return mff.export_forcefield_terms(xyz, perceived=perceived, method='uff')


def _bond(payload, first, second):
    key = (min(first, second), max(first, second))
    for entry in payload['bonds']:
        if (entry['i'], entry['j']) == key:
            return round(entry['r0'], 3), round(entry['k'])
    return None


def test_a_placed_carbon_is_methane():
    """Placing an atom fills its valence, so the structure is a molecule at
    every point rather than a radical waiting to be finished."""
    structure = _empty()
    B.place_atom(structure, 'C', (0.0, 0.0, 0.0))
    assert _formula(structure) == {'C': 1, 'H': 4}

    # And it is a real tetrahedron, not four hydrogens in a heap.
    centre = structure.coords[0]
    lengths = [
        math.dist(centre, structure.coords[j]) for j in range(1, 5)
    ]
    assert all(abs(d - 1.07) < 0.05 for d in lengths), lengths
    angles = []
    for a in range(1, 5):
        for b in range(a + 1, 5):
            angles.append(mff._angle_degrees(
                structure.coords[a], centre, structure.coords[b]))
    assert all(abs(v - 109.47) < 0.5 for v in angles), angles


def test_growing_twice_gives_propane():
    """A hydrogen already in the way is consumed rather than pushed aside,
    which is what makes repeated growing build a chain instead of a thicket."""
    structure = _empty()
    first = B.place_atom(structure, 'C', (0.0, 0.0, 0.0))
    second = B.grow_from(structure, first, 'C')
    B.grow_from(structure, second, 'C')
    assert _formula(structure) == {'C': 3, 'H': 8}
    assert mff.relax_xyz(B.to_xyz(structure), method='uff')['ok']


def test_a_bond_raised_to_double_takes_the_hydrogens_with_it():
    """Setting a bond to double is not a label: ethane becomes ethene, and the
    two hydrogens that no longer fit have to go. Without that the carbons are
    five-valent, the typing molecule refuses the order, and the force field
    hands back a single bond -- measured on benzene, where a drawn double bond
    came out at 1.514 A, exactly the single it had been asked to stop being.

    Ethane, ethene, ethyne and back, with the UFF parameters each time."""
    structure = _empty()
    B.place_atom(structure, 'C', (0.0, 0.0, 0.0))
    B.grow_from(structure, _carbons(structure)[0], 'C')
    assert _formula(structure) == {'C': 2, 'H': 6}
    assert _bond(_terms(structure), *_carbons(structure)) == (1.514, 700)

    expected = {2: ({'C': 2, 'H': 4}, (1.329, 1035)),
                3: ({'C': 2, 'H': 2}, (1.205, 1386)),
                1: ({'C': 2, 'H': 6}, (1.514, 700))}
    for order, (formula, parameters) in expected.items():
        B.set_bond_order(structure, *_carbons(structure), order)
        assert _formula(structure) == formula, order
        assert _bond(_terms(structure), *_carbons(structure)) == parameters, order


def test_deleting_an_atom_takes_its_own_hydrogens_with_it():
    """Leaving them behind turned deleting one carbon of ethane into a carbon
    with seven hydrogens: three orphans, plus the ones grown to fill the
    valence the deletion had just freed."""
    structure = _empty()
    B.place_atom(structure, 'C', (0.0, 0.0, 0.0))
    B.grow_from(structure, _carbons(structure)[0], 'C')
    B.delete_atoms(structure, [_carbons(structure)[1]])
    assert _formula(structure) == {'C': 1, 'H': 4}
    # A hydrogen shared with something that stays is not an orphan.
    assert all(structure.neighbours(j) for j in range(len(structure)))


def test_changing_the_element_re_satisfies_the_valence():
    """Methane, ammonia, water -- and a metal, whose coordination number is
    the user's decision rather than a table's."""
    structure = _empty()
    B.place_atom(structure, 'C', (0.0, 0.0, 0.0))
    for element, formula in (('N', {'N': 1, 'H': 3}),
                             ('O', {'O': 1, 'H': 2}),
                             ('C', {'C': 1, 'H': 4})):
        B.set_element(structure, 0, element)
        assert _formula(structure) == formula, element

    B.set_element(structure, 0, 'Pt')
    assert structure.symbols[0] == 'Pt'
    assert 'H' not in [structure.symbols[j] for j in structure.neighbours(0)][:0]
    assert B.default_valence('Pt') is None


def test_hydrogens_are_only_adjusted_where_the_edit_touched():
    """Filling every unsatisfied valence would quietly hydrogenate a radical
    or a coordinating donor the user put there on purpose."""
    structure = B.Structure(
        ['C', 'H', 'H', 'H', 'C', 'H', 'H', 'H'],
        [(0, 0, 0), (1.07, 0, 0), (-0.5, 0.9, 0), (-0.5, -0.9, 0),
         (5, 0, 0), (6.07, 0, 0), (4.5, 0.9, 0), (4.5, -0.9, 0)],
        {(0, 1): 1, (0, 2): 1, (0, 3): 1, (4, 5): 1, (4, 6): 1, (4, 7): 1},
    )
    # Two methyl radicals. Touching neither must change nothing.
    assert B.adjust_hydrogens(structure, []) == {}
    assert _formula(structure) == {'C': 2, 'H': 6}
    # Touching one fills that one only.
    B.adjust_hydrogens(structure, [0])
    assert _formula(structure) == {'C': 2, 'H': 7}


def test_an_edit_never_leaves_a_stale_index_behind():
    """Removals used to happen as the loop ran, which renumbered the structure
    underneath it: the atoms still queued, and the index the caller was holding
    on to, became something else. Growing with a double bond returned the wrong
    atom and the bond it had just made looked absent."""
    structure = _empty()
    first = B.place_atom(structure, 'C', (0.0, 0.0, 0.0))
    second = B.grow_from(structure, first, 'C', order=2)
    assert structure.symbols[second] == 'C'
    assert structure.order(first, second) == 2
    assert _formula(structure) == {'C': 2, 'H': 4}


def test_elements_offered_are_real_and_metals_take_no_hydrogens():
    for symbol in B.DRAW_ELEMENTS:
        assert B.normalise_element(symbol) == symbol
    assert B.normalise_element('c') == 'C'
    assert B.normalise_element('cL') == 'Cl'
    for bad in ('', 'Xx', '4', 'Carbon', None):
        assert B.normalise_element(bad) is None
    assert B.default_valence('C') == 4
    assert B.default_valence('Fe') is None


def test_a_bond_cannot_be_given_an_order_its_ends_cannot_carry():
    """A tap on a C-H bond asking for a double bond has to be refused rather
    than producing a C=H that nothing downstream can type."""
    structure = _empty()
    B.place_atom(structure, 'C', (0.0, 0.0, 0.0))
    assert B.set_bond_order(structure, 0, 1, 2) is False
    assert _formula(structure) == {'C': 1, 'H': 4}


def test_raising_a_bond_keeps_the_octet():
    """Methanol becomes formaldehyde and methanethiol thioformaldehyde: the
    hydrogens that no longer fit come off both ends, not just the carbon."""
    for element, formula in (('O', {'C': 1, 'H': 2, 'O': 1}),
                             ('S', {'C': 1, 'H': 2, 'S': 1}),
                             ('N', {'C': 1, 'H': 3, 'N': 1})):
        structure = _empty()
        B.place_atom(structure, 'C', (0.0, 0.0, 0.0))
        B.grow_from(structure, 0, element)
        partner = structure.symbols.index(element)
        assert B.set_bond_order(structure, 0, partner, 2) is True, element
        assert _formula(structure) == formula, element


def _angles_at(structure, index):
    import itertools
    neighbours = structure.neighbours(index)
    return sorted(round(mff._angle_degrees(
        structure.coords[a], structure.coords[index], structure.coords[b]), 1)
        for a, b in itertools.combinations(neighbours, 2))


def test_hydrogens_are_placed_in_the_shape_the_bonds_imply():
    """The shape at a centre follows from its steric number -- the number of
    *partners* plus the lone pairs -- not from the sum of its bond orders.

    Counting the orders made ethene's carbons tetrahedral, because the double
    bond was counted twice; and leaving the lone pairs out made ammonia planar
    and water linear."""
    cases = [
        ('C', 109.5), ('N', 109.5), ('O', 109.5), ('S', 109.5), ('P', 109.5),
    ]
    for element, want in cases:
        structure = _empty()
        B.place_atom(structure, element, (0.0, 0.0, 0.0))
        assert abs(_angles_at(structure, 0)[0] - want) < 1.5, element


def test_raising_a_bond_moves_the_hydrogens_to_the_new_shape():
    """Ethene's carbons are trigonal planar, not two tetrahedra sharing an
    edge. Leaving the hydrogens where they were made a molecule that looked
    like a twisted ethane and perceived like one."""
    structure = _empty()
    B.place_atom(structure, 'C', (0.0, 0.0, 0.0))
    B.grow_from(structure, 0, 'C')
    assert _angles_at(structure, _carbons(structure)[0])[0] == 109.5

    B.set_bond_order(structure, *_carbons(structure), 2)
    assert _angles_at(structure, _carbons(structure)[0]) == [120.0, 120.0, 120.0]
    B.set_bond_order(structure, *_carbons(structure), 3)
    assert _angles_at(structure, _carbons(structure)[0]) == [180.0]


def test_a_deleted_hydrogen_stays_deleted():
    """Refilling the valence it freed grew it straight back, so removing one
    was a no-op: methane came back as methane. Deleting a heavy atom is
    different -- there the neighbour is left with a hole nobody asked for."""
    structure = _empty()
    B.place_atom(structure, 'C', (0.0, 0.0, 0.0))
    B.delete_atoms(structure, [1])
    assert _formula(structure) == {'C': 1, 'H': 3}
    B.delete_atoms(structure, [1])
    assert _formula(structure) == {'C': 1, 'H': 2}


def test_joining_two_fragments_brings_them_together_first():
    """Dragging one atom onto another across the viewer can ask for a bond
    between things six Angstrom apart. Handing that to the force field as a
    bond whose equilibrium is 1.5 A is a very large force on two atoms, and it
    tears the structure on the way in.

    So one fragment is moved as a rigid body until the bond is at its
    equilibrium length -- nothing inside either is strained. Two methanes six
    Angstrom apart become ethane at 1.520 A with the moved fragment's C-H
    bonds still at 1.07."""
    structure = _empty()
    B.place_atom(structure, 'C', (0.0, 0.0, 0.0))
    B.place_atom(structure, 'C', (6.0, 0.0, 0.0))
    first, second = _carbons(structure)
    assert B.join(structure, first, second) is True

    first, second = _carbons(structure)
    joined = math.dist(structure.coords[first], structure.coords[second])
    assert abs(joined - 1.52) < 0.05, joined
    assert _formula(structure) == {'C': 2, 'H': 6}
    for hydrogen in structure.hydrogens_on(second):
        length = math.dist(structure.coords[second], structure.coords[hydrogen])
        assert abs(length - 1.07) < 0.02, length
    assert mff.relax_xyz(B.to_xyz(structure), method='uff')['ok']


def test_closing_a_ring_moves_nothing():
    """The two ends are already connected, so there is no fragment to move
    without breaking what holds them. The field closes it instead."""
    structure = _empty()
    B.place_atom(structure, 'C', (0.0, 0.0, 0.0))
    for _ in range(3):
        B.grow_from(structure, _carbons(structure)[-1], 'C')
    ends = (_carbons(structure)[0], _carbons(structure)[-1])
    carbons_before = [structure.coords[i] for i in _carbons(structure)]
    assert B.join(structure, *ends) is True
    carbons_after = [structure.coords[i] for i in _carbons(structure)]
    for before, after in zip(carbons_before, carbons_after):
        assert math.dist(before, after) < 1e-9
    assert _formula(structure) == {'C': 4, 'H': 8}


def test_tapping_an_atom_fills_a_valence_that_is_short():
    """A tap is also how a half-built centre gets finished. Returning early
    when the element was unchanged meant tapping a carbon with C selected did
    nothing at all."""
    structure = B.Structure(['C'], [(0.0, 0.0, 0.0)], {})
    assert B.set_element(structure, 0, 'C') is True
    assert _formula(structure) == {'C': 1, 'H': 4}
    # And a second tap on a full one changes nothing.
    assert B.set_element(structure, 0, 'C') is False


def test_joining_does_not_leave_atoms_on_top_of_each_other():
    """Bringing two fragments together along the bond axis is right for the
    bond and says nothing about anything else: the hydrogens on either side
    arrive wherever the rotation about that axis happens to leave them.

    Two propanes joined end to end left a pair of hydrogens 0.61 A apart, and
    perception then read one of them as bonded to two atoms -- the doubly
    bonded hydrogen that was reported.

    Pushing atoms out of the way is the wrong remedy, because it strains
    geometries that were correct. Turning is the free coordinate: both
    fragments stay rigid and internally perfect, and the only thing that
    changes is the one thing nothing had decided. The closest contact goes
    from 0.61 A to about 1.3."""
    import itertools

    def propane_pair(separation):
        structure = _empty()
        B.place_atom(structure, 'C', (0.0, 0.0, 0.0))
        B.grow_from(structure, 0, 'C')
        B.grow_from(structure, _carbons(structure)[-1], 'C')
        anchor = len(structure)
        B.place_atom(structure, 'C', (separation, 0.3, 0.2))
        B.grow_from(structure, anchor, 'C')
        B.grow_from(structure, _carbons(structure)[-1], 'C')
        return structure

    for separation in (6.0, 3.0, 2.2):
        structure = propane_pair(separation)
        carbons = _carbons(structure)
        assert B.join(structure, carbons[0], carbons[3]) is True

        perceived = mff.perceive_molecule(B.to_xyz(structure))
        adjacency = perceived.neighbours()
        stray = [i for i, s in enumerate(perceived.symbols)
                 if s == 'H' and len(adjacency[i]) != 1]
        assert not stray, (separation, stray)
        closest = min(
            math.dist(structure.coords[i], structure.coords[j])
            for i, j in itertools.combinations(range(len(structure)), 2)
            if not structure.order(i, j)
        )
        assert closest > 1.1, (separation, closest)
        assert mff.relax_xyz(B.to_xyz(structure), method='uff')['ok']
