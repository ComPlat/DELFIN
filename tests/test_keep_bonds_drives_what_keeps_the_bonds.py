"""Under Keep bonds the hand drives what keeps the bonds, never a bond.

From a user's own session (GFN2 in DMSO, 42 atoms, the hand a pull at 0.65
of a bond, Keep bonds on): a hydrogen drawn off its carbon was driven *on its
own bond* -- the perception named the C-H as the coordinate the hand had
changed most, which it had -- the push took it apart at 72 kcal/mol/A, the
wall read the bonding off the answer and refused it, and the picture stood
still through twenty-seven answers under a line saying the bonding was being
kept.  Kept it was.  Nothing else happened, and nothing said why.

A bond is the one coordinate Keep bonds exists to leave alone, so it is not
on offer while the switch is on: what the hand can drive is a torsion or an
angle, which moves an atom without moving it off its neighbours, and a
terminal atom is given the torsion about the bond one further in.  A gesture
that asks for nothing but a bond to stretch is answered with nothing, and the
line says which bond and what to do instead.
"""
from __future__ import annotations

import math

import pytest

from delfin.dashboard import climb as _climb
from delfin.dashboard import gfn_optimize as gfn

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

#: H7 hangs on C1.  Along its bond is (-1.01, 0, -0.4) normalised; across it
#: is the y direction.
_H = 7
_C = 1


def _moved(xyz, index, by):
    from delfin.dashboard.structure_editor import xyz_line
    here = gfn.coordinates_of(xyz)
    rows = []
    for i, line in enumerate(gfn.atom_lines(xyz)):
        x, y, z = here[3 * i], here[3 * i + 1], here[3 * i + 2]
        if i == index:
            x, y, z = x + by[0], y + by[1], z + by[2]
        rows.append(xyz_line(line.split()[0], x, y, z))
    return f'{len(rows)}\nDELFIN drag-follow held={index}\n' + '\n'.join(rows)


def _along_the_bond(scale):
    here = gfn.coordinates_of(_ETHANE)
    axis = [here[3 * _H + n] - here[3 * _C + n] for n in range(3)]
    length = math.sqrt(sum(one * one for one in axis))
    return [scale * one / length for one in axis]


def test_a_radial_pull_has_nothing_to_drive_and_a_sideways_one_a_torsion():
    graph = gfn.bond_graph(_ETHANE)
    outward = _moved(_ETHANE, _H, _along_the_bond(0.6))
    # Without Keep bonds the bond is the coordinate, as it always was.
    plain = gfn.contacts_holding(outward, [_H], was=_ETHANE)
    assert plain and plain[0]['kind'] == 'distance'
    assert sorted(plain[0]['atoms']) == [_C, _H]
    # With it, a bond is not on offer.  What is left to this gesture is the
    # H...C distance across the carbon -- the angle in another currency --
    # which moves the hydrogen without moving it off its carbon.
    kept = gfn.contacts_holding(outward, [_H], was=_ETHANE, keeping=graph)
    assert not any(one['kind'] == 'distance'
                   and tuple(sorted(one['atoms'])) in graph for one in kept)
    # Across the bond the torsion about C0-C1 carries the hydrogen.
    sideways = _moved(_ETHANE, _H, [0.0, 0.6, 0.0])
    turned = gfn.contacts_holding(sideways, [_H], was=_ETHANE, keeping=graph)
    assert turned, 'a sideways pull under Keep bonds drives nothing'
    assert turned[0]['kind'] == 'dihedral'
    assert turned[0]['atoms'][0] == _H and turned[0]['atoms'][1] == _C
    assert all(not (one['kind'] == 'distance'
                    and tuple(sorted(one['atoms'])) in graph) for one in turned)
    # The opening question, asked the same way, gives the same torsion.
    opening = gfn.contacts_holding(sideways, [_H], was=_ETHANE, opening=True,
                                   keeping=graph)
    assert opening and opening[0]['kind'] == 'dihedral'


def test_a_terminal_atom_is_given_its_torsion_only_when_asked():
    from delfin.atom_mapping import cov_radius
    rows = [line.split() for line in gfn.atom_lines(_ETHANE)]
    where = [(float(r[1]), float(r[2]), float(r[3])) for r in rows]
    radius = [cov_radius(r[0]) for r in rows]
    across = [0.0, 1.0, 0.0]
    assert gfn.turn_for(where, radius, [_H], {_H}, pulled=across) == []
    turn = gfn.turn_for(where, radius, [_H], {_H}, pulled=across, terminal=True)
    assert turn and turn[0]['kind'] == 'dihedral'
    assert turn[0]['atoms'][:2] == [_H, _C]


def test_the_stretched_bond_is_named_for_the_line():
    outward = _moved(_ETHANE, _H, _along_the_bond(0.6))
    said = gfn.stretched_bond(_ETHANE, outward, [_H], gfn.bond_graph(_ETHANE))
    assert said is not None
    assert said['atoms'] == (_C, _H)
    assert said['names'] == 'C2-H8'
    assert said['now'] - said['was'] == pytest.approx(0.6, abs=1e-3)
    sideways = _moved(_ETHANE, _H, [0.0, 0.01, 0.0])
    assert gfn.stretched_bond(_ETHANE, sideways, [_H]) is None


def _a_part(structure):
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import test_the_budget_prices_a_relaxed_path as budget
    return budget, budget._a_part(structure)


def _drag(part, budget, held, cursor, answers, reach=1.0):
    """The page's own loop: the last answer with the held atom moved up to a
    reach towards a standing cursor."""
    from delfin.dashboard.structure_editor import xyz_line
    answer = part._current_xyz()
    for n in range(answers):
        here = gfn.coordinates_of(answer)
        rows = []
        for i, line in enumerate(gfn.atom_lines(answer)):
            x, y, z = here[3 * i], here[3 * i + 1], here[3 * i + 2]
            if i == held:
                want = [cursor[0] - x, cursor[1] - y, cursor[2] - z]
                far = math.sqrt(sum(one * one for one in want)) or 1e-9
                step = min(reach, far) / far
                x, y, z = x + want[0] * step, y + want[1] * step, z + want[2] * step
            rows.append(xyz_line(line.split()[0], x, y, z))
        part.submit_manip_sync.value = (
            f'{len(rows)}\nDELFIN drag-follow held={held} f{n}\n' + '\n'.join(rows))
        budget._quiet(part.state)
        answer = part.state.get('thermal_was') or answer
    return answer


def _length(xyz, i, j):
    here = gfn.coordinates_of(xyz)
    return math.dist(here[3 * i:3 * i + 3], here[3 * j:3 * j + 3])


@_needs_xtb
@pytest.mark.parametrize('hand', ['pull', 'move'])
def test_under_keep_bonds_the_hydrogen_swings_and_never_leaves(hand):
    budget, part = _a_part(_ETHANE)
    part.submit_ff_dd.value = 'gfnff'
    part.submit_relax_btn.value = True
    part.submit_hand_dd.value = hand
    part.submit_pull_slider.value = 0.65
    part.submit_topology_btn.value = True
    assert part._begin_gfn_follow()
    here = gfn.coordinates_of(_ETHANE)
    at = here[3 * _H:3 * _H + 3]
    # Straight out along the bond, four angstrom away.
    out = _along_the_bond(4.0)
    answer = _drag(part, budget, _H, [at[n] + out[n] for n in range(3)], 6)
    assert _length(answer, _C, _H) < 1.25, 'the bond was stretched'
    assert gfn.graph_holds(gfn.bond_graph(_ETHANE), answer)[0]
    said = ' '.join(part.state.get('mol_status_lines') or ())
    assert 'bonds kept' in said, said
    # Across it: the torsion turns and every bond stays.
    was = gfn._dihedral([tuple(here[3 * k:3 * k + 3]) for k in range(8)],
                        _H, _C, 0, 2)
    answer = _drag(part, budget, _H, [at[0], at[1] + 4.0, at[2]], 8)
    now_here = gfn.coordinates_of(answer)
    now = gfn._dihedral([tuple(now_here[3 * k:3 * k + 3]) for k in range(8)],
                        _H, _C, 0, 2)
    assert gfn._turned_by(now, was) > 15.0, (was, now)
    assert gfn.graph_holds(gfn.bond_graph(_ETHANE), answer)[0]
    assert _length(answer, _C, _H) < 1.25
    said = ' '.join(part.state.get('mol_status_lines') or ())
    assert 'bonds kept' in said and 'only stretches' not in said, said
    part.submit_relax_btn.value = False


_H2 = """2
hydrogen
H  0.000000  0.000000  0.370000
H  0.000000  0.000000 -0.370000
"""


@_needs_xtb
def test_a_gesture_that_can_only_stretch_a_bond_is_told_so():
    """A diatomic has one coordinate, and under Keep bonds it is not on
    offer.  The line has to say which bond and what to do instead, because
    on screen this is a drag that does nothing."""
    budget, part = _a_part(_H2)
    part.submit_ff_dd.value = 'gfnff'
    part.submit_relax_btn.value = True
    part.submit_hand_dd.value = 'pull'
    part.submit_topology_btn.value = True
    assert part._begin_gfn_follow()
    answer = _drag(part, budget, 1, [0.0, 0.0, -4.0], 4)
    assert _length(answer, 0, 1) < 1.0
    said = ' '.join(part.state.get('mol_status_lines') or ())
    assert 'bonds kept' in said, said
    assert 'only stretches H1-H2' in said, said
    assert 'Keep bonds off' in said
    part.submit_relax_btn.value = False
