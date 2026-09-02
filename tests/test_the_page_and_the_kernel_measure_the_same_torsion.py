"""A torsion the page reports and one the kernel holds have to be one number.

The page measures a dihedral in `dihedralV` and the kernel measures it in
`gfn_optimize._dihedral`.  They disagreed in sign: four points that the kernel
and RDKit both put at -81.870 degrees were shown as +81.870.

That is not only the measure box.  `updateInternalReadout` writes the page's
number into the toolbar's value field, and that field is what Set, Hold and
arming a scan all read -- so "hold this torsion where it is" asked the kernel
for the mirror image of where it was, and a scan armed from a selection walked
away from the structure on screen rather than through it.

Checked by transcribing the page's arithmetic into Python and running it on
the same points, because the two live on opposite sides of the browser and
nothing else compares them.  The transcription is line for line: if it drifts
from the JS this test stops meaning anything, so it asserts the shape of the
JS as well.
"""
from __future__ import annotations

import math
import pathlib
import re

import pytest

from delfin.dashboard import gfn_optimize as gfn

_VIEWER = (pathlib.Path(__file__).resolve().parents[1]
           / 'delfin' / 'dashboard' / 'molecule_viewer.py').read_text()

#: Four points with an unmistakable handedness -- nothing near 0 or 180, where
#: a sign error is invisible.
_FOUR = """4
four points, a clearly handed torsion
C  1.000000  1.000000  0.000000
C  0.000000  0.000000  0.000000
C  0.000000  0.000000  1.500000
C  0.800000 -0.600000  2.000000
"""


def _points(text):
    flat = gfn.coordinates_of(text)
    return [tuple(flat[3 * i:3 * i + 3]) for i in range(len(flat) // 3)]


def _cross(u, v):
    return [u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0]]


def _as_the_page_measures_it(a, b, c, d):
    """dihedralV, transcribed line for line."""
    b1 = [b[i] - a[i] for i in range(3)]
    b2 = [c[i] - b[i] for i in range(3)]
    b3 = [d[i] - c[i] for i in range(3)]
    nb2 = math.sqrt(sum(v * v for v in b2))
    if nb2 < 1e-9:
        return 0.0
    b2n = [v / nb2 for v in b2]
    n1, n2 = _cross(b1, b2), _cross(b2, b3)
    m1 = _cross(b2n, n1)
    x = sum(p * q for p, q in zip(n1, n2))
    y = sum(p * q for p, q in zip(m1, n2))
    return math.degrees(math.atan2(y, x))


def test_the_page_and_the_kernel_agree_about_a_torsion():
    """Including its sign, which is the whole of what was wrong."""
    at = _points(_FOUR)
    page = _as_the_page_measures_it(*at)
    kernel = gfn._dihedral(at, 0, 1, 2, 3)
    assert page == pytest.approx(kernel, abs=1e-6), (page, kernel)
    # And it is a torsion with a side to be on, so the test could see a flip.
    assert abs(page) > 30.0, page


def test_and_so_does_rdkit_which_is_what_the_rest_of_the_world_uses():
    """A convention two of three agree on is still a convention; three is not."""
    Chem = pytest.importorskip('rdkit.Chem')
    from rdkit.Chem import rdMolTransforms

    mol = Chem.MolFromXYZBlock(_FOUR)
    assert mol is not None
    said = rdMolTransforms.GetDihedralDeg(mol.GetConformer(), 0, 1, 2, 3)
    assert said == pytest.approx(gfn._dihedral(_points(_FOUR), 0, 1, 2, 3),
                                 abs=1e-3)


def test_the_transcription_above_still_matches_the_page():
    """This test is worth nothing if the JS moves and the copy does not.

    So the frame the page builds is asserted here in its own words: the cross
    product that decides the sign is b2n x n1, and not the other way round.
    """
    body = _VIEWER.split('function dihedralV(')[1].split('\n    }')[0]
    assert 'var m1x = b2ny*n1z - b2nz*n1y;' in body, body
    assert 'var m1y = b2nz*n1x - b2nx*n1z;' in body
    assert 'var m1z = b2nx*n1y - b2ny*n1x;' in body
    assert 'Math.atan2(y, x)' in body
    # And it is the one the toolbar's value field is filled from, which is
    # what makes the sign a held value rather than a decoration.
    readout = _VIEWER.split('function updateInternalReadout(')[1]
    readout = readout.split('\n    function ')[0]
    assert 'readInternal(scopeKey)' in readout
    assert re.search(r"submit-internal-value", readout)


# ---------------------------------------------------------------------------
# And the same atoms, in the same order
# ---------------------------------------------------------------------------

class _Atom:
    def __init__(self, serial):
        self.serial = serial


def _as_the_page_reports_picks(atoms, serials):
    """ffIndicesOf, transcribed line for line."""
    where = {}
    for i, one in enumerate(atoms):
        if one.serial not in where:
            where[one.serial] = i
    return [where[s] for s in serials if s in where]


def test_the_selection_reaches_the_kernel_in_the_order_it_was_clicked():
    """A value measured over one chain was written onto another.

    Everything the user can see honours the order the atoms were clicked: the
    highlight colours are indexed by it, the measure box prints the chain
    [H5 -> C1 -> C0 -> H2], and the toolbar's value is the torsion over
    exactly that chain.  What was sent to the kernel walked the model from 0
    upwards instead, so Python was told 0,1,2,5.

    Hold then paired that sorted list with the value measured over the clicked
    one: the row read "Holding C0-C1-H2-H5 = -59.9 deg" and xtb's input got
    "dihedral: 1, 2, 3, 6, -59.900000".  For two atoms the order makes no
    difference; for three it moves the vertex of the angle; for four it is a
    different torsion outright.
    """
    atoms = [_Atom(i + 1) for i in range(8)]      # ethane, serials 1..8
    clicked = [6, 2, 1, 3]                        # H5, C1, C0, H2
    assert _as_the_page_reports_picks(atoms, clicked) == [5, 1, 0, 2]
    # And that is not the ascending order, or the test would pass either way.
    assert sorted(_as_the_page_reports_picks(atoms, clicked)) == [0, 1, 2, 5]


def test_a_serial_the_model_does_not_have_is_left_out():
    """Rather than becoming an index into somebody else's atom."""
    atoms = [_Atom(i + 1) for i in range(3)]
    assert _as_the_page_reports_picks(atoms, [3, 99, 1]) == [2, 0]


def test_the_transcription_of_the_picks_still_matches_the_page():
    """Same rule as the torsion above: the copy has to be tied to the JS."""
    body = _VIEWER.split('function ffIndicesOf(')[1].split('\n    }')[0]
    assert 'for (var s = 0; s < serials.length; s++)' in body, body
    assert 'out.push(where[serials[s]])' in body
    # The old shape walked the model and asked whether each atom was picked,
    # which is where the order was lost.
    assert 'serials.indexOf(atoms[i].serial)' not in body
