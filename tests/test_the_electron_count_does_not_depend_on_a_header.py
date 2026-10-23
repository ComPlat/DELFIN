"""The parity guard is handed both shapes of structure, and must read both.

``electron_parity`` skipped two lines before counting, on the assumption that
every structure arrives as a whole XYZ document.  The Optimise button hands
over the coordinate box, which is one; Optimise all hands over the tab's
frames, which are bare atom lines.  For those, the two lines skipped were the
first two *atoms*, so the parity inverted whenever their atomic numbers summed
odd -- and both directions of the guard then failed on the same molecule: a
legal run refused with arithmetic that is wrong, an illegal one run and its
multiplicity asserted on the status line.

Water and methane are the two the finder caught, and they are exactly the two
whose first two atoms sum odd: O + H = 9, C + H = 7.
"""
from __future__ import annotations

import pytest

from delfin.dashboard import gfn_optimize as gfn

_WATER = ('O   0.000000   0.000000   0.000000\n'
          'H   0.960000   0.000000   0.000000\n'
          'H  -0.240000   0.930000   0.000000')
_METHANE = ('C   0.000000   0.000000   0.000000\n'
            'H   0.629000   0.629000   0.629000\n'
            'H  -0.629000  -0.629000   0.629000\n'
            'H  -0.629000   0.629000  -0.629000\n'
            'H   0.629000  -0.629000  -0.629000')
_AMMONIA = ('N   0.000000   0.000000   0.000000\n'
            'H   1.010000   0.000000   0.000000\n'
            'H  -0.340000   0.950000   0.000000\n'
            'H  -0.340000  -0.470000   0.820000')


def _document(body, name):
    return '%d\n%s\n%s\n' % (len(body.splitlines()), name, body)


@pytest.mark.parametrize('name,body,electrons', [
    ('water', _WATER, 10),
    ('methane', _METHANE, 10),
    ('ammonia', _AMMONIA, 10),
])
def test_a_body_and_a_document_have_the_same_electrons(name, body, electrons):
    """One molecule, one parity, however it was written down."""
    want = electrons % 2
    assert gfn.electron_parity(_document(body, name), 0) == want
    assert gfn.electron_parity(body, 0) == want


def test_and_the_charge_still_moves_it():
    """The guard is about the electron count, not about the atoms."""
    assert gfn.electron_parity(_WATER, 0) == 0
    assert gfn.electron_parity(_WATER, 1) == 1
    assert gfn.electron_parity(_document(_WATER, 'water'), 1) == 1


def test_the_two_readings_agree_on_anything_with_four_columns():
    """Because that is the one test the counting loop makes.

    A comment line that happens to have four fields is what the header skip
    was really guarding against, and :func:`gfn_optimize.atom_lines` guards
    against it by asking whether the last three are numbers.
    """
    odd = 'C 0.0 0.0 0.0\nH 1.0 0.0 0.0'          # 6 + 1 = 7, odd
    assert gfn.electron_parity(odd, 0) == 1
    assert gfn.electron_parity(_document(odd, 'a comment 1.0 2.0 3.0'), 0) == 1
