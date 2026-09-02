"""A held value may not be one the atoms cannot have.

A scan is armed to a number somebody typed and has been bounded for a while.
A drag is not: it reads its value off the page's model, and the page will
report whatever the cursor has made -- including a hydrogen pushed most of the
way into the carbon it is bonded to.  The perception passed that on as the
distance to hold, and xtb was asked to meet it.

Measured on one grab of forty-nine answers from a real session: the wish put
H36 at 0.862, 0.778 and 0.776 A from C12, where a C-H is 1.09.  The first
three answers cost +27.4, +64.7 and +65.7 kcal/mol before the hand had gone
anywhere; with the floor they read 0.909 A and +15.2, +16.5 and +16.6.

The bound is on what is *asked for*, never on what the structure may do.  A
relaxation that finds its own way to a short contact is answering a question;
a hand that demands one is asking a question with no answer.

Both floors come out of covalent radii and nothing else, so they hold for any
element and any molecule rather than for the session they were found in.
"""
from __future__ import annotations

import math
import pathlib

import pytest

from delfin.atom_mapping import cov_radius
from delfin.dashboard import gfn_optimize as gfn

_EDITOR = (pathlib.Path(__file__).resolve().parents[1]
           / 'delfin' / 'dashboard' / 'structure_editor.py').read_text()


def test_a_distance_the_hand_asks_for_has_a_floor():
    radius = [cov_radius('C'), cov_radius('H')]
    floor = gfn.NO_CLOSER_THAN * (radius[0] + radius[1])
    where = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.5)]
    assert gfn._no_closer_than(where, radius, 0, 1, 0.5) == pytest.approx(
        floor)
    # A C-H is 1.09 A, so the floor is under it and not on top of it.
    assert floor < 1.09


def test_and_an_ordinary_distance_is_handed_on_untouched():
    radius = [cov_radius('C'), cov_radius('H')]
    for span in (1.09, 1.5, 3.0):
        assert gfn._no_closer_than([(0.0, 0.0, 0.0), (0.0, 0.0, span)],
                                   radius, 0, 1, span) == pytest.approx(span)


@pytest.mark.parametrize('one,two', [
    ('C', 'C'), ('C', 'H'), ('O', 'H'), ('Pd', 'Br'), ('Cs', 'I'),
])
def test_the_floor_is_the_radii_and_nothing_else(one, two):
    """Which is what makes it the same rule for every element."""
    radius = [cov_radius(one), cov_radius(two)]
    got = gfn._no_closer_than([(0.0, 0.0, 0.0), (0.0, 0.0, 0.01)],
                              radius, 0, 1, 0.01)
    assert got == pytest.approx(gfn.NO_CLOSER_THAN * sum(radius))


def test_the_perception_puts_it_on_both_ways_a_distance_is_offered():
    """The opening answer of a drag and every one after it."""
    body = (pathlib.Path(__file__).resolve().parents[1] / 'delfin'
            / 'dashboard' / 'gfn_optimize.py').read_text()
    body = body.split('def contacts_holding(')[1].split('\ndef ')[0]
    assert body.count('_no_closer_than(where, radius, i, j') == 2


def test_an_angle_stops_where_its_outer_two_atoms_stop():
    """Zero was the bound before -- true, and no use.

    An H-C-H closed to five degrees is two hydrogens inside one another, and
    the walk went there and priced it.  Closing a bend drives the outer pair
    together, and they stop where any other pair stops; with the arms at *a*
    and *b* and that pair no closer than *d*, the law of cosines gives the
    angle exactly.
    """
    body = _EDITOR.split('def _scan_floor_for(')[1].split('\n    def ')[0]
    assert "if leg['kind'] == 'angle' and len(leg['atoms']) == 3:" in body
    assert 'math.degrees(math.acos(' in body
    assert 'return 0.0\n        if' not in body, 'zero is not a floor'
    # The same share the distance uses, so one rule and not two.
    assert body.count('_SCAN_NO_CLOSER') == 2

    # And the arithmetic it stands on, checked against the shapes it is for.
    for tag, arm, outer, want in (('methane H-C-H', 1.09, ('H', 'H'), 28),
                                  ('propane C-C-C', 1.53, ('C', 'C'), 50),
                                  ('water H-O-H', 0.96, ('H', 'H'), 32)):
        closest = 0.85 * sum(cov_radius(s) for s in outer)
        cosine = (2 * arm * arm - closest * closest) / (2 * arm * arm)
        got = math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
        assert round(got) == want, (tag, got)


def test_a_torsion_is_left_periodic():
    """289 degrees is a place a structure can be in.

    What a torsion can drive two atoms into is real and is not a bound on its
    own value -- see the grab that reached +6302 kcal/mol on one, which no
    floor on a number would have caught.
    """
    body = _EDITOR.split('def _scan_floor_for(')[1].split('\n    def ')[0]
    assert 'periodic' in body
    assert "leg['kind'] != 'distance'" in body
