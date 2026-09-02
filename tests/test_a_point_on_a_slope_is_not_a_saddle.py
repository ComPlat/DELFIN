"""What a Hessian says is a verdict only where the structure stopped moving.

``climb_to_saddle`` knows whether the gradient vanished -- it says so itself,
"The climb did not converge in N steps" -- and ``reach_the_reaction`` replaced
that sentence with one built from the mode count alone.  All three of its
branches assert arrival: "The climb settled somewhere with no mode going the
wrong way, so it is a minimum", "The climb reached a saddle at X cm-1", "it
converged onto a minimum".  A curvature spectrum says nothing about being at a
stationary point, so a geometry the step ceiling cut off was described as
having settled -- and, when its modes happened to come out right, written into
the box and filed as a transition state with a frequency quoted to the cm-1.

Measured on a cyclohexene whose C0-C1 was driven to 1.90 A and then climbed:
the press said "it converged onto a different saddle, at -205 cm-1" while the
gradient at that geometry was 3.3e-03 Hartree per Bohr per coordinate -- and
the editor's own slope test, handed the same structure, answered "this
structure is neither a minimum nor a saddle".  One codebase, two answers.

The sentences are checked here rather than the climb, because the climb is
several seconds of xtb and the wording is the whole of what was wrong.
"""
from __future__ import annotations

import pathlib

import pytest

from delfin.dashboard import climb

_CLIMB = (pathlib.Path(__file__).resolve().parents[1]
          / 'delfin' / 'dashboard' / 'climb.py').read_text()
_EDITOR = (pathlib.Path(__file__).resolve().parents[1]
           / 'delfin' / 'dashboard' / 'structure_editor.py').read_text()

#: A landing with one mode going the wrong way, at a point on a slope.
_ON_A_SLOPE = {'count': 1, 'modes': [-205.0], 'lowest': -205.0, 'ok': True,
               'share': 0.9, 'reaction': True, 'still': False,
               'gmax': 2.3e-2, 'grms': 3.3e-3}
#: The same landing, but standing still.
_STILL = dict(_ON_A_SLOPE, still=True, gmax=1.0e-4, grms=2.0e-5)

_CLAIMS = ('converged', 'settled', 'reached a saddle', 'is a minimum')


def test_a_landing_on_a_slope_is_never_called_converged():
    """Whichever of the two phrases the sentence is built from."""
    for said in (climb._named_landing(_ON_A_SLOPE),
                 climb._why_not(_ON_A_SLOPE, 'C0-C1')):
        assert not any(one in said for one in _CLAIMS), said
        assert 'out of tries' in said, said
        # The modes are still reported -- they are what the walk was
        # following, and hiding them leaves nothing to read.
        assert '-205' in said, said


def test_and_the_gradient_it_stopped_at_is_in_the_sentence():
    """A number the user can act on, against the one it is measured against."""
    said = climb._why_not(_ON_A_SLOPE, 'C0-C1')
    assert '2.3e-02' in said, said
    assert f'{climb.GRADIENT_MAX:.0e}' in said, said


def test_and_the_second_place_it_appears_in_says_it_shorter():
    """The ladder prefixes every later line with why the first rung missed.

    Said at length twice, the press ends up as one line saying one thing
    twice -- which is what it did the first time this was fixed.
    """
    long_form = climb._why_not(_ON_A_SLOPE, 'C0-C1')
    short = climb._named_landing(_ON_A_SLOPE)
    assert len(short) < len(long_form) / 2, (short, long_form)
    assert '-205' in short, short


def test_a_climb_that_ran_out_of_tries_says_so_and_says_what_to_do():
    said = climb._why_not(_ON_A_SLOPE, 'C0-C1')
    assert 'did not converge' in said, said
    assert 'Pressing again' in said, said
    # And not when the ladder has rungs left, or every later line repeats it.
    quiet = climb._why_not(_ON_A_SLOPE, 'C0-C1', advise=False)
    assert 'Pressing again' not in quiet, quiet
    assert not any(one in said
                   for one in ('settled', 'reached a saddle')), said


def test_a_landing_that_did_stand_still_reads_exactly_as_it_always_did():
    """The wording is only touched where the claim was false."""
    assert climb._named_landing(_STILL) == ('it converged onto a different '
                                            'saddle, at -205 cm-1')
    assert climb._why_not(dict(_STILL, count=0), 'C0-C1') == (
        'The climb settled somewhere with no mode going the wrong way, so it '
        'is a minimum and not a transition state.')


def test_and_so_does_one_from_a_producer_that_does_not_measure_it():
    """``still`` absent, or None, is unknown -- which is not the same as
    moving, and must not turn a working sentence into a refusal."""
    for shape in ({'count': 0, 'modes': []},
                  {'count': 0, 'modes': [], 'still': None}):
        assert climb._named_landing(shape) == 'it converged onto a minimum'


def test_the_verdict_measures_whether_anything_was_standing_still():
    """Off the gradient the walk already has, and against its own threshold."""
    body = _CLIMB.split('    def verdict(')[1].split('\n    def ')[0]
    assert "'still': still" in body
    assert 'gmax < GRADIENT_MAX and grms < GRADIENT_RMS' in body
    # None where there is nothing to read: a verdict asked of a structure this
    # walk has not stepped on.
    assert 'still = gmax = grms = None' in body
    assert 'if self.gradient is not None:' in body


def test_the_climb_keeps_its_own_answer_where_the_ladder_cannot_overwrite_it():
    """``ok`` is the ladder's to rewrite; whether it converged is not."""
    body = _CLIMB.split('def climb_to_saddle(')[1].split('\ndef ')[0]
    assert "'ok': arrived, 'converged': arrived" in body


@pytest.mark.parametrize('rung', [
    "landed = reached.get('still') is not False",
    "softest.get('still') is not False and softest.get('reaction')",
    "orca.get('still') is not False and orca.get('reaction')",
])
def test_no_rung_of_the_ladder_accepts_a_landing_on_a_slope(rung):
    """Rung one, rung two and ORCA, each asked the same question."""
    body = _CLIMB.split('def reach_the_reaction(')[1]
    assert rung in body, rung


def test_the_editor_does_not_file_a_point_its_own_search_refused():
    """Two refusals, both of them the search's own.

    :func:`saddle.verdict` applies autodE's -40 cm-1 floor and renames
    anything shallower; the press printed that refusal, wrote it into the box,
    and filed the same geometry as a transition state anyway -- which is the
    one of the three the user meets again later.  And a Hessian taken where
    nothing was standing still is not a verdict about a saddle at all.
    """
    body = _EDITOR.split('def _note_the_saddle(')[1].split('\n    def ')[0]
    assert "if first_order is False or shape.get('still') is False:" in body
    assert 'order = 0' in body
    # And the two presses that have the verdict in hand pass it.
    assert _EDITOR.count("first_order=said['first_order'])") == 2
