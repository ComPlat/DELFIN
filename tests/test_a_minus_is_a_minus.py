"""Five correct answers scored 0/5, sigma 0.9 -- the instrument's signature.

GLM writes the typographic minus before every negative energy:
"\u2212613.4171289012 Eh". The audit task's rubric asks for the figure as
"-613\\.417", with the ASCII hyphen-minus, and the number tokeniser
reads signs the same way. So the model that had said everything the
task wanted, with the right number, missed expected[1] five times out
of five (2026-09-11). A minus is a minus: every haystack and every
number is read with the dashes and minus signs a model might write
folded to "-" first.
"""

from __future__ import annotations

from delfin.agent import benchmark as B


def _traj(text: str) -> B.Trajectory:
    return B.Trajectory(text=text)


def test_a_typographic_minus_matches_an_ascii_pattern():
    sig = B.Signal(pattern=r"-613\.417", against="text")
    assert B._signal_matches(sig, _traj("Finaler Single-Point: **\u2212613.4171289012 Eh**"))
    assert B._signal_matches(sig, _traj("SPE = -613.417 Eh"))
    assert not B._signal_matches(sig, _traj("SPE = 613.417 Eh"))


def test_the_match_is_found_for_evidence_too():
    sig = B.Signal(pattern=r"-613\.417", against="text")
    assert B._signal_match(sig, _traj("E = \u2212613.4171 Eh")) is not None


def test_a_number_with_a_typographic_minus_is_read_as_negative():
    assert -113.562 in B.numbers_in("das Minimum liegt bei \u2212113.562 Eh")
    assert -113.562 in B.numbers_in("Minimum: \u2013113.562")      # an en dash, as pasted text has


def test_every_dash_a_model_writes_folds_to_one():
    for sign in ("\u2212", "\u2010", "\u2011", "\u2012", "\u2013", "\u2014", "\u2015", "\ufe63", "\uff0d"):
        assert B.ascii_minus(f"{sign}1") == "-1"
    assert B.ascii_minus("plain-text") == "plain-text"
