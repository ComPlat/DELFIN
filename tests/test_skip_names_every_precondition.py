"""A skip names every missing precondition, on any machine.

Several test modules need two independent things at once: the optional
``opi`` package and an ORCA binary. A module-level
``pytest.importorskip("opi")`` fires at collection, so the whole file came
out as ONE skip naming the package, and the ORCA reason -- the only reason
there is on a machine that has ``opi`` -- was never shown.

That masking is a decision, not an observed outcome, and a file that can
only skip on this host proves nothing about it. This file tests the
decision instead: it drives ``tests/_skip_reasons.skip_listing_all``
directly with each combination of preconditions and reads the reason out of
the ``Skipped`` exception.

Input: one reason per precondition, ``None`` for the ones in place -- all
four combinations.
Output: whether a skip is raised and which reasons its message names.

Both directions are simulated, so nothing here depends on what this host
has installed: a reason string stands for an absent precondition, ``None``
for a present one. The helper takes strings and does the joining, so no
package is imported to decide the question.

Semantics pinned: a missing precondition is named, EVERY missing one is
named, and nothing is skipped when all are in place. The rejected
alternative -- checking in order and returning at the first -- is pinned
too: it passes the first three cases and fails the fourth.
"""
from __future__ import annotations

import pytest

from _skip_reasons import absent_reason, skip_listing_all

OPI = "opi is not importable (ModuleNotFoundError: No module named 'opi')"
ORCA = "no ORCA binary for the parser probe ($ORCA_BINARY unset, ...)"


def _reason_for(opi: bool, orca: bool) -> str | None:
    """The skip reason for one combination, or None when nothing skips.

    pytest.skip raises through BaseException, not Exception: caught as
    ``except Exception`` it propagates and SKIPS this test instead of being
    read -- and a skipped test asserts nothing.
    """
    try:
        skip_listing_all([OPI if opi else None, ORCA if orca else None])
    except pytest.skip.Exception as exc:
        return str(getattr(exc, "msg", exc))
    return None


def test_nothing_skips_when_both_are_in_place():
    assert _reason_for(opi=False, orca=False) is None


def test_a_missing_package_is_named():
    why = _reason_for(opi=True, orca=False)
    assert why and "opi" in why
    assert "ORCA" not in why, "it is there; do not say it is not"


def test_a_missing_binary_is_named():
    """The case the module-level skip hid on a machine that has opi."""
    why = _reason_for(opi=False, orca=True)
    assert why and "ORCA" in why
    assert "opi" not in why, "it is importable; do not say it is not"


def test_both_missing_are_named_together():
    """The one that fails if the two are checked in order instead.

    Ordering passes the three cases above and fails here: whichever is
    checked first returns its own reason and hides the other.
    """
    why = _reason_for(opi=True, orca=True)
    assert why
    assert "opi" in why
    assert "ORCA" in why


def test_the_order_of_the_reasons_does_not_decide_the_answer():
    """Neither precondition may be an early return: swap and read again."""
    forward = _reason_for(opi=True, orca=True)
    try:
        skip_listing_all([ORCA, OPI])
    except pytest.skip.Exception as exc:
        reverse = str(getattr(exc, "msg", exc))
    assert forward and reverse
    assert "opi" in reverse and "ORCA" in reverse


def test_a_package_that_is_really_absent_is_named_by_its_own_error():
    """``absent_reason`` reports the failure, not a guess about it.

    A name no distribution provides, so the answer is the same on every
    machine: the reason has to carry the exception's text, because "not
    installed" and "installed but broken" share this symptom.
    """
    reason = absent_reason("delfin_no_such_optional_package_zz9")
    assert reason is not None
    assert "delfin_no_such_optional_package_zz9" in reason
    assert "ModuleNotFoundError" in reason


def test_a_package_that_is_present_has_no_reason():
    assert absent_reason("json") is None
