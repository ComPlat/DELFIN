"""A skip names every missing precondition, on any machine.

Several test modules need two independent things at once: the optional
``opi`` package and an ORCA binary. A module-level
``pytest.importorskip("opi")`` fires at collection, so the whole file came
out as ONE skip naming the package, and the ORCA reason -- the only reason
there is on a machine that has ``opi`` -- was never shown.

That masking is a decision, not an observed outcome, and a file that can
only skip on this host proves nothing about it. This file tests the
decision instead: it drives ``tests/_skip_reasons.missing`` directly with
each combination of preconditions and reads what it returns.

Input: one reason per precondition, ``None`` for the ones in place -- all
four combinations.
Output: ``None`` when nothing is missing, else a message naming every
missing one.

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

from _skip_reasons import absent_reason, missing

OPI = "opi is not importable (ModuleNotFoundError: No module named 'opi')"
ORCA = "no ORCA binary for the parser probe ($ORCA_BINARY unset, ...)"


def test_nothing_is_missing_when_both_are_in_place():
    assert missing(None, None) is None


def test_a_missing_package_is_named():
    why = missing(OPI, None)
    assert why and "opi" in why
    assert "ORCA" not in why, "it is there; do not say it is not"


def test_a_missing_binary_is_named():
    """The case the module-level skip hid on a machine that has opi."""
    why = missing(None, ORCA)
    assert why and "ORCA" in why
    assert "opi" not in why, "it is importable; do not say it is not"


def test_both_missing_are_named_together():
    """The one that fails if the two are checked in order instead.

    Ordering passes the three cases above and fails here: whichever is
    checked first returns its own reason and hides the other.
    """
    why = missing(OPI, ORCA)
    assert why
    assert "opi" in why
    assert "ORCA" in why


def test_the_order_of_the_reasons_does_not_decide_the_answer():
    """Neither precondition may be an early return: swap and read again."""
    forward = missing(OPI, ORCA)
    reverse = missing(ORCA, OPI)
    assert forward and reverse
    assert "opi" in reverse and "ORCA" in reverse


def test_the_skip_is_raised_where_the_reader_can_see_the_test():
    """The reason is returned, and the test raises the skip itself.

    Measured with the skip raised from inside a helper: pytest reported
    ``SKIPPED [1] tests/_skip_reasons.py:71``, naming the helper instead of
    the test. A skip reason that cannot be attributed to a test is half of
    the problem this fixes, so the helper returns and the test raises.

    pytest.skip raises through BaseException, not Exception: caught as
    ``except Exception`` it propagates and SKIPS this test instead of being
    read -- and a skipped test asserts nothing.
    """
    why = missing(OPI, ORCA)
    with pytest.raises(pytest.skip.Exception) as caught:
        if why:
            pytest.skip(why)
    assert str(caught.value.msg) == why


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


def test_a_package_this_host_may_lack_is_measured_not_assumed():
    """``absent_reason`` answers about the package, never about the host.

    Both directions are read here -- whichever one this host happens to be
    in, the answer has to agree with a fresh import of the same name. A
    hard-coded "opi is absent" would pass on CI and lie on a machine that
    has it, which is the mistake this whole change is about.
    """
    import importlib.util

    installed = importlib.util.find_spec("opi") is not None
    assert (absent_reason("opi") is None) is installed
