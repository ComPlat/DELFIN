"""Skip reasons that name every missing precondition, not only the first.

A module-level ``pytest.importorskip`` fires at collection, so a file with
two independent preconditions reports only the one that was checked first.
The per-test reason is then never shown, on any machine -- including one
that satisfies the module-level condition, where it is the only reason
there is.

Input: a reason per precondition, ``None`` for the ones that are in place.
Output: nothing when none is missing, else a ``pytest.skip`` whose reason
names every missing one, joined with ``"; "``.

Semantics pinned here: a missing precondition is named, *every* missing one
is named, and nothing is skipped when all are in place.

Rejected alternative: checking them in order with an early return. That
gives the same answer for every case except the one that matters -- two
missing -- and there it hides one behind the other, which is the masking
this module removes. ``tests/test_skip_names_every_precondition.py`` pins
that difference as a failing case.

Kept in one place rather than copied into each file: eleven test modules
need this construct, and an eleventh copy is an eleventh place to drift.
"""
from __future__ import annotations

from typing import Iterable, Optional

import pytest


def absent_reason(module: str) -> Optional[str]:
    """``None`` if ``module`` imports, else a reason naming the failure.

    The reason carries the exception's own text: "not installed" and
    "installed but broken" are different problems with the same symptom,
    and the reader needs to tell them apart.
    """
    try:
        __import__(module)
    except Exception as exc:                                    # noqa: BLE001
        return f"{module} is not importable ({type(exc).__name__}: {exc})"
    return None


def orca_reason() -> Optional[str]:
    """``None`` when an ORCA binary will answer a probe, else the reason.

    Same lookup the probe uses (``find_orca``: ``$ORCA_BINARY``, then
    ``orca`` on ``PATH``, then the newest ``/opt/orca*/orca``), so this
    skip and an ``unavailable`` verdict cannot disagree about what counts
    as an ORCA.
    """
    from delfin.orca_knowledge.probe import find_orca

    if not find_orca():
        return ("no ORCA binary for the parser probe ($ORCA_BINARY unset, no "
                "'orca' on PATH, no /opt/orca*/orca)")
    return None


def skip_listing_all(reasons: Iterable[Optional[str]]) -> None:
    """Skip naming every reason that is not ``None``; return if none is.

    Raises ``pytest.skip``, which goes through ``BaseException``: a caller
    that means to READ the reason must catch ``pytest.skip.Exception``, not
    ``Exception``, or it is skipped instead of asserting.
    """
    missing = [reason for reason in reasons if reason]
    if missing:
        pytest.skip("; ".join(missing))
