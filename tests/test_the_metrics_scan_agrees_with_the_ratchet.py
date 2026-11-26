"""The metrics row's reachability scan and the ratchet are one measurement.

``delfin/agent/commit_metrics.py`` carries a copy of the scan in
``tests/test_every_capability_has_a_caller.py``, for the three reasons its
docstring states (``tests/`` is not an importable package, the product must
not depend on the suite, and importing a pytest module would drag collection
into every command that records a row). A copy is only honest while it
answers the same thing as the original, and until this file existed nothing
said so: the duplicate could have drifted in either direction with both
suites green -- the ratchet reporting one unmounted set and the row a
different one, and each agreeing with itself.

What this pins, and the input that would make each fail:

* the **definitions** both scans work from: a differing ``definitions()``
  (an extra skip path, a changed suffix list) gives different totals;
* the **unmounted set**, not just its size: a scan that missed one name and
  invented another would keep the count and still measure a different
  population, so the symmetric difference is asserted empty;
* the **tier** of each name, because ``nowhere`` and ``test-only`` are
  different findings and a copy that reported the wrong one would disagree
  with the fixture the ratchet enforces.

The counts are not compared against a number written here. A hardcoded 8
would have to be edited by whoever pays the ratchet down, which is the
moment the pin should be silent rather than in the way.
"""

from __future__ import annotations

import pytest

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _ratchet():
    """The ratchet module itself, loaded from ``tests/`` by path.

    ``tests/`` is not a package, so there is no import statement to write.
    Loading the real file rather than restating its functions is the point:
    a pin against a second copy of the original would pin nothing.
    """
    path = ROOT / "tests" / "test_every_capability_has_a_caller.py"
    spec = importlib.util.spec_from_file_location("_the_ratchet", path)
    assert spec is not None and spec.loader is not None, path
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Scans the whole repository (30-100 s): in the slow suite, so a fast
# CI shard stays within its 30-minute limit (2026-11-26).
@pytest.mark.slow
def test_both_scans_find_the_same_definitions():
    """The failing input: a scan that skips a directory the other reads."""
    from delfin.agent import commit_metrics as cm

    mine = cm.definitions(ROOT)
    theirs = _ratchet()._definitions()
    assert len(mine) == len(theirs), (len(mine), len(theirs))
    assert set(mine) == set(theirs), sorted(set(mine) ^ set(theirs))[:20]


# Scans the whole repository (30-100 s): in the slow suite, so a fast
# CI shard stays within its 30-minute limit (2026-11-26).
@pytest.mark.slow
def test_both_scans_agree_on_every_unmounted_name():
    """The failing input: one name swapped for another, count unchanged."""
    from delfin.agent import commit_metrics as cm

    ratchet = _ratchet()
    defs = ratchet._definitions()
    reached = ratchet._mentions(ROOT / "delfin", defs, True)
    theirs = {n for n in defs if reached[n] == 0}

    row = cm.reachability(ROOT)
    mine = set(row["unmounted_names"])
    assert mine == theirs, {
        "only in the row": sorted(mine - theirs)[:20],
        "only in the ratchet": sorted(theirs - mine)[:20],
    }
    # And the row's own count is the size of the set it printed, not a
    # second tally kept beside it.
    assert row["unmounted"] == len(mine)


# Scans the whole repository (30-100 s): in the slow suite, so a fast
# CI shard stays within its 30-minute limit (2026-11-26).
@pytest.mark.slow
def test_both_scans_agree_on_the_tier_of_each_name():
    """The failing input: a copy that calls a test-only name dead code."""
    from delfin.agent import commit_metrics as cm

    ratchet = _ratchet()
    defs = ratchet._definitions()
    reached = ratchet._mentions(ROOT / "delfin", defs, True)
    suite = ratchet._mentions(ROOT / "tests", defs, False)
    theirs = {n: ("test-only" if suite.get(n) else "nowhere")
              for n in defs if reached[n] == 0}

    row = cm.reachability(ROOT)
    expected = {"test-only": 0, "nowhere": 0}
    for tier in theirs.values():
        expected[tier] += 1
    assert row["tiers"] == expected, (row["tiers"], expected)
