"""A path that appears in the checkout is reported with who left it.

Input: the checkout's top-level names, sampled as each test starts and
once more at the end of the run. Output: the session guard's message,
with "(after <nodeid>)" beside a path whose arrival it can place.

Why. The guard reported three paths and no name, and finding the single
wrong argument behind them took a bisect of the suite. The paths were
created by the LAST test of its file, which is why sampling only at the
start of the next test was not enough on its own.

Shallow on purpose, and the number is the reason: the full walk the
session guard uses costs 7.0 ms, which is about 250 s over a suite this
size once per test; listing the checkout's direct children costs 0.016 ms,
or 0.6 s. A path created deeper is still reported, just not placed.

Attribution only. The session guard decides, on the whole tree; this
answers the question it could not.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import conftest as C  # noqa: E402


def test_a_name_that_appears_between_two_starts_is_placed(monkeypatch):
    monkeypatch.setattr(C, "_CHECKOUT_TOP_BLAME", {})
    monkeypatch.setattr(C, "_CHECKOUT_TOP_LAST", None)
    monkeypatch.setattr(C, "_CHECKOUT_TOP_LAST_NODE", "(before)")
    seen = [frozenset({"a"}), frozenset({"a"}), frozenset({"a", "leaked"})]
    monkeypatch.setattr(C, "_checkout_top", lambda: seen.pop(0))

    class _Item:
        def __init__(self, nodeid):
            self.nodeid = nodeid

    C.pytest_runtest_setup(_Item("t::one"))
    C.pytest_runtest_setup(_Item("t::two"))
    C.pytest_runtest_setup(_Item("t::three"))
    assert C._CHECKOUT_TOP_BLAME == {"leaked": "t::two"}, (
        "the name belongs to the test that ran between the two samples")


def test_the_first_sample_blames_nobody(monkeypatch):
    """Whatever is in the checkout before the first test is not a leak."""
    monkeypatch.setattr(C, "_CHECKOUT_TOP_BLAME", {})
    monkeypatch.setattr(C, "_CHECKOUT_TOP_LAST", None)
    monkeypatch.setattr(C, "_checkout_top", lambda: frozenset({"a", "b"}))

    class _Item:
        nodeid = "t::one"

    C.pytest_runtest_setup(_Item())
    assert C._CHECKOUT_TOP_BLAME == {}


def test_the_reader_is_told_which_test(monkeypatch):
    """The message is what a person acts on, so the name has to be in it."""
    monkeypatch.setattr(C, "_CHECKOUT_TOP_BLAME",
                        {"child_home": "tests/t.py::test_x"})
    import os

    def _named(path):
        who = C._CHECKOUT_TOP_BLAME.get(os.path.basename(str(path)))
        return f"{path} (after {who})" if who else str(path)

    assert _named("/repo/child_home") == (
        "/repo/child_home (after tests/t.py::test_x)")
    assert _named("/repo/other") == "/repo/other"


def test_the_guard_places_what_it_can_and_prints_the_rest():
    """A nested path keeps its place in the report even unplaced -- the
    shallow sample cannot name it, and dropping it would hide a leak."""
    import inspect

    src = inspect.getsource(C._the_suite_does_not_write_into_the_checkout)
    assert "_named(p)" in src
    assert "_CHECKOUT_TOP_LAST" in src, "the last test's leak is placed too"


def test_sampling_is_shallow_and_says_why():
    import inspect

    src = inspect.getsource(C._checkout_top)
    assert "listdir" in src
    assert "0.016" in src and "7.0" in src, "the trade-off carries its numbers"
