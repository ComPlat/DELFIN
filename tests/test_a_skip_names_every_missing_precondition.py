"""A skip names every missing precondition, on any machine.

``tests/test_ds1_real_run.py`` needs two independent things: the optional
``morfeus`` package and the practice geometry
``.gate/practice/znh2o2cl2.xyz``, which was never committed. A module-level
``pytest.importorskip("morfeus")`` fired at collection, so the file came out
as ONE skip naming the package and the geometry reason was never shown --
including on a machine that has morfeus, where it is the only reason there
is.

That file can only ever skip here, so its own run proves nothing: it holds
two tests, and measured on delfin-ci -- where neither precondition is
present -- both skip with both reasons named, 2 skipped and nothing passed.
Two is also the ceiling: a file of two tests cannot report more. This file
tests the decision
instead of the outcome. It calls ``_needs_geometry`` directly with each
combination of the two preconditions put in place, and reads the reason out
of the ``Skipped`` exception. Both directions are simulated, so the result
does not depend on what this host happens to have installed: ``morfeus`` is
put into ``sys.modules`` to make it present and set to ``None`` there to make
it absent (which is what makes ``import`` raise), and the geometry path is
pointed at a file that does or does not exist.

Input: the two preconditions, in all four combinations. Output: whether a
skip is raised and what its reason names.

Semantics pinned here: a precondition that is missing is named, every
missing one is named, and nothing is skipped when both are in place. The
rejected alternative -- checking them in order -- is pinned too: it passes
the first three cases and fails the fourth, because ordering only moves the
masking from one reason to the other.
"""
from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

import pytest

MODULE_PATH = Path(__file__).resolve().parent / "test_ds1_real_run.py"


@pytest.fixture
def ds1():
    """The module under test, loaded from its path.

    Importable at all only because the module-level skip is gone: with it,
    importing the module raised ``Skipped`` at import time, which is why the
    reason below could never be reached from outside.
    """
    spec = importlib.util.spec_from_file_location("_ds1_under_test", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def _reason(ds1, monkeypatch, *, morfeus: bool, geometry: bool,
            tmp_path: Path) -> str | None:
    """The skip's reason for one combination, or None when nothing skips."""
    if morfeus:
        stand_in = types.ModuleType("morfeus")
        stand_in.__version__ = "0.0.0-stand-in"
        monkeypatch.setitem(sys.modules, "morfeus", stand_in)
    else:
        # None in sys.modules is what makes `import morfeus` raise, on a
        # host that has the package as well as on one that does not.
        monkeypatch.setitem(sys.modules, "morfeus", None)

    target = tmp_path / "znh2o2cl2.xyz"
    if geometry:
        target.write_text("1\n\nZn 0.0 0.0 0.0\n")
    monkeypatch.setattr(ds1, "GEOMETRY", target)

    # pytest.skip raises through BaseException, not Exception: caught as
    # `except Exception` it propagates and SKIPS this test instead of being
    # read -- three of the four cases below came out as skips that way, and a
    # skipped test asserts nothing.
    try:
        ds1._needs_geometry()
    except pytest.skip.Exception as exc:
        return str(getattr(exc, "msg", exc))
    return None


def test_nothing_skips_when_both_are_in_place(ds1, monkeypatch, tmp_path):
    assert _reason(ds1, monkeypatch, morfeus=True, geometry=True,
                   tmp_path=tmp_path) is None


def test_a_missing_package_is_named(ds1, monkeypatch, tmp_path):
    why = _reason(ds1, monkeypatch, morfeus=False, geometry=True,
                  tmp_path=tmp_path)
    assert why and "morfeus" in why
    assert "scratch geometry" not in why, "it is there; do not say it is not"


def test_a_missing_geometry_is_named(ds1, monkeypatch, tmp_path):
    """The case the module-level skip hid on a machine that has morfeus."""
    why = _reason(ds1, monkeypatch, morfeus=True, geometry=False,
                  tmp_path=tmp_path)
    assert why and "scratch geometry" in why
    assert "morfeus" not in why, "it is importable; do not say it is not"


def test_both_missing_are_named_together(ds1, monkeypatch, tmp_path):
    """The one that fails if the two are checked in order instead.

    Ordering passes the three cases above and fails here: whichever is
    checked first returns its own reason and hides the other.
    """
    why = _reason(ds1, monkeypatch, morfeus=False, geometry=False,
                  tmp_path=tmp_path)
    assert why
    assert "morfeus" in why
    assert "scratch geometry" in why


def test_the_file_has_no_module_level_skip_left():
    """Loading it must not skip. A module-level ``importorskip`` -- or an
    autouse fixture doing the same thing -- brings the masking back.

    It reads the source rather than taking the ``ds1`` fixture on purpose:
    with the masking in place the fixture raises ``Skipped`` at import, so
    every test that takes it would come out SKIPPED, and a skipped test
    reports nothing. This one fails there instead.
    """
    source = MODULE_PATH.read_text()
    assert "importorskip" not in source.split("def _needs_geometry")[0], (
        "a module-level importorskip masks the per-test reason again")
    assert "autouse" not in source, (
        "an autouse fixture that skips the module is the same masking")


def test_the_skip_count_in_this_docstring_cannot_exceed_the_file():
    """Pin the number the docstring above states, to the file it is about.

    The docstring said 3 skipped where the file holds two tests, so the
    figure was not merely stale: no run of a two-test file can report it.
    An uncheckable number drifts silently, and this one was read as
    evidence that the real chain had been exercised.

    Input: the test functions defined in ``test_ds1_real_run.py``, counted
    from its source. Output: the count, compared against the figure this
    module's own docstring quotes. Counting from the source rather than
    from a run keeps the pin host-independent -- a run here skips, and a
    skipped test reports nothing.
    """
    import re

    defined = re.findall(r"^def (test_\w+)", MODULE_PATH.read_text(),
                         flags=re.MULTILINE)
    quoted = re.search(r"(\d+) skipped and nothing passed", __doc__ or "")
    assert quoted, "the docstring no longer states a skip count"
    assert int(quoted.group(1)) == len(defined), (
        f"the docstring says {quoted.group(1)} skipped; "
        f"{MODULE_PATH.name} defines {len(defined)} tests: {defined}")
