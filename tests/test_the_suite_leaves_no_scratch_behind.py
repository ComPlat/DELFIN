"""The suite's own scratch goes where pytest can collect it.

105 calls across 65 test files make a directory with ``tempfile.mkdtemp``
instead of the ``tmp_path`` fixture, and not one of them removes it.
Measured on the login node: 20836 entries in /tmp, among them 1827
``ws_*``, 427 ``office_*``, 295 ``planmode_*``, 233 ``askuser_*`` -- each
named after the test that made it.

Rewriting 105 call sites would fix the ones that exist and none of the
ones written next week, so ``tempfile.tempdir`` points at the directory
pytest already manages. pytest keeps the last three runs and removes what
is older: bounded by the suite rather than by the calendar.

Control, four of those files (65 tests) with a freshly emptied TMPDIR:

    origin/main   29 directories left behind
    with this      0
"""

from __future__ import annotations

import pathlib
import tempfile


def test_a_bare_mkdtemp_lands_under_pytests_own_root(tmp_path_factory):
    folder = pathlib.Path(tempfile.mkdtemp())
    base = tmp_path_factory.getbasetemp().resolve()
    assert base in folder.resolve().parents, (
        f"{folder} is outside {base}, so nothing will ever collect it")


def test_a_named_temporary_file_lands_there_too(tmp_path_factory):
    with tempfile.NamedTemporaryFile(suffix=".xyz") as handle:
        where = pathlib.Path(handle.name).resolve()
    base = tmp_path_factory.getbasetemp().resolve()
    assert base in where.parents


def test_the_answer_is_the_same_one_the_product_would_get(tmp_path_factory):
    """Product code asks ``gettempdir()``; it must see the same place, or
    the redirect covers the tests and misses what they drive."""
    base = tmp_path_factory.getbasetemp().resolve()
    assert base in pathlib.Path(tempfile.gettempdir()).resolve().parents or \
        base == pathlib.Path(tempfile.gettempdir()).resolve().parent
