"""The calc tools read ~/calc and ~/archive, whatever the user configured.

``_calc_dirs`` is set by exactly one caller in the whole codebase --
delfin/dashboard/tab_agent.py. In the CLI, in a headless run and in every
benchmark attempt it is empty, and ``_ensure_calc_loaded`` then fell back
to the literals ``~/calc`` and ``~/archive``.

For a user who set ``paths.calculations_dir`` to anything else, that is
an agent indexing a folder they do not use, silently. ``calc_summary``
reports a corpus that is not theirs and ``search_calcs`` answers "no runs
match" about a directory it never opened -- a negative finding from a
search of the wrong place, which is the same failure as the missing-
argument sweep one layer up.

``mcp_isolation._roots_for`` has read ``paths.calculations_dir`` and
``paths.archive_dir`` since the isolation work. There should be one
convention for where a user's calculations live, not two, so the calc
index now reads the same keys with the same fallbacks.

And the corpus names itself. "777 calculations" is a number an answer
will quote; until it says which directories it counted, nobody can tell a
complete answer from one about the wrong folder. Provenance is the first
of the evidence rules.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

import delfin.agent.api_client as A


@pytest.fixture
def calc_fixture(tmp_path):
    """A tiny calculation tree, in a directory that is nobody's default."""
    calc = tmp_path / "my-calcs"
    (calc / "run_a").mkdir(parents=True)
    (calc / "run_a" / "run_a.out").write_text(
        "ORCA\n! PBE0 def2-SVP\nFINAL SINGLE POINT ENERGY  -76.4\n",
        encoding="utf-8")
    return calc


def _fresh_executor():
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    ex._calc_engine = None
    ex._calc_dirs = {}
    ex._calc_roots = {}
    return ex


def test_the_configured_directory_is_the_one_that_is_read(
        monkeypatch, calc_fixture, tmp_path):
    archive = tmp_path / "my-archive"
    archive.mkdir()
    monkeypatch.setattr(
        A._DocToolExecutor, "_configured_calc_dirs",
        staticmethod(lambda: {"calc": str(calc_fixture),
                              "archive": str(archive)}))
    ex = _fresh_executor()
    assert ex._ensure_calc_loaded()
    assert ex._calc_roots["calc"] == str(calc_fixture)
    assert ex._calc_roots["archive"] == str(archive)


def test_an_explicit_override_still_wins(monkeypatch, calc_fixture):
    """The dashboard sets _calc_dirs for a reason; settings must not
    quietly take that back."""
    monkeypatch.setattr(
        A._DocToolExecutor, "_configured_calc_dirs",
        staticmethod(lambda: {"calc": "/nowhere/at/all"}))
    ex = _fresh_executor()
    ex._calc_dirs = {"calc": str(calc_fixture)}
    assert ex._ensure_calc_loaded()
    assert ex._calc_roots["calc"] == str(calc_fixture)


def test_with_nothing_configured_the_old_default_stands(monkeypatch):
    """An addition, not a change: a user who configured nothing gets
    exactly what they got before."""
    monkeypatch.setattr(
        A._DocToolExecutor, "_configured_calc_dirs", staticmethod(dict))
    ex = _fresh_executor()
    ex._ensure_calc_loaded()
    home = str(Path("~/calc").expanduser())
    assert ex._calc_roots.get("calc", "") in ("", home)


def test_the_summary_names_the_directories_it_counted(
        monkeypatch, calc_fixture, tmp_path):
    monkeypatch.setattr(
        A._DocToolExecutor, "_configured_calc_dirs",
        staticmethod(lambda: {"calc": str(calc_fixture)}))
    ex = _fresh_executor()
    with tempfile.TemporaryDirectory() as ws:
        perms = A.KitToolPermissions(mode="default", workspace=ws)
        out = json.loads(ex._execute_calc("calc_summary", {}))
    assert out["indexed_from"]["calc"] == str(calc_fixture)
    assert perms is not None


def test_a_missing_directory_is_said_rather_than_implied(monkeypatch):
    """"0 calculations" and "there is no calculation directory" are
    different answers, and only one of them is about the user's science."""
    monkeypatch.setattr(
        A._DocToolExecutor, "_configured_calc_dirs",
        staticmethod(lambda: {"calc": "/definitely/not/here",
                              "archive": "/nor/here"}))
    ex = _fresh_executor()
    ex._ensure_calc_loaded()
    out = json.loads(ex._execute_calc("calc_summary", {}))
    assert out["indexed_from"] == "no calculation directory found"


def test_settings_are_read_the_same_way_mcp_isolation_reads_them():
    """One convention, not two. If either side moves to another key, this
    is where the divergence shows up."""
    import inspect

    from delfin.agent import mcp_isolation

    theirs = inspect.getsource(mcp_isolation)
    mine = inspect.getsource(A._DocToolExecutor._configured_calc_dirs)
    for key in ("calculations_dir", "archive_dir"):
        assert key in theirs, key
        assert key in mine, key
