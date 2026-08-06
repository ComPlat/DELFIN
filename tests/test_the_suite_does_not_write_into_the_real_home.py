"""The test suite wrote into the user's real state, and into the checkout.

Two user-scoped stores resolve through ``Path.home()``, so the ``tmp_path``
a test passes as its repo root does not bound them:

  ``~/.delfin/memory``               the user-wide memory store
  ``~/.delfin/office_workspaces.json``  the office-folder registry

After an ordinary full run the registry held 178 entries. 176 were reprs of
mock objects that tests had passed where a path was expected, and the
directories named in them had really been created -- 752 kB of
``MagicMock/mock._permissions.workspace/<id>/.delfin/`` inside the DELFIN
checkout. They stayed invisible because a ``.delfin/`` ignore rule covers
them, so they were never at risk of being committed; they simply
accumulated.

That registry is one of the carriers of the folder lock: what it contains
decides which directory counts as an office workspace. The repository root
was among the 178 entries, which made a full-suite run treat DELFIN's own
checkout as an office folder -- and one memory test that passes in
isolation failed in the full run for that reason alone. It looked exactly
like flakiness and was not: it was a test writing into the state a later
test reads.

These tests pin the redirect. They assert about the running process, so
they hold for every test in the suite, not just this file.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.agent import memory_store as ms


_REAL_HOME = pathlib.Path.home() / ".delfin"


def _under_real_home(p: pathlib.Path) -> bool:
    try:
        p.resolve().relative_to(_REAL_HOME.resolve())
        return True
    except (OSError, ValueError):
        return False


# ---------------------------------------------------------------------------
# The two user-scoped stores are redirected
# ---------------------------------------------------------------------------

def test_the_user_wide_memory_store_is_not_the_real_one():
    assert not _under_real_home(ms._delfin_global_memory_dir())


def test_the_office_registry_is_not_the_real_one():
    assert not _under_real_home(ms._office_ws_file())


def test_a_global_memory_does_not_land_in_the_users_home(tmp_path):
    p, _, _ = ms.save_typed_memory(
        "user: prefers metric units", repo_root=tmp_path, scope="user")
    assert not _under_real_home(p)


def test_registering_an_office_folder_does_not_touch_the_real_registry(
    tmp_path,
):
    before = _REAL_HOME / "office_workspaces.json"
    stamp = before.stat().st_mtime if before.exists() else None
    ms.register_office_workspace(tmp_path)
    after = before.stat().st_mtime if before.exists() else None
    assert after == stamp
    assert str(tmp_path.resolve()) in ms._load_office_workspaces()


# ---------------------------------------------------------------------------
# ...and the redirect does not leak between tests
# ---------------------------------------------------------------------------

def test_the_registry_cache_does_not_survive_the_previous_test(tmp_path):
    """A module-global cache would outlive the per-test redirect."""
    assert str(tmp_path.resolve()) not in ms._load_office_workspaces()


@pytest.mark.parametrize("run", [1, 2])
def test_two_tests_do_not_see_each_others_office_folders(tmp_path, run):
    known_before = set(ms._load_office_workspaces())
    ms.register_office_workspace(tmp_path)
    assert known_before == set() or tmp_path not in known_before


# ---------------------------------------------------------------------------
# The thing that produced the mock directories in the first place
# ---------------------------------------------------------------------------

def test_a_mock_is_not_accepted_as_an_office_folder():
    """Registration requires a real directory, whatever it was handed."""
    from unittest.mock import MagicMock
    before = set(ms._load_office_workspaces())
    ms.register_office_workspace(MagicMock().workspace)
    assert set(ms._load_office_workspaces()) == before


def test_a_path_that_does_not_exist_is_not_registered(tmp_path):
    before = set(ms._load_office_workspaces())
    ms.register_office_workspace(tmp_path / "never_created")
    assert set(ms._load_office_workspaces()) == before


def test_a_file_is_not_an_office_folder(tmp_path):
    f = tmp_path / "notes.txt"
    f.write_text("x", encoding="utf-8")
    before = set(ms._load_office_workspaces())
    ms.register_office_workspace(f)
    assert set(ms._load_office_workspaces()) == before
