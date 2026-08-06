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


# ---------------------------------------------------------------------------
# ...and neither does anything else the suite writes to
# ---------------------------------------------------------------------------

def test_every_writable_user_state_sink_is_redirected(_user_state_targets):
    """The sink table is resolved once per session; this asserts the running
    process actually points every one of them somewhere disposable."""
    assert _user_state_targets, "the sink table resolved to nothing"
    leaked = [
        f"{mod.__name__}.{attr}"
        for mod, attr, _rel in _user_state_targets
        if _under_real_home(getattr(mod, attr))
    ]
    assert not leaked, f"still writing into the user's home: {leaked}"


def test_the_locator_index_is_redirected():
    """203 job records had accumulated here, 141 naming pytest tmp dirs."""
    from delfin.agent import bash_jobs as bj
    assert not _under_real_home(bj._INDEX_PATH)


def test_a_started_job_does_not_reach_the_real_index(tmp_path):
    from delfin.agent import bash_jobs as bj
    real = _REAL_HOME / "bash_jobs_index.json"
    stamp = real.stat().st_mtime if real.exists() else None
    bj._note_job_workspace("test-job", str(tmp_path))
    after = real.stat().st_mtime if real.exists() else None
    assert after == stamp
    assert bj._lookup_job_workspace("test-job") == str(tmp_path)


# ---------------------------------------------------------------------------
# The sinks resolved per call, including the one that matters most
# ---------------------------------------------------------------------------

def test_every_per_call_sink_is_redirected(_user_state_resolvers, tmp_path):
    assert _user_state_resolvers, "the resolver table resolved to nothing"
    leaked = []
    for mod, attr, _rel in _user_state_resolvers:
        fn = getattr(mod, attr)
        try:
            p = fn(tmp_path)
        except TypeError:
            p = fn()
        if _under_real_home(p):
            leaked.append(f"{mod.__name__}.{attr}")
    assert not leaked, f"still writing into the user's home: {leaked}"


def test_the_project_memory_store_is_not_in_the_users_home(tmp_path):
    """Its path is ~/.delfin/projects/<slug>/memory, and the slug comes from
    the repo root -- so passing tmp_path only changed the NAME of the
    directory it created in the user's home, it did not move it."""
    p, _, _ = ms.save_typed_memory("project: x", repo_root=tmp_path)
    assert not _under_real_home(p)


def test_the_users_permission_settings_are_not_the_real_file():
    """A test persisting an allow-rule edited the user's real settings: the
    file was found holding allow_patterns ["^.*$"] and default_mode
    "bypassPermissions" -- every command auto-approved, from a fixture."""
    from delfin.agent import kit_settings as ks
    assert not _under_real_home(ks.USER_SETTINGS_PATH)


def test_persisting_an_allow_rule_does_not_reach_the_real_settings():
    from delfin.agent import hooks_editor as he
    real = _REAL_HOME / "settings.json"
    stamp = real.stat().st_mtime if real.exists() else None
    he._write_settings({"kit": {"allow_patterns": ["^.*$"]}})
    after = real.stat().st_mtime if real.exists() else None
    assert after == stamp
