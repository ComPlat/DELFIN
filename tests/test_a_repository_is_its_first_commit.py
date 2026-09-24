"""Opt-in: key the memory store by the history, so clones share one.

Input: a workspace path and the ``agent.memory_key`` setting. Output: the
slug naming its store.

    "path"  (default)  the main worktree's path — clones stay apart
    "repo"             the first commit of the history — clones share

Measured before this existed: a clone of a repository got its own store,
so a fresh checkout after a disaster started empty and several clones of
one project each learned separately.

Why the first commit and not the origin URL: it survives a remote being
renamed or moved, exists in a repository that has no remote at all, and
is identical in every clone of the same history. A repository with more
than one root commit takes the first in rev-list order, which is stable
for a given history.

Why "path" stays the default. A key that is too narrow costs recall; one
that is too wide puts one project's notes into another's prompt. Widening
silently would merge stores that are separate today, which nobody asked
for. Switching is a decision, and after it the notes under the old key
are still on disk — untouched, not migrated.

Every failure answers the path: no git, no commit, a timeout, an
unreadable repository. Splitting a store is the cheaper error.
"""

from __future__ import annotations

import subprocess

import pytest

from delfin.agent import memory_store as M


def _repo(path, body=None):
    """A repository whose first commit carries *body*.

    The body is a parameter because a commit object is its tree, author,
    message and timestamp: two repositories built with identical content
    in the same second get the same first commit, and git then calls them
    one history. Real projects differ; a fixture has to be told to.
    """
    path.mkdir(parents=True, exist_ok=True)
    subprocess.run(["git", "init", "-q"], cwd=str(path), check=True)
    for k, v in (("user.email", "a@b.c"), ("user.name", "t")):
        subprocess.run(["git", "config", k, v], cwd=str(path), check=True)
    (path / "f.txt").write_text(body or (path.name + "\n"), encoding="utf-8")
    subprocess.run(["git", "add", "-A"], cwd=str(path), check=True)
    subprocess.run(["git", "commit", "-q", "-m", "first"], cwd=str(path),
                   check=True)
    return path


@pytest.fixture()
def by_history(monkeypatch):
    monkeypatch.setattr(M, "_memory_key_mode", lambda: "repo")
    M._repo_root_cache.clear()
    yield
    M._repo_root_cache.clear()


@pytest.fixture(autouse=True)
def _clear_cache():
    M._repo_root_cache.clear()
    yield
    M._repo_root_cache.clear()


# -- what the setting changes ----------------------------------------------

def test_by_default_two_clones_stay_apart(tmp_path):
    src = _repo(tmp_path / "orig")
    clone = tmp_path / "clone"
    subprocess.run(["git", "clone", "-q", str(src), str(clone)], check=True)
    assert M._project_slug(src) != M._project_slug(clone)


def test_keyed_by_history_two_clones_share_one_store(tmp_path, by_history):
    src = _repo(tmp_path / "orig")
    clone = tmp_path / "clone"
    subprocess.run(["git", "clone", "-q", str(src), str(clone)], check=True)
    assert M._project_slug(src) == M._project_slug(clone)


def test_keyed_by_history_two_repositories_still_stay_apart(tmp_path,
                                                            by_history):
    """Different histories, different first commits: the costly error is
    one project's notes in another's prompt."""
    a = _repo(tmp_path / "a", "one project\n")
    b = _repo(tmp_path / "b", "another project\n")
    assert M._project_slug(a) != M._project_slug(b)


def test_two_byte_identical_beginnings_are_one_history(tmp_path, by_history):
    """The limit of the key, written down rather than left to be found.

    A commit object is its tree, author, message and timestamp. Two
    repositories created with all four identical -- a scaffolding script
    running twice within one second -- produce the same first commit, and
    share a store under this setting. Measured, not argued: this test
    fails the day git stops doing it.

    Not defended against, because every defence costs the property the
    setting exists for. Mixing the path back in would separate a clone
    from its original; mixing in the creation time would separate a
    repository from its own backup. The exposure is bounded: the setting
    is opt-in, the default key is the path, and two projects a user keeps
    apart in their work differ in their first commit.
    """
    a = _repo(tmp_path / "a", "same\n")
    b = _repo(tmp_path / "b", "same\n")
    if M._first_commit(a) != M._first_commit(b):
        pytest.skip("the two commits landed in different seconds")
    assert M._project_slug(a) == M._project_slug(b)


def test_keyed_by_history_a_worktree_still_shares(tmp_path, by_history):
    repo = _repo(tmp_path / "repo")
    wt = tmp_path / "wt"
    subprocess.run(["git", "worktree", "add", "-q", str(wt), "-b", "side"],
                   cwd=str(repo), check=True)
    assert M._project_slug(wt) == M._project_slug(repo)


def test_the_slug_is_one_segment_and_carries_the_hash(tmp_path, by_history):
    """No directory name in it: that is the one part that differs between
    a clone and its original, and putting it in would make the setting do
    nothing."""
    repo = _repo(tmp_path / "repo")
    slug = M._project_slug(repo)
    assert "/" not in slug and slug.startswith("-")
    first = M._first_commit(repo)
    assert first[:12] in slug
    assert "tmp" not in slug, "no path fragment may leak into the key"


# -- falling back ----------------------------------------------------------

def test_a_directory_outside_any_repository_keeps_its_path(tmp_path,
                                                           by_history):
    plain = tmp_path / "plain"
    plain.mkdir()
    assert M._project_slug(plain) == "-" + str(plain.resolve()).replace(
        "/", "-").lstrip("-")


def test_a_repository_with_no_commit_yet_keeps_its_path(tmp_path, by_history):
    empty = tmp_path / "empty"
    empty.mkdir()
    subprocess.run(["git", "init", "-q"], cwd=str(empty), check=True)
    assert M._project_slug(empty) == "-" + str(empty.resolve()).replace(
        "/", "-").lstrip("-")


def test_git_failing_falls_back_rather_than_raising(tmp_path, by_history,
                                                    monkeypatch):
    def _boom(*a, **k):
        raise OSError("no git here")
    monkeypatch.setattr(M.subprocess, "run", _boom)
    plain = tmp_path / "anything"
    plain.mkdir()
    assert M._project_slug(plain) == "-" + str(plain.resolve()).replace(
        "/", "-").lstrip("-")


# -- the setting itself ----------------------------------------------------

def test_the_default_is_the_path():
    assert M._memory_key_mode() == "path"


def test_an_unknown_value_reads_as_the_default(monkeypatch):
    from delfin import user_settings

    monkeypatch.setattr(user_settings, "load_settings",
                        lambda *a, **k: {"agent": {"memory_key": "wat"}})
    assert M._memory_key_mode() == "path"


def test_the_setting_is_shipped_with_its_default():
    from delfin.user_settings import DEFAULT_SETTINGS

    assert (DEFAULT_SETTINGS.get("agent") or {}).get("memory_key") == "path"


# -- the control ----------------------------------------------------------
#
# Source-level, like the other settings-UI tests here: building the tab
# needs a dashboard context and probes the machine for runtimes, which is
# not what these assertions are about. What they do cover is the drift
# that costs something -- a control writing a key nothing reads.

def _tab_source():
    import pathlib as _pl

    return (_pl.Path(M.__file__).resolve().parent.parent / "dashboard"
            / "tab_settings.py").read_text(encoding="utf-8")


def test_the_dashboard_offers_the_choice():
    """Two sensible values, so a dropdown: a typed third would silently
    point the agent at an empty store."""
    src = _tab_source()
    assert "memory_key_input = widgets.Dropdown(" in src
    assert "('this checkout (default)', 'path')" in src
    assert "('this repository (all clones)', 'repo')" in src


def test_both_save_paths_write_the_key():
    """The panel has its own save button and the tab has a save-all; a key
    written by only one of them changes back the next time the other runs."""
    src = _tab_source()
    assert "payload['agent']['memory_key'] = str(" in src
    assert "settings_payload['agent']['memory_key'] = str(" in src


def test_the_load_path_reads_it_back():
    src = _tab_source()
    assert "agent_payload.get('memory_key')" in src
    assert "memory_key_input.value = 'repo' if _key == 'repo' else 'path'" in src


def test_the_control_and_the_store_agree_on_the_key(monkeypatch):
    """The round trip that matters: the value the dropdown persists is the
    value the store switches on."""
    from delfin import user_settings

    for written, expected in (("repo", "repo"), ("path", "path")):
        monkeypatch.setattr(user_settings, "load_settings",
                            lambda *a, _w=written, **k: {
                                "agent": {"memory_key": _w}})
        assert M._memory_key_mode() == expected


def test_the_panel_says_switching_moves_nothing():
    """Notes under the old key stay on disk. A user who is not told that
    reads a suddenly empty store as data loss."""
    src = _tab_source()
    assert "Switching does not move or delete anything" in src
