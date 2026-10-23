"""One shared directory, and it is one directory.

Three sessions worked on a repository, each in its own worktree, and
handed each other FILENAMES. Every handover failed with "path is outside
the allowed workspace roots" -- the file existed, in the author's
worktree, which is not a path the reader may open. The containment was
right: a session that can read another session's worktree can read its
whole workspace.

So this is not a wider boundary. It is one named directory inside the
repository's own `.delfin`, shared by the repository and every worktree
of it. Most of what follows asserts the LIMITS of that door rather than
the door, because a shared root is the one place two sessions can both
write, which makes it the one place to plant a symlink.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

from delfin.agent import exchange as H


def _git(where, *args):
    return subprocess.run(["git", "-C", str(where), *args],
                          capture_output=True, text=True, timeout=30)


@pytest.fixture
def repo(tmp_path, monkeypatch):
    """A repository with two worktrees, as the three sessions had."""
    from delfin.agent import session_presence as P
    monkeypatch.setattr(P, "_git_cache", {})

    root = tmp_path / "project"
    root.mkdir()
    _git(root, "init", "-q")
    _git(root, "config", "user.email", "t@example.invalid")
    _git(root, "config", "user.name", "t")
    (root / "README.md").write_text("x\n")
    _git(root, "add", "README.md")
    _git(root, "commit", "-qm", "first")
    trees = root / ".delfin" / "worktrees"
    trees.mkdir(parents=True)
    for name in ("wt-a", "wt-b"):
        _git(root, "worktree", "add", "-q", "-b", name, str(trees / name))
    return root


class TestOnePlaceForOneRepository:
    def test_every_worktree_agrees_on_the_path(self, repo):
        """The point. The checkout ROOT differs per worktree -- that is
        what a worktree is -- so the common directory is the anchor."""
        a = H.directory_for(repo / ".delfin" / "worktrees" / "wt-a")
        b = H.directory_for(repo / ".delfin" / "worktrees" / "wt-b")
        main = H.directory_for(repo)
        assert a == b == main
        assert a is not None

    def test_it_is_inside_the_repository_that_owns_it(self, repo):
        room = H.directory_for(repo)
        assert room == (repo / ".delfin" / "handover").resolve()

    def test_another_repository_gets_another_place(self, repo, tmp_path,
                                                   monkeypatch):
        from delfin.agent import session_presence as P
        monkeypatch.setattr(P, "_git_cache", {})
        other = tmp_path / "other"
        other.mkdir()
        _git(other, "init", "-q")
        assert H.directory_for(other) != H.directory_for(repo)

    def test_outside_a_repository_there_is_none(self, tmp_path, monkeypatch):
        from delfin.agent import session_presence as P
        monkeypatch.setattr(P, "_git_cache", {})
        loose = tmp_path / "loose"
        loose.mkdir()
        assert H.directory_for(loose) is None

    def test_it_is_owner_only(self, repo):
        room = H.directory_for(repo)
        assert room.is_dir()
        assert (os.stat(room).st_mode & 0o777) == 0o700, oct(
            os.stat(room).st_mode)

    def test_asking_without_creating_touches_nothing(self, repo):
        room = H.directory_for(repo, create=False)
        assert room is not None
        assert not room.exists()


class TestItIsOneDoorAndNotAParent:
    def test_it_is_never_the_repository(self, repo):
        room = H.directory_for(repo)
        assert room != repo.resolve()
        assert repo.resolve() in room.parents

    def test_it_is_never_dot_delfin_itself(self, repo):
        """`.delfin` holds the worktrees, the session store and the
        settings. Handing it over would hand over everything this keeps
        apart."""
        room = H.directory_for(repo)
        assert room != (repo / ".delfin").resolve()

    def test_it_does_not_contain_the_worktrees(self, repo):
        room = H.directory_for(repo)
        trees = (repo / ".delfin" / "worktrees").resolve()
        assert room not in trees.parents and room != trees

    def test_a_symlink_inside_it_is_not_a_way_out(self, repo, tmp_path):
        """Both sessions can write here, so both can plant one."""
        room = H.directory_for(repo)
        outside = tmp_path / "beyond"
        outside.mkdir()
        (outside / "secret.txt").write_text("no\n")
        (room / "doorway").symlink_to(outside)
        assert not H.contains(room / "doorway" / "secret.txt", repo)
        assert H.contains(room / "report.md", repo)

    def test_a_dotdot_path_is_not_inside_it(self, repo):
        room = H.directory_for(repo)
        assert not H.contains(room / ".." / "settings.json", repo)
        assert not H.contains(room / ".." / "worktrees" / "wt-a", repo)


class TestItReachesTheSession:
    def test_granting_does_not_create_the_directory(self, repo):
        """Granting a root runs whenever a session's permissions are
        built -- thousands of times in one suite run. Creating there put
        a directory into the real checkout as a side effect of
        constructing an object: it appeared in the repository, the guard
        against writing into the checkout saw it, and 14 unrelated tests
        failed. `write_file` creates its own parents, so the first
        session that writes a report makes the directory.
        """
        room = H.directory_for(repo, create=False)
        assert room is not None
        assert H.with_exchange(repo, ()) == (room,)
        assert not room.exists(), (
            "granting a root must not write to the repository")

    def test_a_write_into_it_creates_it(self, repo):
        """The other half: granted-but-absent has to be usable."""
        from delfin.agent import api_client as A

        mine = repo / ".delfin" / "worktrees" / "wt-a"
        room = H.directory_for(repo, create=False)
        perms = A.KitToolPermissions(
            workspace=mine, extra_workspace_dirs=H.with_exchange(mine, ()))
        perms.mode = "bypassPermissions"
        eng = A._DocToolExecutor.__new__(A._DocToolExecutor)
        eng._permissions = perms
        eng._execute_write_file(
            {"path": str(room / "schema.md"), "content": "the contract\n"},
            perms)
        assert (room / "schema.md").read_text() == "the contract\n"

    def test_the_root_is_granted_once(self, repo):
        room = H.directory_for(repo, create=False)
        first = H.with_exchange(repo, ())
        assert first == (room,)
        assert H.with_exchange(repo, first) == first, "granted twice"

    def test_an_existing_grant_is_kept(self, repo, tmp_path):
        other = tmp_path / "granted"
        other.mkdir()
        out = H.with_exchange(repo, (other,))
        assert other in out and H.directory_for(repo) in out

    def test_the_setting_turns_it_off(self, repo, monkeypatch):
        import delfin.user_settings as us
        monkeypatch.setattr(us, "load_settings",
                            lambda *a, **k: {"agent": {"handover_dir": False}})
        assert H.with_exchange(repo, ()) == ()

    def test_it_is_on_without_a_setting(self, repo, monkeypatch):
        import delfin.user_settings as us
        monkeypatch.setattr(us, "load_settings", lambda *a, **k: {})
        assert H.with_exchange(repo, ()) == (H.directory_for(repo),)

    def test_an_unreadable_settings_file_does_not_switch_it_off(self, repo,
                                                               monkeypatch):
        import delfin.user_settings as us
        monkeypatch.setattr(us, "load_settings",
                            lambda *a, **k: (_ for _ in ()).throw(OSError("x")))
        assert H.enabled() is True

    def test_outside_a_repository_nothing_is_granted(self, tmp_path,
                                                     monkeypatch):
        from delfin.agent import session_presence as P
        monkeypatch.setattr(P, "_git_cache", {})
        loose = tmp_path / "loose2"
        loose.mkdir()
        assert H.with_exchange(loose, ()) == ()


class TestTheRefusalNamesIt:
    def test_a_read_of_another_worktree_says_where_to_put_it(self, repo,
                                                             monkeypatch):
        """Driven through the real resolver: this message is where the
        need showed up, and a hint nobody is handed is not a hint."""
        from delfin.agent import api_client as A
        from delfin.agent import session_presence as P
        monkeypatch.setattr(P, "_git_cache", {})

        mine = repo / ".delfin" / "worktrees" / "wt-a"
        theirs = repo / ".delfin" / "worktrees" / "wt-b"
        (theirs / "report.md").write_text("their work\n")

        perms = A.KitToolPermissions(workspace=mine)
        eng = A._DocToolExecutor.__new__(A._DocToolExecutor)
        eng._permissions = perms
        _resolved, err = eng._resolve_in_workspace(
            str(theirs / "report.md"), perms, for_read=True)
        assert err, "the other worktree must still be refused"
        assert "another session's worktree" in err, err
        assert str(H.directory_for(repo)) in err, err

    def test_an_unrelated_path_gets_no_handover_hint(self, repo, tmp_path,
                                                     monkeypatch):
        """A hint that fires everywhere is one nobody reads."""
        from delfin.agent import api_client as A
        from delfin.agent import session_presence as P
        monkeypatch.setattr(P, "_git_cache", {})

        mine = repo / ".delfin" / "worktrees" / "wt-a"
        perms = A.KitToolPermissions(workspace=mine)
        eng = A._DocToolExecutor.__new__(A._DocToolExecutor)
        eng._permissions = perms
        _resolved, err = eng._resolve_in_workspace(
            str(tmp_path / "elsewhere" / "x.md"), perms, for_read=True)
        assert err
        # On the hint's own words, not on "worktree": the roots list in
        # every such refusal names the session's own worktree path.
        assert "another session's worktree" not in err, err
        assert "hand work over" not in err, err

    def test_the_shared_directory_itself_is_readable(self, repo, monkeypatch):
        """The grant has to actually work, not only be announced."""
        from delfin.agent import api_client as A
        from delfin.agent import session_presence as P
        monkeypatch.setattr(P, "_git_cache", {})

        mine = repo / ".delfin" / "worktrees" / "wt-a"
        room = H.directory_for(repo)
        (room / "schema.md").write_text("the frozen contract\n")

        perms = A.KitToolPermissions(
            workspace=mine, extra_workspace_dirs=H.with_exchange(mine, ()))
        eng = A._DocToolExecutor.__new__(A._DocToolExecutor)
        eng._permissions = perms
        resolved, err = eng._resolve_in_workspace(
            str(room / "schema.md"), perms, for_read=True)
        assert err is None, err
        assert resolved == (room / "schema.md").resolve()


class TestWhatTheModelIsTold:
    def test_the_line_names_the_path_and_the_limit(self, repo):
        said = H.describe(repo)
        assert str(H.directory_for(repo, create=False)) in said
        assert "other worktrees" in said
        assert "not readable" in said, (
            "a line that offers the shared directory without saying the "
            "worktrees stay closed invites the next request for one")

    def test_outside_a_repository_it_says_nothing(self, tmp_path, monkeypatch):
        from delfin.agent import session_presence as P
        monkeypatch.setattr(P, "_git_cache", {})
        loose = tmp_path / "loose3"
        loose.mkdir()
        assert H.describe(loose) == ""
