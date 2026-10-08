"""A missing path that says what it looks like saves a round.

Both shapes are from the field reports of 2026-10-08, three sessions
working together on a cluster.

An ELISION: `/clusterfs/clusterfs/groups/team/u12345/.../ChemDarwin_ABC_TEAM/README.md`
reached the file layer twice. The model had abbreviated a long path in
its own prose and then sent the abbreviation. "File not found" is true
and useless -- the path is not a path.

A HOME DIRECTORY BUILT BY HAND: `/home/u12345/agent_workspace/...`,
where the account's home is `/clusterfs/groups/team/u12345`. The guess is
reasonable and wrong, and nothing in the refusal said where home
actually is, so the answer cost a round every time.

Driven through the real tool executor rather than the helper alone: a
message nobody is handed is a message nobody reads.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture
def executor(tmp_path):
    perms = KitToolPermissions(workspace=tmp_path)
    perms.mode = "bypassPermissions"      # nothing here is about asking
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    eng._permissions = perms
    return eng, perms


def _error(out) -> str:
    return str(json.loads(out).get("error") or "")


class TestTheHelper:
    def test_an_elided_path_is_named_as_one(self):
        said = A._why_that_path_is_not_there("/clusterfs/.../README.md")
        assert "'...' segment" in said
        assert "abbreviation" in said

    def test_a_relative_elision_counts_too(self):
        assert "'...' segment" in A._why_that_path_is_not_there("a/.../b.md")

    def test_a_parent_directory_is_not_an_elision(self):
        """`..` is a real path component and means something."""
        assert A._why_that_path_is_not_there("../sibling/file.md") == ""
        assert A._why_that_path_is_not_there("a/../b.md") == ""

    def test_a_home_built_by_hand_is_told_where_home_is(self, monkeypatch,
                                                        tmp_path):
        home = tmp_path / "home" / "ka" / "ka_ibcs" / "u12345"
        home.mkdir(parents=True)
        monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
        said = A._why_that_path_is_not_there("/home/u12345/workspace/x.md")
        assert str(home) in said, said
        assert "u12345" in said

    def test_a_path_under_the_real_home_gets_no_lecture(self, monkeypatch,
                                                        tmp_path):
        home = tmp_path / "home" / "u12345"
        home.mkdir(parents=True)
        monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
        assert A._why_that_path_is_not_there(str(home / "gone.md")) == ""

    def test_an_ordinary_typo_says_nothing(self, monkeypatch, tmp_path):
        home = tmp_path / "home" / "u12345"
        home.mkdir(parents=True)
        monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
        assert A._why_that_path_is_not_there("/etc/passwdd") == ""
        assert A._why_that_path_is_not_there("notes/readme.md") == ""

    def test_a_relative_path_is_not_read_as_a_home_guess(self, monkeypatch,
                                                         tmp_path):
        home = tmp_path / "home" / "u12345"
        home.mkdir(parents=True)
        monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
        assert A._why_that_path_is_not_there("u12345/notes.md") == ""

    def test_it_never_raises(self):
        for odd in ("", None, "\x00", "/" * 300, "a" * 5000):
            assert isinstance(A._why_that_path_is_not_there(odd), str)


class TestItReachesTheTool:
    def test_read_file_hands_the_advice_back(self, executor, monkeypatch,
                                             tmp_path):
        eng, perms = executor
        home = tmp_path / "home" / "ka" / "ka_ibcs" / "u12345"
        home.mkdir(parents=True)
        monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
        # Inside the workspace, so the read gate is not what answers.
        missing = tmp_path / "u12345" / "gone.md"
        out = eng._execute_read_file({"path": str(missing)}, perms)
        said = _error(out)
        assert "File not found" in said
        assert str(home) in said, said

    def test_an_elision_inside_the_workspace_is_named(self, executor):
        eng, perms = executor
        out = eng._execute_read_file({"path": "reports/.../plan.md"}, perms)
        said = _error(out)
        assert "not found" in said.lower()
        assert "abbreviation" in said, said

    def test_a_plain_miss_is_still_plain(self, executor):
        eng, perms = executor
        out = eng._execute_read_file({"path": "reports/plan.md"}, perms)
        said = _error(out)
        assert "not found" in said.lower()
        assert "abbreviation" not in said and "home directory" not in said

    def test_the_advice_does_not_reach_a_refusal(self, executor, monkeypatch,
                                                 tmp_path):
        """A path outside the roots is refused by the read gate, and that
        refusal says what to do about the gate. Advice about how the path
        was spelled would be a second answer to a different question --
        and it would confirm whether the file exists out there."""
        eng, perms = executor
        home = tmp_path / "home" / "ka" / "ka_ibcs" / "u12345"
        home.mkdir(parents=True)
        monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
        out = eng._execute_read_file(
            {"path": "/home/u12345/secret/notes.md"}, perms)
        said = _error(out)
        assert said, out
        assert "home directory" not in said, (
            "the gate's refusal must not be turned into a hint about a "
            f"path outside the workspace: {said}")
