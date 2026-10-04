"""A repository can narrow where the file tools may write.

``.delfin/settings.json`` {"kit": {"write_allow_globs": [...]}} names the
paths a session in that repo may write with write_file/edit_file. A write
elsewhere is refused with a reason. The list only narrows: it is taken
from the repo file alone, never copied into the user's file, and it
cannot reach outside the workspace.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A
from delfin.agent import kit_settings as ks


def _repo(tmp_path, globs):
    repo = tmp_path / "repo"
    (repo / ".delfin").mkdir(parents=True)
    (repo / ".delfin" / "settings.json").write_text(
        json.dumps({"kit": {"write_allow_globs": globs}}), encoding="utf-8")
    return repo


def test_the_repo_scope_is_loaded(tmp_path):
    repo = _repo(tmp_path, ["delfin/agent/foo.py", "tests/test_foo*.py"])
    got = ks.load(repo, user_path=tmp_path / "user.json")
    assert got.write_allow_globs == ["delfin/agent/foo.py", "tests/test_foo*.py"]


def test_a_user_file_cannot_set_or_widen_it(tmp_path):
    user = tmp_path / "user.json"
    user.write_text(json.dumps({"kit": {"write_allow_globs": ["**"]}}),
                    encoding="utf-8")
    assert ks.load(None, user_path=user).write_allow_globs == []
    repo = _repo(tmp_path, ["calc/**"])
    assert ks.load(repo, user_path=user).write_allow_globs == ["calc/**"]


def test_saving_the_merged_view_never_copies_it_to_the_user_file(tmp_path):
    repo = _repo(tmp_path, ["calc/**"])
    assert "write_allow_globs" not in ks.load(
        repo, user_path=tmp_path / "u.json").to_dict()


@pytest.mark.parametrize("rel, ok", [
    ("delfin/agent/foo.py", True),
    ("tests/test_foo_x.py", True),
    ("tests/sub/test_foo_x.py", False),   # * does not cross a directory
    ("delfin/agent/engine.py", False),
    ("calc/a/b/c.inp", True),
    ("calc", False),
])
def test_globs_match_repo_relative_paths(tmp_path, rel, ok):
    perms = A.KitToolPermissions(
        workspace=tmp_path,
        write_allow_globs=("delfin/agent/foo.py", "tests/test_foo*.py",
                           "calc/**"))
    assert A._write_in_scope(tmp_path / rel, perms) is ok


def test_no_scope_means_no_restriction(tmp_path):
    perms = A.KitToolPermissions(workspace=tmp_path)
    assert A._write_in_scope(tmp_path / "anything.py", perms) is True


def test_a_path_outside_the_workspace_matches_no_glob(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    perms = A.KitToolPermissions(workspace=ws, write_allow_globs=("**",))
    assert A._write_in_scope(tmp_path / "other" / "x.py", perms) is False


def test_the_write_gate_refuses_with_a_reason(tmp_path):
    ws = tmp_path / "ws"
    (ws / "src").mkdir(parents=True)
    perms = A.KitToolPermissions(
        workspace=ws, mode="acceptEdits", write_allow_globs=("src/ok.py",))
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    err = ex._gate_write_path(str(ws / "src" / "other.py"), perms,
                              "write_file", {"content": "x"})
    assert err and "write scope" in err and "src/ok.py" in err
    assert ex._gate_write_path(str(ws / "src" / "ok.py"), perms,
                               "write_file", {"content": "x"}) is None
