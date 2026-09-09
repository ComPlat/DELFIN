"""The path refusals for the mistake nobody was making on purpose.

A workspace is usually a subdirectory of something, and whoever wrote the
request was standing outside it. The agent IS that directory, so the path
it copies out of the request resolves one level of the same name deeper,
or arrives as an absolute cwd naming the directory above.

Measured over ~/.delfin/audit.log, 2026-09-09: 27 refusals reading `path
is outside the allowed workspace roots` and 9 reading `cwd is not a
directory`, all one wording mistake, and 10 of the 13 in the first group
passed the project root as an absolute cwd. Second only to the auto-allow
list, and unlike that one it costs nothing to fix -- the refusal already
had the two paths it needed to name the mistake, and named neither.

The hint is a hint on purpose. Silently accepting
`tests/fixtures/ws/run.py` from inside `tests/fixtures/ws` would hide a
genuine nested directory of the same name, and would teach the model
nothing about where it is.

Companion to test_the_one_repo_path_that_means_here.py, which stops the
mistake happening; this one is about what the agent is told when it does.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor
from delfin.agent.api_client import _DocToolExecutor as E


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "tests" / "fixtures" / "user_project_workspace"
    d.mkdir(parents=True)
    (d / "run.py").write_text("print('hi')\n", encoding="utf-8")
    return d


# ---------------------------------------------------------------------------
# The path the caller meant
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("given,meant", [
    ("tests/fixtures/user_project_workspace/run.py", "run.py"),
    ("tests/fixtures/user_project_workspace", "."),
    ("fixtures/user_project_workspace/a/b.py", "a/b.py"),
    ("user_project_workspace/run.py", "run.py"),
])
def test_the_duplicated_head_is_recognised(ws, given, meant):
    assert E._names_the_workspace_from_outside(given, ws) == meant


@pytest.mark.parametrize("given", [
    "run.py",
    "src/run.py",
    "delfin/agent/cli.py",
    "../outside.py",
    "/absolute/path.py",
    "",
    "user_project_workspace_other/run.py",   # a different directory
])
def test_an_ordinary_path_gets_no_hint(ws, given):
    """A hint that fires for everyone is noise, and noise is what people
    learn to skip."""
    assert E._names_the_workspace_from_outside(given, ws) is None
    assert E._hint_if_it_names_the_workspace(given, ws) == ""


def test_a_directory_containing_the_workspace_is_recognised(ws):
    root = ws.parents[2]
    assert E._is_ancestor_of_workspace(str(root), ws) is True
    assert E._is_ancestor_of_workspace(str(ws), ws) is False
    assert E._is_ancestor_of_workspace(str(ws / "sub"), ws) is False


# ---------------------------------------------------------------------------
# ...where the agent actually meets it
# ---------------------------------------------------------------------------

def _err(args, ws) -> str:
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    out = _doc_executor.execute("bash", args, perms)
    return str(json.loads(out).get("error", ""))


def test_the_cwd_refusal_names_the_path_that_works(ws):
    err = _err({"command": "ls",
                "cwd": "tests/fixtures/user_project_workspace"}, ws)
    assert "cwd is not a directory" in err
    assert "already IS 'user_project_workspace'" in err
    assert "'.'" in err


def test_the_root_as_cwd_is_told_what_it_named(ws):
    err = _err({"command": "ls", "cwd": str(ws.parents[2])}, ws)
    assert "CONTAINS your workspace" in err
    assert str(ws) in err


def test_an_unrelated_bad_cwd_stays_plain(ws):
    err = _err({"command": "ls", "cwd": "nope"}, ws)
    assert "cwd is not a directory" in err
    assert "already IS" not in err
    assert "CONTAINS your workspace" not in err


def test_the_grant_advice_is_not_lost(ws):
    """The outside-roots refusal still has to say the thing it is for --
    that a path can be granted -- or the hint has replaced it."""
    err = _err({"command": "ls", "cwd": str(ws.parents[2])}, ws)
    assert "--add-dir" in err or "extra_workspace_dirs" in err


def test_the_path_that_works_actually_works(ws):
    """A hint pointing at a command that is also refused would be worse
    than no hint."""
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    out = json.loads(_doc_executor.execute(
        "bash", {"command": "ls", "cwd": "."}, perms))
    assert out.get("exit_code") == 0, out
    assert "run.py" in str(out.get("stdout", ""))
