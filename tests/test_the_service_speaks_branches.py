"""The agent service speaks branches: ``branch.*`` service methods.

Input: a private agent service (its own HOME, tests/tui_service.py's rig)
and a fake git checkout under it. Output: the new service methods

* ``branch.list`` -- the branch rows of a view (``view``: mine/all/hidden/
  merged/remote), of the session's workspace;
* ``branch.pull`` -- fast-forward-only update of main, never a switch;
* ``branch.test`` -- a worktree at the branch, its changed tests run there,
  pass/fail counts and the first failure, worktree removed unless kept;
* ``branch.cleanup`` -- what may be deleted, and the deletion itself with
  ``confirm=true`` (nothing is deleted without it).

Every method answers even when the workspace is not a git checkout, with an
``ok: false`` and a reason -- no crash, no silent success.
"""

from __future__ import annotations

import os
import subprocess

import pytest

from tests.tui_service import ServiceProcess


def _run(cwd, *args, check=True):
    out = subprocess.run(["git", "-c", "user.email=t@t", "-c", "user.name=t",
                          *args], cwd=str(cwd), capture_output=True, text=True,
                         env={**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
                              "GIT_CONFIG_SYSTEM": os.devnull})
    if check and out.returncode != 0:
        raise AssertionError(f"git {' '.join(args)}: {out.stderr}")
    return out.stdout.strip()


def _commit(cwd, name):
    (cwd / name).write_text(name + "\n")
    _run(cwd, "add", name)
    _run(cwd, "commit", "-q", "-m", name)


@pytest.fixture
def git_service(tmp_path):
    """A private service and a git checkout as one of its workspaces."""
    svc = ServiceProcess(tmp_path, start=False)
    ws = svc.workspace("gitws")
    origin = tmp_path / "origin.git"
    _run(tmp_path, "init", "-q", "--bare", "-b", "main", str(origin))
    _run(tmp_path, "clone", "-q", str(origin), str(ws))
    _commit(ws, "a.txt")
    _run(ws, "push", "-q", "-u", "origin", "main")
    _run(ws, "branch", "old/merged")
    _run(ws, "checkout", "-q", "-b", "feat")
    _commit(ws, "f.txt")
    _run(ws, "checkout", "-q", "main")
    svc.start()
    try:
        yield svc, ws
    finally:
        svc.stop()


def test_branch_list_views(git_service):
    svc, ws = git_service
    with svc.client() as c:
        c.request("hello", {"protocolVersion": 1})
        reply = c.request("branch.list", {"workspace": str(ws), "view": "mine"},
                          timeout=60)
        assert reply["ok"] is True
        refs = [r["ref"] for r in reply["rows"]]
        assert refs[0] == "main"
        assert "feat" in refs
        merged = c.request("branch.list", {"workspace": str(ws),
                                           "view": "merged"}, timeout=60)
        assert "old/merged" in [r["ref"] for r in merged["rows"]]


def test_branch_pull(git_service):
    svc, ws = git_service
    other = tmp_other(ws)
    _commit(other, "b.txt")
    _run(other, "push", "-q", "origin", "main")
    with svc.client() as c:
        c.request("hello", {"protocolVersion": 1})
        reply = c.request("branch.pull", {"workspace": str(ws)}, timeout=60)
    assert reply["ok"] is True
    assert _run(ws, "rev-parse", "HEAD") == _run(other, "rev-parse", "HEAD")


def tmp_other(ws):
    other = ws.parent / "pull-origin-clone"
    _run(ws.parent, "clone", "-q",
         str(ws.parent / "origin.git"), str(other))
    return other


def test_branch_pull_refuses_off_main(git_service):
    svc, ws = git_service
    _run(ws, "checkout", "-q", "-b", "side")
    with svc.client() as c:
        c.request("hello", {"protocolVersion": 1})
        reply = c.request("branch.pull", {"workspace": str(ws)}, timeout=60)
    assert reply["ok"] is False
    assert "main" in reply["text"].lower()
    # the checkout is still where it was: a pull never switches
    assert _run(ws, "branch", "--show-current") == "side"


def test_branch_test_runs_and_cleans_up(git_service, monkeypatch):
    svc, ws = git_service
    _run(ws, "checkout", "-q", "-b", "fix-it")
    (ws / "tests").mkdir()
    (ws / "tests" / "test_the_fix.py").write_text("def test_fix():\n"
                                                  "    assert 1 == 1\n")
    _run(ws, "add", "tests/test_the_fix.py")
    _run(ws, "commit", "-q", "-m", "fix")
    _run(ws, "checkout", "-q", "main")
    with svc.client() as c:
        c.request("hello", {"protocolVersion": 1})
        reply = c.request("branch.test", {"workspace": str(ws),
                                          "ref": "fix-it"}, timeout=120)
    assert reply["ok"] is True, reply
    assert reply["passed"] >= 1 and reply["failed"] == 0
    # the worktree is gone again, the checkout still on main
    assert not (ws.parent / "delfin-preview" / "fix-it").exists()
    assert _run(ws, "branch", "--show-current") == "main"


def test_branch_cleanup_needs_confirm_and_then_deletes(git_service):
    svc, ws = git_service
    with svc.client() as c:
        c.request("hello", {"protocolVersion": 1})
        listing = c.request("branch.cleanup", {"workspace": str(ws)},
                            timeout=60)
        assert listing["ok"] is True
        assert "old/merged" in [r["name"] for r in listing["candidates"]]
        refused = c.request("branch.cleanup", {"workspace": str(ws),
                                               "delete": ["old/merged"]},
                            timeout=60)
        assert refused["ok"] is False
        assert "confirm" in refused["text"].lower()
        done = c.request("branch.cleanup", {"workspace": str(ws),
                                            "delete": ["old/merged"],
                                            "confirm": True}, timeout=60)
        assert done["ok"] is True
    assert _run(ws, "branch", "--list", "old/merged") == ""


def test_branch_methods_answer_a_non_git_workspace(git_service):
    svc, _ws = git_service
    empty = svc.workspace("empty")
    with svc.client() as c:
        c.request("hello", {"protocolVersion": 1})
        for method, params in (("branch.list", {"workspace": str(empty)}),
                               ("branch.pull", {"workspace": str(empty)}),
                               ("branch.test", {"workspace": str(empty),
                                                "ref": "x"}),
                               ("branch.cleanup", {"workspace": str(empty)})):
            reply = c.request(method, params, timeout=60)
            assert reply.get("ok") is False, (method, reply)
            assert reply.get("text"), method
