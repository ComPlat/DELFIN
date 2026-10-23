"""A contributor's change reaches main through a pull request.

Asked for on 2026-09-15: only the maintainer pushes to main; everyone else
pushes a branch and opens a pull request the maintainer accepts. GitHub's
branch rule already refuses the push for anyone without bypass rights, but
an agent learned that only from the rejected push. The gate now refuses a
contributor's push to the default branch before it leaves the machine --
whatever the user asked -- and a branch push hands over the link that opens
the pull request.
"""

from __future__ import annotations

import json
import subprocess

import pytest

from delfin.agent import api_client as A
from delfin.agent import job_monitor as jm
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture(autouse=True)
def _watch_index(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    subprocess.run(["git", "-C", str(root), "checkout", "-q", "-b", "main"],
                   check=True)
    (root / "f.txt").write_text("x\n")
    subprocess.run(["git", "-C", str(root), "add", "f.txt"], check=True)
    subprocess.run(["git", "-C", str(root), "-c", "user.email=a@b",
                    "-c", "user.name=t", "commit", "-qm", "init"], check=True)
    return root


def _granted(root):
    perms = KitToolPermissions(workspace=root, mode="bypassPermissions")
    A._grant_push_from(perms, "commit and push it", new_request=True)
    return perms


def _gate(perms, cmd):
    return _DocToolExecutor()._run_permission_gate(
        "bash", {"command": cmd}, perms)


@pytest.mark.parametrize("cmd", [
    "git push origin main",
    "git push origin HEAD:main",
    "git push",
    "git push -u origin",
    "git add f.txt && git commit -m fix && git push origin main",
    "git push --force-with-lease origin main",
])
def test_a_contributor_cannot_push_to_main_even_when_asked(repo, monkeypatch, cmd):
    monkeypatch.setattr(A, "_git_role", lambda: "contributor")
    err = _gate(_granted(repo), cmd)
    assert err and "pull request" in err, err


def test_a_contributor_pushes_a_branch(repo, monkeypatch):
    monkeypatch.setattr(A, "_git_role", lambda: "contributor")
    assert _gate(_granted(repo), "git push -u origin mh/viewer-notice") is None


def test_the_maintainer_pushes_to_main(repo, monkeypatch):
    monkeypatch.setattr(A, "_git_role", lambda: "maintainer")
    assert _gate(_granted(repo), "git push origin main") is None


def test_a_branch_push_hands_over_the_pull_request_link(repo):
    perms = KitToolPermissions(workspace=repo, mode="bypassPermissions")
    note = A._after_push(
        {"command": "git push -u origin mh/viewer-notice"},
        json.dumps({"exit_code": 0, "stdout": "", "stderr": (
            "To github.com:ComPlat/DELFIN.git\n"
            " * [new branch]      mh/viewer-notice -> mh/viewer-notice\n")}),
        perms)
    assert ("https://github.com/ComPlat/DELFIN/compare/main...mh/viewer-notice"
            "?expand=1") in note


def test_contributor_is_the_default_role(monkeypatch):
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: {})
    assert A._git_role() == "contributor"
    monkeypatch.setattr(us, "load_settings",
                        lambda *a, **k: {"agent": {"git_role": "Maintainer"}})
    assert A._git_role() == "maintainer"
    monkeypatch.setattr(us, "load_settings",
                        lambda *a, **k: {"agent": {"git_role": "owner"}})
    assert A._git_role() == "contributor"


def test_the_targets_a_push_names():
    assert A._push_targets("git push origin main", ".") == {"main"}
    assert A._push_targets("git push origin HEAD:refs/heads/main", ".") == {"main"}
    assert A._push_targets("git push -u origin mh/topic", ".") == {"mh/topic"}
    assert A._push_targets("git push origin a b:c", ".") == {"a", "c"}
