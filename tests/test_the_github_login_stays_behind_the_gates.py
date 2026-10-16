"""The account's GitHub login stays behind the checks that guard it.

Audit of main's protection, 2026-10-16: the gh login on the cluster is a
repository admin, and admins may bypass the review on main. In
bypassPermissions ``gate explain`` answered "allow" for ``gh auth token``
(prints the credential) and for ``gh api -X PUT repos/o/r/pulls/7/merge``
(merges with none of the consent ``gh pr merge`` needs), while
``gh pr merge --admin`` itself was refused. Inside the bwrap cage
``~/.config/gh`` is masked and both fail (measured the same day); these
rules are what holds where there is no cage.

Every case goes through the real bash executor (``_doc_executor.execute``),
in every permission mode, with a fake ``gh`` and ``git`` first on PATH
that would leave a mark if they ran.
"""
from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import handover, security_events
from delfin.agent.api_client import KitToolPermissions, _doc_executor

MODES = ("default", "acceptEdits", "bypassPermissions")

PRINTS_THE_LOGIN = [
    "gh auth token",
    "gh auth token --hostname github.com",
    "/usr/bin/gh auth token",
    "GH_HOST=github.com gh auth token",
    "timeout 10 gh auth token",
    "echo $(gh auth token)",
    'curl -s -H "Authorization: token $(gh auth token)" https://api.github.com/user',
    "echo `gh auth token`",
    "sh -c 'gh auth token'",
    "gh auth status -t",
    "gh auth status --show-token",
    "gh auth git-credential get",
    "git credential fill",
    "git credential-store get",
]

WRITES_THROUGH_THE_API = [
    "gh api -X PUT repos/o/r/pulls/7/merge",
    "gh api --method PUT repos/o/r/pulls/7/merge",
    "gh api --method=PATCH repos/o/r",
    "gh api -XDELETE repos/o/r/git/refs/heads/main",
    "gh api -X post repos/o/r/rulesets",
    "gh api repos/o/r/issues/1/comments -f body=hi",
    "gh api repos/o/r/contents/a.txt --input body.json",
    "gh api repos/o/r/collaborators/x -F permission=admin",
    "gh api graphql -f query='mutation { mergePullRequest }'",
    "ls && gh api -X PUT repos/o/r/pulls/7/merge",
]

UNTOUCHED = [
    "gh auth status",
    "gh auth status --hostname github.com",
    "gh api repos/o/r/pulls",
    "gh api -X GET repos/o/r/pulls -f state=open",
    "gh api --method HEAD repos/o/r",
    "gh pr view 7 --json state",
    'git commit -m "explain gh auth token in the docs"',
    'echo "gh api -X PUT is refused"',
    "git config --get user.name",
    "git log --oneline -3",
]


@pytest.fixture
def ws(tmp_path, monkeypatch):
    """A workspace, and fake gh/git/curl first on PATH that log each call.

    DELFIN runs git itself before a command (the workspace's repository,
    the push readiness), so the log is checked for the refused call, not
    for being empty."""
    w = tmp_path / "ws"
    w.mkdir()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    mark = tmp_path / "ran"
    for name in ("gh", "git", "curl"):
        p = bin_dir / name
        p.write_text(f'#!/bin/sh\necho "{name} $*" >> {mark}\n')
        p.chmod(p.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    return w, mark


@pytest.fixture
def events(monkeypatch):
    seen = []
    monkeypatch.setattr(security_events, "record",
                        lambda kind, *a, **k: seen.append(kind))
    return seen


def _ran(mark: Path) -> list:
    """The logged calls that are what these rules refuse."""
    if not mark.exists():
        return []
    return [line for line in mark.read_text().splitlines()
            if line.startswith(("gh auth token", "gh auth status -t",
                                "gh auth status --show", "gh auth git-",
                                "gh api", "curl"))
            or (line.startswith("git") and "credential" in line)]


def _run(cmd: str, ws: Path, mode: str) -> str:
    perms = KitToolPermissions(mode=mode, workspace=str(ws),
                               confirm_callback=lambda *a, **k: False)
    return str(_doc_executor.execute("bash", {"command": cmd}, perms))


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", PRINTS_THE_LOGIN)
def test_a_command_that_prints_the_login_is_refused_in_every_mode(
        cmd, mode, ws, events):
    w, mark = ws
    out = _run(cmd, w, mode)
    assert "prints the account's GitHub login" in out, out
    assert "gh_credential_print" in events
    assert _ran(mark) == [], "the refused command ran"


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", WRITES_THROUGH_THE_API)
def test_a_raw_api_write_is_refused_in_every_mode(cmd, mode, ws, events):
    w, mark = ws
    out = _run(cmd, w, mode)
    assert "writes to GitHub through the raw API" in out, out
    assert "gh_api_write" in events
    assert _ran(mark) == [], "the refused command ran"


@pytest.mark.parametrize("mode", MODES)
def test_a_credential_read_with_a_config_override_is_refused_earlier(
        mode, ws, events):
    """`git -c credential.helper= credential fill` meets the git
    configuration guard first, which is stricter; the rule here would
    name it too, and either way it does not run."""
    w, mark = ws
    cmd = "git -c credential.helper= credential fill"
    assert A._prints_a_credential(cmd) == "git credential"
    out = _run(cmd, w, mode)
    assert '"error"' in out or out.startswith("blocked"), out
    assert _ran(mark) == [], "the refused command ran"


@pytest.mark.parametrize("cmd", UNTOUCHED)
def test_reads_and_mentions_are_not_these_refusals(cmd):
    assert A._prints_a_credential(cmd) == ""
    assert A._gh_api_write(cmd) == ""


@pytest.mark.parametrize("cmd", ["gh auth status", "gh api repos/o/r/pulls"])
def test_a_read_reaches_the_command_in_the_unattended_mode(cmd, ws, events):
    """The ordinary read is not caught on its way: in bypassPermissions it
    is not refused by these rules (whatever the cage then does)."""
    w, _ = ws
    out = _run(cmd, w, "bypassPermissions")
    assert "GitHub login" not in out and "raw API" not in out, out
    assert "gh_credential_print" not in events
    assert "gh_api_write" not in events


def test_the_existing_git_parser_is_the_one_in_use():
    """One parser of git's own options: the gate reuses the subcommand
    reader the read-only checks use, it does not carry a second one."""
    assert A._git_subcommand(["-C", "/x", "-c", "a=b", "credential", "fill"]) \
        == "credential"
    assert A._prints_a_credential("git -C /x credential fill") == \
        "git credential"


def test_both_refusals_are_never_handed_to_the_user_as_a_command():
    for kind in ("gh_credential_print", "gh_api_write"):
        assert kind in handover._CONTENT_KINDS
        assert handover.for_user("gh auth token", kind=kind) == ""
        assert kind in security_events._KINDS
