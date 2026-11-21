# SPDX-License-Identifier: MIT
"""Everything the uncaged push route must refuse, as a table.

``_github_push_plan`` decides whether a ``git push`` runs OUTSIDE the
process cage: ``_run_github_push`` then calls ``subprocess.run`` on the
plan's argv directly, with no cage around it. The route exists because a
push needs a credential that the cage masks (field report 2026-10-14), and
it is safe only to the extent that the parse is: there is no sandbox behind
it to catch a mistake.

Every refusal it relies on was readable in the code and **none of them was
measured**. The three existing tests cover the gh route, the fallback when
``gh`` is missing, a checkout that configures the transport, and that no
caged command sees the token. None covers the flag allowlist, the URL
allowlist, a delete, a force, or a shell construct -- so an edit to
``_PLAIN_PUSH_FLAGS`` or to the pattern would widen the uncaged route with
CI green.

Input: one command string, or one remote URL, against a throwaway repo.
Output: whether a plan comes back. A plan means the command runs uncaged,
so each row asserts the direction the design claims, in both directions --
a table that only refused things would pass if the parser refused
everything.

Found by writing this table: a newline separated tokens like a space,
because the pattern used ``\\s``. ``git push origin topic\\nmain`` was
accepted and ran as ``git push origin topic main``, pushing a branch nobody
asked for, while the second line was neither run nor refused. No arbitrary
command came through -- the flag allowlist refuses ``rm -rf /tmp/x`` -- but
an unrequested publish did.
"""
from __future__ import annotations

import subprocess

import pytest

from delfin.agent import api_client as A

_GITHUB_SSH = "git@github.com:ComPlat/DELFIN.git"


def _git(cwd, *args):
    return subprocess.run(["git", *args], cwd=str(cwd), capture_output=True,
                          text=True, check=False)


def _repo(root, url: str):
    """A repo on branch ``topic`` whose ``origin`` push URL is ``url``."""
    ws = root / "ws"
    ws.mkdir(parents=True)
    _git(ws, "init", "-q", "-b", "main")
    _git(ws, "config", "user.email", "fixture@example.invalid")
    _git(ws, "config", "user.name", "fixture")
    (ws / "f.txt").write_text("x\n")
    _git(ws, "add", "f.txt")
    _git(ws, "commit", "-qm", "first")
    _git(ws, "remote", "add", "origin", url)
    _git(ws, "checkout", "-q", "-b", "topic")
    return ws


#: (command, may_run_uncaged). The False rows are the security claim; the
#: True rows keep the parser from passing this file by refusing everything.
_COMMANDS = [
    ("git push origin topic", True),
    ("git push -u origin topic", True),
    ("git push --dry-run origin topic", True),
    ("git push --quiet --tags origin topic", True),
    ("git push origin topic 2>&1", True),
    ("git push origin topic 2>&1 | tail -n 20", True),
    # Flags that make git run a command, or move where it pushes. None is in
    # _PLAIN_PUSH_FLAGS, and none contains a character the pattern allows.
    ("git push --receive-pack=/tmp/x origin topic", False),
    ("git push --exec=/tmp/x origin topic", False),
    ("git push --upload-pack=/tmp/x origin topic", False),
    ("git push --repo=https://evil.invalid/a/b origin topic", False),
    ("git push --push-option=x origin topic", False),
    # Irreversible or history-rewriting forms: the cage is where they belong,
    # because the gate that asks about them runs there.
    ("git push origin :topic", False),
    ("git push --delete origin topic", False),
    ("git push origin +topic", False),
    ("git push --force origin topic", False),
    ("git push --force-with-lease origin topic", False),
    ("git push --mirror origin", False),
    ("git push --all origin", False),
    # Shell constructs. The route runs argv with no shell, so anything that
    # needs a shell must not reach it -- it would otherwise be dropped
    # silently rather than executed.
    ("git push origin topic && rm -rf /tmp/x", False),
    ("git push origin topic; id", False),
    ("git push origin $(id)", False),
    ("git push origin `id`", False),
    ("git push origin topic > /tmp/x", False),
    ("git push origin topic | sh", False),
    # A second LINE: the defect this table found.
    ("git push origin topic\nid", False),
    ("git push origin topic\nmain", False),
    ("git push origin topic\rmain", False),
    ("git push\norigin topic", False),
    # A remote that is not configured here.
    ("git push elsewhere topic", False),
]


@pytest.mark.parametrize("cmd,may_run_uncaged", _COMMANDS,
                         ids=[c.replace("\n", "\\n").replace("\r", "\\r")
                              for c, _ in _COMMANDS])
def test_only_a_plain_push_leaves_the_cage(tmp_path, cmd, may_run_uncaged):
    ws = _repo(tmp_path, _GITHUB_SSH)
    plan = A._github_push_plan(cmd, ws)
    if may_run_uncaged:
        assert plan is not None, (
            "a plain push must still take the credentialed route, or every "
            "push fails in the cage with the login the user just made")
        assert plan["argv"][0] == "git"
    else:
        assert plan is None, (
            f"{cmd!r} would run OUTSIDE the cage as {plan and plan['argv']}")


#: (url, may_run_uncaged). Only github.com over ssh or https may.
_URLS = [
    (_GITHUB_SSH, True),
    ("https://github.com/ComPlat/DELFIN.git", True),
    ("git@github.com.evil.invalid:a/b.git", False),
    ("https://github.com.evil.invalid/a/b.git", False),
    ("https://evil.invalid/github.com/a/b.git", False),
    ("git://github.com/a/b.git", False),
    ("file:///tmp/other", False),
    ("/tmp/other", False),
    ("ext::sh -c id", False),
    ("git@gitlab.example.invalid:a/b.git", False),
]


@pytest.mark.parametrize("url,may_run_uncaged", _URLS)
def test_only_a_github_remote_leaves_the_cage(tmp_path, url, may_run_uncaged):
    """The URL decides, not the remote's name: 'origin' proves nothing."""
    ws = _repo(tmp_path, url)
    plan = A._github_push_plan("git push origin topic", ws)
    assert (plan is not None) is may_run_uncaged, (
        f"remote url {url!r} -> {'uncaged' if plan else 'cage'}")


def test_the_second_line_is_not_folded_into_the_refspecs(tmp_path):
    """The defect, pinned as its own case rather than only as a row.

    A two-line command is two commands. Folding the second line into the
    push turns it into a refspec -- so a push of 'topic' also published
    'main' -- and the line itself then ran nowhere. Both halves are wrong,
    and the refusal has to be visible as one claim.
    """
    ws = _repo(tmp_path, _GITHUB_SSH)
    assert A._github_push_plan("git push origin topic\nmain", ws) is None
    # And the single-line form it must not be confused with still works.
    plan = A._github_push_plan("git push origin topic main", ws)
    assert plan is not None and plan["argv"][-2:] == ["topic", "main"]


def test_the_pattern_separates_on_spaces_and_tabs_only(tmp_path):
    """Why the fix is in the separator and not in a blocklist of characters.

    A tab is a separator a shell accepts, so it stays one here. Every other
    whitespace character ends the command, which is what a shell does too.
    """
    ws = _repo(tmp_path, _GITHUB_SSH)
    assert A._github_push_plan("git\tpush\torigin\ttopic", ws) is not None
    for ws_char in ("\n", "\r", "\v", "\f"):
        assert A._github_push_plan(f"git push origin{ws_char}topic", ws) is None, (
            f"{ws_char!r} separated two tokens")
