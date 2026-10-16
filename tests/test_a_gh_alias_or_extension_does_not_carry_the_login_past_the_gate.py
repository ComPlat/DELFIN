"""A gh alias or extension does not carry the GitHub login past the gate.

The rules of tests/test_the_github_login_stays_behind_the_gates.py read
``gh auth token`` and ``gh api -X PUT ...`` in the command text. An alias
moves that text out of the command: after ``gh alias set t 'auth token'``
the call is ``gh t``; gh appends the call's own arguments to a plain
alias, so an alias for ``api``, ``auth`` or ``config`` turns ``gh a -X PUT
...`` or ``gh s --show-token`` into the refused command. An extension
installs a program that runs with the login and that no rule reads.

Input: ``gh alias set``, ``gh alias import`` (with a file under the
workspace) and ``gh extension install|upgrade`` commands, through the
real bash executor in every mode, with a stand-in ``gh`` first on PATH
that logs each call. Output: the refusal, the security event, and the
log.

Semantics: a plain alias whose expansion starts with a refused group
(auth, api, config, alias, extension) is refused; a shell alias
(``--shell`` or ``!``) whose body prints the login or writes through the
API, or that builds its command at run time (``$``, backticks), is
refused; a body read from stdin, or an import file that cannot be read,
is refused; ``gh extension install`` and ``upgrade`` are refused. Other
aliases and listing stay allowed.
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

ALIAS_REFUSED = [
    "gh alias set tok 'auth token'",
    "gh alias set st 'auth status'",
    "gh alias set a api",
    "gh alias set merge 'api -X PUT repos/o/r/pulls/7/merge'",
    "gh alias set cfg 'config get oauth_token'",
    "gh alias set --clobber tok 'auth token'",
    "gh alias set --shell tok 'gh auth token'",
    "gh alias set tok '!gh auth token'",
    "gh alias set -s w 'gh api -X PUT repos/o/r/pulls/7/merge'",
    "gh alias set --shell run 'echo $1 | sh'",
    "gh alias set --shell t 'eval `cat f`'",
    "echo 'auth token' | gh alias set tok -",
    "gh alias set tok - < body.txt",
    "bash -c \"gh alias set tok 'auth token'\"",
    "gh alias import -",
    "gh alias import bad-aliases.yml",
    "gh alias import missing.yml",
    "ls && gh alias set a api",
]

EXTENSION_REFUSED = [
    "gh extension install owner/gh-thing",
    "gh ext install https://github.com/owner/gh-thing",
    "gh extension install .",
    "gh extension upgrade --all",
    "gh extensions upgrade gh-thing",
]

UNTOUCHED = [
    "gh alias set pv 'pr view'",
    "gh alias set co 'pr checkout'",
    "gh alias set --shell igrep 'gh issue list --label=bug | grep flaky'",
    "gh alias list",
    "gh alias delete pv",
    "gh alias import ok-aliases.yml",
    "gh extension list",
    "gh extension remove gh-thing",
    'git commit -m "gh alias set tok auth token"',
    "echo gh extension install is refused",
]


@pytest.fixture
def ws(tmp_path, monkeypatch):
    w = tmp_path / "ws"
    w.mkdir()
    (w / "bad-aliases.yml").write_text("pv: pr view\ntok: auth token\n")
    (w / "ok-aliases.yml").write_text("pv: pr view\nco: pr checkout\n")
    (w / "body.txt").write_text("auth token\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    mark = tmp_path / "ran"
    p = bin_dir / "gh"
    p.write_text(f'#!/bin/sh\necho "gh $*" >> {mark}\n')
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
    if not mark.exists():
        return []
    return [line for line in mark.read_text().splitlines()
            if line.startswith(("gh alias set", "gh alias import",
                                "gh extension", "gh ext", "gh extensions"))]


def _run(cmd: str, ws: Path, mode: str) -> str:
    perms = KitToolPermissions(mode=mode, workspace=str(ws),
                               confirm_callback=lambda *a, **k: True)
    return str(_doc_executor.execute("bash", {"command": cmd}, perms))


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", ALIAS_REFUSED)
def test_an_alias_that_carries_the_login_is_refused(cmd, mode, ws, events):
    w, mark = ws
    out = _run(cmd, w, mode)
    assert "blocked:" in out and "alias" in out, out
    assert "gh_alias_refused" in events
    assert _ran(mark) == [], "the refused command ran"


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", EXTENSION_REFUSED)
def test_an_extension_install_is_refused(cmd, mode, ws, events):
    w, mark = ws
    out = _run(cmd, w, mode)
    assert "blocked:" in out and "extension" in out, out
    assert "gh_extension_refused" in events
    assert _ran(mark) == [], "the refused command ran"


@pytest.mark.parametrize("cmd", UNTOUCHED)
def test_other_aliases_and_listing_are_not_refused(cmd, ws):
    w, _ = ws
    assert A._gh_alias_or_extension(A._prose_blanked(cmd), str(w)) is None


def test_the_refusals_are_never_handed_to_the_user_as_a_command():
    for kind in ("gh_alias_refused", "gh_extension_refused"):
        assert kind in handover._CONTENT_KINDS
        assert kind in security_events._KINDS
