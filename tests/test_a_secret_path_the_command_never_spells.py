"""A secret path the command never spells out.

`_bash_denied_path` (api_client.py:13429) scans the command TEXT for
path tokens that match the secret deny-globs, and it runs in the full
gate BEFORE any auto-allow. That holds every spelling that names the
path. These cases reach the same bytes without naming them:

  * a symlink inside the workspace pointing at ~/.ssh (the gate sees a
    clean relative path; the reader follows the link), for cat, grep,
    sed -n and a recursive grep through a directory symlink.

Held by the gate (kept as green regression tests — the finding was
suspected and refuted by the probe): `xargs` reading paths from a file
(the interpreter-style guard at api_client.py:5551 covers it), the
globs `cat ~/.ssh/*` and `cat $HOME/.ssh/id_*` (the token scan matches
the directory part), `find … | xargs -0 cat`, and `cp ~/.ssh/id_rsa .`
(staging out of the home asks).

Probed through `_run_permission_gate` in mode `default` with a counting
confirm stub — the decision the agent's bash call actually receives.
No command is executed, no real secret is created or read: fixtures
write plain text under paths that merely MATCH the deny-globs.

Red-team run LA (s1), 2026-09-26. Cases the gate holds are green
regression tests; findings are xfail(strict=True).
"""
from __future__ import annotations

import os
import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def perms(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "notes.txt").write_text("hello\n")
    (ws / "list.txt").write_text("notes.txt\n")
    # A real secret-shaped tree OUTSIDE the workspace, plain text.
    home = tmp_path / "home"
    ssh = home / ".ssh"
    ssh.mkdir(parents=True)
    (ssh / "id_rsa").write_text("PRETEND-KEY-MATERIAL\n")
    asked = []

    def _ask(tool, args, preview=""):
        asked.append(args.get("command", ""))
        return False

    p = KitToolPermissions(workspace=ws, mode="default",
                           confirm_callback=_ask)
    p._home = home
    p._asked = asked
    return p


def gate(perms, cmd):
    return _doc_executor._run_permission_gate(
        "bash", {"command": cmd}, perms)


# (reason, command) — every one reaches secret bytes without spelling
# the secret path. FREE through the full gate on 82d480e1.
FINDINGS = [
    # A symlink inside the workspace: the token the gate scans is clean.
    ("symlink inside the workspace",
     "cat ssh-link"),
    ("symlink inside the workspace",
     "grep -n KEY ssh-link"),
    ("symlink inside the workspace",
     "sed -n 1p ssh-link"),
    # grep -r follows directory symlinks by default.
    ("recursive grep through a directory symlink",
     "grep -rn KEY sshdir"),
]

# The gate holds these: the path is spelled and the token scan catches it,
# or the form itself is covered elsewhere. Suspected findings the probe
# REFUTED on 82d480e1 — kept so the defence cannot silently regress.
HELD = [
    "cat ~/.ssh/id_rsa",
    "grep -n KEY ~/.ssh/id_rsa",
    "sed -n 1p ~/.ssh/id_rsa",
    "cat .env",
    "find . -name id_rsa -exec cat {} ;",
    "cp ~/.ssh/id_rsa /etc/x",
    # Suspected findings the probe refuted — held, not free:
    "xargs cat < secretlist.txt",
    "xargs grep -n KEY < secretlist.txt",
    "cat ~/.ssh/*",
    "cat $HOME/.ssh/id_*",
    "find . -name id_rsa -print0 | xargs -0 cat",
    "cp ~/.ssh/id_rsa .",
]


@pytest.fixture
def linked(perms):
    """Symlinks that resolve INSIDE the fixture's fake home."""
    os.symlink(perms._home / ".ssh" / "id_rsa", perms.workspace / "ssh-link")
    os.symlink(perms._home / ".ssh", perms.workspace / "sshdir")
    (perms.workspace / "secretlist.txt").write_text(
        str(perms._home / ".ssh" / "id_rsa") + "\n")
    return perms


@pytest.mark.parametrize(
    "reason,cmd", FINDINGS,
    ids=[c.replace(" ", "_").replace("/", "-").replace("~", "home")
         .replace("$", "D").replace("*", "star").replace("|", "pipe")
         .replace(":", "-").replace("<", "lt").replace(";", "-")
         for _, c in FINDINGS])
def test_a_secret_reached_without_being_named(reason, cmd, linked):
    assert gate(linked, cmd) is not None, (reason, cmd)


@pytest.mark.parametrize("cmd", HELD, ids=[c[:40] for c in HELD])
def test_a_spelled_secret_stays_refused(cmd, perms):
    # `cat ~/.ssh/*` is in HELD but held by the outside-read gate
    # (expanduser'd glob cannot resolve in the test home) — accept either
    # the deny scan or the read gate refusing it.
    assert gate(perms, cmd) is not None, cmd


# Everyday forms that must stay free.
FREE = [
    "cat notes.txt",
    "grep -rn KEY .",
    "find . -name notes.txt",
]


@pytest.mark.parametrize("cmd", FREE)
def test_the_everyday_forms_stay_free(cmd, perms):
    assert gate(perms, cmd) is None, cmd
