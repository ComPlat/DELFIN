"""A push with no usable credential fails at once instead of hanging.

Measured on the shared bug archive, 2026-09-16 .. 10-14: 63 push and
`gh pr create` calls, 9 of them timeouts (37 min), and two sessions that
spent 10 and 8 attempts on transports, proxies and temp directories.
The captured output named the cause each time:

    mktemp: failed to create file via template '/scratch/tmp.XXXX':
        Read-only file system
    Missing or invalid credentials.
    Error: connect EACCES /run/user/<uid>/vscode-git-<id>.sock

Reproduced on main 2026-10-14 through the real executor (bubblewrap
cage, the editor's GIT_ASKPASS inherited, TMPDIR outside the writable
roots): `git push --dry-run` to an HTTPS remote ran into the 60 s limit.
The editor's askpass script makes a temp file with mktemp, which fails on
the read-only TMPDIR, and then runs `cat` on the empty name -- reading
the remote helper's stdin, a pipe that stays open: it waits forever.

  inside the cage                 GIT_TERMINAL_PROMPT=0, no editor askpass
                                  (its socket is refused there anyway),
                                  TMPDIR writable
  a push that still waits         ended at PUSH_TIMEOUT_S with the cause
                                  and who can fix it
  a push that fails               the push_diagnosis line beside it
"""

from __future__ import annotations

import json
import os
import stat
import subprocess
import time

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture
def repo(tmp_path, monkeypatch):
    ws = tmp_path / "ws"
    ws.mkdir()
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    subprocess.run(["git", "init", "-q", str(ws)], check=True, env=env)
    subprocess.run(["git", "-C", str(ws), "-c", "user.name=t", "-c",
                    "user.email=t@example.org", "commit", "-q",
                    "--allow-empty", "-m", "c"], check=True, env=env)
    bindir = ws / "bin"
    bindir.mkdir()
    # The editor's askpass, as shipped: mktemp, the helper, then cat.
    askpass = bindir / "askpass.sh"
    askpass.write_text(
        "#!/bin/sh\nP=`mktemp`\n\"$VSCODE_GIT_ASKPASS_NODE\" x $*\n"
        "cat $P\nrm -f $P\n", encoding="utf-8")
    askpass.chmod(askpass.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("GIT_ASKPASS", str(askpass))
    monkeypatch.setenv("VSCODE_GIT_ASKPASS_NODE", "/bin/false")
    monkeypatch.setenv("VSCODE_GIT_IPC_HANDLE",
                       "/run/user/0/vscode-git-0000000000.sock")
    monkeypatch.delenv("GIT_TERMINAL_PROMPT", raising=False)
    # Writable on the host, outside every writable root of the cage.
    outside = tmp_path / "host-tmp"
    outside.mkdir()
    monkeypatch.setenv("TMPDIR", str(outside))
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    perms = KitToolPermissions(workspace=ws, mode="bypassPermissions")
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    eng._permissions = perms

    def _run(cmd, **extra):
        return json.loads(eng._execute_bash({"command": cmd, **extra}, perms))
    return ws, bindir, _run


def test_the_cage_has_no_prompt_and_no_editor_askpass(repo):
    ws, bindir, run = repo
    out = run("printenv GIT_TERMINAL_PROMPT; echo \"[${GIT_ASKPASS:-none}]\";"
              " echo \"[${VSCODE_GIT_IPC_HANDLE:-none}]\"")
    lines = out["stdout"].split()
    assert lines[0] == "0"
    assert lines[1] == "[none]"
    assert lines[2] == "[none]"


def test_a_temp_file_can_be_made(repo):
    ws, bindir, run = repo
    out = run("f=$(mktemp) && echo made && rm -f \"$f\"")
    assert out["exit_code"] == 0, out
    assert "made" in out["stdout"]


def test_a_push_that_waits_is_ended_and_says_why(repo, monkeypatch):
    ws, bindir, run = repo
    helper = bindir / "git-remote-hang"
    helper.write_text("#!/bin/sh\nsleep 60\n", encoding="utf-8")
    helper.chmod(helper.stat().st_mode | stat.S_IXUSR)
    subprocess.run(["git", "-C", str(ws), "remote", "add", "far", "hang::x"],
                   check=True)
    monkeypatch.setattr(A, "PUSH_TIMEOUT_S", 2)
    t0 = time.monotonic()
    out = run("git push far HEAD:refs/heads/x", timeout_s=300)
    assert time.monotonic() - t0 < 20
    err = out.get("error", "")
    assert "credential" in err
    assert "user" in err


def test_a_push_that_fails_gets_the_diagnosis(repo):
    ws, bindir, run = repo
    helper = bindir / "git-remote-nocred"
    helper.write_text(
        "#!/bin/sh\necho \"fatal: could not read Username for "
        "'https://example.org': terminal prompts disabled\" >&2\nexit 128\n",
        encoding="utf-8")
    helper.chmod(helper.stat().st_mode | stat.S_IXUSR)
    subprocess.run(["git", "-C", str(ws), "remote", "add", "far",
                    "nocred::x"], check=True)
    out = run("git push far HEAD:refs/heads/x")
    assert out["exit_code"] != 0
    assert "credential" in out.get("hint", "")
