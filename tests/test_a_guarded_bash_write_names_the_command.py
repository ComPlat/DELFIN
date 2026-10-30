"""A guarded bash write names the command and the file it targets.

Input: a shell command (``cp scratch/f.py delfin/agent/api_client.py``)
whose write target is a protected core file. Output: the confirmation
preview the gate builds. Semantics: the preview shows the command
(credentials masked, cut to 400 characters with a marker when cut) and
the file it would write, so whoever answers sees what the yes lets run.
Before the fix the preview was only the guard warning -- no command, no
file -- and an approval one cannot read is one one cannot give.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

from delfin.agent import api_client as A


def _guarded_bash(cmd: str, workspace: Path, *, previews: list):
    (workspace / "delfin" / "agent").mkdir(parents=True, exist_ok=True)
    target = workspace / "delfin" / "agent" / "api_client.py"
    target.write_text("# stub\n")
    src = workspace / "scratch" / "new_api_client.py"
    src.parent.mkdir(exist_ok=True)
    src.write_text("# new\n")

    def cb(name, args, preview):
        previews.append(preview)
        return True

    perms = A.KitToolPermissions(
        workspace=workspace, mode="default",
        path_protected_globs=("delfin/agent/api_client.py",),
        confirm_callback=cb,
    )
    # A real session reads the file before editing it; the shell-overwrite
    # guard then lets the self-mod guard ask instead of tripping the
    # read-baseline gate.
    perms.read_tracker[str(target.resolve())] = target.stat().st_mtime
    err = A._DocToolExecutor()._gate_bash_write_targets(cmd, {}, perms)
    assert err is None  # the callback approved
    return previews


def test_guarded_bash_preview_names_the_command_and_file(tmp_path):
    previews = []
    _guarded_bash(
        "cp scratch/new_api_client.py delfin/agent/api_client.py",
        tmp_path, previews=previews)
    assert previews
    preview = previews[0]
    assert "[SELF-MODIFICATION GUARD]" in preview
    assert "cp scratch/new_api_client.py delfin/agent/api_client.py" in preview
    # The file it targets is named in the preview.
    assert "api_client.py" in preview


def test_a_credential_in_the_guarded_command_is_masked(tmp_path):
    previews = []
    key = "sk-" + "A1b2C3d4E5f6G7h8I9j0K1l2M3n4O5p6"
    _guarded_bash(
        f"OPENAI_API_KEY={key} cp scratch/new_api_client.py "
        "delfin/agent/api_client.py",
        tmp_path, previews=previews)
    assert previews
    assert key not in previews[0]
    assert "cp scratch/new_api_client.py" in previews[0]


def test_a_long_guarded_command_is_cut(tmp_path):
    previews = []
    pad = "x" * 2000
    _guarded_bash(
        f"cp scratch/new_api_client.py delfin/agent/api_client.py # {pad}",
        tmp_path, previews=previews)
    assert previews
    assert "x" * 401 not in previews[0]


def test_approvals_show_of_a_guarded_bash_request_names_command_and_file(
        tmp_path, monkeypatch, capsys):
    """`delfin-agent approvals show <id>` prints the command and the file
    for a guarded bash request, through the redirected approvals room."""
    from delfin.agent import file_confirm as _fc
    room = tmp_path / "approvals"
    room.mkdir()
    monkeypatch.setattr(_fc, "requests_dir", lambda: room)
    record = {
        "id": "g1", "session_id": "sess-x", "kind": "confirm",
        "tool": "bash",
        "command": "cp scratch/new_api_client.py delfin/agent/api_client.py",
        "preview": ("[SELF-MODIFICATION GUARD]\n"
                    "This file is part of the agent's own safety layer.\n\n"
                    "It would be run as:\n"
                    "  cp scratch/new_api_client.py delfin/agent/api_client.py"
                    "\n\n"
                    "It writes:\n  delfin/agent/api_client.py\n\n"),
        "asked_at": time.time(), "pid": os.getpid(), "host": "here",
    }
    (room / "g1.request.json").write_text(
        json.dumps(record), encoding="utf-8")
    os.chmod(room / "g1.request.json", 0o600)

    import argparse
    from delfin.agent import cli
    ns = argparse.Namespace(approvals_action="show", request_id="g1")
    rc = cli.cmd_approvals(ns)
    out = capsys.readouterr().out
    assert rc == 0
    # The command and the file it targets are both readable in `show`.
    assert "cp scratch/new_api_client.py delfin/agent/api_client.py" in out
    assert "delfin/agent/api_client.py" in out
    assert "[SELF-MODIFICATION GUARD]" in out
