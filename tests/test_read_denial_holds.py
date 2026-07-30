"""A refused read must hold against every tool, not just the one that asked.

Field case 2026-07-30: the user denied three read_file requests for files in
another project; seconds later the agent obtained all three with `cat`. Our
own prompt tells it not to work around a block — it did anyway, so the
enforcement has to be mechanical.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from delfin.agent.api_client import (
    _DocToolExecutor,
    KitToolPermissions,
    _bash_outside_reads,
    _bash_reads_denied_path,
)


def _setup():
    ws = Path(tempfile.mkdtemp()) / "ws"
    ws.mkdir(parents=True)
    outside = Path(tempfile.mkdtemp()) / "other_project"
    outside.mkdir(parents=True)
    (outside / "parser.py").write_text("SOURCE\n", encoding="utf-8")
    ex = _DocToolExecutor()
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = "bypassPermissions"      # the guard must hold even here
    return ex, perms, ws, outside


# --- pure helpers ----------------------------------------------------------

def test_denied_path_is_recognised_in_any_command():
    denied = {"/data/secret/parser.py"}
    for cmd in ("cat /data/secret/parser.py",
                "python3 -c \"print(open('/data/secret/parser.py').read())\"",
                "cp /data/secret/parser.py ./here.py"):
        assert _bash_reads_denied_path(cmd, denied), cmd


def test_denied_file_also_blocks_reaching_its_directory():
    denied = {"/data/secret/parser.py"}
    assert _bash_reads_denied_path("ls -la /data/secret", denied)
    assert _bash_reads_denied_path("grep -rn x /data/secret/", denied)


def test_unrelated_commands_are_untouched():
    denied = {"/data/secret/parser.py"}
    for cmd in ("ls -la", "pytest -q", "cat README.md",
                "cat /data/public/notes.txt"):
        assert _bash_reads_denied_path(cmd, denied) == "", cmd


def test_content_readers_with_absolute_paths_are_detected():
    assert _bash_outside_reads("cat /etc/hosts") == ["/etc/hosts"]
    assert _bash_outside_reads("head -20 /var/log/x.log") == ["/var/log/x.log"]
    assert _bash_outside_reads("base64 /home/u/key") == ["/home/u/key"]


def test_relative_and_non_reader_commands_are_not_detected():
    assert _bash_outside_reads("cat README.md") == []
    assert _bash_outside_reads("ls -la /etc") == []          # not a dump
    assert _bash_outside_reads("python3 /usr/bin/tool.py") == []
    assert _bash_outside_reads("cat /etc/*.conf") == []      # glob, not literal


# --- executor --------------------------------------------------------------

def test_refused_read_cannot_be_repeated_with_cat():
    ex, perms, _, outside = _setup()
    target = outside / "parser.py"
    perms.denied_paths.add(str(target))
    out = json.loads(ex.execute("bash", {"command": f"cat {target}"}, perms))
    assert "blocked" in out["error"]
    assert "refusal covers the data" in out["error"]
    assert "SOURCE" not in out["error"]


def test_refusal_survives_a_different_reader_and_background_bash():
    ex, perms, _, outside = _setup()
    perms.denied_paths.add(str(outside / "parser.py"))
    for tool, cmd in (("bash", f"tail -3 {outside}/parser.py"),
                      ("bash_background", f"cat {outside}/parser.py")):
        out = json.loads(ex.execute(tool, {"command": cmd}, perms))
        assert "blocked" in out["error"], cmd


def test_outside_dump_without_a_callback_is_refused_not_silently_run():
    """No confirm callback (headless): a dump outside the workspace must be
    refused with the reason, exactly as read_file would be."""
    ex, perms, _, outside = _setup()
    out = json.loads(ex.execute(
        "bash", {"command": f"cat {outside}/parser.py"}, perms))
    assert "blocked" in out["error"] and "outside" in out["error"]


def test_work_inside_the_workspace_is_unaffected():
    ex, perms, ws, _ = _setup()
    (ws / "own.txt").write_text("mine\n", encoding="utf-8")
    out = json.loads(ex.execute("bash", {"command": "cat own.txt"}, perms))
    assert out.get("exit_code") == 0 and "mine" in out.get("stdout", "")


def test_denial_ledger_records_a_real_refusal_only(monkeypatch):
    """A timeout is absence, not refusal — it must NOT poison the path."""
    ex, perms, _, outside = _setup()
    target = outside / "parser.py"

    class _Broker:
        last_timed_out = False

        def callback(self, tool, args, preview):
            return False

    broker = _Broker()
    perms.confirm_callback = broker.callback
    res = ex.execute("read_file", {"path": str(target)}, perms)
    assert "declined" in res and str(target) in perms.denied_paths

    perms2 = KitToolPermissions(workspace=str(perms.workspace))
    broker.last_timed_out = True
    perms2.confirm_callback = broker.callback
    res2 = ex.execute("read_file", {"path": str(target)}, perms2)
    assert "TIMED OUT" in res2 and not perms2.denied_paths
