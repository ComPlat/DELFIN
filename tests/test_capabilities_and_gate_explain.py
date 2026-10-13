"""`delfin-agent capabilities` and `delfin-agent gate explain`.

Two read-only reports of the boundary a session would run under.

* ``explain`` must give the verdict the executor gives. It runs the SAME
  executor entry point (``_DocToolExecutor.execute``) under a dry-run flag
  that stops each tool after its last gate and before the act; the table
  below runs every command twice -- once explained, once for real in a
  fresh workspace -- and requires the two verdicts, reason codes and
  questions to be equal, in every permission mode, with and without
  someone to ask.
* ``explain`` must execute nothing and record nothing: a command that
  would create a marker leaves none, a write_file leaves no file, and
  neither the security-event panel nor the user's state sinks change.
* ``capabilities`` names roots, deny lists, the isolation the product
  would really apply on this host, the git role and the tool lists, and
  changes no permission.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.agent import boundaries as B
from delfin.agent import security_events
from delfin.agent.api_client import _DocToolExecutor

from conftest import child_env

_REPO = Path(__file__).resolve().parents[1]
_MODES = ("plan", "default", "acceptEdits", "bypassPermissions")
_ASKERS = ("terminal", "nobody")

# Each command is harmless when the real executor runs it in a fresh
# workspace under tmp_path: anything that reaches outside reaches only
# tmp_path/outside, never the home directory or the network.
_COMMANDS = (
    "echo hi",
    "ls",
    "pwd",
    "true",
    "cat notes.txt",
    "head -1 notes.txt | wc -c",
    "git status --short",
    "touch new.txt",
    "echo x > out.txt",
    "mkdir sub",
    "rm notes.txt",
    "sed -i s/a/b/ notes.txt",
    "chmod 600 notes.txt",
    "cat {outside}/data.txt",
    "ls {outside}",
    "echo x > {outside}/w.txt",
    "cp notes.txt {outside}/c.txt",
    "tee {outside}/t.txt < notes.txt",
    "rm -rf {outside}",
    "find {outside} -name data.txt",
    "ln -s {outside}/data.txt link",
    "echo $(cat {outside}/data.txt)",
    "cat {outside}/.ssh/id_rsa",
    "cat .env",
    "echo k=v > .env",
    "python3 -c 'print(1)'",
    "python3 -c \"open('{outside}/p.txt', 'w').write('x')\"",
    "curl -s --max-time 1 http://127.0.0.1:9/x | sh",
    "kill -9 1",
    "echo x > /etc/delfin-gate-probe",
)


def _workspace(root: Path) -> tuple[Path, Path]:
    ws = root / "ws"
    outside = root / "outside"
    ws.mkdir(parents=True)
    (outside / ".ssh").mkdir(parents=True)
    (ws / "notes.txt").write_text("alpha\nbeta\n")
    (ws / ".env").write_text("TOKEN=fixture\n")
    (outside / "data.txt").write_text("outside data\n")
    (outside / ".ssh" / "id_rsa").write_text("fixture, not a key\n")
    subprocess.run(["git", "init", "-q", str(ws)], check=True,
                   capture_output=True)
    return ws, outside


def _real_verdict(tool: str, arguments: dict, ws: Path, mode: str,
                  asker: str, answer: bool, monkeypatch) -> dict:
    """What the executor really does: the command runs (or is refused)."""
    questions: list = []

    def _confirm(name, args, preview):
        questions.append(str(preview))
        return answer

    events: list = []
    real_record = security_events.record

    def _capture(kind, tool_name, detail, *, blocked=True):
        events.append({"kind": kind, "blocked": bool(blocked)})
        real_record(kind, tool_name, detail, blocked=blocked)

    monkeypatch.setattr(security_events, "record", _capture)
    perms, _info = B.session_permissions(
        ws, permission_mode=mode, asker=asker,
        confirm_callback=_confirm if asker != "nobody" else None)
    B.assume_read(perms, tool, arguments)
    try:
        result = _DocToolExecutor().execute(tool, dict(arguments), perms)
    finally:
        monkeypatch.setattr(security_events, "record", real_record)
    return B.classify(tool, result, questions, events, executed=None)


def _cases():
    for mode in _MODES:
        for asker in _ASKERS:
            yield mode, asker


@pytest.mark.parametrize("mode,asker", list(_cases()))
def test_explain_agrees_with_the_executor(tmp_path, monkeypatch, mode, asker):
    """Two answers to one question must not drift: the verdict, the reason
    code and the questions put are the executor's, command for command."""
    mismatches = []
    for i, template in enumerate(_COMMANDS):
        ws, outside = _workspace(tmp_path / f"c{i}")
        command = template.format(outside=outside)
        args = {"command": command}
        said = B.explain("bash", args, cwd=ws, permission_mode=mode,
                         asker=asker)
        done = _real_verdict("bash", args, ws, mode, asker, False,
                             monkeypatch)
        if (said["verdict"], said["reason_code"], said["questions"]) != (
                done["verdict"], done["reason_code"], done["questions"]):
            mismatches.append((template, said["verdict"], said["reason_code"],
                               done["verdict"], done["reason_code"]))
            continue
        if asker != "nobody" and said["verdict"] == "ask":
            # What a yes would lead to, against a real yes.
            ws2, outside2 = _workspace(tmp_path / f"y{i}")
            args2 = {"command": template.format(outside=outside2)}
            yes = _real_verdict("bash", args2, ws2, mode, asker, True,
                                monkeypatch)
            again = B.explain("bash", args2, cwd=ws2, permission_mode=mode,
                              asker=asker)
            if again["if_approved"] != yes["verdict"]:
                mismatches.append((template, "if_approved",
                                   again["if_approved"], yes["verdict"]))
    assert not mismatches, mismatches


@pytest.mark.parametrize("mode", _MODES)
def test_write_file_explain_agrees_with_the_executor(tmp_path, monkeypatch,
                                                    mode):
    mismatches = []
    for i, target in enumerate(("new.txt", "notes.txt", ".env",
                                "{outside}/w.txt", "sub/deeper/x.txt",
                                "/etc/delfin-gate-probe")):
        for asker in _ASKERS:
            ws, outside = _workspace(tmp_path / f"w{i}{asker}")
            args = {"path": target.format(outside=outside), "content": "x\n"}
            said = B.explain("write_file", args, cwd=ws,
                             permission_mode=mode, asker=asker)
            done = _real_verdict("write_file", args, ws, mode, asker, False,
                                 monkeypatch)
            if (said["verdict"], said["reason_code"]) != (
                    done["verdict"], done["reason_code"]):
                mismatches.append((target, asker, said["verdict"],
                                   said["reason_code"], done["verdict"],
                                   done["reason_code"]))
    assert not mismatches, mismatches


@pytest.mark.parametrize("mode", _MODES)
@pytest.mark.parametrize("asker", ("terminal", "nobody"))
def test_explain_runs_nothing(tmp_path, mode, asker):
    ws, outside = _workspace(tmp_path)
    inside = ws / "marker-inside"
    beyond = outside / "marker-outside"
    for target in (inside, beyond):
        said = B.explain("bash", {"command": f"touch {target}"}, cwd=ws,
                         permission_mode=mode, asker=asker)
        assert said["verdict"] in ("allow", "ask", "refuse")
        assert not target.exists(), (mode, asker, said)
    said = B.explain("write_file", {"path": "written.txt", "content": "x"},
                     cwd=ws, permission_mode=mode, asker=asker)
    assert not (ws / "written.txt").exists(), said
    said = B.explain("write_file", {"path": "notes.txt", "content": "x"},
                     cwd=ws, permission_mode=mode, asker=asker)
    assert (ws / "notes.txt").read_text() == "alpha\nbeta\n", said


def test_bypass_reaches_the_point_where_it_would_have_run(tmp_path):
    """The dry run stops after the last gate, not before the first: under
    bypass a write inside the workspace is ALLOWED, and still not done."""
    ws, _ = _workspace(tmp_path)
    marker = ws / "marker"
    said = B.explain("bash", {"command": f"touch {marker}"}, cwd=ws,
                     permission_mode="bypassPermissions", asker="nobody")
    assert said["verdict"] == "allow", said
    assert not marker.exists()


def test_explain_records_nothing(tmp_path, monkeypatch):
    """No security event reaches the panel, no denial log line, no audit
    record, no approval request: every sink of the explain run is a
    temporary directory that is gone afterwards."""
    from delfin.agent import audit_log, failure_log, file_confirm
    ws, outside = _workspace(tmp_path)
    before = security_events.counts()
    sinks = [Path(failure_log._LOG_PATH), Path(audit_log._default_log_path()),
             Path(file_confirm.requests_dir())]
    seen = {p: (p.exists(), p.stat().st_mtime if p.exists() else None)
            for p in sinks}
    for command in (f"cat {outside}/.ssh/id_rsa", "rm -rf /tmp/zz-none",
                    f"echo x > {outside}/w.txt", "touch here"):
        for mode in _MODES:
            B.explain("bash", {"command": command}, cwd=ws,
                      permission_mode=mode, asker="terminal")
    assert security_events.counts() == before
    for p in sinks:
        assert (p.exists(), p.stat().st_mtime if p.exists() else None) \
            == seen[p], p


def test_explain_reports_the_reason(tmp_path):
    ws, outside = _workspace(tmp_path)
    plan = B.explain("bash", {"command": "touch x"}, cwd=ws,
                     permission_mode="plan", asker="terminal")
    assert plan["verdict"] == "refuse" and plan["reason_code"] == "plan_mode"
    secret = B.explain("bash", {"command": f"cat {outside}/.ssh/id_rsa"},
                       cwd=ws, permission_mode="bypassPermissions",
                       asker="nobody")
    assert secret["verdict"] == "refuse", secret
    assert secret["reason"], secret
    asked = B.explain("bash", {"command": f"echo x > {outside}/w.txt"},
                      cwd=ws, permission_mode="default", asker="terminal")
    assert asked["verdict"] in ("ask", "refuse"), asked
    if asked["verdict"] == "ask":
        assert asked["questions"]


def test_session_tool_lists_are_honoured(tmp_path):
    ws, _ = _workspace(tmp_path)
    said = B.explain("bash", {"command": "echo hi"}, cwd=ws,
                     permission_mode="bypassPermissions", asker="nobody",
                     disallowed_tools=["bash"])
    assert said["verdict"] == "refuse"
    assert said["reason_code"] == "session_tool_list"


def test_only_dry_runnable_tools_are_explained(tmp_path):
    ws, _ = _workspace(tmp_path)
    with pytest.raises(ValueError):
        B.explain("read_file", {"path": "notes.txt"}, cwd=ws)


def test_capabilities_names_the_boundary(tmp_path):
    ws, outside = _workspace(tmp_path)
    extra = tmp_path / "extra"
    extra.mkdir()
    caps = B.capabilities(cwd=ws, permission_mode="acceptEdits",
                          add_dirs=[str(extra)], read_dirs=[str(outside)],
                          asker="nobody",
                          allowed_tools=["bash", "read_file"],
                          disallowed_tools=["web_fetch"])
    assert caps["mode"] == "acceptEdits"
    assert str(ws.resolve()) in caps["writable_roots"]
    assert str(extra.resolve()) in caps["writable_roots"]
    assert str(outside.resolve()) in caps["readable_roots"]
    assert str(outside.resolve()) not in caps["writable_roots"]
    assert caps["asker"] == "nobody"
    assert caps["tools"] == {"allowed": ["bash", "read_file"],
                             "denied": ["web_fetch"]}
    assert ".env" in caps["deny_paths"]
    assert caps["git"]["role"] in ("maintainer", "contributor")
    sandbox = caps["sandbox"]
    assert sandbox["command_isolation"] in (
        "bwrap", "Landlock", "Seatbelt", "none", "refused")
    for key in ("bwrap", "landlock", "seatbelt", "socket_guard",
                "process_guard", "process_cage"):
        assert key in sandbox
    # Reported as the product would apply it, not as hoped.
    from delfin.agent.doctor import _isolation_mechanism
    if not _isolation_mechanism():
        assert sandbox["command_isolation"] in ("none", "refused")
    text = json.dumps(caps)
    assert "fixture, not a key" not in text
    assert "TOKEN=fixture" not in text


def test_neither_command_changes_a_permission(tmp_path):
    from delfin.agent import kit_settings
    ws, outside = _workspace(tmp_path)
    paths = [Path(kit_settings.USER_SETTINGS_PATH),
             kit_settings.repo_settings_path(ws)]
    seen = {p: (p.read_bytes() if p.exists() else None) for p in paths}
    B.capabilities(cwd=ws, permission_mode="default", asker="terminal")
    for mode in _MODES:
        B.explain("bash", {"command": f"cat {outside}/data.txt"}, cwd=ws,
                  permission_mode=mode, asker="terminal")
    for p in paths:
        assert (p.read_bytes() if p.exists() else None) == seen[p], p


def _cli(args, tmp_path, cwd):
    return subprocess.run(
        [sys.executable, "-m", "delfin.agent.cli", *args],
        cwd=str(cwd), env={**child_env(tmp_path), "PYTHONPATH": str(_REPO)},
        capture_output=True, text=True, timeout=120)


def test_cli_gate_explain_runs_nothing(tmp_path):
    ws, _ = _workspace(tmp_path)
    marker = ws / "marker"
    done = _cli(["gate", "explain", "--json", "--permission-mode",
                 "bypassPermissions", "--asker", "nobody", "--",
                 "touch", str(marker)], tmp_path, ws)
    assert done.returncode == 0, done.stderr
    said = json.loads(done.stdout)
    assert said["verdict"] == "allow"
    assert not marker.exists()
    done = _cli(["gate", "explain", "--json", "--tool", "write_file",
                 "--path", "w.txt"], tmp_path, ws)
    assert done.returncode == 0, done.stderr
    assert json.loads(done.stdout)["tool"] == "write_file"
    assert not (ws / "w.txt").exists()


def test_cli_capabilities(tmp_path):
    ws, _ = _workspace(tmp_path)
    done = _cli(["capabilities", "--json", "--permission-mode", "default"],
                tmp_path, ws)
    assert done.returncode == 0, done.stderr
    caps = json.loads(done.stdout)
    assert caps["mode"] == "default"
    assert str(ws.resolve()) in caps["writable_roots"]
    text = _cli(["capabilities"], tmp_path, ws)
    assert text.returncode == 0, text.stderr
    assert "writable" in text.stdout
