"""A pending approval is answered by the human, never by an agent.

Input: a request waiting in each answer room -- ``file_confirm``'s
approvals directory, ``terminal_confirm``'s published-question directory
-- and the CLI backend's gate state file (``cli_gate``), all under a
private HOME. Each route an agent session has to answer one is driven
through the real gates: the shell command (``delfin-agent approvals
approve|deny|answer`` and its other spellings, the Python functions that
write an answer, a direct write of the answer file), the write tools,
and the CLI backend's PreToolUse hook. Every permission mode is tried,
with shell isolation off and on.

Output: every route is refused before anything is written, with a
``blocked:`` message that tells the model to show the request to the
user; no dialog is raised for it (the person would be asked to approve
the agent approving itself). The second layer, ``cmd_approvals``, refuses
when it runs with the marker every agent child process carries.

Controls: the same command in a plain shell (no marker) still answers;
``approvals ls`` and ``show`` from the agent are not refused by this
guard; a write elsewhere in the workspace still works.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import cli_gate
from delfin.agent import file_confirm as fc
from delfin.agent import terminal_confirm as tc


MARKER = "DELFIN_AGENT_SESSION"


class Yes:
    """A dialog that approves everything: if the gate asks, it passes."""

    last_timed_out = False

    def __init__(self):
        self.shown = []

    def ask(self, name, args, preview):
        self.shown.append((name, preview))
        return True


@pytest.fixture()
def home(tmp_path, monkeypatch):
    """A private HOME whose rooms are the ones every module and every
    child process resolves."""
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setenv("HOME", str(h))
    monkeypatch.delenv("DELFIN_SCRATCH_STATE", raising=False)
    monkeypatch.delenv(MARKER, raising=False)
    approvals = h / ".delfin" / "approvals"
    terminal = h / ".delfin" / "terminal_confirmations"
    gate = h / ".delfin" / "cli_gate"
    monkeypatch.setattr(fc, "requests_dir", lambda: approvals)
    monkeypatch.setattr(tc, "_PENDING_DIR", terminal)
    monkeypatch.setattr(cli_gate, "state_dir", lambda: gate)
    monkeypatch.setattr(A, "_classifier_settings", lambda: None)
    return h


@pytest.fixture()
def ws(tmp_path):
    w = tmp_path / "ws"
    w.mkdir()
    return w.resolve()


def _pending(home: Path) -> tuple[str, Path]:
    broker = fc.FileConfirmBroker(session_id="waiting-session", timeout_s=30)
    rid, path = broker._write_request(
        "write_file", {"path": "x"}, "[SELF-MODIFICATION GUARD] ...")
    return rid, path


MODES = [
    pytest.param("plan", True, id="plan"),
    pytest.param("default", True, id="default"),
    pytest.param("acceptEdits", True, id="acceptEdits"),
    pytest.param("bypassPermissions", True, id="bypass-attended"),
    pytest.param("bypassPermissions", False, id="bypass-unattended"),
]


def _perms(ws, mode, attended, **kw):
    dialog = Yes() if attended else None
    perms = A.KitToolPermissions(
        workspace=ws, mode=mode,
        confirm_callback=dialog.ask if dialog else None, **kw)
    return perms, dialog


def _spellings(rid: str, room: Path) -> dict:
    py = sys.executable
    answer = room / f"{rid}.answer.json"
    body = json.dumps({"id": rid, "decision": "approve"})
    return {
        "delfin-agent approve": f"delfin-agent approvals approve {rid}",
        "delfin-agent deny": f"delfin-agent approvals deny {rid} --reason ok",
        "delfin-agent answer": f"delfin-agent approvals answer {rid} 1",
        "python -m": f"python -m delfin.agent.cli approvals approve {rid}",
        "interpreter path": f"{py} -m delfin.agent.cli approvals approve {rid}",
        "after cd": f"cd /tmp && delfin-agent approvals approve {rid}",
        "quoted pieces": f"del'fin'-agent \"approvals\" appro\\ve {rid}",
        "nested shell": f"bash -c 'delfin-agent approvals approve {rid}'",
        "xargs": f"echo {rid} | xargs delfin-agent approvals approve",
        "env -u marker": f"env -u {MARKER} delfin-agent approvals approve {rid}",
        "file_confirm.answer": (
            f"python3 -c \"from delfin.agent import file_confirm as f; "
            f"f.answer('{rid}', True)\""),
        "answer_waiting": (
            f"python3 -c \"from delfin.agent.terminal_confirm import "
            f"answer_waiting; answer_waiting('{rid}', True)\""),
        "write_choice_answer": (
            "python3 -c \"from delfin.agent import approval_answers as a; "
            f"a.write_choice_answer({{'id': '{rid}'}}, ['1'])\""),
        "interpreter open": (
            f"python3 -c \"open('{answer}','w').write('{{}}')\""),
        "redirect": f"printf '%s' '{body}' > {answer}",
        "redirect via ~": (
            f"printf x > ~/.delfin/approvals/{rid}.answer.json"),
    }


# -- (a) the shell route, at the gate ----------------------------------------

@pytest.mark.parametrize("mode,attended", MODES)
def test_every_spelling_of_answering_is_refused_at_the_gate(
        home, ws, mode, attended):
    rid, _ = _pending(home)
    perms, dialog = _perms(ws, mode, attended)
    ex = A._DocToolExecutor()
    open_routes = []
    for label, cmd in _spellings(rid, fc.requests_dir()).items():
        err = ex._run_permission_gate("bash", {"command": cmd}, perms)
        if err is None or "blocked:" not in err:
            open_routes.append(f"{label}: {str(err)[:90]!r}")
    assert not open_routes, (
        f"{mode} (attended={attended}) leaves these routes to answering "
        "an approval open:\n  " + "\n  ".join(open_routes))
    assert not (dialog and dialog.shown), (
        "the user was asked to approve the agent answering an approval: "
        f"{dialog.shown[:2]}")


def test_the_refusal_says_to_show_the_request_to_the_user(home, ws):
    rid, _ = _pending(home)
    perms, _ = _perms(ws, "bypassPermissions", False)
    err = A._DocToolExecutor()._run_permission_gate(
        "bash", {"command": f"delfin-agent approvals approve {rid}"}, perms)
    assert err and err.startswith("blocked:"), err
    assert "user" in err and "approvals ls" in err, err


def test_bash_background_takes_the_same_refusal(home, ws):
    rid, _ = _pending(home)
    perms, _ = _perms(ws, "bypassPermissions", False)
    out = json.loads(A._doc_executor.execute(
        "bash_background",
        {"command": f"delfin-agent approvals approve {rid}"}, perms))
    assert "blocked:" in str(out.get("error", "")), out


# -- (a) the shell route, executed, isolation off and on ----------------------

@pytest.mark.parametrize("isolation", ["off", "auto"])
@pytest.mark.parametrize("mode,attended", [
    pytest.param("acceptEdits", True, id="acceptEdits"),
    pytest.param("bypassPermissions", False, id="bypass-unattended"),
])
def test_an_executed_answer_command_writes_no_answer(
        home, ws, monkeypatch, isolation, mode, attended):
    monkeypatch.setattr(A, "_BASH_ISOLATION_OVERRIDE", isolation)
    rid, req = _pending(home)
    perms, _ = _perms(ws, mode, attended)
    answer = req.with_name(f"{rid}.answer.json")
    cmd = f"{sys.executable} -m delfin.agent.cli approvals approve {rid}"
    A._doc_executor.execute("bash", {"command": cmd}, perms)
    assert not answer.exists(), (
        f"isolation={isolation}, {mode}: the agent's shell approved "
        f"{rid}: {answer.read_text()[:200]}")


# -- (b) the write tools and shell writes into the rooms ----------------------

def _room_targets(home: Path, rid: str) -> dict:
    return {
        "approvals answer": home / ".delfin" / "approvals" / f"{rid}.answer.json",
        "terminal answer": (home / ".delfin" / "terminal_confirmations"
                            / f"{rid}.answer.json"),
        "cli_gate state": home / ".delfin" / "cli_gate" / "tok.json",
    }


@pytest.mark.parametrize("isolation", ["off", "auto"])
@pytest.mark.parametrize("granted", ["none", "the .delfin dir"])
@pytest.mark.parametrize("mode,attended", MODES)
def test_no_write_route_reaches_an_answer_room(
        home, ws, monkeypatch, mode, attended, granted, isolation):
    """Driven through ``execute``: the bash write-target gate runs in the
    executor, not in ``_run_permission_gate``."""
    monkeypatch.setattr(A, "_BASH_ISOLATION_OVERRIDE", isolation)
    rid, _ = _pending(home)
    for d in ("approvals", "terminal_confirmations", "cli_gate"):
        (home / ".delfin" / d).mkdir(parents=True, exist_ok=True)
    extra = ((home / ".delfin",) if granted != "none" else ())
    perms, dialog = _perms(ws, mode, attended, extra_workspace_dirs=extra)
    open_routes, asked = [], []
    for label, target in _room_targets(home, rid).items():
        for tool, args, existing in (
                ("write_file", {"path": str(target), "content": "{}"}, False),
                ("write_file", {"path": str(target), "content": "{}"}, True),
                ("edit_file", {"path": str(target), "old_string": "old",
                               "new_string": "{}"}, True),
                ("bash", {"command": f"echo '{{}}' > {target}"}, False),
                ("bash", {"command": f"python3 -c \"open('{target}','w')"
                                     f".write('{{}}')\""}, False)):
            route = f"{tool} -> {label}{' (existing)' if existing else ''}"
            if existing:
                target.write_text("old", encoding="utf-8")
                # Read first, so the read-before-write rule is not what
                # holds the write.
                A._doc_executor.execute("read_file", {"path": str(target)},
                                        perms)
            elif target.exists():
                target.unlink()
            before = len(dialog.shown) if dialog else 0
            A._doc_executor.execute(tool, args, perms)
            if dialog and len(dialog.shown) > before:
                asked.append(route)
            changed = (target.read_text(encoding="utf-8") != "old"
                       if existing else target.exists())
            if changed:
                open_routes.append(route)
    assert not open_routes and not asked, (
        f"{mode} (attended={attended}, granted={granted}, "
        f"isolation={isolation}):\n  written: {open_routes}\n"
        f"  asked the user: {asked}")


@pytest.mark.parametrize("mode,attended", MODES)
def test_the_write_tool_refusal_is_blocked_not_asked(home, ws, mode, attended):
    """Granted or not, the refusal names the room and asks nobody."""
    rid, _ = _pending(home)
    (home / ".delfin" / "approvals").mkdir(parents=True, exist_ok=True)
    perms, _ = _perms(ws, mode, attended,
                      extra_workspace_dirs=(home / ".delfin",))
    target = home / ".delfin" / "approvals" / f"{rid}.answer.json"
    out = A._doc_executor.execute(
        "write_file", {"path": str(target), "content": "{}"}, perms)
    assert not target.exists()
    if mode != "plan":
        assert "blocked:" in str(out), out


# -- (c) the CLI backend's hook --------------------------------------------

@pytest.mark.parametrize("mode", ["default", "acceptEdits", "bypassPermissions"])
def test_the_cli_backend_hook_refuses_answering(home, ws, mode):
    rid, _ = _pending(home)
    cli_gate.write_state("tok", workspace=str(ws), mode=mode, can_ask=False,
                         extra_dirs=[str(home / ".delfin")])
    for cmd in (f"delfin-agent approvals approve {rid}",
                f"python -m delfin.agent.cli approvals deny {rid}"):
        why = cli_gate.decide("tok", {"tool_name": "Bash",
                                      "tool_input": {"command": cmd},
                                      "cwd": str(ws)})
        assert why and "blocked:" in why, (cmd, why)
    state = cli_gate.state_dir() / "tok.json"
    why = cli_gate.decide("tok", {"tool_name": "Write",
                                  "tool_input": {"file_path": str(state),
                                                 "content": "{}"},
                                  "cwd": str(ws)})
    assert why, "the hook let the session rewrite its own gate state"


# -- the second layer: cmd_approvals ----------------------------------------

def _ns(action, rid, **kw):
    return argparse.Namespace(approvals_action=action, request_id=rid,
                              reason=kw.get("reason", ""),
                              choice=kw.get("choice", "1"))


@pytest.mark.parametrize("action", ["approve", "deny", "answer"])
def test_cmd_approvals_refuses_inside_an_agent_session(
        home, monkeypatch, capsys, action):
    from delfin.agent import cli
    rid, req = _pending(home)
    monkeypatch.setenv(MARKER, "1")
    rc = cli.cmd_approvals(_ns(action, rid))
    assert rc != 0
    assert not req.with_name(f"{rid}.answer.json").exists()
    assert "user" in capsys.readouterr().err


def test_the_answer_functions_refuse_inside_an_agent_session(
        home, monkeypatch):
    rid, req = _pending(home)
    monkeypatch.setenv(MARKER, "1")
    assert fc.answer(rid, True) is False
    assert tc.answer_waiting(rid, True) is False
    assert not req.with_name(f"{rid}.answer.json").exists()


def test_every_agent_child_carries_the_marker(home, ws):
    assert A._scrubbed_bash_env().get(MARKER)
    perms, _ = _perms(ws, "bypassPermissions", False)
    out = json.loads(A._doc_executor.execute(
        "bash", {"command": f"printenv {MARKER}"}, perms))
    assert out.get("exit_code") == 0 and out.get("stdout", "").strip(), out


def test_no_room_path_is_widened_inside_the_bwrap_cage(home, ws):
    """Granted ~/.delfin is writable inside the cage, except the rooms."""
    (home / ".delfin" / "approvals").mkdir(parents=True)
    perms, _ = _perms(ws, "bypassPermissions", False,
                      extra_workspace_dirs=(home / ".delfin",))
    argv = A._bash_isolation_argv("true", ws, perms, mode="bwrap")
    room = str((home / ".delfin" / "approvals").resolve())
    pairs = list(zip(argv, argv[1:]))
    assert any(a in ("--ro-bind", "--ro-bind-try") and b == room
               for a, b in pairs), argv


# -- controls ---------------------------------------------------------------

def test_a_human_in_a_plain_shell_still_answers(home, tmp_path):
    rid, req = _pending(home)
    env = {k: v for k, v in os.environ.items() if k != MARKER}
    env["HOME"] = str(home)
    proc = subprocess.run(
        [sys.executable, "-m", "delfin.agent.cli", "approvals", "approve", rid],
        env=env, capture_output=True, text=True, timeout=120, cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr[-800:]
    record = json.loads(req.with_name(f"{rid}.answer.json").read_text())
    assert record["decision"] == "approve"


def test_cmd_approvals_without_the_marker_answers(home, capsys):
    from delfin.agent import cli
    rid, req = _pending(home)
    assert cli.cmd_approvals(_ns("deny", rid, reason="not now")) == 0
    assert json.loads(req.with_name(f"{rid}.answer.json").read_text())[
        "decision"] == "deny"


@pytest.mark.parametrize("mode,attended", MODES)
def test_reading_the_queue_is_not_refused_by_this_guard(home, ws, mode, attended):
    rid, _ = _pending(home)
    perms, _ = _perms(ws, mode, attended)
    ex = A._DocToolExecutor()
    for cmd in ("delfin-agent approvals ls", f"delfin-agent approvals show {rid}",
                "delfin-agent approvals"):
        err = ex._run_permission_gate("bash", {"command": cmd}, perms)
        assert err is None or "blocked:" not in err, (cmd, err)


def test_a_commit_message_naming_the_command_is_not_refused(home, ws):
    perms, _ = _perms(ws, "bypassPermissions", False)
    err = A._DocToolExecutor()._run_permission_gate(
        "bash", {"command": 'git commit -m "refuse delfin-agent approvals '
                            'approve from an agent"'}, perms)
    assert err is None or "approval" not in err, err


def test_a_write_elsewhere_in_the_workspace_still_works(home, ws):
    perms, _ = _perms(ws, "acceptEdits", True)
    out = A._doc_executor.execute(
        "write_file", {"path": "notes/answer.json", "content": "{}"}, perms)
    assert (ws / "notes" / "answer.json").exists(), out
