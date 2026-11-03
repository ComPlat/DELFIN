"""With the CLI backend, a changed command means what it means everywhere.

The CLI backend asks through its hook: the question goes into the
approvals directory, DELFIN's Relay puts it to the dialog of this
process (terminal or dashboard broker), and the answer goes back as a
file. A change must survive that trip:

- the request carries the command, so the dialog offers to change it;
- the answer carries the changed command, so the hook's refusal tells
  the model which command to run instead, word for word;
- the hook's next question for exactly that command is answered by the
  grant in this process's broker: it runs once without a new dialog.
"""

from __future__ import annotations

import threading
import time

import pytest

from delfin.agent import cli_gate, kit_settings
from delfin.agent import terminal_confirm as tc


@pytest.fixture
def scene(tmp_path, monkeypatch):
    from delfin.agent import audit_log
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH",
                        home / ".delfin" / "settings.json")
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
    monkeypatch.setattr(audit_log, "append", lambda *a, **kw: None)
    ws = tmp_path / "ws"
    ws.mkdir()
    token = cli_gate.new_token()
    cli_gate.write_state(token, workspace=str(ws), mode="default", can_ask=True)
    return token, ws


def _ev(command, cwd):
    return {"tool_name": "Bash", "tool_input": {"command": command}, "cwd": str(cwd)}


def _decide(token, event, broker, answer=None):
    relay = cli_gate.Relay(token, broker.callback, poll_s=0.05).start()
    seen = {}

    def answer_it():
        for _ in range(400):
            req = broker.take()
            if req is not None:
                seen["req"] = req
                try:
                    answer(req)
                except ValueError as exc:   # not changeable: refuse, never hang
                    seen["error"] = str(exc)
                    broker.resolve(req, False)
                return
            time.sleep(0.01)

    t = threading.Thread(target=answer_it, daemon=True) if answer else None
    if t:
        t.start()
    try:
        return cli_gate.decide(token, event), seen
    finally:
        relay.stop()
        if t:
            t.join(5)


def test_a_change_travels_through_the_hook_and_runs_once(scene):
    token, ws = scene
    broker = tc.TerminalConfirmBroker(timeout_s=10.0)
    reason, seen = _decide(token, _ev("some-unlisted-command -x", ws), broker,
                           answer=lambda req: broker.change(req, "other-unlisted-command -y"))
    assert seen["req"].command == "some-unlisted-command -x"
    assert "c" in {o.key for o in tc.options_for(seen["req"])}
    assert reason and "changed the command" in reason
    assert "other-unlisted-command -y" in reason
    # The model runs the changed command: no new dialog, allowed once.
    reason, seen = _decide(token, _ev("other-unlisted-command -y", ws), broker)
    assert reason is None and not seen
    # Once: asked again, and refused here.
    reason, seen = _decide(token, _ev("other-unlisted-command -y", ws), broker,
                           answer=lambda req: broker.resolve(req, False))
    assert seen and reason
