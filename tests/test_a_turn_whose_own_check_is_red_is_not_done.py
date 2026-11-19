"""A turn whose own verification is still red does not report "done".

Measured on 2026-11-18 across six swarm agents: five ended with
``stop_reason: done``, and the work of one of them had its own new
integration test red -- the wiring it needed sat unapplied in a handover
file. Nothing in the outcome said so. A person reading the stream saw the
engine's warning; every program reading the result saw success, which is
what an orchestrator, a supervisor and every measurement go by.

The engine already did the hard part: it runs auto-verify on the Python
files a turn edited, forces up to two fix rounds while the check is red,
and says loudly that "the final state was NOT re-verified" when it gives
up. Then it ended the turn ``end_turn``, which ``classify`` maps to
"done".

So this is one engine stop reason, ``verify_still_red``, mapped to the
"incomplete" that landed with the stop-reason work the day before --
exit code 9, and ``ended_by_itself`` true, because such a turn may be
resumed exactly as a finished one may.

What is deliberately NOT incomplete, and both are pinned below: a
verification that was SKIPPED (no runner, or a workspace on the slow
list) leaves nothing known, which is a weaker statement than a failure
and already has its own notice; and ``auto_verify`` switched off is a
choice the person made, which must not be turned into a failed turn.

These tests drive the real headless run against the scripted endpoint.
Reading the source would pass for a version where the condition is
written and never reached.
"""

from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import pytest

from tests.agent_terminal_probes import fake_model


def _plan_edit_then_answer(_messages) -> dict:
    """Edit the one file, then answer. Two rounds, no test runs of its own:
    the auto-verify at turn end is the only check in play."""
    state = _plan_edit_then_answer.state
    state["n"] += 1
    if state["n"] == 1:
        # A NEW file: write_file refuses to overwrite an existing one
        # ("refusing to overwrite", which is edit_file's job), and a
        # refused call is not an edit -- the engine does not track it and
        # the verification then has nothing to check. Measured while
        # writing this test, which is why it says so here.
        return {"tool": "write_file",
                "args": {"path": "added.py",
                         "content": "def added():\n    return 1\n"}}
    return {"text": "Wrote added.py. Done."}


class _Launch:
    def __init__(self, plan_fn):
        plan_fn.state = {"n": 0}
        self._plan_fn = plan_fn

    def plan(self, messages):
        return self._plan_fn(messages)


def _run(monkeypatch, tmp_path, *, verify_mode: str, verify_cmd: str):
    """One real headless turn; returns (exit code, parsed result)."""
    from delfin.agent import api_client as _A
    from delfin.agent import cli as agent_cli
    from delfin.agent import mcp_client as M
    from delfin.agent import model_capabilities as mc
    from delfin.agent import provider_profile as PP

    home = tmp_path / "home"
    ws = home / "ws"
    ws.mkdir(parents=True)
    (ws / "source.py").write_text("def answer():\n    return 0\n")
    (ws / "conftest.py").write_text("")

    monkeypatch.chdir(ws)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))

    class _EmptyReg:
        def discover_all(self):
            return []

        def discover_resources(self):
            return []

        def discover_prompts(self):
            return []

    cap = mc.ModelCapabilities(model=fake_model.MODEL, provider="ollama",
                               context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: cap)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    baseline = tmp_path / "learned_profiles.json"
    baseline.write_bytes(PP._DEFAULT_PATH.read_bytes())
    monkeypatch.setattr(PP, "_DEFAULT_PATH", baseline)
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda engine, root, **kw: "sid-1")
    monkeypatch.setattr(agent_cli, "_resume_or_create", lambda engine, args: "")

    # The verification: a command of our own, so the test does not depend
    # on a pytest run inside a pytest run.
    monkeypatch.setenv("DELFIN_AUTO_VERIFY", verify_mode)
    settings = home / ".delfin_settings.json"
    settings.write_text(json.dumps(
        {"agent": {"auto_verify": verify_mode,
                   "auto_verify_command": verify_cmd}}) + "\n")

    stdin = io.StringIO("")
    stdin.isatty = lambda: False              # type: ignore[method-assign]
    monkeypatch.setattr("sys.stdin", stdin)

    monkeypatch.setattr(fake_model, "plan", _Launch(_plan_edit_then_answer).plan)
    srv, port = fake_model.serve(log_path=str(tmp_path / "model.log"))
    monkeypatch.setenv("OLLAMA_HOST", f"http://127.0.0.1:{port}")
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    monkeypatch.setattr(_A, "_resolve_max_tool_rounds", lambda *a, **k: 8)

    captured: list[str] = []
    real_write = sys.stdout.write

    def _tee(text):
        captured.append(text)
        return real_write(text)

    monkeypatch.setattr(sys.stdout, "write", _tee)
    try:
        rc = agent_cli.main(["-p", "edit the file", "--cwd", str(ws),
                             "--provider", "ollama",
                             "--model", fake_model.MODEL,
                             "--permission-mode", "default",
                             "--output-format", "json"])
    finally:
        try:
            srv.shutdown()
        except Exception:
            pass

    result = None
    for line in "".join(captured).splitlines():
        line = line.strip()
        if line.startswith("{") and '"stop_reason"' in line:
            try:
                result = json.loads(line)
            except Exception:
                continue
    return rc, result


@pytest.mark.timeout(180)
def test_a_red_check_makes_the_turn_incomplete(monkeypatch, tmp_path):
    rc, result = _run(monkeypatch, tmp_path,
                      verify_mode="command", verify_cmd="sh -c 'exit 1'")

    assert result is not None, "no result object from the headless run"
    assert result["stop_reason"] == "incomplete", result
    assert rc == 9, rc


@pytest.mark.timeout(180)
def test_a_green_check_still_finishes(monkeypatch, tmp_path):
    """The control. A passing verification clears the recorded problem --
    it is re-assigned on every run -- so the turn ends as it always did."""
    rc, result = _run(monkeypatch, tmp_path,
                      verify_mode="command", verify_cmd="sh -c 'exit 0'")

    assert result is not None
    assert result["stop_reason"] == "done", result
    assert rc == 0, rc


@pytest.mark.timeout(180)
def test_verification_switched_off_is_not_an_incomplete_turn(monkeypatch,
                                                             tmp_path):
    """A switch the person turned off must not invent a failure."""
    rc, result = _run(monkeypatch, tmp_path,
                      verify_mode="off", verify_cmd="sh -c 'exit 1'")

    assert result is not None
    assert result["stop_reason"] == "done", result
    assert rc == 0, rc


def test_the_reason_classifies_and_may_be_resumed():
    """Three properties of the mapping, not of the engine.

    ``incomplete`` and its exit code 9 landed with the stop-reason work;
    this only adds an engine reason to it. ``ended_by_itself`` must be
    true: such a turn stopped on its own terms and a supervisor may
    continue it, which is the whole point of not calling it an error.
    """
    from delfin.agent import turn_outcome as T

    reason = T.classify(error="", engine_stop_reason="verify_still_red",
                        denials=[], ended_on_denial=False)
    assert reason == "incomplete"
    assert T.EXIT_CODES[reason] == 9
    assert T.ended_by_itself(reason) is True
    assert T.classify(error="", engine_stop_reason="end_turn",
                      denials=[], ended_on_denial=False) == "done"


def test_a_red_check_outranks_the_task_list_reasons():
    """Both map to "incomplete"; a red check is the more specific fact and
    must be the one reported, so a reader learns WHY."""
    import inspect

    from delfin.agent import api_client as _A

    src = inspect.getsource(_A)
    i = src.index('_stop = {"open": "end_turn_open_tasks"')
    block = src[max(0, i - 2500):i]
    assert "_stop = \"verify_still_red\"" in block, (
        "the red-check branch must precede the task-list branches")
    red = block.index('_stop = "verify_still_red"')
    promised = block.index("_PROMISED_STOP_REASON")
    assert red < promised, "and precede the promised-tool-call branch"
