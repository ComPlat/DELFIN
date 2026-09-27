"""Tests for delfin/agent/session_end.py (package 2, phase 4).

The shared session-end: auto-memory distillation AND skill learning
run from one place, used by BOTH the dashboard and the CLI chat. Each
part fails silently on its own; the setting ``agent.skill_learning.
enabled`` (default on) gates the learning. Tests never write to the
real ~/.delfin: the tool-trace dir is monkeypatched to tmp_path and
the LLM is injected.
"""

import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import delfin.agent.tool_trace as tt  # noqa: E402
from delfin.agent import session_end  # noqa: E402

GOOD_SKILL = """---
name: fix-failing-gate-tests
description: Drive a failing gate test red-green with a broken control
---

# Fix failing gate tests

1. Reproduce the failure first.
"""


def _bash_msg(cmd: str, out: str = "", exit_code: int = 0) -> dict:
    res = f'{{"exit_code": {exit_code}, "stdout": "{out}"}}'
    return {"role": "tool", "content":
            f'<span class="tool-name">$</span> {cmd}'
            f'<details><summary> &rarr; {res}</details>'}


def _runner_msg(target: str, passed: int = 7) -> dict:
    return {"role": "tool", "content":
            '<span class="tool-name">TestRunner</span>  '
            f'<span class="tool-param">{target}</span>'
            '<details><summary> &rarr; '
            f'{{"summary": {{"passed": {passed}, "failed": 0, "errors": 0}}}}'
            '</details>'}


QUALIFYING_SESSION = [
    {"role": "user", "content": "Fix the failing test in module X."},
    _bash_msg("grep -n KEYWORD delfin/x.py"),
    _bash_msg("sed -n 100,140p delfin/x.py"),
    _bash_msg("ruff check delfin/x.py"),
    _bash_msg("git diff --stat"),
    _runner_msg("tests/test_x.py"),
]


class _StubPropose:
    def __init__(self):
        self.calls = []

    def __call__(self, name, text, *, evidence, source, base_version=""):
        self.calls.append(dict(name=name, evidence=evidence, source=source))
        return ("proposal", name)


class TestLearnAtSessionEnd:
    def test_qualifying_session_with_trace_evidence_yields_one_proposal(
            self, monkeypatch, tmp_path):
        monkeypatch.setattr(tt, "_DIR", tmp_path / "traces")
        sid = "sess-abc"
        for i in range(4):
            tt.record(sid, tool="Bash", tool_input="grep x", output="ok")
        tt.record(sid, tool="TestRunner", tool_input="tests/test_x.py",
                  output='{"summary": {"passed": 7, "failed": 0}}')
        stub = _StubPropose()
        out = session_end.learn_at_session_end(
            QUALIFYING_SESSION, session_id=sid,
            llm=lambda *a: GOOD_SKILL, _propose=stub)
        assert out is not None
        assert len(stub.calls) == 1
        assert stub.calls[0]["evidence"][0].kind == "test"

    def test_no_evidence_in_messages_nor_trace_no_proposal(
            self, monkeypatch, tmp_path):
        monkeypatch.setattr(tt, "_DIR", tmp_path / "traces")
        msgs = [m for m in QUALIFYING_SESSION if "TestRunner" not in
                m.get("content", "")]
        stub = _StubPropose()
        out = session_end.learn_at_session_end(
            msgs, session_id="sess-none", llm=lambda *a: GOOD_SKILL,
            _propose=stub)
        assert out is None
        assert stub.calls == []

    def test_setting_off_disables_learning(self, monkeypatch, tmp_path):
        monkeypatch.setattr(tt, "_DIR", tmp_path / "traces")
        stub = _StubPropose()
        out = session_end.learn_at_session_end(
            QUALIFYING_SESSION, session_id="sess-off",
            settings={"agent": {"skill_learning": {"enabled": False}}},
            llm=lambda *a: GOOD_SKILL, _propose=stub)
        assert out is None
        assert stub.calls == []

    def test_learning_failure_never_raises(self, monkeypatch, tmp_path):
        monkeypatch.setattr(tt, "_DIR", tmp_path / "traces")

        def boom(*a, **k):
            raise RuntimeError("learning exploded")

        out = session_end.learn_at_session_end(
            QUALIFYING_SESSION, session_id="sess-boom", llm=boom)
        assert out is None

    def test_trace_supplies_evidence_for_a_cli_shaped_session(
            self, monkeypatch, tmp_path):
        """The CLI's saved chat messages carry no tool results at all
        (cli._display_messages keeps user/assistant only). The session's
        tool trace must supply BOTH the qualification basis and the
        evidence — this is the whole reason _synthetic_tool_messages
        exists."""
        monkeypatch.setattr(tt, "_DIR", tmp_path / "traces")
        sid = "sess-cli-shape"
        for i in range(5):
            tt.record(sid, tool="Bash", tool_input="grep x",
                      output="match")
        tt.record(sid, tool="TestRunner", tool_input="tests/test_x.py",
                  output='{"summary": {"passed": 4, "failed": 0}}')
        cli_msgs = [  # what the CLI actually saves: plain text only
            {"role": "user", "content": "Fix the failing test in module X."},
            {"role": "assistant", "content": "Done, tests green now."},
        ]
        stub = _StubPropose()
        out = session_end.learn_at_session_end(
            cli_msgs, session_id=sid,
            llm=lambda *a: GOOD_SKILL, _propose=stub)
        assert out is not None, \
            "the tool trace must qualify + evidence a CLI session"
        assert len(stub.calls) == 1
        assert stub.calls[0]["evidence"][0].kind == "test"


class TestSessionEndSettings:
    def test_default_is_enabled(self):
        cfg = session_end.skill_learning_settings(None)
        assert cfg["enabled"] is True

    def test_explicit_off_is_respected(self):
        cfg = session_end.skill_learning_settings(
            {"agent": {"skill_learning": {"enabled": False}}})
        assert cfg["enabled"] is False


# --- Integration: the PUBLIC session-end paths call the shared hook ----

class _StubEngine:
    session_id = "sess-cli-1"
    token_usage = {"input": 1, "output": 1}

    def __init__(self):
        self.messages = [{"role": "user", "content": "do it"}]

    def stream_response(self, user_message="", **kwargs):
        self.messages.append({"role": "assistant", "content": "ANSWER"})
        return "ANSWER"

    def export_state(self):
        return {"engine_state": {}}


_RESULT = {"text": "ANSWER", "tool_calls": [], "input_tokens": 1,
           "output_tokens": 1, "error": "", "session_id": "sess-cli-1"}


class TestPublicSessionEndHooks:
    def test_cmd_run_route_calls_learn_at_session_end(
            self, monkeypatch, tmp_path, capsys):
        """The public CLI path (``delfin-agent -p``) runs the shared
        session-end hook after the answer — the integration the module
        tests cannot prove."""
        import argparse
        from delfin.agent import cli as agent_cli

        monkeypatch.setattr(Path, "home",
                            classmethod(lambda cls: tmp_path))
        monkeypatch.setattr(agent_cli, "_build_engine",
                            lambda args: _StubEngine())
        monkeypatch.setattr(agent_cli, "_run_once",
                            lambda engine, prompt, **kw: dict(_RESULT))
        monkeypatch.setattr(agent_cli, "_save_session",
                            lambda *a, **k: "sess-cli-1")
        monkeypatch.setattr("sys.stdin", io.TextIOWrapper(
            io.BytesIO(b"")))  # not a tty -> no extra read

        called: list[dict] = []

        def fake_learn(msgs, *, session_id="", **kw):
            called.append({"n_msgs": len(msgs), "sid": session_id})
            return None

        monkeypatch.setattr(
            "delfin.agent.session_end.learn_at_session_end", fake_learn)
        rc = agent_cli.main(["-p", "do it"])
        assert rc == 0
        assert called == [{"n_msgs": 1, "sid": "sess-cli-1"}]

    def test_dashboard_distill_thread_calls_learn_at_session_end(
            self, monkeypatch, tmp_path):
        """The dashboard's _maybe_distill_session runs the shared hook
        alongside auto-memory — via source inspection of the built
        closure, because building the whole dashboard UI is out of
        scope here. The call must be inside the thread's try/except so
        a learning failure cannot break the session save."""
        import inspect
        import delfin.dashboard.tab_agent as tab

        src = inspect.getsource(tab)
        assert "from delfin.agent.session_end import learn_at_session_end" \
            in src, "dashboard must use the shared session_end hook"
        # The hook call sits after the memory distill, both inside _run.
        distill_at = src.index("distill_and_save(msgs, repo_root=")
        hook_at = src.index("learn_at_session_end(msgs, session_id=")
        assert 0 < distill_at < hook_at
