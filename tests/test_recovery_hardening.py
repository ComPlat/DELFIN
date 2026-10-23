"""Recovery hardening: mid-stream retry, degenerate-loop detection,
FAIL-outcome recording, repeat-failure heads-up, and event-driven job
completion wiring.

Previously: a transient error after ANY partial output killed the whole
turn; a model alternating two error shapes (or repeating the identical
successful call) looped to the round cap; errored turns recorded no
outcome so profiles only learned from successes; the repeat-failure
detector was collected-but-never-consumed.
"""

from __future__ import annotations

import json
import textwrap
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import mcp_client as M
from delfin.agent import model_capabilities as mc


# --- Fake OpenAI streaming primitives (mirrors test_loop_hardening) ---------

class _FnDelta:
    def __init__(self, name=None, arguments=None):
        self.name = name
        self.arguments = arguments


class _TCDelta:
    def __init__(self, index, id=None, name=None, arguments=None):
        self.index = index
        self.id = id
        self.function = _FnDelta(name, arguments)


class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _PartialThenRaise:
    """Stream that yields one text chunk, then dies mid-stream."""

    def __init__(self, exc):
        self._exc = exc

    def __iter__(self):
        yield _Chunk([_Choice(_Delta(content="partial answer "))])
        raise self._exc

    def close(self):
        pass


class APITimeoutError(Exception):
    """Name matches the transient classifier's 'timeout' hint."""


def _tool_round(tool_name, args_json, tc_id="c1"):
    return _Stream([
        _Chunk([_Choice(_Delta(tool_calls=[
            _TCDelta(0, id=tc_id, name=tool_name, arguments=args_json)]))]),
        _Chunk([_Choice(_Delta(), finish="tool_calls")], usage=_Usage()),
    ])


def _final(txt="done"):
    return _Stream([
        _Chunk([_Choice(_Delta(content=txt), finish="stop")],
               usage=_Usage()),
    ])


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


@pytest.fixture(autouse=True)
def _isolate_failure_log(monkeypatch, tmp_path):
    """Keep the cross-session failure log out of the user's real home so
    repeat-failure heads-ups are deterministic in tests."""
    from delfin.agent import failure_log as _fl
    monkeypatch.setattr(_fl, "_LOG_PATH", tmp_path / "failure_log.jsonl")


@pytest.fixture
def client(monkeypatch, tmp_path):
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    import delfin.agent.api_client as A
    monkeypatch.setattr(A.time, "sleep", lambda *a, **k: None)
    return A.create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(tmp_path))


def _drive_seq(client, factory):
    """factory(call_index) -> stream or raises; returns (events,)"""
    calls = {"n": 0}

    def _create(**kwargs):
        idx = calls["n"]
        calls["n"] += 1
        return factory(idx)

    client.client.chat.completions.create = _create
    events = list(client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100))
    return events, calls["n"]


# ---------------------------------------------------------------------------
# Mid-stream retry after partial output
# ---------------------------------------------------------------------------

def test_transient_error_after_partial_output_retries_round(client):
    def _factory(idx):
        if idx == 0:
            return _PartialThenRaise(APITimeoutError("upstream reset"))
        return _final("full answer")

    events, n_calls = _drive_seq(client, _factory)
    said = "".join(e.text for e in events
                   if e.type == "text_delta" and e.text)
    told = "".join(e.text for e in events
                   if e.type == "notice" and e.text)
    assert n_calls == 2
    assert "full answer" in said
    # The retry banner is the harness talking about itself. It still
    # reaches the reader, on the channel that is not the model's answer —
    # it used to be scored as one.
    assert "retrying" in told.lower()
    assert "restarts below" in told              # partial-output marker
    assert "retrying" not in said.lower()
    # The turn ended normally, not via an exception.
    assert any(e.type == "message_delta" for e in events)


# ---------------------------------------------------------------------------
# Degenerate-loop detection
# ---------------------------------------------------------------------------

def test_alternating_error_signatures_abort(client):
    def _factory(idx):
        if idx % 2 == 0:
            return _tool_round("zzqqxx_alpha", "{}", tc_id=f"a{idx}")
        return _tool_round("zzqqxx_beta", "{}", tc_id=f"b{idx}")

    events, n_calls = _drive_seq(client, _factory)
    stops = [getattr(e, "stop_reason", "") for e in events
             if e.type == "message_delta"]
    assert "consecutive_identical_errors" in stops
    assert n_calls <= 6                          # far below the round cap


def test_identical_successful_rounds_steer_then_abort(client):
    args = json.dumps({"status": "approve", "evidence": "looks fine"})

    def _factory(idx):
        return _tool_round("report_verdict", args, tc_id=f"r{idx}")

    events, n_calls = _drive_seq(client, _factory)
    stops = [getattr(e, "stop_reason", "") for e in events
             if e.type == "message_delta"]
    assert "no_progress_loop" in stops
    assert n_calls <= 7


# ---------------------------------------------------------------------------
# FAIL outcome recorded on errored turns
# ---------------------------------------------------------------------------

@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    # Named `solo` because every retired mode name now migrates onto
    # it; a fixture built around `quick` describes a manifest the
    # loader no longer has. The multi-role route is kept so the
    # engine's role advancement stays under test.
    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)
    return tmp_path


def test_errored_turn_records_fail_outcome(agent_tree):
    from delfin.agent.engine import AgentEngine

    def _boom(*a, **k):
        raise RuntimeError("backend exploded")
        yield  # pragma: no cover — makes this a generator function

    fake = MagicMock()
    fake.stream_message = MagicMock(side_effect=_boom)
    with patch("delfin.agent.engine.create_client", return_value=fake):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli",
                             mode="quick", pack_dir=agent_tree)
    recorded = []
    engine.record_cycle_outcome = (
        lambda verdict, task, **kw: recorded.append((verdict, kw)))
    with pytest.raises(RuntimeError):
        engine.stream_response("please do the thing")
    assert recorded and recorded[0][0] == "FAIL"
    assert recorded[0][1].get("error_type")


# ---------------------------------------------------------------------------
# Repeat-failure heads-up is consumed, not just collected
# ---------------------------------------------------------------------------

def test_repeat_failure_heads_up_appended(monkeypatch, tmp_path):
    from delfin.agent import api_client as A
    from delfin.agent import failure_log as fl
    monkeypatch.setattr(fl, "_LOG_PATH", tmp_path / "failure_log.jsonl")
    out = ""
    for _ in range(3):
        out = A._doc_executor.execute("zzqqxx", {})
    obj = json.loads(out)                 # payload stays machine-parseable
    assert "error" in obj
    assert "Change approach" in obj.get("heads_up", "")


# ---------------------------------------------------------------------------
# Event-driven job completion wiring
# ---------------------------------------------------------------------------

def test_finished_jobs_block_renders_drained_events(agent_tree, monkeypatch):
    from unittest.mock import MagicMock as MM
    from delfin.agent.engine import AgentEngine
    fake = MM()
    fake.stream_message = MM(side_effect=lambda *a, **k: iter(()))
    with patch("delfin.agent.engine.create_client", return_value=fake):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli",
                             mode="quick", pack_dir=agent_tree)

    class _Perms:
        workspace = agent_tree
    engine.client._permissions = _Perms()

    import delfin.agent.bash_jobs as bj
    monkeypatch.setattr(bj, "drain_finished_events", lambda ws: [{
        "job_id": "bg_1", "command": "orca emitter.inp", "description": "",
        "exit_code": 0, "runtime_s": 7200.0,
        "stdout_tail": "ORCA TERMINATED NORMALLY", "stderr_tail": "",
    }])
    import delfin.agent.job_monitor as jm
    monkeypatch.setattr(jm, "check_agent_jobs", lambda ws: [])
    block = engine._build_finished_jobs_block()
    assert "bg_1" in block
    assert "ok" in block
    assert "ORCA TERMINATED NORMALLY" in block


def test_finished_jobs_block_empty_without_events(agent_tree, monkeypatch):
    from unittest.mock import MagicMock as MM
    from delfin.agent.engine import AgentEngine
    fake = MM()
    fake.stream_message = MM(side_effect=lambda *a, **k: iter(()))
    with patch("delfin.agent.engine.create_client", return_value=fake):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli",
                             mode="quick", pack_dir=agent_tree)

    class _Perms:
        workspace = agent_tree
    engine.client._permissions = _Perms()

    import delfin.agent.bash_jobs as bj
    import delfin.agent.job_monitor as jm
    monkeypatch.setattr(bj, "drain_finished_events", lambda ws: [])
    monkeypatch.setattr(jm, "check_agent_jobs", lambda ws: [])
    assert engine._build_finished_jobs_block() == ""


def test_watch_job_registers_and_requires_job_id(tmp_path):
    from delfin.agent import api_client as A
    from delfin.agent.api_client import KitToolPermissions
    perms = KitToolPermissions(workspace=tmp_path, mode="bypassPermissions")
    out = json.loads(A._doc_executor.execute(
        "watch_job", {"job_id": "12345", "description": "ORCA opt"},
        permissions=perms))
    assert out["status"] == "watching"
    assert out["kind"] == "slurm"
    from delfin.agent.job_monitor import load_watched
    reg = tmp_path / ".delfin" / "agent_watched_jobs.json"
    assert reg.is_file()
    missing = json.loads(A._doc_executor.execute(
        "watch_job", {}, permissions=perms))
    assert "error" in missing
