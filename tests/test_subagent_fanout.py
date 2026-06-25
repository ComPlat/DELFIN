"""Tests for parallel subagent fan-out (_fan_out_subagents).

When the model emits >=2 `subagent` tool-calls in one turn they must run
concurrently; with <2 the sequential path must stay untouched.  We stub
the module-level _doc_executor.execute with a sleeping fake to prove the
calls overlap, and assert the gating + id-mapping contract.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import api_client as ac


def _tc(i: int, name: str = "subagent") -> dict:
    return {
        "id": f"id{i}",
        "function": {
            "name": name,
            "arguments": json.dumps({
                "subagent_type": "explore",
                "description": f"task{i}",
                "prompt": "do x",
            }),
        },
    }


@pytest.fixture
def sleeping_executor(monkeypatch):
    """Replace _doc_executor.execute with a 0.3 s sleeper that echoes
    the description, so wall-clock reveals parallelism."""
    seen: list[str] = []

    def fake_exec(name, args, permissions=None):
        seen.append(args.get("description", ""))
        time.sleep(0.3)
        return json.dumps({"ran": args.get("description")})

    monkeypatch.setattr(ac._doc_executor, "execute", fake_exec)
    return seen


def test_three_subagents_run_in_parallel(sleeping_executor):
    tcs = [_tc(1), _tc(2), _tc(3)]
    t0 = time.time()
    futures, executor = ac._fan_out_subagents(tcs, permissions=None)
    results = {k: f.result() for k, f in futures.items()}
    elapsed = time.time() - t0
    if executor is not None:
        executor.shutdown(wait=True)

    assert len(futures) == 3
    # Sequential would be ~0.9 s; parallel must be well under that.
    assert elapsed < 0.6, f"fan-out not parallel (took {elapsed:.2f}s)"
    # Results map back to the right tool_call id / description.
    assert json.loads(results["id2"])["ran"] == "task2"


def test_single_subagent_keeps_sequential_path(sleeping_executor):
    futures, executor = ac._fan_out_subagents([_tc(1)], permissions=None)
    assert futures == {}
    assert executor is None


def test_only_subagent_calls_are_fanned_out(sleeping_executor):
    tcs = [_tc(1), _tc(2), _tc(9, name="search_docs")]
    futures, executor = ac._fan_out_subagents(tcs, permissions=None)
    try:
        assert set(futures) == {"id1", "id2"}        # search_docs excluded
    finally:
        if executor is not None:
            executor.shutdown(wait=True)


def test_malformed_arguments_default_to_empty(monkeypatch):
    captured: list[dict] = []

    def fake_exec(name, args, permissions=None):
        captured.append(args)
        return "{}"

    monkeypatch.setattr(ac._doc_executor, "execute", fake_exec)
    bad = {"id": "idX", "function": {"name": "subagent", "arguments": "{not json"}}
    futures, executor = ac._fan_out_subagents([bad, _tc(2)], permissions=None)
    for f in futures.values():
        f.result()
    if executor is not None:
        executor.shutdown(wait=True)
    assert {} in captured        # the unparseable args fell back to {}


# ---------------------------------------------------------------------------
# Live-panel registry (Claude-Code-style subagent monitoring)
# ---------------------------------------------------------------------------

def test_running_registry_roundtrip(tmp_path, monkeypatch):
    from delfin.agent import subagents as sa
    monkeypatch.setattr(sa, "_RUNNING_DIR", tmp_path / "running")
    sa._running_update("abc", {"type": "explore", "description": "map hooks",
                                "started_at": 123.0,
                                "actions": ["read_file core.py"],
                                "last_action": "read_file core.py"})
    reg = sa.read_running()
    assert reg["abc"]["type"] == "explore"
    assert reg["abc"]["last_action"] == "read_file core.py"   # live drill-down
    sa._running_update("abc", None)
    assert "abc" not in sa.read_running()


def test_running_registry_is_per_subagent_no_race(tmp_path, monkeypatch):
    """Two parallel subagents each own their file — updating one never drops
    the other (the shared-dict design race-dropped entries)."""
    from delfin.agent import subagents as sa
    monkeypatch.setattr(sa, "_RUNNING_DIR", tmp_path / "running")
    sa._running_update("aaa", {"type": "code", "description": "type.py",
                               "started_at": 1.0, "actions": [], "last_action": ""})
    sa._running_update("bbb", {"type": "code", "description": "range.py",
                               "started_at": 1.0, "actions": [], "last_action": ""})
    # Rapid updates to bbb must not lose aaa.
    for i in range(5):
        sa._running_update("bbb", {"type": "code", "description": "range.py",
                                   "started_at": 1.0,
                                   "actions": [f"step{i}"], "last_action": f"step{i}"})
    reg = sa.read_running()
    assert set(reg) == {"aaa", "bbb"}
    assert reg["bbb"]["last_action"] == "step4"


def test_format_action_summaries():
    from delfin.agent.subagents import _format_action
    assert _format_action("write_file", '{"path": "/a/b/type.py"}') == "write_file type.py"
    assert _format_action("read_file", {"path": "/x/core.py"}) == "read_file core.py"
    assert _format_action("bash", '{"command": "pytest -q"}') == "bash: pytest -q"


def test_background_flag_in_subagent_schema():
    from pathlib import Path
    src = (Path(__file__).resolve().parent.parent
           / "delfin" / "agent" / "api_client.py").read_text(encoding="utf-8")
    i = src.find('"name": "subagent"')
    assert '"background"' in src[i:i + 3000], "background param missing"
    assert "started_in_background" in src


# ---------------------------------------------------------------------------
# Finished-subagent sessions + resume (Claude-Code SendMessage analog)
# ---------------------------------------------------------------------------

def _ev(**kw):
    from types import SimpleNamespace
    base = dict(type="", text="", tool_name="", tool_input="",
                tool_output="", input_tokens=0, output_tokens=0)
    base.update(kw)
    return SimpleNamespace(**base)


def _mk_client(text: str = "ok", tool=None):
    """Fake client. ``tool=(name, input, output)`` emits a tool_use +
    tool_result pair before the final text, so resume-fidelity can be
    exercised."""
    from unittest.mock import MagicMock
    events = []
    if tool:
        name, inp, out = tool
        events.append(_ev(type="tool_use", tool_name=name, tool_input=inp))
        events.append(_ev(type="tool_result", tool_name=name, tool_output=out))
    events.append(_ev(type="text_delta", text=text))
    events.append(_ev(type="message_delta", input_tokens=10, output_tokens=5))
    client = MagicMock()
    client.stream_message = MagicMock(return_value=iter(events))
    client._permissions = None
    client.set_permissions = MagicMock()
    return client


def test_session_persist_and_list(tmp_path, monkeypatch):
    from delfin.agent import subagents as sa
    monkeypatch.setattr(sa, "_SESSIONS_DIR", tmp_path / "sess")
    res = sa.run_subagent(
        subagent_type="explore", description="find answer",
        prompt="please find the answer to everything",
        parent_client=_mk_client("the answer is 42"), parent_perms=None,
    )
    assert res.sa_id
    rec = sa.load_subagent_session(res.sa_id)
    assert rec["subagent_type"] == "explore"
    assert [m["role"] for m in rec["messages"]] == ["user", "assistant"]
    assert "the answer is 42" in rec["messages"][-1]["content"]
    listed = sa.list_finished()
    assert listed and listed[0]["sa_id"] == res.sa_id


def test_resume_replays_prior_turns(tmp_path, monkeypatch):
    from delfin.agent import subagents as sa
    monkeypatch.setattr(sa, "_SESSIONS_DIR", tmp_path / "sess")
    res1 = sa.run_subagent(
        subagent_type="explore", description="find X",
        prompt="where does X live? look around",
        parent_client=_mk_client("initial finding: X lives in y.py"),
        parent_perms=None,
    )
    second = _mk_client("follow-up: X moved to z.py")
    res2 = sa.run_subagent(
        subagent_type="explore", description="",
        prompt="now check whether X moved recently",
        parent_client=second, parent_perms=None,
        resume_from=res1.sa_id,
    )
    assert res2.sa_id == res1.sa_id          # same session accumulates
    assert res2.error == ""
    sent = second.stream_message.call_args.kwargs["messages"]
    contents = " | ".join(str(m.get("content")) for m in sent)
    assert "initial finding: X lives in y.py" in contents   # replayed
    assert "now check whether X moved recently" in contents


def test_resume_preserves_tool_outputs(tmp_path, monkeypatch):
    """The fidelity fix: a resumed subagent sees the actual tool OUTPUTS
    it produced earlier, not just its own text conclusions."""
    from delfin.agent import subagents as sa
    monkeypatch.setattr(sa, "_SESSIONS_DIR", tmp_path / "sess")
    res1 = sa.run_subagent(
        subagent_type="explore", description="grep",
        prompt="search for the audit log writer",
        parent_client=_mk_client(
            "found it",
            tool=("grep_file", '{"q":"audit"}',
                  "src/audit.py:42: def write_audit(record):")),
        parent_perms=None,
    )
    rec = sa.load_subagent_session(res1.sa_id)
    # The tool output was captured with its result, not dropped.
    assert rec["interactions"]
    assert rec["interactions"][0]["name"] == "grep_file"
    assert "src/audit.py:42" in rec["interactions"][0]["output"]

    # On resume, that concrete output is replayed into the model's context.
    second = _mk_client("continuing")
    sa.run_subagent(
        subagent_type="explore", description="",
        prompt="now open that file and summarise write_audit",
        parent_client=second, parent_perms=None, resume_from=res1.sa_id,
    )
    sent = second.stream_message.call_args.kwargs["messages"]
    contents = " | ".join(str(m.get("content")) for m in sent)
    assert "src/audit.py:42: def write_audit" in contents   # output, not just text
    assert "grep_file" in contents


def test_resume_accumulates_interactions_across_rounds(tmp_path, monkeypatch):
    from delfin.agent import subagents as sa
    monkeypatch.setattr(sa, "_SESSIONS_DIR", tmp_path / "sess")
    r1 = sa.run_subagent(
        subagent_type="explore", description="d",
        prompt="first step, look around here",
        parent_client=_mk_client("step1", tool=("read_file", "a.py", "AAA")),
        parent_perms=None,
    )
    sa.run_subagent(
        subagent_type="explore", description="",
        prompt="second step, keep going please",
        parent_client=_mk_client("step2", tool=("read_file", "b.py", "BBB")),
        parent_perms=None, resume_from=r1.sa_id,
    )
    rec = sa.load_subagent_session(r1.sa_id)
    outputs = [it["output"] for it in rec["interactions"]]
    assert "AAA" in outputs and "BBB" in outputs       # both rounds kept
    roles = [m["role"] for m in rec["messages"]]
    assert roles == ["user", "assistant", "user", "assistant"]


def test_resume_unknown_id_is_contained(tmp_path, monkeypatch):
    from delfin.agent import subagents as sa
    monkeypatch.setattr(sa, "_SESSIONS_DIR", tmp_path / "sess")
    res = sa.run_subagent(
        subagent_type="explore", description="d",
        prompt="some long enough prompt here",
        parent_client=_mk_client(), parent_perms=None,
        resume_from="nope1234",
    )
    assert "unknown resume id" in res.error
    assert res.final_text == ""


def test_resume_id_in_subagent_schema():
    from pathlib import Path
    src = (Path(__file__).resolve().parent.parent
           / "delfin" / "agent" / "api_client.py").read_text(encoding="utf-8")
    i = src.find('"name": "subagent"')
    assert '"resume_id"' in src[i:i + 4000], "resume_id param missing"
