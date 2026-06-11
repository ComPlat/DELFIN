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
    monkeypatch.setattr(sa, "_RUNNING_PATH", tmp_path / "running.json")
    sa._running_update("abc", {"type": "explore", "description": "map hooks",
                                "started_at": 123.0})
    reg = sa.read_running()
    assert reg["abc"]["type"] == "explore"
    sa._running_update("abc", None)
    assert "abc" not in sa.read_running()


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

def _mk_client(text: str = "ok"):
    from types import SimpleNamespace
    from unittest.mock import MagicMock
    events = [
        SimpleNamespace(type="text_delta", text=text, tool_name="",
                        tool_input="", input_tokens=0, output_tokens=0),
        SimpleNamespace(type="message_delta", text="", tool_name="",
                        tool_input="", input_tokens=10, output_tokens=5),
    ]
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
    rec = sa.load_subagent_session(res1.sa_id)
    assert [m["role"] for m in rec["messages"]] == [
        "user", "assistant", "user", "assistant",
    ]


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
