"""The nudge and the cap, driven through the public turn path.

Package 7, phase 4. The unit tests pin maybe_nudge and the hard cap in
isolation; these tests drive them the way DELFIN does in production:

- a LONG turn — a stream_message call whose model keeps calling tools
  for round after round — accumulates enough work (tool calls AND
  result text) to fire exactly ONE nudge, and the nudge text lands in
  the request the model sees next;
- a SHORT turn (one tool round, small results) fires none;
- the injected memory block stays under the hard cap with a 500-fact
  store behind it (the phase-1 measurement store shape).

The client is the real OpenAIClient.stream_message tool loop with a
stub backend (chat.completions.create), the same pattern the
empty-stream tests use; the tool calls execute real read_file calls
against a temp workspace via the doc executor, so the round accounting
counts real work. No network, no real ~/.delfin.
"""

from __future__ import annotations

import threading
import types
from pathlib import Path

import pytest


def _stub_client(script):
    """An OpenAIClient whose backend follows ``script``.

    ``script`` is a list of rounds; each round is either a text answer
    (str, ends the turn) or a list of (tool_name, args_json) pairs the
    model "calls" in that round. The stub records the MESSAGES list of
    every request so tests can see what the model was shown.
    """
    from delfin.agent import api_client as ac

    state = {"round": 0, "requests": []}

    class _Stream:
        def __init__(self, events):
            self._events = events

        def __iter__(self):
            for ev in self._events:
                yield ev

        def close(self):
            pass

    def _events_for(round_idx):
        step = script[round_idx]
        if isinstance(step, str):
            delta = types.SimpleNamespace(
                content=step, tool_calls=None, reasoning_content=None)
            finish = "stop"
            events = [types.SimpleNamespace(
                usage=None, choices=[types.SimpleNamespace(
                    delta=delta, finish_reason=None)]),
                types.SimpleNamespace(
                    usage=types.SimpleNamespace(prompt_tokens=5,
                                                completion_tokens=3),
                    choices=[types.SimpleNamespace(
                        delta=types.SimpleNamespace(
                            content=None, tool_calls=None,
                            reasoning_content=None),
                        finish_reason=finish)])]
            return events
        # A tool-call round: one tool call per scripted pair.
        tcs = []
        for i, (name, args) in enumerate(step):
            tcs.append(types.SimpleNamespace(
                index=i, id="c{}".format(i),
                function=types.SimpleNamespace(name=name, arguments=args)))
        delta = types.SimpleNamespace(
            content=None, tool_calls=tcs, reasoning_content=None)
        return [types.SimpleNamespace(
            usage=None, choices=[types.SimpleNamespace(
                delta=delta, finish_reason=None)]),
            types.SimpleNamespace(
                usage=types.SimpleNamespace(prompt_tokens=5,
                                            completion_tokens=3),
                choices=[types.SimpleNamespace(
                    delta=types.SimpleNamespace(
                        content=None, tool_calls=None,
                        reasoning_content=None),
                    finish_reason="tool_calls")])]

    class _Stub:
        class chat:
            class completions:
                @staticmethod
                def create(**kw):
                    state["requests"].append(kw.get("messages"))
                    idx = min(state["round"], len(script) - 1)
                    evs = _events_for(idx)
                    state["round"] += 1
                    return _Stream(evs)

    c = ac.OpenAIClient.__new__(ac.OpenAIClient)
    c.client = _Stub()
    c.model = "kit.glm-5.3"
    c._provider = "kit"
    c._base_url = "https://ki-toolbox.scc.kit.edu/api/v1"
    c._api_key = "x"
    c.effort = ""
    c._permissions = None
    c.on_model_switched = None
    c._steer_lock = threading.Lock()
    c._steer_msgs = []
    c._run_notes = []
    c._test_evidence = []
    c._observed_files_session = set()
    c._red_test_files = []
    c._thrash = None
    c.max_tool_rounds = 60
    return c, state


def _read_args(path):
    return json.dumps({"path": str(path)})


import json  # noqa: E402  (used by _read_args)


@pytest.fixture
def workspace(tmp_path):
    (tmp_path / "pack" / "shared").mkdir(parents=True)
    (tmp_path / "pack" / "agents").mkdir()
    big = tmp_path / "big.txt"
    big.write_text("finding " * 600)
    return tmp_path, big


def _run_turn(client, message="go"):
    out = []
    for ev in client.stream_message(
            messages=[{"role": "user", "content": message}],
            system="You are a test.", max_tokens=64):
        out.append(ev)
    return out


@pytest.mark.xfail(
    reason="nudge hook in api_client is wired by the coordinator's "
    "bundle commit (package 7 proposal); until then a long turn "
    "cannot fire a nudge. Flips to XPASS once wired — then remove "
    "this marker.", strict=False)
def test_long_turn_fires_exactly_one_nudge(workspace):
    _, big = workspace
    # 20 rounds of a read_file whose result is large (real executor,
    # real result sizes), then a final answer.
    script = [[("read_file", _read_args(big))] for _ in range(20)]
    script.append("done")
    client, state = _stub_client(script)
    # The wiring (coordinator's api_client commit) attaches the state;
    # for the integration path the client carries it after the hook.
    # Before the hook lands this attribute is absent -> the nudge
    # cannot fire -> this test is RED against the unwired tree.
    _run_turn(client)
    nudges = [
        m.get("content", "") for req in state["requests"]
        for m in (req or [])
        if isinstance(m, dict) and m.get("role") == "user"
        and "worth keeping" in str(m.get("content", ""))
    ]
    assert len(nudges) == 1, (
        "a 20-round tool turn must fire exactly one keep-what-you-learned "
        "nudge; got {}".format(len(nudges)))


def test_short_turn_fires_no_nudge(workspace):
    _, big = workspace
    script = [[("read_file", _read_args(big))], "done"]
    client, state = _stub_client(script)
    _run_turn(client)
    nudges = [
        m for req in state["requests"] for m in (req or [])
        if isinstance(m, dict) and m.get("role") == "user"
        and "worth keeping" in str(m.get("content", ""))]
    assert nudges == []


def test_memory_block_stays_under_cap_with_500_facts(
        workspace, monkeypatch):
    """Phase-1 store shape (500 facts), public recall path, hard cap.

    The end-to-end pairing of package 7: the same session that fires
    nudges must never see a memory block past the hard cap — the two
    halves (curating input, bounding injection) only work together.
    """
    ws, _ = workspace
    from delfin.agent import memory_store as ms
    from delfin.agent import prompt_loader as pl
    monkeypatch.setattr(Path, "home", lambda: ws)
    for i in range(500):
        ms.save_typed_memory(
            "feedback: when drawing molecule {} use wedge bonds for "
            "the stereo centre of ligand {}".format(i, i),
            repo_root=ws)
    out = pl.PromptLoader(ws)._load_external_memory_context(
        task_text="draw the molecule")
    assert out
    assert len(out) <= pl.MEMORY_CONTEXT_HARD_CAP, (
        "block of {} chars exceeds the hard cap of {}".format(
            len(out), pl.MEMORY_CONTEXT_HARD_CAP))
