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


def _stub_client(script, workspace_root=None):
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
                    msgs = kw.get("messages") or []
                    state["requests"].append(msgs)
                    # First request of a turn: system+user only (no
                    # assistant yet) -> restart the script so every
                    # turn plays the full tool-round sequence.
                    if not any(isinstance(m, dict)
                               and m.get("role") == "assistant"
                               for m in msgs):
                        state["round"] = 0
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
    # Real permissions over the temp workspace so read_file actually
    # reads (with _permissions=None every tool round returns a 199-char
    # error string and the work thresholds are never met).
    c._permissions = (
        ac.KitToolPermissions(workspace=workspace_root, mode="default")
        if workspace_root is not None else None)
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
    # One DISTINCT file per round: the loop's no-progress detector
    # (identical round signatures) ends a turn that re-reads the same
    # file over and over, before the nudge hook has enough rounds.
    for i in range(24):
        (tmp_path / "big{}.txt".format(i)).write_text(
            "finding {} ".format(i) * 600)
    big = tmp_path / "big0.txt"
    return tmp_path, big


def _round_files(ws, n):
    return [ws / "big{}.txt".format(i) for i in range(n)]


def _run_turn(client, message="go"):
    out = []
    for ev in client.stream_message(
            messages=[{"role": "user", "content": message}],
            system="You are a test.", max_tokens=64):
        out.append(ev)
    return out


def _count_new_nudges(requests):
    """Nudges that APPEARED between two consecutive requests.

    A nudge stays in api_messages for the rest of the turn, so every
    later request carries it again — counting message occurrences
    would count one nudge 20 times. What the hook actually does is
    append it once; the visible event is a request containing a nudge
    its predecessor did not.
    """
    new = 0
    prev = set()
    for req in requests:
        cur = {
            id(m) for m in (req or [])
            if isinstance(m, dict) and m.get("role") == "user"
            and "worth keeping" in str(m.get("content", ""))}
        new += len(cur - prev)
        prev |= cur
    return new


def test_long_turn_fires_exactly_one_nudge(workspace):
    ws, big = workspace
    # 20 rounds, each reading a DISTINCT file (real executor, real
    # result sizes), then a final answer.
    files = _round_files(ws, 20)
    script = [[("read_file", _read_args(f))] for f in files]
    script.append("done")
    client, state = _stub_client(script, workspace_root=ws)
    # The wiring (coordinator's api_client commit) attaches the state;
    # for the integration path the client carries it after the hook.
    # Before the hook lands this attribute is absent -> the nudge
    # cannot fire -> this test is RED against the unwired tree.
    _run_turn(client)
    assert _count_new_nudges(state["requests"]) == 1, (
        "a 20-round tool turn must fire exactly one keep-what-you-learned "
        "nudge; got {}".format(
            _count_new_nudges(state["requests"])))


def test_short_turn_fires_no_nudge(workspace):
    ws, big = workspace
    script = [[("read_file", _read_args(big))], "done"]
    client, state = _stub_client(script, workspace_root=ws)
    _run_turn(client)
    nudges = [
        m for req in state["requests"] for m in (req or [])
        if isinstance(m, dict) and m.get("role") == "user"
        and "worth keeping" in str(m.get("content", ""))]
    assert nudges == []


def test_two_consecutive_long_turns_one_nudge_each(workspace):
    """The cap resets at the ENTRY of the public turn method, so a
    second user turn with fresh work nudges again — once each, never
    twice in one turn."""
    ws, big = workspace
    files = _round_files(ws, 20)
    script = [[("read_file", _read_args(f))] for f in files]
    script.append("done")
    client, state = _stub_client(script, workspace_root=ws)
    _run_turn(client)
    _run_turn(client)
    # Nudges appear in the requests of BOTH turns, exactly once each.
    assert _count_new_nudges(state["requests"]) == 2, (
        "two consecutive long turns must fire exactly one nudge each; "
        "got {}".format(_count_new_nudges(state["requests"])))


def test_the_nudge_says_it_is_an_automatic_note(workspace):
    """The injected user-role text must not read as the human speaking."""
    ws, big = workspace
    files = _round_files(ws, 20)
    script = [[("read_file", _read_args(f))] for f in files]
    script.append("done")
    client, state = _stub_client(script, workspace_root=ws)
    _run_turn(client)
    nudges = [
        str(m.get("content", "")) for req in state["requests"]
        for m in (req or [])
        if isinstance(m, dict) and m.get("role") == "user"
        and "worth keeping" in str(m.get("content", ""))]
    assert nudges
    for text in nudges:
        assert text.startswith("[DELFIN note]"), (
            "the nudge must announce itself as a system note, not the "
            "user: {}".format(text[:80]))


def test_setting_disables_the_nudge_on_the_turn_path(
        workspace, monkeypatch):
    ws, big = workspace
    files = _round_files(ws, 20)
    script = [[("read_file", _read_args(f))] for f in files]
    script.append("done")
    client, state = _stub_client(script, workspace_root=ws)
    from delfin.agent import memory_nudge as mn
    monkeypatch.setattr(
        mn, "_enabled",
        staticmethod(lambda settings=None: False))
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
