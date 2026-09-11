"""A turn's tool rounds stay in the history for the turns after it.

Two field cases on 2026-09-11 share one root: the client kept its
rounds in a request-local list and the engine's history received only
the final text, so every correction turn and every later question was
answered by a model that could not see what it had read or run. GLM,
told two cited paths did not exist, wrote "ich habe die Dateien nicht
gelesen" after seven tool calls; DeepSeek, asked to say its answer in
German, added "eigentlich hätte ich die Dateien zunächst lesen sollen"
after reading them.
"""
import inspect
import textwrap
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import turn_history as th
from delfin.agent.api_client import StreamEvent


def _call(cid, name="read_file", args='{"path": "a.py"}'):
    return {"id": cid, "type": "function",
            "function": {"name": name, "arguments": args}}


def _round(cid, result="1  x = 1\n", name="read_file"):
    return [{"role": "assistant", "content": "", "tool_calls": [_call(cid, name)]},
            {"role": "tool", "tool_call_id": cid, "content": result}]


# --- compact_turn_rows ----------------------------------------------------

def test_a_round_is_kept_with_its_pairing():
    rows = th.compact_turn_rows(_round("c1"))
    assert [r["role"] for r in rows] == ["assistant", "tool"]
    assert rows[0]["tool_calls"][0]["id"] == "c1"
    assert rows[1]["tool_call_id"] == "c1"


def test_a_long_result_is_cut_to_a_head():
    rows = th.compact_turn_rows(_round("c1", result="x" * 5000))
    assert len(rows[1]["content"]) < 1400
    assert "not kept in the history" in rows[1]["content"]


def test_long_arguments_stay_json():
    import json
    rows = th.compact_turn_rows([
        {"role": "assistant", "content": None,
         "tool_calls": [_call("c1", "write_file", json.dumps({"content": "y" * 5000}))]},
        {"role": "tool", "tool_call_id": "c1", "content": "ok"}])
    args = json.loads(rows[0]["tool_calls"][0]["function"]["arguments"])
    assert "_not_kept_in_history" in args
    assert rows[0]["content"] == ""


def test_an_unanswered_call_is_dropped():
    rows = th.compact_turn_rows(
        [{"role": "assistant", "content": "", "tool_calls": [_call("c1")]}]
        + _round("c2"))
    assert [r.get("tool_call_id") for r in rows if r["role"] == "tool"] == ["c2"]
    assert all(c["id"] == "c2" for r in rows if r["role"] == "assistant"
               for c in r["tool_calls"])


def test_an_answer_to_nothing_is_dropped():
    rows = th.compact_turn_rows(
        [{"role": "tool", "tool_call_id": "ghost", "content": "?"}] + _round("c1"))
    assert len(rows) == 2


def test_steers_and_interim_text_travel_too():
    rows = th.compact_turn_rows(
        _round("c1") + [{"role": "user", "content": "nur mod1 bitte"},
                        {"role": "assistant", "content": "Verstanden."}]
        + _round("c2"))
    assert [r["role"] for r in rows] == [
        "assistant", "tool", "user", "assistant", "assistant", "tool"]


def test_a_long_turn_keeps_its_last_rounds_and_says_so():
    rows = []
    for i in range(20):
        rows += _round(f"c{i}", name="grep_file" if i % 2 else "read_file")
    kept = th.compact_turn_rows(rows, keep_rounds=4)
    ids = [r["tool_call_id"] for r in kept if r["role"] == "tool"]
    assert ids == ["c16", "c17", "c18", "c19"]
    first = next(r for r in kept if r["role"] == "assistant")
    assert "16 earlier tool round(s)" in first["content"]
    assert "read_file ×8" in first["content"] and "grep_file ×8" in first["content"]


def test_system_rows_never_come_along():
    rows = th.compact_turn_rows(
        [{"role": "system", "content": "prompt"}] + _round("c1"))
    assert rows[0]["role"] == "assistant"


# --- repair_tool_pairing ---------------------------------------------------

def test_a_cut_through_a_round_is_repaired():
    history = [{"role": "user", "content": "q"},
               {"role": "tool", "tool_call_id": "old", "content": "orphan"},
               {"role": "assistant", "content": "", "tool_calls": [_call("c1")]},
               {"role": "assistant", "content": "answer"}]
    out = th.repair_tool_pairing(history)
    assert [m["role"] for m in out] == ["user", "assistant"]
    assert out[0] is history[0] and out[1] is history[3]


def test_a_whole_history_passes_untouched():
    history = [{"role": "user", "content": "q"}] + _round("c1") + [
        {"role": "assistant", "content": "answer"}]
    assert th.repair_tool_pairing(history) == history


# --- through the engine ----------------------------------------------------

@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    (lite_dir / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


def _client(replies, rounds, keeps=True):
    """A client that, like OpenAIClient, leaves its request rows behind."""
    fake = MagicMock()
    fake._observed_files_session = set()
    if keeps:
        fake.KEEPS_TOOL_HISTORY = True
    else:
        del fake.KEEPS_TOOL_HISTORY
    calls = {"n": 0}
    seen = []

    def _stream(*a, **k):
        i = calls["n"]
        calls["n"] += 1
        msgs = list(k.get("messages") or (a[0] if a else []))
        seen.append(msgs)
        base = [{"role": m["role"], "content": m["content"]} for m in msgs]
        fake._turn_rows = base + list(rounds[min(i, len(rounds) - 1)])
        fake._turn_rows_base = len(base)
        yield StreamEvent(type="text_delta",
                          text=replies[min(i, len(replies) - 1)])
        yield StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0)

    fake.stream_message = MagicMock(side_effect=_stream)
    fake._seen = seen
    return fake


def _engine(agent_tree, client):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        return AgentEngine(repo_dir=agent_tree, backend="cli",
                           mode="quick", pack_dir=agent_tree)


def test_the_rounds_go_into_the_history_before_the_answer(agent_tree):
    client = _client(["Die Datei setzt x auf 1."], [_round("c1")])
    engine = _engine(agent_tree, client)
    engine.stream_response("Was steht in a.py?")
    roles = [m["role"] for m in engine.messages]
    assert roles[-4:] == ["user", "assistant", "tool", "assistant"]
    assert engine.messages[-3]["tool_calls"][0]["id"] == "c1"
    assert engine.messages[-1]["content"].startswith("Die Datei setzt")


def test_the_next_turn_sees_what_the_last_one_ran(agent_tree):
    client = _client(["Die Datei setzt x auf 1.", "Ja, x = 1, ich habe sie gelesen."],
                     [_round("c1"), []])
    engine = _engine(agent_tree, client)
    engine.stream_response("Was steht in a.py?")
    engine.stream_response("Hast du sie wirklich gelesen?")
    second = client._seen[1]
    assert any(m.get("role") == "tool" and m.get("tool_call_id") == "c1"
               for m in second)
    assert any(m.get("tool_calls") for m in second if m.get("role") == "assistant")


def test_a_client_that_does_not_declare_it_changes_nothing(agent_tree):
    client = _client(["Antwort."], [_round("c1")], keeps=False)
    engine = _engine(agent_tree, client)
    engine.stream_response("Was steht in a.py?")
    assert [m["role"] for m in engine.messages[-2:]] == ["user", "assistant"]


def test_the_sanitiser_repairs_a_broken_round(agent_tree):
    client = _client(["Antwort."], [[]])
    engine = _engine(agent_tree, client)
    engine.messages.extend([
        {"role": "user", "content": "q"},
        {"role": "tool", "tool_call_id": "gone", "content": "orphan"},
        {"role": "assistant", "content": "a"}])
    engine._sanitize_messages()
    assert not any(m.get("role") == "tool" for m in engine.messages)


def test_the_client_passes_the_pairing_keys_through():
    from delfin.agent.api_client import OpenAIClient
    src = inspect.getsource(OpenAIClient.stream_message)
    assert 'for key in ("tool_calls", "tool_call_id", "name")' in src
    assert "self._turn_rows_base = len(api_messages)" in src
    assert OpenAIClient.KEEPS_TOOL_HISTORY is True
