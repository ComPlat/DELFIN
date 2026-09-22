"""Six remembers and a fragment where the answer should have been.

The memory addendum asks every role to persist durable facts as it works,
and one model takes that as the work. Measured 2026-09-08 through the
engine, workflow_verify_after_modify in acceptEdits mode, kit.glm-5.3:

    event kinds: thinking 64, tool 6, result 6, text 6
    tools: remember x 6
    answer: "Functional setzen und danach anzeigen:"

Six consecutive `remember` calls and then a fragment ending in a colon,
where the two ACTION lines belonged. DeepSeek called it zero times on
every dashboard task of the same run and solved that task in two calls at
quality 88 — so this is one model over-applying a shared rule, and the fix
belongs to the model, not to the rule.

The existing no-progress guard cannot catch it: it keys on tool name AND
arguments, and six remembers with different content read as progress.
"""

import json
import threading
import types

import pytest

from delfin.agent import api_client as ac
from delfin.agent.model_profiles import (
    ModelProfile, _COERCE, get_profile)


def test_the_cap_is_universal_not_per_model():
    """The loop-shape is not a model quirk: any model can take the
    memory addendum as the work. The cap is the default for every
    profile — kit.glm-5.3 merely measured it first (six calls, a
    fragment for an answer, 2026-09-08).

    Threshold from the trace archive (2026-09-18, 33 sessions, 99
    ts-gap-grouped turns): the highest number of memory writes any
    recorded turn legitimately made is ONE (`remember` once in the
    whole archive, `forget` zero). Every tool that legitimately
    repeats — bash to 31/turn, task_create to 5 — is outside the
    memory-write class. Two allowed, refused from the third on, has
    zero false positives in the entire archive.
    """
    assert get_profile("kit.glm-5.3").max_memory_writes_per_turn == 2
    # Not a GLM special case: every model carries the same default.
    for model in ("kit.deepseek-v4-flash", "kit.qwen3.5-397b-A17b",
                  "sonnet", "azure.gpt-5.4", "totally-unknown-model"):
        assert get_profile(model).max_memory_writes_per_turn == 2, model
    assert ModelProfile().max_memory_writes_per_turn == 2
    # 0 stays a valid opt-OUT for a user who wants no cap at all.
    assert _COERCE["max_memory_writes_per_turn"] is int


def test_only_writes_are_counted():
    """A memory READ is not the loop this breaks."""
    assert ac._MEMORY_WRITE_TOOLS == {"remember", "forget"}
    for read_tool in ("history_search", "history_get", "read_file"):
        assert read_tool not in ac._MEMORY_WRITE_TOOLS


def _client(model: str, calls: int):
    """A client whose model asks to remember ``calls`` times, then stops."""
    state = {"round": 0, "executed": []}

    class _Stream:
        def __iter__(self):
            state["round"] += 1
            if state["round"] > calls:
                delta = types.SimpleNamespace(
                    content="fertig", tool_calls=None, reasoning_content=None)
                return iter((types.SimpleNamespace(
                    usage=None,
                    choices=[types.SimpleNamespace(
                        delta=delta, finish_reason="stop")]),))
            delta = types.SimpleNamespace(
                content=None, reasoning_content=None,
                tool_calls=[types.SimpleNamespace(
                    index=0, id=f"c{state['round']}",
                    function=types.SimpleNamespace(
                        name="remember",
                        arguments=json.dumps(
                            {"text": f"fact {state['round']}"})))])
            return iter((types.SimpleNamespace(
                usage=None,
                choices=[types.SimpleNamespace(
                    delta=delta, finish_reason="tool_calls")]),))

        def close(self):
            pass

    class _Stub:
        class chat:
            class completions:
                @staticmethod
                def create(**kw):
                    return _Stream()

    c = ac.OpenAIClient.__new__(ac.OpenAIClient)
    c.client = _Stub()
    c.model = model
    c._provider = "kit"
    c._base_url = "https://ki-toolbox.scc.kit.edu/api/v1"
    c._api_key = "x"
    c.effort = ""
    c._permissions = None
    c.on_model_switched = None
    c._steer_lock = threading.Lock()
    c._steer_queue = []
    c._run_notes = []
    c._steer_msgs = []
    c._stop_flag = False
    return c, state


def _run(client, rounds: int = 8):
    results = []
    try:
        for ev in client.stream_message(
                messages=[{"role": "user", "content": "setz das functional"}],
                system="t", max_tokens=64):
            if getattr(ev, "type", "") == "tool_result":
                results.append(ev.tool_output or "")
            if len(results) >= rounds:
                break
    except Exception:
        pass
    return results


def test_a_capped_model_is_refused_after_its_budget():
    client, _ = _client("kit.glm-5.3", calls=6)
    results = _run(client)
    refusals = [r for r in results if "memory write refused" in r]
    assert len(results) >= 3, results
    # Two go through, the rest are refused.
    assert not any("refused" in r for r in results[:2]), results[:2]
    assert refusals, "the third remember was not held back"


def test_the_refusal_says_what_to_do_instead():
    """A silent success is a loop; the model has to hear that it is doing
    the wrong job."""
    client, _ = _client("kit.glm-5.3", calls=6)
    results = _run(client)
    refusal = next(r for r in results if "memory write refused" in r)
    assert "finish the task" in refusal
    assert "2 fact" in refusal


def test_every_other_model_is_capped_the_same_way():
    """Not a model-specific special case: the detector is keyed on the
    tool NAME alone (api_client._MEMORY_WRITE_TOOLS = remember/forget),
    so a model that never had the problem carries the same guard.
    Two go through, the third is refused, and the turn is NOT aborted
    — the user still gets their answer.
    """
    client, _ = _client("kit.deepseek-v4-flash", calls=6)
    results = _run(client)
    assert results, "no tool ran at all"
    assert not any("refused" in r for r in results[:2]), results[:2]
    assert any("memory write refused" in r for r in results), results
    # The turn continues past the refusal rather than dying: a later
    # round still produced a tool result.
    assert len(results) >= 3, results


def test_a_legitimate_repeater_is_untouched():
    """bash_status polling a running job is the same tool, called over
    and over with different round-trip timings — exactly the shape the
    round-signature guard is blind to and this class must NOT cover.
    Only memory WRITES are capped; nothing else in _MEMORY_WRITE_TOOLS
    outside remember/forget is affected.
    """
    assert "bash_status" not in ac._MEMORY_WRITE_TOOLS
    assert "bash" not in ac._MEMORY_WRITE_TOOLS
    assert ac._MEMORY_WRITE_TOOLS == frozenset({"remember", "forget"})
