"""Tests for tool-result-aware context editing in long tool-call loops.

Contract: when accumulated ``role=="tool"`` output exceeds the budget,
the OLDEST tool results are elided to a placeholder while the most recent
ones stay verbatim; user/system/assistant messages and every
``tool_call_id`` are preserved (protocol-safe). Under budget -> no-op.
"""

from __future__ import annotations

from delfin.agent.api_client import _elide_old_tool_results, _ELIDED_PREFIX


def _tool(i, n):
    return {"role": "tool", "tool_call_id": f"c{i}", "content": "X" * n}


def test_noop_under_budget():
    msgs = [{"role": "system", "content": "s"},
            {"role": "user", "content": "u"},
            _tool(1, 100), _tool(2, 100)]
    before = [dict(m) for m in msgs]
    assert _elide_old_tool_results(msgs, char_budget=10000) == 0
    assert msgs == before                       # untouched


def test_oldest_elided_recent_kept():
    msgs = [{"role": "system", "content": "s"}]
    # 10 tool results of 1000 chars each = 10000; budget 4000, keep 3.
    for i in range(10):
        msgs.append({"role": "assistant", "content": f"reasoning {i}",
                     "tool_calls": [{"id": f"c{i}"}]})
        msgs.append(_tool(i, 1000))
    n = _elide_old_tool_results(msgs, char_budget=4000, keep_recent=3)
    assert n >= 1
    tools = [m for m in msgs if m.get("role") == "tool"]
    # the three most recent are kept verbatim ...
    assert all(not t["content"].startswith(_ELIDED_PREFIX) for t in tools[-3:])
    assert all(len(t["content"]) == 1000 for t in tools[-3:])
    # ... the oldest were elided ...
    assert tools[0]["content"].startswith(_ELIDED_PREFIX)
    # ... but every tool_call_id survives (protocol stays valid) ...
    assert [t["tool_call_id"] for t in tools] == [f"c{i}" for i in range(10)]
    # ... and reasoning/system/user are never touched.
    assert all(m["content"].startswith("reasoning")
               for m in msgs if m.get("role") == "assistant")


def test_assistant_and_user_never_elided():
    big = "Y" * 50000
    msgs = [{"role": "user", "content": big},
            {"role": "assistant", "content": big},
            _tool(1, 50000)]
    _elide_old_tool_results(msgs, char_budget=1000, keep_recent=0)
    assert msgs[0]["content"] == big            # user huge -> untouched
    assert msgs[1]["content"] == big            # assistant huge -> untouched


def test_already_elided_not_doubled():
    msgs = [_tool(1, 2000), _tool(2, 2000), _tool(3, 2000)]
    _elide_old_tool_results(msgs, char_budget=1500, keep_recent=0)
    # run again — placeholders must not be re-wrapped
    _elide_old_tool_results(msgs, char_budget=1500, keep_recent=0)
    elided = [m for m in msgs if m["content"].startswith(_ELIDED_PREFIX)]
    for m in elided:
        assert m["content"].count(_ELIDED_PREFIX) == 1


# --- Minimal fake OpenAI stream, enough to drive N tool rounds -------------

class _D:
    def __init__(self, content=None, tool_calls=None):
        self.content, self.tool_calls = content, tool_calls


class _C:
    def __init__(self, delta, finish=None):
        self.delta, self.finish_reason = delta, finish


class _U:
    prompt_tokens, completion_tokens = 5, 3


class _Ch:
    def __init__(self, choices, usage=None):
        self.choices, self.usage = choices, usage


class _S:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _TC:
    """One streamed tool_call delta — a read_file the gate will refuse,
    which is enough: the round happened and its result is in the list."""
    def __init__(self, i):
        self.index, self.id, self.type = 0, f"c{i}", "function"
        self.function = type("F", (), {"name": "read_file",
                                       "arguments": '{"path": "nope.txt"}'})()


def _driven_client(monkeypatch, *, tool_rounds):
    from delfin.agent import mcp_client as M
    from delfin.agent import model_capabilities as mc
    from delfin.agent.api_client import KitToolPermissions, create_client

    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=262_144, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: type(
        "R", (), {"discover_all": lambda s: [],
                  "discover_resources": lambda s: [],
                  "discover_prompts": lambda s: []})())

    import tempfile
    ws = tempfile.mkdtemp()
    client = create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=ws)
    client.set_permissions(KitToolPermissions(workspace=ws, mode="default"))

    state = {"n": 0}

    def _create(**kw):
        state["n"] += 1
        if state["n"] <= tool_rounds:
            return _S([_Ch([_C(_D(tool_calls=[_TC(state["n"])]),
                               finish="tool_calls")], usage=_U())])
        return _S([_Ch([_C(_D(content="done"), finish="stop")], usage=_U())])

    client.client.chat.completions.create = _create
    return client


def test_called_every_round_with_the_context_scaled_budget(monkeypatch):
    """Elision runs once per tool round, on the live message list.

    This was a source-text assertion: it read api_client.py and required
    the call to appear within 2500 characters of the loop head. Anything
    added at the top of that loop pushed it out of the window and the
    test failed while the behaviour was untouched — a failure that names
    the wrong thing is worse than no test. What actually matters is that
    the call happens each round, with the budget scaled to the model's
    context and not the fixed default.
    """
    import delfin.agent.api_client as A

    seen: list[int] = []
    real = A._elide_old_tool_results

    def _spy(messages, *, char_budget, **kw):
        seen.append(char_budget)
        return real(messages, char_budget=char_budget, **kw)

    monkeypatch.setattr(A, "_elide_old_tool_results", _spy)

    client = _driven_client(monkeypatch, tool_rounds=2)
    list(client.stream_message("sys", [{"role": "user", "content": "go"}],
                               max_tokens=64))

    assert len(seen) >= 3, "elision did not run on every round"
    assert all(b == seen[0] for b in seen), "the budget changed mid-turn"
    # A 262k-context model must keep far more tool output than the fixed
    # default, or it spends the turn re-paging files it already read.
    assert seen[0] > A._TOOL_CONTEXT_CHAR_BUDGET, (
        "the fixed default was used instead of the model's context")


# ---------------------------------------------------------------------------
# Context-scaled elision budget (bug 172455: a 262k model must keep far more
# tool output than the fixed 60k-char default, so it stops re-paging files)
# ---------------------------------------------------------------------------

from delfin.agent.api_client import (
    _tool_context_char_budget, _TOOL_CONTEXT_CHAR_BUDGET,
)


class _Caps:
    def __init__(self, ctx):
        self.context_window = ctx


def test_budget_floor_when_caps_missing_or_small():
    # None caps, zero/None window, and small models all keep the 60k floor.
    assert _tool_context_char_budget(None) == _TOOL_CONTEXT_CHAR_BUDGET
    assert _tool_context_char_budget(_Caps(0)) == _TOOL_CONTEXT_CHAR_BUDGET
    assert _tool_context_char_budget(_Caps(None)) == _TOOL_CONTEXT_CHAR_BUDGET
    # 8k model: 8192*0.45*4 = 14745 < 60000 → floor, never smaller than today.
    assert _tool_context_char_budget(_Caps(8192)) == _TOOL_CONTEXT_CHAR_BUDGET


def test_budget_scales_up_for_big_context():
    big = _tool_context_char_budget(_Caps(262_144))
    assert big == int(262_144 * 0.45 * 4)        # ~471k chars
    assert big > _TOOL_CONTEXT_CHAR_BUDGET * 5   # far more than the fixed default


def test_big_budget_keeps_what_small_budget_elides():
    # 12 tool results of 8000 chars = 96000 total.
    def _mk():
        m = [{"role": "system", "content": "s"}]
        for i in range(12):
            m.append({"role": "assistant", "content": f"r{i}",
                      "tool_calls": [{"id": f"c{i}"}]})
            m.append({"role": "tool", "tool_call_id": f"c{i}",
                      "content": "X" * 8000})
        return m
    small = _mk()
    big = _mk()
    # Default (small) budget elides; a 262k-scaled budget (~471k) keeps all.
    n_small = _elide_old_tool_results(
        small, char_budget=_tool_context_char_budget(_Caps(8192)))
    n_big = _elide_old_tool_results(
        big, char_budget=_tool_context_char_budget(_Caps(262_144)))
    assert n_small >= 1            # 96k > 60k floor → some elided
    assert n_big == 0              # 96k < 471k → nothing elided, no re-paging
