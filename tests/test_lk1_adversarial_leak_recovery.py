"""Adversarial tests for LK1 (reviewer nacht-s11).

The defect under review: a leaked tool call written as text is re-requested
at most ONCE per user turn (engine.py ``_dsml_leak_spent``).  A long
autonomous turn that leaks call A (re-requested, then makes a real tool
call), then leaks call B, ends the turn with B never running and only a
notice -- the 4-hour idle stall of 06.11.

These tests drive the real path: the engine's DSML check and, once the
builder ships it, the policy object ``leak_recovery``.  Until the fix
exists, the progress cases are RED (expected findings) and the
no-progress cases are GREEN (behaviour that must not change).
"""
from __future__ import annotations

import pytest

from delfin.agent.text_sanitize import leaked_tool_call_dominates


# ---------------------------------------------------------------- helpers

def _leak(name: str = "bash", command: str = "git status") -> str:
    """A leaked <invoke> block as deepseek-v4-flash writes it."""
    return (
        f'<invoke name="{name}">\n'
        f'  <parameter name="command">{command}</parameter>\n'
        f"</invoke>"
    )


def _announced_leak(name: str = "bash", command: str = "git status") -> str:
    """The field shape: announce the step, then write the call as text."""
    return f"Let me run the next step.\n{_leak(name, command)}"


# ------------------------------------------------- no-progress (GREEN now)

class TestNoProgressStaysStrict:
    """Two leaks in a row without a real tool call between them: today's
    behaviour is annotate-and-stop.  The LK1 fix must keep exactly this."""

    def test_two_leaks_in_a_row_dominates_each(self):
        # Both blocks dominate the answer individually; the POLICY decides
        # not to re-request the second one because no real call ran
        # between the leaks.  The detector itself must see each one.
        first = _announced_leak("bash", "git add -A")
        second = _announced_leak("bash", "git commit -m x")
        assert leaked_tool_call_dominates(first) == (
            "bash", {"command": "git add -A"})
        assert leaked_tool_call_dominates(second) == (
            "bash", {"command": "git commit -m x"})

    def test_two_leaks_back_to_back_first_dominatees(self):
        both = _announced_leak("bash", "git add -A") + "\n" + _leak(
            "bash", "git commit -m x")
        # The detector sees only the FIRST block as dominant (nothing
        # after the last one, calls back to back) -- the engine latches
        # and a second re-request never happens without a real call
        # between.  That must stay true after LK1.
        result = leaked_tool_call_dominates(both)
        assert result is not None
        assert result[0] == "bash"


# ----------------------------------------------------- progress (RED now)

class TestProgressAllowsSecondRerequest:
    """The LK1 target: after one re-request AND one real tool call
    completed, a SECOND leak in the same user turn is re-requested too
    (bounded by max_per_turn, default 3).  Red on unpatched code: the
    latch spent, the second leak only annotated."""

    def test_policy_allows_after_real_progress(self):
        from delfin.agent.leak_recovery import DSMLReplayPolicy

        p = DSMLReplayPolicy()
        p.spend()                       # leak A re-requested
        p.note_real_tool_call()         # a REAL tool call completed
        assert p.allowed_for(dominates=True) is True


# ------------------------------------------- policy edge cases (adversarial)

class TestPolicyEdgeCases:
    """Attack the budget object itself (builder commit 61c9d01c6)."""

    def test_zero_budget_grants_nothing(self):
        """max_per_turn=0 says NO re-requests this turn -- yet the
        _used==0 branch of allowed_for grants one.  Self-inconsistent
        with exhausted==True (FINDING F1, reviewer nacht-s11)."""
        from delfin.agent.leak_recovery import DSMLReplayPolicy

        p = DSMLReplayPolicy(max_per_turn=0)
        assert p.exhausted is True
        assert p.allowed_for(dominates=True) is False, (
            "exhausted says the turn is spent; allowed_for must agree")

    def test_negative_budget_grants_nothing(self):
        from delfin.agent.leak_recovery import DSMLReplayPolicy

        p = DSMLReplayPolicy(max_per_turn=-1)
        assert p.exhausted is True
        assert p.allowed_for(dominates=True) is False

    def test_falsy_dominates_is_a_refusal_that_spends_nothing(self):
        from delfin.agent.leak_recovery import DSMLReplayPolicy

        for bad in (None, "", 0, [], {}):
            p = DSMLReplayPolicy()
            assert p.allowed_for(bad) is False
            assert p.used == 0, f"refused re-request spent budget: {bad!r}"

# ------------------------------------------------ engine path (real caller)

ENGINE_LEAK = (
    '<invoke name="bash">'
    '<parameter name="command">git status</parameter>'
    "</invoke>"
)
ENGINE_NORMAL = "Done."


def _stream_client(script):
    """client.stream_message follows a script of stream functions."""
    import textwrap
    from unittest.mock import MagicMock, patch

    from delfin.agent.api_client import StreamEvent

    client = MagicMock()
    client.model = "deepseek-v4-flash"
    it = iter(script)

    def stream(system, messages, max_tokens=4096, session_id="",
               thinking_budget=0, **kw):
        fn = next(it)
        for text in fn():
            yield StreamEvent(type="text_delta", text=text)

    client.stream_message = MagicMock(side_effect=stream)

    def _engine(agent_tree):
        from delfin.agent.engine import AgentEngine
        with patch("delfin.agent.engine.create_client",
                   return_value=client):
            eng = AgentEngine(repo_dir=agent_tree, backend="api",
                              provider="deepseek", model=client.model,
                              mode="quick", pack_dir=agent_tree)
        eng.client = client
        return eng

    return client, _engine


@pytest.fixture
def agent_tree(tmp_path):
    """A minimal pack tree the real AgentEngine constructor expects."""
    import textwrap
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# solo mode")
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


class TestEngineNeverLoops:
    """Engine-path adversarial: bounded stream calls, no execution of
    text, no re-request for prose.  All valid on unpatched main (calls==2
    today) AND after the fix (calls<=1+max_per_turn)."""

    def test_leak_forever_is_bounded(self, agent_tree):
        """A model that leaks on EVERY stream call: the number of stream
        calls is bounded by 1 + max_per_turn (default 3), never a loop."""
        import importlib.util
        spec = importlib.util.find_spec("delfin.agent.leak_recovery")
        bound = 3
        if spec is not None:
            lr = importlib.import_module("delfin.agent.leak_recovery")
            bound = getattr(
                getattr(lr, "DSMLReplayPolicy", None), "MAX_PER_TURN", 3)
        client, make = _stream_client(
            [lambda: iter([ENGINE_LEAK])] * 20)
        eng = make(agent_tree)
        eng.stream_response("list commits")
        n = client.stream_message.call_count
        assert n <= 1 + bound, (
            f"stream_message called {n} times; must stay <= 1+{bound}")

    def test_fenced_leak_in_code_example_is_not_rerequested(self, agent_tree):
        """A complete <invoke> block INSIDE a code fence is documentation,
        not a leak: no re-request, the answer is returned as written."""
        fenced = (
            "Here is the markup you must never emit:\n"
            "```\n" + ENGINE_LEAK + "\n```"
        )
        client, make = _stream_client([lambda: iter([fenced])])
        eng = make(agent_tree)
        out = eng.stream_response("show me the format")
        assert client.stream_message.call_count == 1
        assert fenced in out

    def test_unknown_tool_leak_bounded_and_never_executed(self, agent_tree):
        """A leak naming a tool the session does not have: bounded stream
        calls, the text is returned verbatim (never parsed and run as a
        tool call).  The 'did not run' notice itself is emitted inside
        the real stream_message (api_client.py:34478) which this fake
        client replaces -- notice wording is pinned by the U5 tests, not
        by this engine-level fake."""
        unknown = (
            '<invoke name="definitely_not_a_real_tool">'
            '<parameter name="x">1</parameter></invoke>'
        )
        client, make = _stream_client([lambda: iter([unknown])])
        eng = make(agent_tree)
        out = eng.stream_response("do it")
        n = client.stream_message.call_count
        # Bounded either way; and the text was never executed.
        assert n <= 4
        assert unknown in out, "leak text must be returned as written"


# ------------------------------------------------- fresh-turn budget (wiring)

class TestFreshTurnBudgetContract:
    """The policy holds a per-TURN budget but has no reset API.  The
    protected fresh-turn re-arm (engine.py ~:3157) must therefore build a
    FRESH DSMLReplayPolicy() on every real user turn; reusing one instance
    across turns silently bleeds _used and eventually starves every turn.
    These pin that contract."""

    def test_reusing_one_instance_across_turns_bleeds_the_budget(self):
        """If the wiring reuses an instance (no reset), a later turn's FIRST
        re-request is wrongly refused after earlier turns spent the budget.
        max_per_turn=1: turn1 spends 1 -> exhausted; the SAME instance can
        never again grant a fresh turn's first leak."""
        from delfin.agent.leak_recovery import DSMLReplayPolicy
        p = DSMLReplayPolicy(max_per_turn=1)
        # turn 1: first re-request allowed, spent.
        assert p.allowed_for(dominates=True) is True
        p.spend()
        assert p.exhausted is True
        # turn 2, SAME reused instance (the wiring mistake): a fresh turn's
        # first leak must NOT be refused, but a non-reset instance refuses it.
        assert p.allowed_for(dominates=True) is False, (
            "a reused instance cannot honour a fresh user turn - the wiring "
            "must build a new DSMLReplayPolicy() per real user turn")

    def test_fresh_instance_grants_each_turn_first_rerequest(self):
        """The correct wiring builds a fresh instance per real user turn:
        every new turn starts with allowed first re-request."""
        from delfin.agent.leak_recovery import DSMLReplayPolicy
        for _ in range(3):  # three consecutive fresh turns
            p = DSMLReplayPolicy()
            assert p.allowed_for(dominates=True) is True
            assert p.used == 0
            p.spend()
