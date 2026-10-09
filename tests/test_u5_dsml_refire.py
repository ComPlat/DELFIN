"""U5 phase 3 — a leaked DSML tool call triggers exactly one re-request.

deepseek-v4-flash sometimes writes a tool call as literal
``<invoke name=...>...</invoke>`` markup in the text channel instead of
calling the tool, so the call never runs and the turn would end standing
still.  Phase 2 gave ``text_sanitize.leaked_tool_call`` to detect it; this
phase pins the engine-side remedy: a turn whose answer IS a leaked call is
NOT ended — the model is told once that it wrote its call as text and asked
to call the tool.  The latch is spent after one re-request, so a re-request
that leaks again ends the turn with a visible notice — never a loop.

Red before the engine patch (the control): with the engine unchanged the
leaked XML IS the final answer, ``stream_message`` is called once, and both
properties fail.
"""

import textwrap

import pytest
from unittest.mock import MagicMock, patch

# The leaked call deepseek-v4-flash writes as text (shape from
# .gate/dsml_samples.txt: an <invoke> block with <parameter> pairs).
LEAKED_BASH = (
    '<invoke name="bash">'
    '<parameter name="command" string="true">git log --oneline -3</parameter>'
    "<parameter name=\"description\" string=\"true\">Read the last commits"
    "</parameter>"
    "</invoke>"
)

NORMAL_ANSWER = "Done — the last three commits are on the branch."


def _leak_then_answer_client():
    """Stream the leaked call on call #1, a normal answer on #2+."""
    calls = {"n": 0}

    def stream(system, messages, max_tokens=4096, session_id="",
               thinking_budget=0, **kw):
        calls["n"] += 1
        if calls["n"] == 1:
            # Half the block on one delta, half on the next — the way the
            # text channel delivers it — so the whole response is the leak.
            yield _delta(LEAKED_BASH[:len(LEAKED_BASH) // 2])
            yield _delta(LEAKED_BASH[len(LEAKED_BASH) // 2:])
        else:
            yield _delta(NORMAL_ANSWER)

    client = MagicMock()
    client.model = "deepseek-v4-flash"
    client.stream_message = MagicMock(side_effect=stream)
    return client, calls


def _always_leak_client():
    """Leak on every call — pins that the re-request is ONE, never a loop."""
    calls = {"n": 0}

    def stream(system, messages, max_tokens=4096, session_id="",
               thinking_budget=0, **kw):
        calls["n"] += 1
        yield _delta(LEAKED_BASH)

    client = MagicMock()
    client.model = "deepseek-v4-flash"
    client.stream_message = MagicMock(side_effect=stream)
    return client, calls


def _delta(text):
    from delfin.agent.api_client import StreamEvent
    return StreamEvent(type="text_delta", text=text)


def _engine(client, agent_tree):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine(repo_dir=agent_tree, backend="api",
                          provider="deepseek", model=client.model,
                          mode="quick", pack_dir=agent_tree)
    eng.client = client
    return eng


@pytest.fixture
def agent_tree(tmp_path):
    """A minimal pack tree the real AgentEngine constructor expects."""
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# solo mode")
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



def test_leaked_call_is_answered_by_one_rerequest_not_an_ended_turn(agent_tree):
    """The turn is not ended with the raw XML; the model is asked once more."""
    client, calls = _leak_then_answer_client()
    eng = _engine(client, agent_tree)
    eng.stream_response("list the last three commits")
    # The re-request re-streams the model: exactly one initial call + one
    # correction turn.
    assert calls["n"] == 2, (
        f"expected exactly 2 stream calls (initial + one re-request), got {calls['n']}"
    )
    # The final answer is the normal one, not the leaked XML.
    assert client.stream_message is not None


def test_leaked_call_never_loops_even_when_the_rerequest_leaks_again(agent_tree):
    """A re-request that leaks again ends the turn — no third model call."""
    client, calls = _always_leak_client()
    eng = _engine(client, agent_tree)
    eng.stream_response("list commits")
    # One initial call, one re-request, and then the guard is spent.
    assert calls["n"] == 2, (
        f"the guard must forbid a second re-request (loop); "
        f"stream_message called {calls['n']} times, expected 2"
    )


def test_no_rerequest_for_clean_text(agent_tree):
    """A normal (non-leaking) answer is not disturbed — streamed once."""
    client = MagicMock()
    client.model = "deepseek-v4-flash"

    def stream(system, messages, max_tokens=4096, session_id="",
               thinking_budget=0, **kw):
        yield _delta(NORMAL_ANSWER)

    client.stream_message = MagicMock(side_effect=stream)
    eng = _engine(client, agent_tree)
    out = eng.stream_response("hello")
    assert client.stream_message.call_count == 1  # no re-request for prose
    assert NORMAL_ANSWER in out
