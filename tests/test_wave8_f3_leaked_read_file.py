"""Wave8/f3 — pin the exact field-reported DSML leak shape.

deepseek-v4-flash sometimes writes a tool call as literal markup in the text
channel instead of calling the tool, so the call never runs. The field report
this task pins is a ``read_file`` call whose arguments include an integer-like
``<parameter name="limit">201</parameter>`` and a worktree-relative
``<parameter name="path">scratch/x.txt</parameter>``, written with whitespace
between parameters:

    <invoke name="read_file"> <parameter name="limit">201</parameter>
    <parameter name="path">scratch/x.txt</parameter> </invoke>

These tests pin that exact shape, in the two layers that matter:

1. the parser recognises it and recovers both arguments (the integer value is
   carried as its textual form, exactly as a real call would place it);
2. the engine remedy for a round whose whole answer is that leak is to ask
   once more — it NEVER executes the leaked call, because executing text would
   bypass the argument validation and the permission gate a real call goes
   through. The parser is wired to a re-request, never to a dispatcher.
"""
import textwrap

import pytest
from unittest.mock import MagicMock, patch

from delfin.agent.text_sanitize import (
    leaked_tool_call,
    leaked_tool_call_dominates,
)

# The field-report shape: read_file, integer limit, scratch/ path, whitespace
# between the <parameter> blocks.
LEAKED_READ_FILE = (
    '<invoke name="read_file">'
    ' <parameter name="limit">201</parameter>'
    ' <parameter name="path">scratch/x.txt</parameter>'
    ' </invoke>'
)
NORMAL_ANSWER = "scratch/x.txt holds 201 lines of scratch notes."


def test_parser_recovers_the_field_reported_read_file():
    """The exact leaked read_file markup is a call with both args recovered.
    The integer limit is carried as its textual form ('201'), matching how a
    real call serialises an integer argument."""
    got = leaked_tool_call(LEAKED_READ_FILE)
    assert got is not None
    assert got[0] == "read_file"
    assert got[1] == {"limit": "201", "path": "scratch/x.txt"}


def test_dominates_accepts_the_bare_field_reported_read_file():
    """A round whose FINAL text is exactly this leak (no real tool call, no
    surrounding prose) is the leak, so the re-request fires."""
    assert leaked_tool_call_dominates(LEAKED_READ_FILE) == (
        "read_file", {"limit": "201", "path": "scratch/x.txt"})
    # Whitespace alone around it changes nothing: still the bare call.
    assert leaked_tool_call_dominates(
        "  \n" + LEAKED_READ_FILE + "\n  ") is not None


def _delta(text):
    from delfin.agent.api_client import StreamEvent
    return StreamEvent(type="text_delta", text=text)


# --- engine level -------------------------------------------------------------


def test_engine_re_requests_once_for_a_leaked_read_file(agent_tree):
    """A turn that is only the leaked read_file markup (no real tool call) is
    NOT executed and NOT left standing still: the model is asked exactly once
    to call the tool, and the leak is never dispatched as a call."""
    calls = {"n": 0}

    def stream(system, messages, max_tokens=4096, session_id="",
               thinking_budget=0, **kw):
        calls["n"] += 1
        yield _delta(LEAKED_READ_FILE) if calls["n"] == 1 else _delta(NORMAL_ANSWER)

    client = MagicMock()
    client.model = "deepseek-v4-flash"
    client.stream_message = MagicMock(side_effect=stream)
    eng = _engine(client, agent_tree)
    out = eng.stream_response("read scratch/x.txt")
    # one initial call + exactly one re-request; the re-request takes the
    # correction string naming the leaked tool.
    assert calls["n"] == 2, (
        f"expected exactly one re-request (never a loop), got {calls['n']}")
    # The final answer is the real one, not the leaked markup.
    assert NORMAL_ANSWER in out


def test_an_answer_that_shows_the_markup_is_not_a_leaked_call(agent_tree):
    """A real answer that merely quotes the markup survives unchanged — the
    re-request fires only when the round produced no real tool call and the
    final text IS the leaked call."""
    prose = (
        "To read a file you would call " + LEAKED_READ_FILE
        + " which opens it. Here is what it contains instead."
    )
    client = MagicMock()
    client.model = "deepseek-v4-flash"

    def stream(system, messages, max_tokens=4096, session_id="",
               thinking_budget=0, **kw):
        yield _delta(prose)

    client.stream_message = MagicMock(side_effect=stream)
    eng = _engine(client, agent_tree)
    out = eng.stream_response("how do I read a file?")
    assert client.stream_message.call_count == 1  # no re-request for prose
    assert prose in out


# ---- shared engine fixture ---------------------------------------------------

def _engine(client, agent_tree):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine(repo_dir=agent_tree, backend="api",
                          provider="deepseek", model=client.model,
                          mode="quick", pack_dir=agent_tree)
    eng.client = client
    return eng


# A pytest-style agent_tree fixture (the real AgentEngine constructor needs a
# minimal pack tree). Replicated here so this file is self-contained.
@pytest.fixture
def agent_tree(tmp_path):
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
