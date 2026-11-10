"""LK1 — engine-level regression: a second leaked DSML call after REAL
progress IS re-requested within the same user turn.

The 06.11 defect (RED control, current engine): within one user turn the
model leaks call A (re-requested, runs), makes one REAL tool call, then
leaks call B.  The one-shot latch (``_dsml_leak_spent``) refuses the second
re-request, so B never runs and the turn ends with only a "did not run"
notice — a session sat idle ~4 hours with uncommitted work because of it.

This test drives the REAL ``AgentEngine.stream_response`` with a fake model
stream that (call 1) leaks A, (call 2) makes a real tool call then leaks B.
The DESIRED behaviour is that B is re-requested (a third stream call).
On the unpatched engine only TWO stream calls happen and B is never
re-requested, so this test is RED.  After the LK1 engine hook (the policy
in ``leak_recovery.py`` + the re-request loop) it is GREEN.

No tool is actually executed here: the fake streams ``tool_use`` /
``tool_result`` events like the real api_client would, and the engine's
``tool_result`` handling (which, after the fix, records real progress on
the re-request policy) is exactly what lets the second leak legitimately be
re-requested.  ``leaked_tool_call_dominates`` must not see the real-call
prose as making the answer non-dominant, so the answer is written so the
leak dominates at the end.
"""
from __future__ import annotations

import pytest
from unittest.mock import MagicMock, patch

from delfin.agent.api_client import StreamEvent


LEAK_A = (
    '<invoke name="bash">'
    '<parameter name="command">git commit -m x</parameter>'
    "</invoke>"
)
LEAK_B = (
    '<invoke name="bash">'
    '<parameter name="command">git push</parameter>'
    "</invoke>"
)


def _progress_client(calls):
    """client.stream_message follows a script; yields raw StreamEvents."""
    client = MagicMock()
    client.model = "deepseek-v4-flash"
    it = iter(calls)

    def stream(system, messages, max_tokens=4096, session_id="",
               thinking_budget=0, **kw):
        fn = next(it)
        yield from fn()

    client.stream_message = MagicMock(side_effect=stream)
    return client


def _engine(client, agent_tree):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client",
               return_value=client):
        eng = AgentEngine(repo_dir=agent_tree, backend="api",
                          provider="deepseek", model=client.model,
                          mode="quick", pack_dir=agent_tree)
    eng.client = client
    return eng


def _text(t):
    return StreamEvent(type="text_delta", text=t)


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


def test_a_tool_result_records_real_progress_on_the_rerequest_policy(
        agent_tree):
    """The progress signal works: a tool_result event is what lets a later
    leak be re-requested again (guard against the hook being a no-op)."""
    from delfin.agent.leak_recovery import DSMLReplayPolicy
    p = DSMLReplayPolicy()
    assert p.allowed_for(True) is True
    p.spend()                       # leak A re-requested
    assert p.allowed_for(True) is False        # no progress yet
    p.note_real_tool_call()         # a real tool call completed
    assert p.allowed_for(True) is True


def test_second_leak_after_a_real_tool_call_is_rerequested(agent_tree):
    """Target 1 of the brief.  call1 leaks A; call2 makes a real tool call
    (tool_use + tool_result) and then leaks B.  The engine must re-request
    B (call 3).  RED on the current one-shot engine, GREEN after the hook.
    """
    def run_only():
        return [  # call 1: leak A (dominant)
            _text(LEAK_A),
        ]

    def correct_with_progress_then_leak_b():
        return [  # call 2: real tool call A, then leak B
            _text("I ran it.\n"),
            StreamEvent(type="tool_use", tool_name="Bash",
                        tool_input='{"command": "git commit -m x"}'),
            StreamEvent(type="tool_result", tool_name="Bash",
                        tool_output="[main abc123d] x\n"),
            _text(LEAK_B),
        ]

    def second_correction():
        return [  # call 3 (after the fix): re-request for B
            _text("Done running the push."),
        ]

    client = _progress_client([run_only, correct_with_progress_then_leak_b,
                               second_correction])
    eng = _engine(client, agent_tree)
    out = eng.stream_response("commit and push")
    # The re-request for B must happen: call 3 is consumed.  On the
    # unpatched one-shot engine only call 1 and call 2 run, call 3 is never
    # consumed, and the assertion below fails => RED.
    assert client.stream_message.call_count == 3, (
        f"stream_message called {client.stream_message.call_count} times; "
        f"the second leaked call B must be re-requested (expected 3)")
    assert LEAK_B not in out.replace("\n", ""), (
        "a re-requested leak B must not be returned as the answer")


def test_two_leaks_in_a_row_are_not_re_requested(agent_tree):
    """Never-loops core: the same leak twice (no real call between) is
    annotated and stopped — no second re-request, one stream call after the
    first leak, exactly as today.  Must stay green after the fix."""
    def run_only():
        return [_text(LEAK_A)]

    def leak_again():
        return [_text(LEAK_A)]

    client = _progress_client([run_only, leak_again])
    eng = _engine(client, agent_tree)
    out = eng.stream_response("commit")
    # Leak A is re-requested once (call 2); the correction leaks A again,
    # and without a real tool call between there is NO second re-request —
    # the turn makes exactly 2 stream calls.  This is the never-loops
    # invariant; the "did not run" caveat only surfaces on the follow-up
    # turn, not within the same turn.
    assert client.stream_message.call_count == 2
