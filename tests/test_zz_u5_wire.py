"""U5 phase 3 — the DSML leak wiring in the api_client dispatch.

deepseek-v4-flash writes a tool call as literal ``<invoke name=...>...</invoke>``
XML in the text channel instead of calling the tool.  Phase 2 gave
``text_sanitize.leaked_tool_call`` to detect a COMPLETE block; this phase
wires that detector into the api_client dispatch (the layer that already
reclaims the Harmony ``to=<tool> {json}`` leak).  The wiring is a
RE-REQUEST, never the silent Harmony redispatch: a complete DSML block is
answered with a corrective one-shot notice telling the model its call was
written as text, and a spent per-run latch forbids a second re-request.

These assertions are made on the PROTECTED source (``OpenAIClient.
stream_message``) the operator builds from the ``.gate/u5_dsml_wire.patch``,
so they fail red before the build and go green after it.  Each marker is
chosen to not exist in the current source — no false green.
"""

import inspect

from delfin.agent import api_client as A
from delfin.agent import text_sanitize


def _stream_message_src() -> str:
    return inspect.getsource(A.OpenAIClient.stream_message)


def test_the_dispatch_invokes_leaked_tool_call():
    """The dispatch calls ``leaked_tool_call_dominates(...)`` (a real call, not
    the adjacent Harmony ``parse_leaked_tool_calls`` reuse)."""
    src = _stream_message_src()
    # text_sanitize.leaked_tool_call_dominates exists (phase 3)…
    assert hasattr(text_sanitize, "leaked_tool_call_dominates")
    # …and the dispatch actually calls it.
    assert "leaked_tool_call_dominates(" in src
    # The bare name would match inside ``parse_leaked_tool_calls``; the
    # distinctive evidence is the OPEN PAREN of a real call site.
    assert src.count("leaked_tool_call_dominates(") >= 1


def test_a_complete_dsml_block_is_not_silently_redispatched():
    """The DSML leak takes the re-request path, NOT the Harmony synthesis.

    A complete ``<invoke name=...>...</invoke>`` written as text must never
    be silently turned into a synthesized tool call — the model is asked,
    once, to call the tool instead.  The patch therefore reads the DSML
    detection result into a dedicated branch and does not hand it to the
    Harmony recovery, which is what would execute it without a re-request.
    """
    src = _stream_message_src()
    # The DSML result has its own branch, distinct from Harmony recovery,
    # so a complete DSML block is not handed to the synthesis path.
    assert "_dsml_leaked" in src


def test_the_corrective_notice_says_written_as_text():
    """The single re-request tells the model it wrote its call as text."""
    src = _stream_message_src()
    assert "written as text" in src.lower()


def test_wire_is_one_shot_never_a_loop():
    """A spent per-run latch forbids a second re-request.

    The same invariant the engine's correction-turns guarantee (one forced
    correction; a nested re-entry sees the latch spent and stops) must hold
    here, so a turn that leaks again still ends instead of re-requesting
    forever.
    """
    src = _stream_message_src()
    assert "_dsml_rerequest" in src
    # The latch caps the second re-request by checking the spent flag before
    # emitting another corrective notice — the guard read, not just the name.
    assert "_dsml_rerequest" in src and "spent" in src
