"""U5 phase 3 — adversarial wiring guard: instructive DSML prose must NOT re-request.

The builder's red control pins the happy path (a real leaked call → one
re-request, never a loop) and the plain-prose negative control. This guard
covers the adversarial sibling the wiring is most likely to get wrong: an
answer that merely *shows* ``<invoke>`` markup as an example during normal
conversation (teaching the user the DSL) is NOT a leaked call, so the engine
must not issue a re-request — a wiring that re-requests on any
``leaked_tool_call``-shaped content would loop the model on instructive prose.

The detection seam already guarantees ``text_sanitize.leaked_tool_call``
returns ``None`` for instructive-prose input (phase 2 corpus:
test_prose_mention_is_not_a_call, test_code_snippet_is_not_a_call,
test_xml_doc_is_not_a_call). This test pins the *wiring* to honour that
``None``: stream_message is called exactly once and the prose text is
returned verbatim.

Green now (no wiring yet → no re-request) and green after a correct patch;
it flips red only if a future wiring re-requests on instructive prose — which
is precisely the false-positive-into-loop we must stop. Committed green as a
regression guard; the red-control burden is carried by test_u5_dsml_refire.py.
"""
import textwrap

import pytest
from unittest.mock import MagicMock, patch

# Prose that shows a COMPLETE, valid-looking <invoke> block as an example,
# but is clearly explanatory, not a call to run now. leaked_tool_call must
# classify this as not-a-call (and the wiring must not re-request).
INSTRUCTIVE_DSML = (
    "To run a command, the model writes markup like "
    '<invoke name="bash">'
    '<parameter name="command" string="true">git status</parameter>'
    '<parameter name="description" string="true">show working tree</parameter>'
    "</invoke>"
    " — that is the DSL shape, but I am just explaining it to you, "
    "not scheduling anything now."
)


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


def test_instructive_dsml_prose_does_not_trigger_rerequest(agent_tree):
    """An answer that merely shows <invoke> markup must be streamed once."""
    client = MagicMock()
    client.model = "deepseek-v4-flash"
    client.stream_message = MagicMock(
        side_effect=lambda system, messages, max_tokens=4096,
        session_id="", thinking_budget=0, **kw: iter([_delta(INSTRUCTIVE_DSML)]))
    eng = _engine(client, agent_tree)
    out = eng.stream_response("explain the markup")
    assert client.stream_message.call_count == 1, (
        f"instructive DSML prose must not trigger a re-request; "
        f"stream_message called {client.stream_message.call_count} times")
    assert INSTRUCTIVE_DSML in out, "the instructive text must be returned verbatim"
