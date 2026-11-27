"""A sub-agent's request must carry its prompt in ONE system row, at index 0.

The symptom was measured, not theorised: three `subagent` calls in one real
session each died 0.1 s after being dispatched, every one reporting

    Error code: 400 - System message must be at the beginning.

so no specialist preset could run at all, and the delegation table, the
per-preset mode clamp and the routing cues were unreachable in practice even
though the tool itself is offered.

Two places hand the child its prompt: ``run_subagent`` seeds ``messages``
with a ``system`` row AND passes the same text as the ``system=`` argument of
``stream_message``. An OpenAI-compatible client cannot honour both -- it puts
``system=`` at index 0 and then copies every ``messages`` entry verbatim, so
the request carries a second system row at index 1 and the endpoint refuses
the whole call. The Anthropic path takes ``system`` as its own parameter and
the CLI backend runs its own prompt, which is why only the OpenAI-shaped
providers (KIT, blablador, openai, ollama/vLLM) broke, and why nothing in the
suite noticed.

Semantics pinned here: the child's prompt reaches the wire exactly once, as
the first row, and the delegate's instructions are still there -- dropping the
duplicate must not drop the prompt. The delegated task must still arrive as a
user turn. Nothing about the delegate's tool surface changes.
"""

from __future__ import annotations

import tempfile
from types import SimpleNamespace

import delfin.agent.api_client as A
from delfin.agent import subagents as SA


class _RecordingCompletions:
    """Stands in for ``client.chat.completions`` and keeps each payload.

    ``close()`` is not decoration: the client closes the stream it is handed,
    so a bare iterator makes the run die on an AttributeError before the
    payload can be read back.
    """

    def __init__(self, sink: list[dict]):
        self._sink = sink

    def create(self, **kwargs):
        self._sink.append(kwargs)

        class _Stream:
            def __iter__(self):
                return iter([])

            def close(self):
                pass

        return _Stream()


def _openai_shaped_client(sink: list[dict]):
    """The product's own ``OpenAIClient`` with a recording transport.

    Built with ``__new__`` so no key or network is touched; the request
    assembly under test is the shipped code, not a re-implementation.
    """
    client = A.OpenAIClient.__new__(A.OpenAIClient)
    client.model = "gpt-4o"
    client._model = "gpt-4o"
    client._provider = "kit"
    client._base_url = "https://endpoint.invalid/v1"
    client._api_key = "k"
    client._permissions = None
    client._turn_rows = None
    client._turn_rows_base = 0
    client._turn_thinking_budget = 0
    client.client = SimpleNamespace(
        chat=SimpleNamespace(completions=_RecordingCompletions(sink)))
    return client


_TASK = "Report which level of theory each of the two runs used."


def _child_payload(tmp_path=None):
    """Dispatch a real ``run_subagent`` and return the child's wire payload.

    The child client is a shallow copy of the parent, so it inherits the
    recording transport and the REAL ``stream_message``. Only the things
    beside the request are replaced: the permission derivation, the session
    and telemetry writes, and the live-registry updates.
    """
    tmp = str(tmp_path) if tmp_path is not None else tempfile.mkdtemp()
    sink: list[dict] = []
    parent = _openai_shaped_client(sink)

    saved = {}
    replacements = {
        "_isolate_sub_client": lambda _c: None,
        "_derive_perms": lambda *_a, **_k: SimpleNamespace(
            mode="default", workspace=tmp, subagent_depth=1,
            session_allowed_tools=None, session_denied_tools=None),
        "_save_subagent_session": lambda *_a, **_k: None,
        "_write_telemetry": lambda *_a, **_k: None,
        "_running_update": lambda *_a, **_k: None,
        "drop_notes": lambda *_a, **_k: None,
        "read_running": lambda **_k: {},
    }
    for key, value in replacements.items():
        saved[key] = getattr(SA, key)
        setattr(SA, key, value)
    try:
        SA.run_subagent(
            subagent_type="explore",
            description="find the level of theory",
            prompt=_TASK,
            parent_client=parent,
            parent_perms=SimpleNamespace(mode="default", workspace=tmp,
                                         subagent_depth=0),
        )
    finally:
        for key, value in saved.items():
            setattr(SA, key, value)

    assert sink, "the child request never reached the transport"
    return sink[0]


def test_the_child_request_has_exactly_one_system_row():
    rows = _child_payload().get("messages", [])
    system_rows = [i for i, r in enumerate(rows)
                   if r.get("role") in ("system", "developer")]
    assert len(system_rows) == 1, (
        f"the child request carries {len(system_rows)} system rows at "
        f"{system_rows}; an OpenAI-compatible endpoint refuses the request "
        "with 400 'System message must be at the beginning', so the "
        "delegate never runs")
    assert system_rows[0] == 0, "the system row must be the first row"


def test_the_delegate_still_gets_its_preset_prompt():
    """Silencing the 400 must not cost the child its instructions.

    A fix that dropped BOTH rows would stop the error and hand the delegate
    no role, no caps and no read-only instruction -- a quieter and worse
    failure, and the reason this assertion sits beside the one above.
    """
    rows = _child_payload().get("messages", [])
    system_text = " ".join(str(r.get("content") or "") for r in rows
                           if r.get("role") in ("system", "developer"))
    assert system_text.strip(), "the child was sent with no system prompt"
    assert "explore" in system_text.lower(), (
        "the preset's own instructions are missing from the child request")


def test_the_task_reaches_the_child_as_a_user_row():
    rows = _child_payload().get("messages", [])
    users = [r for r in rows if r.get("role") == "user"]
    assert any(_TASK in str(r.get("content") or "") for r in users), (
        "the delegated task never reached the child as a user turn")


def test_the_deduplication_grants_the_child_no_extra_tools():
    """The fix is about the message list, not about what the child may call.

    The child's surface is its preset's own; a request assembled after the
    change must not be wider than one assembled before it.
    """
    payload = _child_payload()
    names = {str((t.get("function") or {}).get("name") or "")
             for t in (payload.get("tools") or [])}
    assert "subagent" not in names, (
        "a depth-1 delegate must not be offered the spawn tool: the nesting "
        "cap keeps a child from fanning out its own children")
