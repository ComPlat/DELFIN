"""Driven: a stated wait we will not sit out ends the round.

The helpers are checked in test_the_server_says_when_to_come_back.py.
This drives the real ``stream_message`` loop, because the defect lives in
the retry DECISION and a helper that answers correctly while the loop
ignores it would leave it exactly in place.

What makes it concrete (AI@JSC, Blablador API access, read 2026-11-19):
exceeding a model's per-minute request limit answers 429 with
Retry-After and blocks that model for 60 s; every further violation
doubles the block, to a maximum of 3600 s. Our own transient ladder
retries up to seven times with steps to 240 s. So a 30-minute block used
to be asked again inside the block, four times over -- each one a fresh
violation that doubles it. The retry made the outage longer.
"""

from __future__ import annotations

import pytest

import delfin.agent.api_client as A
import delfin.agent.model_capabilities as mc
import delfin.agent.mcp_client as M


class _Delta:
    def __init__(self, content=None):
        self.content = content
        self.tool_calls = None


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish
        self.index = 0


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


def _rate_limited(retry_after: str):
    """A genuine openai.RateLimitError carrying the header.

    Not a stand-in: the retry path asks _is_transient_api_error about the
    exception's type and reads the header off the transport's own
    case-insensitive mapping. A hand-made object would satisfy neither.
    """
    openai = pytest.importorskip("openai")
    try:
        import httpx2 as transport
    except ImportError:                      # pragma: no cover - older SDK
        transport = pytest.importorskip("httpx")
    request = transport.Request(
        "POST", "https://api.blablador.fz-juelich.de/v1/chat/completions")
    response = transport.Response(429, request=request,
                                  headers={"Retry-After": retry_after})
    return openai.RateLimitError("rate limited", response=response, body=None)


@pytest.fixture
def client(monkeypatch, tmp_path):
    # No sleeping in the suite, and a ladder long enough that a retry
    # would be unmistakable if one happened.
    monkeypatch.setenv("DELFIN_STREAM_RETRY_DELAY_SCALE", "0")
    # Driven through an OpenAI-compatible provider that needs no key:
    # the retry decision is in the shared stream loop and knows nothing
    # about which gateway answered, so the guard is checked where it
    # lives rather than on one provider's branch.
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    cl = A.create_client(backend="api", provider="ollama", model="m",
                         cwd=str(tmp_path))
    cl.set_permissions(A.KitToolPermissions(
        workspace=tmp_path, mode="default", task_session_id="sess-block"))
    return cl


def _drain(client, **kw):
    """Events the round produced, and the exception it ended with."""
    events, raised = [], None
    try:
        for event in client.stream_message(
                "sys", [{"role": "user", "content": "hi"}], max_tokens=50,
                **kw):
            events.append(event)
    except Exception as exc:                 # the round ends with it
        raised = exc
    return events, raised


def _notices(events):
    return "".join(getattr(e, "text", "") or "" for e in events
                   if getattr(e, "type", "") == "notice")


def test_a_half_hour_block_is_asked_once_and_not_again(client):
    calls = []

    def _create(**kw):
        calls.append(kw.get("model"))
        raise _rate_limited("1800")

    client.client.chat.completions.create = _create
    events, raised = _drain(client)

    assert len(calls) == 1, (
        f"the round asked again inside a 30-minute block ({len(calls)} "
        f"requests); each retry is a fresh violation that doubles it")
    assert raised is not None, "and the round does end rather than hang"


def test_the_round_says_how_long_and_why_it_stopped(client):
    def _create(**kw):
        raise _rate_limited("1800")

    client.client.chat.completions.create = _create
    note = _notices(_drain(client)[0])

    assert "30 min" in note, f"the wait is named: {note!r}"
    assert "1800" in note, "in seconds too, for a log"
    assert "lengthens the block" in note, (
        "and why asking again is not the answer")


def test_a_short_block_is_still_retried_and_then_succeeds(client):
    """The guard must not have turned every 429 into a dead round. A wait
    we will sit out is exactly what the header reading was built for."""
    calls = []

    def _create(**kw):
        calls.append(kw.get("model"))
        if len(calls) == 1:
            raise _rate_limited("30")
        return _Stream([_Chunk([_Choice(_Delta("Hallo."), finish="stop")],
                               usage=_Usage())])

    client.client.chat.completions.create = _create
    events, raised = _drain(client)
    text = "".join(getattr(e, "text", "") or "" for e in events
                   if getattr(e, "type", "") == "text_delta")

    assert len(calls) == 2, "a 30 s wait is sat out and the round retried"
    assert text == "Hallo.", "and the answer arrives"
    assert raised is None


def test_an_error_without_a_header_still_uses_the_ladder(client):
    """Absent is not long. A transient error nobody stated a wait for is
    the case the ladder was measured for, and it must still retry."""
    calls = []

    def _create(**kw):
        calls.append(kw.get("model"))
        if len(calls) == 1:
            openai = pytest.importorskip("openai")
            try:
                import httpx2 as transport
            except ImportError:              # pragma: no cover
                transport = pytest.importorskip("httpx")
            request = transport.Request("POST", "https://example.invalid/v1")
            raise openai.APIStatusError(
                "bad gateway",
                response=transport.Response(502, request=request), body=None)
        return _Stream([_Chunk([_Choice(_Delta("Hallo."), finish="stop")],
                               usage=_Usage())])

    client.client.chat.completions.create = _create
    events, raised = _drain(client)

    assert len(calls) == 2, "the ladder still rides out an unstated outage"
    assert raised is None
