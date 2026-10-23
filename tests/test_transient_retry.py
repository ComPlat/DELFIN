"""Transient-API-error retry on the streaming call.

A flaky shared-proxy hiccup (timeout / 5xx / rate-limit) on a long KIT run
should be retried with backoff instead of killing the whole turn — and the
retry must NOT duplicate already-streamed output, nor fire on a deterministic
error (400/401) that would just fail again. No live model: create() is faked.
"""

from __future__ import annotations

import pytest

from delfin.agent import mcp_client as M
from delfin.agent import model_capabilities as mc


# --- minimal OpenAI-stream fakes -------------------------------------------
class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


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


def _final_text(txt="done"):
    return _Stream([_Chunk([_Choice(_Delta(content=txt), finish="stop")],
                           usage=_Usage())])


class _FakeRegistry:
    def discover_all(self):
        return []

    def discover_resources(self):
        return []

    def discover_prompts(self):
        return []


class APITimeoutError(Exception):
    """Class name matches the transient classifier's 'timeout' hint."""


class BadRequestError(Exception):
    def __init__(self, msg):
        super().__init__(msg)
        self.status_code = 400


@pytest.fixture
def client(monkeypatch, tmp_path):
    # provider="ollama" needs no API key (the retry path is provider-agnostic),
    # so this runs in CI's key-free, fully-mocked environment. Using "kit" here
    # would error at construction when KIT_TOOLBOX_API_KEY is absent.
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _FakeRegistry())
    import delfin.agent.api_client as A
    monkeypatch.setattr(A.time, "sleep", lambda *a, **k: None)  # skip real backoff
    return A.create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(tmp_path))


def _drive(client):
    """Everything the USER sees. The harness's own sentences -- a retry
    banner, a stop, a cost ceiling -- carry their own event type so they
    stay out of the model's answer, but they are still shown, so a test
    about what the user is told has to collect both."""
    return [ev.text for ev in client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100)
        if ev.type in ("text_delta", "notice") and ev.text]


def test_transient_error_is_retried_then_succeeds(client):
    calls = {"n": 0}

    def _create(**kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise APITimeoutError("upstream timed out")
        return _final_text("done")

    client.client.chat.completions.create = _create
    joined = "".join(_drive(client))
    assert calls["n"] == 2                  # retried exactly once
    assert joined.count("done") == 1        # answer present once — no duplication
    assert "retrying" in joined.lower()     # the retry was surfaced to the user


def test_nontransient_error_is_not_retried(client):
    calls = {"n": 0}

    def _create(**kwargs):
        calls["n"] += 1
        raise BadRequestError("bad request")

    client.client.chat.completions.create = _create
    with pytest.raises(BadRequestError):
        _drive(client)
    assert calls["n"] == 1                   # a deterministic 400 is not retried


# --- KIT/litellm/vLLM proxy 'Extra data' hiccup: a 400 that IS transient -----

class _Status400(Exception):
    def __init__(self, msg):
        super().__init__(msg)
        self.status_code = 400


def test_proxy_extra_data_400_is_transient():
    from delfin.agent.api_client import _is_transient_api_error
    e = _Status400(
        "Error code: 400 - {'detail': 'litellm.BadRequestError: "
        'Hosted_vllmException - {"error":{"message":"Extra data: line 1 '
        'column 505 (char 504)","type":"BadRequestError","code":400}}\'}')
    assert _is_transient_api_error(e) is True


def test_genuine_bad_request_400_is_not_transient():
    from delfin.agent.api_client import _is_transient_api_error
    # real client errors must NOT be retried (they'd just fail again)
    assert _is_transient_api_error(_Status400("model not found")) is False
    assert _is_transient_api_error(_Status400("context length exceeded")) is False
    # 'Extra data' WITHOUT the proxy markers is too generic to retry
    assert _is_transient_api_error(_Status400("Extra data in user field")) is False


# ---------------------------------------------------------------------------
# A gateway that runs out of its own resources reports a 400
# ---------------------------------------------------------------------------


def _err(msg, status=400):
    exc = RuntimeError(msg)
    exc.status_code = status
    return exc


def test_gateway_database_exhaustion_is_transient():
    """Field case 2026-07-30: a run died mid-task on a 400 whose body was a
    Postgres connection-slot error from the gateway — its own capacity, not
    a malformed request."""
    from delfin.agent.api_client import _is_transient_api_error
    body = ("Error code: 400 - {'detail': '(psycopg.OperationalError) "
            "connection failed: FATAL: remaining connection slots are "
            "reserved for roles with the SUPERUSER attribute "
            "(Background on this error at: https://sqlalche.me/e/20/e3q8)'}")
    assert _is_transient_api_error(_err(body))


def test_pool_and_availability_messages_are_transient():
    from delfin.agent.api_client import _is_transient_api_error
    for body in ("Error code: 400 - too many connections for role",
                 "Error code: 400 - connection pool exhausted",
                 "Error code: 400 - service temporarily unavailable"):
        assert _is_transient_api_error(_err(body)), body


def test_real_client_errors_are_still_not_retried():
    """The markers must never rescue a genuinely bad request."""
    from delfin.agent.api_client import _is_transient_api_error
    for body in (
        "Error code: 400 - {'error': {'message': 'model not found'}}",
        "Error code: 400 - maximum context length is 131072 tokens",
        "Error code: 400 - Invalid value for 'temperature'",
        "Error code: 401 - invalid api key",
        "Error code: 404 - no such model",
    ):
        assert not _is_transient_api_error(_err(body)), body


def test_the_retry_banner_is_not_the_model_speaking(client):
    """The live incident, at its source. Two office tasks were scored as a
    model that answered badly, at quality 35, and the entire recorded
    answer was three of these banners -- because the harness's own
    sentence arrived on the same event type as the model's words.

    Asserted on the event stream rather than on the module's source: what
    matters is which type the banner actually carries when a request
    fails, not which literal appears near which literal.
    """
    calls = {"n": 0}

    def _create(**kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise APITimeoutError("upstream timed out")
        return _final_text("done")

    client.client.chat.completions.create = _create
    events = list(client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100))

    banners = [e for e in events if "retrying" in (e.text or "").lower()]
    assert banners, "the retry was not surfaced at all"
    assert all(e.type == "notice" for e in banners), \
        [e.type for e in banners]
    # ... and the model's own word still comes through as the answer.
    answer = "".join(e.text for e in events if e.type == "text_delta")
    assert answer.strip() == "done"
