"""Tests for _augment_api_error (server error body recording).

Input:  a litellm-style 4xx/5xx exception (status_code attribute,
optional response.text / message holding the server's error body).
Output: the same exception object with args rewritten to carry the
server's body cut to 1000 characters and scrubbed of credentials; for
a context-length 400 additionally a plain-words explanation with the
measured prompt token count. Non-4xx/5xx exceptions pass through
unchanged. Rejected alternative: raising a new exception type -- every
existing classifier reads status_code and str(exc), so the class and
the idempotency guard (exc._delfin_augmented) must stay stable.
"""

import json

import pytest

from delfin.agent.api_client import _augment_api_error

# A fake key long enough (>=16 chars, >=2 char classes) to pass the
# _plausible_secret heuristic used by _redact_tool_result. Built from
# parts so no secret-shaped literal sits in this file.
FAKE_KEY = "-".join(["sk", "SUPERSECRET"]) + "-1234567890ABCDEF"


class FakeApiError(Exception):
    """Mimics the litellm/APIConnectionError shape the client raises."""

    def __init__(self, status_code, body_text):
        super().__init__(f"Error code: {status_code} - {body_text}")
        self.status_code = status_code
        self.response = type("R", (), {"text": body_text})()
        self.message = body_text


def _body(msg="The model is overloaded", **extra):
    payload = {"error": {"message": msg}}
    payload.update(extra)
    return json.dumps(payload)


def test_400_body_is_recorded_and_key_masked():
    out = _augment_api_error(FakeApiError(400, _body("bad request",
                                                     api_key=FAKE_KEY)))
    text = str(out)
    assert "bad request" in text
    assert "[HTTP 400 from the model endpoint:" in text
    assert FAKE_KEY not in text
    assert "[redacted:" in text
    assert out.status_code == 400


def test_key_masked_when_body_already_in_str_exc():
    # Regression: str(exc) already embeds the raw body, so the dedupe rule
    # left args untouched and the raw key reached the log. The final text
    # must be the scrubbed one in every path.
    body = _body("bad request", api_key=FAKE_KEY)
    out = _augment_api_error(FakeApiError(400, body))
    assert FAKE_KEY not in str(out)
    assert FAKE_KEY not in "".join(str(a) for a in out.args)


def test_body_cut_to_1000_chars():
    out = _augment_api_error(FakeApiError(500, _body("x" * 5000)))
    text = str(out)
    marker = "[HTTP 500 from the model endpoint: "
    idx = text.find(marker)
    assert idx != -1
    tail = text[idx + len(marker):].rstrip("]")
    assert len(tail) <= 1000
    assert "x" * 900 in tail


def test_idempotent_second_call_unchanged():
    exc = FakeApiError(400, _body())
    first = str(_augment_api_error(exc))
    second = str(_augment_api_error(exc))
    assert first == second
    assert "The model is overloaded" in second


def test_non_4xx_5xx_and_bodyless_unchanged():
    ok = FakeApiError(200, _body())
    assert _augment_api_error(ok) is ok
    assert "[HTTP" not in str(ok)

    no_status = RuntimeError("boom")
    assert _augment_api_error(no_status) is no_status
    assert "[HTTP" not in str(no_status)


def test_context_length_400_gets_wording_and_token_count():
    body = json.dumps(
        {"error": {"message": "This model's maximum context length is "
                              "4096 tokens, however the messages are "
                              "5000 tokens long."}})
    msgs = [{"role": "user", "content": "hello world " * 200}]
    out = _augment_api_error(FakeApiError(400, body),
                             request_kwargs={"messages": msgs})
    text = str(out)
    assert "context is too long" in text
    assert "measured about" in text
    assert "tokens" in text


def test_context_length_400_without_kwargs_still_worded():
    body = json.dumps({"error": {"message": "prompt is too long"}})
    out = _augment_api_error(FakeApiError(400, body))
    text = str(out)
    assert "context is too long" in text
    assert "measured about 0 tokens" in text


def test_404_and_502_body_recorded():
    for status, msg in ((404, "model not found"), (502, "bad gateway")):
        out = _augment_api_error(FakeApiError(status, _body(msg)))
        text = str(out)
        assert f"[HTTP {status} from the model endpoint:" in text
        assert msg in text


def test_base_error_text_is_bounded():
    # The provider's own str(exc) may already embed a huge body (litellm
    # repeats the request body in the 400 message). The base must be capped
    # too, or that embedded copy floods the log even though the standalone
    # body is cut to 1000.
    huge = "h" * 60000
    out = _augment_api_error(FakeApiError(400, huge))
    text = str(out)
    assert len(text) <= 2000 + 1000 + 200  # base cap + body cap + marker slack
    assert "h" * 4000 not in text


def _kit_400(detail: str) -> FakeApiError:
    """A KIT-proxy-shaped 400: its body is a JSON "detail" string."""
    return FakeApiError(400, json.dumps({"detail": detail}))


def test_run_once_json_error_carries_body_without_key():
    # cli._run_once builds the JSON result's "error" field as str(exc) and
    # the session log records the same string. Drive that contract: a stub
    # engine whose stream_response raises the augmented 400 - the same
    # exception shape OpenAIClient.stream_message produces at the bare
    # raise - and assert result["error"] keeps the server body, loses the
    # key.
    from delfin.agent.cli import _run_once

    body = _kit_400("Model 'glm-5.3-flash' is temporarily unavailable "
                    + "(api_key=" + FAKE_KEY + ")")
    exc = _augment_api_error(body)

    class _StubEngine:
        last_turn_stop_reason = ""
        token_usage = {}

        def stream_response(self, **kw):
            raise exc

    result = _run_once(_StubEngine(), "hi", max_tokens=8)
    assert "glm-5.3-flash" in result["error"]
    assert "temporarily unavailable" in result["error"]
    assert FAKE_KEY not in result["error"]


def test_stream_400_body_reaches_the_raised_error():
    # End-to-end through the REAL raise site: OpenAIClient.stream_message
    # with the one HTTP call scripted to return KIT's 400 body. The
    # exception the generator raises must carry the server message and not
    # the key - this is exactly what cli._run_once turns into the JSON
    # result's "error".
    import threading

    from delfin.agent import api_client as ac

    body = _kit_400("Model 'glm-5.3-flash' is temporarily unavailable "
                    + "(api_key=" + FAKE_KEY + ")")

    class _RaisingStub:
        class chat:
            class completions:
                @staticmethod
                def create(**kw):
                    # _kit_400 already returns a FakeApiError whose body is
                    # json.dumps({"detail": ...}); raise it verbatim so the
                    # augment path sees the real KIT 400 shape.
                    raise body

    c = ac.OpenAIClient.__new__(ac.OpenAIClient)
    c.client = _RaisingStub()
    c.model = "kit.glm-5.3"
    c._provider = "kit"
    c._base_url = "https://ki-toolbox.scc.kit.edu/api/v1"
    c._api_key = "x"
    c.effort = ""
    c._permissions = None
    c.on_model_switched = None
    c._steer_lock = threading.Lock()
    c._steer_queue = []
    c._run_notes = []
    c._stop_flag = False
    c.metrics_session = None

    raised = None
    try:
        for _ev in c.stream_message(
                messages=[{"role": "user", "content": "hi"}],
                system="t", max_tokens=32):
            pass
    except Exception as exc:  # noqa: BLE001 - we assert on this exact raise
        raised = exc
    assert raised is not None, "stream_message should have raised"
    text = str(raised)
    assert "glm-5.3-flash" in text
    assert "temporarily unavailable" in text
    assert FAKE_KEY not in text


def test_fallback_non_streaming_400_body_recorded():
    # The stream-unsupported 400 is retried ONCE without streaming; if that
    # non-streaming retry hits a 400 of its own, the body must still reach
    # the log. The first create (stream=True) fails with litellm's
    # __aiter__/ModelResponse shape; the fallback create (stream=False)
    # returns KIT's 400 body.
    import threading

    from delfin.agent import api_client as ac

    body = json.dumps({"detail": "Model 'glm-5.3-flash' is temporarily "
                       + "unavailable (api_key=" + FAKE_KEY + ")"})

    class _FakeStreamUnsupported(Exception):
        def __init__(self):
            super().__init__(
                "'async for' requires an object with __aiter__ method, "
                "got ModelResponse")

    raw_kw = {}

    class _FallbackStub:
        class chat:
            class completions:
                @staticmethod
                def create(**kw):
                    raw_kw.update(kw)
                    if kw.get("stream", True) is True:
                        raise _FakeStreamUnsupported()
                    raise FakeApiError(400, body)

    c = ac.OpenAIClient.__new__(ac.OpenAIClient)
    c.client = _FallbackStub()
    c.model = "kit.glm-5.3"
    c._provider = "kit"
    c._base_url = "https://ki-toolbox.scc.kit.edu/api/v1"
    c._api_key = "x"
    c.effort = ""
    c._permissions = None
    c.on_model_switched = None
    c._steer_lock = threading.Lock()
    c._steer_queue = []
    c._run_notes = []
    c._stop_flag = False
    c.metrics_session = None

    raised = None
    try:
        for _ev in c.stream_message(
                messages=[{"role": "user", "content": "hi"}],
                system="t", max_tokens=32):
            pass
    except Exception as exc:  # noqa: BLE001 - assert on this raise
        raised = exc
    assert raised is not None, "fallback create's 400 should have raised"
    # The fallback must actually have been taken (stream=False issued).
    assert any(kw.get("stream") is False for kw in [raw_kw])
    text = str(raised)
    assert "glm-5.3-flash" in text
    assert "temporarily unavailable" in text
    assert FAKE_KEY not in text


if __name__ == "__main__":
    pytest.main([__file__])
