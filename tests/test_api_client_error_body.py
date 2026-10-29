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


def test_augmented_text_survives_json_result_round_trip():
    # cli.py builds the JSON "error" field as str(exc) and the session log
    # records the same string; the augmented text must survive both
    # unchanged.
    out = _augment_api_error(FakeApiError(400, _body("bad request",
                                                     api_key=FAKE_KEY)))
    text = str(out)
    payload = json.loads(json.dumps({"error": text}))
    assert "bad request" in payload["error"]
    assert FAKE_KEY not in payload["error"]


if __name__ == "__main__":
    pytest.main([__file__])
