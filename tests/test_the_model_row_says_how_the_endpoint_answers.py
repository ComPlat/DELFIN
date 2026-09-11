"""A model picked in the dashboard was a name in a dropdown until the
first turn, and the first turn on a queued endpoint was a spinner for
ten minutes. Measured on the KIT deployment on 2026-09-11: kit.glm-5.3
took 459 s to its first token for a one-word question, kit.deepseek-v4-
flash 0.6 s. The dashboard now sends one word as soon as the engine
exists and says how long the endpoint took -- before the user waits.
"""

from __future__ import annotations

import pathlib

from delfin.agent.endpoint_probe import ProbeResult, describe, probe_first_token


class _Stream:
    def __init__(self, chunks, delay_s=0.0):
        self._chunks, self._delay = list(chunks), delay_s
        self.closed = False

    def __iter__(self):
        import time
        for c in self._chunks:
            if self._delay:
                time.sleep(self._delay)
            yield c

    def close(self):
        self.closed = True


class _SDK:
    """The slice of openai.OpenAI the probe touches, recording the call."""

    def __init__(self, stream=None, raise_=None):
        self.calls = []
        outer = self

        class _Completions:
            def create(self, **kwargs):
                outer.calls.append(kwargs)
                if raise_ is not None:
                    raise raise_
                return stream

        class _Chat:
            completions = _Completions()

        self.chat = _Chat()


def test_the_probe_stops_at_the_first_chunk_and_reports_the_wait():
    stream = _Stream(["a", "b", "c"], delay_s=0.05)
    sdk = _SDK(stream)
    res = probe_first_token(sdk, "kit.glm-5.3", timeout_s=5.0, reasoning_effort="low")
    assert res.answered and 0.04 <= res.seconds < 2.0
    assert stream.closed, "the probe must not read the answer, one chunk is the measurement"
    call = sdk.calls[0]
    assert call["model"] == "kit.glm-5.3" and call["stream"] is True
    assert call["max_tokens"] <= 8 and call["timeout"] == 5.0
    assert call["reasoning_effort"] == "low"
    assert call["messages"] == [{"role": "user", "content": "hi"}]


def test_a_model_that_does_not_reason_is_not_sent_the_parameter():
    sdk = _SDK(_Stream(["a"]))
    probe_first_token(sdk, "kit.deepseek-v4-flash", timeout_s=5.0)
    assert "reasoning_effort" not in sdk.calls[0]


def test_a_timeout_is_an_answer_too():
    class _Timeout(Exception):
        pass

    res = probe_first_token(_SDK(raise_=_Timeout("no bytes")), "kit.glm-5.3", timeout_s=1.0)
    assert not res.answered and res.error == "_Timeout"


def test_an_empty_stream_is_not_called_an_answer():
    res = probe_first_token(_SDK(_Stream([])), "m", timeout_s=1.0)
    assert not res.answered and "without a chunk" in res.error


def test_the_line_says_what_the_number_means():
    fast = describe(ProbeResult("kit.deepseek-v4-flash", 0.6, True))
    assert fast == "kit.deepseek-v4-flash answered a one-word probe in 0.6 s."
    slow = describe(ProbeResult("kit.glm-5.3", 47.0, True))
    assert slow.startswith("kit.glm-5.3 answered a one-word probe in 47 s.")
    assert "queueing" in slow
    none = describe(ProbeResult("kit.glm-5.3", 60.0, False, "APITimeoutError"),
                    first_token_budget_s=600.0)
    assert none.startswith("kit.glm-5.3: no answer to a one-word probe after 60 s (APITimeoutError).")
    assert "queueing" in none and "10 min" in none and "another model" in none


def test_the_dashboard_probes_when_the_engine_is_built():
    text = (pathlib.Path(__file__).resolve().parents[1]
            / "delfin" / "dashboard" / "tab_agent.py").read_text(encoding="utf-8")
    i = text.index("def _ensure_engine(")
    ensure = text[i:i + 20000]
    assert "_probe_endpoint_in_background(engine, provider, model)" in ensure
    j = text.index("def _probe_endpoint_in_background(")
    body = text[j:j + 3000]
    assert "probe_first_token(sdk, model" in body
    assert 'get("probe_endpoint", True)' in body, "off with agent.probe_endpoint: false"
    assert 'if state.get("engine") is not engine:' in body, "a stale number must not be posted"
    assert "daemon=True" in body
