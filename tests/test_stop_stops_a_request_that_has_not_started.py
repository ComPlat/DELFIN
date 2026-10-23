"""A user pressed Stop while the first turn on kit.glm-5.3 waited for
its first byte, watched nothing happen, and typed into a turn that was
still running (field report 2026-09-11: "ich hab stop gesendet jetzt
kann ich nicht mehr stoppen"). The stop flag is read between stream
events; a request that has not started delivers none. The OpenAI
client's stop now cuts the transport, the blocked read raises at once,
and the stream loop ends the turn as stopped instead of retrying.
"""

from __future__ import annotations

import inspect
import threading
import time

import pytest

from delfin.agent import api_client


class _BlockingSDK:
    """The slice of openai.OpenAI the client touches: create() blocks
    until close() is called, then raises like a cut connection."""

    def __init__(self):
        self.closed = threading.Event()
        outer = self

        class _Completions:
            def create(self, **kwargs):
                outer.closed.wait(timeout=20)
                raise ConnectionError("connection closed")

        class _Chat:
            completions = _Completions()

        self.chat = _Chat()
        self.base_url = "https://ki-toolbox.scc.kit.edu/api/v1"

    def close(self):
        self.closed.set()


def test_a_stop_cuts_the_transport_and_the_turn_ends_stopped(monkeypatch):
    sdk = _BlockingSDK()
    fresh = _BlockingSDK()
    import openai
    monkeypatch.setattr(openai, "OpenAI", lambda **kw: fresh)
    client = api_client.OpenAIClient.__new__(api_client.OpenAIClient)
    client.client = sdk
    client._sdk_kwargs = {"api_key": "x"}
    stopped = {"flag": False}
    client.should_stop = lambda: stopped["flag"]

    assert not sdk.closed.is_set()
    stopped["flag"] = True
    client.signal_stop()
    assert sdk.closed.is_set(), "the blocked read must be released at once"
    assert client.client is fresh, "the next turn needs a transport that is open"


def test_the_stream_loop_ends_stopped_before_it_would_retry():
    from pathlib import Path
    src = Path(api_client.__file__).read_text(encoding="utf-8")
    i = src.index("except Exception as _stream_exc:")
    window = src[i:i + 2200]
    assert "if self._stop_was_requested():" in window
    assert window.index("if self._stop_was_requested():") < window.index("_is_transient_api_error(")
    assert 'stop_reason="stopped"' in window


def test_the_base_client_stop_is_still_a_no_op():
    """Backends without a transport of their own keep the cooperative stop."""
    assert api_client._BaseClient.signal_stop(object()) is None


def test_the_sdk_does_not_retry_on_its_own():
    """Three silent 600 s waits before the stream loop could say a word."""
    src = inspect.getsource(api_client.OpenAIClient.__init__)
    assert 'kwargs["max_retries"] = 0' in src


def test_a_notice_is_not_the_models_output_in_the_dashboard():
    from pathlib import Path
    text = (Path(__file__).resolve().parents[1] / "delfin" / "dashboard" / "tab_agent.py").read_text(encoding="utf-8")
    i = text.index("def _on_notice(text):")
    body = text[i:text.index("def _on_token(text):", i)]
    assert "_mark_output()" not in body, "a notice must not stamp the first token"
    assert "_append_system_message(text.strip())" in body
    assert "on_notice=_on_notice," in text
