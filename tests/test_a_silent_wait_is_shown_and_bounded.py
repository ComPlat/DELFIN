"""A turn that hangs on the endpoint's first byte is the loudest silence
in DELFIN: the user sees nothing for minutes while the request is still
inside its (600 s) timeout. Field report 2026-09-27: `delfin-agent run
--provider kit --model kit.glm-5.3` hung in about half the turns, once a
full 700 s, the main thread parked in httpcore's
_receive_response_headers inside stream_message.

These tests pin the contract for that wait:

- the deadline really bounds the wait for the FIRST header (not just
  between chunks), and a stall is abandoned on that schedule;
- the wait is not silent: a "waiting for the model … Ns" indicator is
  emitted while no byte has arrived, so the turn looks alive;
- after the deadline the user reads a clear message and the round is
  retried (the existing transient-retry path), not dropped silently.
"""
import threading
import time
import types

import pytest

from delfin.agent import api_client as ac


def _client(stub):
    c = ac.OpenAIClient.__new__(ac.OpenAIClient)
    c.client = stub
    c.model = "kit.glm-5.3"
    c._provider = "kit"
    c._base_url = "https://ki-toolbox.scc.kit.edu/api/v1"
    c._api_key = "x"
    c.effort = ""
    c._permissions = None
    c.on_model_switched = None
    c._steer_lock = threading.Lock()
    c._steer_queue = []
    c._steer_msgs = []
    c._run_notes = []
    c._stop_flag = False
    return c


class _HangThenAnswer:
    """Stub provider: the first N create() calls accept the connection and
    never send a header until ``stall_s`` has passed; then a normal answer.
    That is exactly the observed KIT stall, compressed into test time."""

    def __init__(self, stall_s: float, hangs: int):
        self.stall_s = stall_s
        self.hangs = hangs
        self.calls = 0
        self.chat = types.SimpleNamespace(
            completions=types.SimpleNamespace(create=self._create))

    def _create(self, **kw):
        self.calls += 1
        if self.calls <= self.hangs:
            # Accept the connection, emit no header for stall_s, then raise
            # like httpx does when the read timeout fires mid-headers.
            time.sleep(self.stall_s)
            raise TimeoutError("timed out waiting for response headers")
        delta = types.SimpleNamespace(content="done", tool_calls=None)
        chunk = types.SimpleNamespace(
            usage=types.SimpleNamespace(prompt_tokens=3, completion_tokens=2),
            choices=[types.SimpleNamespace(delta=delta,
                                           finish_reason="stop")])

        class _Stream:
            def __iter__(self):
                return iter((chunk,))

            def close(self):
                pass

        return _Stream()


def test_a_header_stall_is_abandoned_on_the_deadline_and_says_so(monkeypatch):
    """The stall must end at DELFIN_REQUEST_TIMEOUT_S (compressed here to
    1 s), not hang forever; the user must read a retry notice; the round
    must then be retried and answer."""
    monkeypatch.setenv("DELFIN_REQUEST_TIMEOUT_S", "1")
    stub = _HangThenAnswer(stall_s=30.0, hangs=1)
    client = _client(stub)

    events = list(client.stream_message(
        messages=[{"role": "user", "content": "hi"}],
        system="You are a test.", max_tokens=64))

    texts = [getattr(e, "text", "") for e in events
             if getattr(e, "type", "") in ("text_delta", "notice")]
    notices = [t for t in texts if "retry" in t.lower() or "no byte" in t.lower()]
    assert any(t == "done" or "done" in t for t in texts), (
        "the retried round never produced its answer")
    assert notices, (
        "a stalled first header ended without any message the user could read")


def test_the_wait_for_the_first_byte_is_not_silent(monkeypatch):
    """While no byte has arrived, a waiting indicator must be emitted —
    'waiting for the model … Ns' — so a legitimately slow endpoint is told
    apart from a dead one by the person watching the screen."""
    monkeypatch.setenv("DELFIN_REQUEST_TIMEOUT_S", "3")
    stub = _HangThenAnswer(stall_s=2.0, hangs=1)
    client = _client(stub)

    saw_waiting = []
    for ev in client.stream_message(
            messages=[{"role": "user", "content": "hi"}],
            system="You are a test.", max_tokens=64):
        t = getattr(ev, "type", "")
        if t == "waiting":
            saw_waiting.append(getattr(ev, "elapsed_s", None))
        if t == "text_delta" and saw_waiting:
            break
    assert saw_waiting, (
        "the turn waited seconds for the first byte and said nothing")
