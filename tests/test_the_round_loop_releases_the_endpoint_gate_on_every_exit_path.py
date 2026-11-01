"""The round loop releases the endpoint gate on EVERY exit path.

The gate (``delfin.agent.endpoint_gate``) limits concurrent model requests
host-wide over flock'd slot files. A slot held longer than one request is a
host-wide stall: every other session of the user queues behind it. So the
wiring in ``OpenAIClient.stream_message`` must release on success, on retry,
on cancel, on exception, and when the consumer stops iterating early -- and
never hold two rounds' slots at once.

No live endpoint: ``chat.completions.create`` is faked. The gate itself is
real, over a private runtime dir (the gate tests' fixture shape); the limit
is 1 so a leak shows up immediately as the next request blocking.
"""

from __future__ import annotations

import os
import threading
import time

import pytest

from delfin.agent import endpoint_gate as G
from delfin.agent import model_capabilities as mc


URL = "https://fake-endpoint.example/v1"


class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls
        self.reasoning_content = None


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
    """Iterable with a real close(); a consumer that stops iterating early
    leaves the generator suspended in the for-loop -- the path
    GeneratorExit exercises."""

    def __init__(self, chunks):
        self._chunks = list(chunks)
        self.closed = False

    def __iter__(self):
        return self

    def __next__(self):
        if not self._chunks:
            raise StopIteration
        return self._chunks.pop(0)

    def close(self):
        self.closed = True


class _FakeRegistry:
    def discover_all(self):
        return []

    def discover_resources(self):
        return []

    def discover_prompts(self):
        return []


class APITimeoutError(Exception):
    """Class-name transient (matches _is_transient_api_error's hints)."""


class RateLimited(Exception):
    status_code = 429


def _make_client(monkeypatch, tmp_path, model="qwen2.5-coder:7b"):
    """An OpenAIClient against the fake endpoint, gate ENABLED over a
    private runtime dir (shared by every client built in one test)."""
    rt = tmp_path / "rt"
    rt.mkdir(exist_ok=True)
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(rt))
    monkeypatch.setenv("DELFIN_ENDPOINT_MAX_CONCURRENCY", "1")
    monkeypatch.delenv("DELFIN_ENDPOINT_PER_MODEL_LIMIT", raising=False)
    monkeypatch.setattr("delfin.user_settings.load_settings", lambda: {})
    import delfin.agent.service.rundir as rundir
    gate_root = rundir.ensure_runtime_dir() / "endpoint_gates"
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(
        "delfin.agent.mcp_client.get_registry",
        lambda *a, **k: _FakeRegistry())
    import delfin.agent.api_client as A
    monkeypatch.setattr(A.time, "sleep", lambda *a, **k: None)
    return A, A.create_client(
        backend="api", provider="ollama", model=model,
        api_key=URL, cwd=str(tmp_path)), gate_root


def _drive(client):
    return [ev.text for ev in client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100)
        if ev.type in ("text_delta", "notice") and ev.text]


def _held_slots(gate_dir):
    """Number of slot files in ONE endpoint's gate dir currently flocked."""
    count = 0
    if not gate_dir.is_dir():
        return 0
    for f in gate_dir.iterdir():
        if not f.name.startswith("slot-") or not f.is_file():
            continue
        fd = os.open(str(f), os.O_RDWR | os.O_CREAT, 0o600)
        try:
            import fcntl
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            count += 1     # locked: someone holds it
        else:
            import fcntl
            fcntl.flock(fd, fcntl.LOCK_UN)
        finally:
            os.close(fd)
    return count


def _held_slots_for_url(gate_root, base_url, model="qwen2.5-coder:7b"):
    """Slots files carrying a flock for one model's round: the per-model
    adaptive gate AND the endpoint-wide total-cap gate (both are held while
    a request runs once per-model limits are the default)."""
    g = G.EndpointGate(base_url, model)
    held = _held_slots(gate_root / g._model_key)
    if g._cap_key:
        held += _held_slots(gate_root / g._cap_key)
    return held


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_slot_released_after_success(monkeypatch, tmp_path):
    A, client, gate_root = _make_client(monkeypatch, tmp_path)
    stream = _Stream([_Chunk([_Choice(_Delta(content="hello"), finish="stop")],
                             usage=_Usage())])
    held_during = {}

    def _create(**kwargs):
        # Positive control: while the request runs, the slots ARE held --
        # two of them: one per-model slot and one endpoint-wide cap slot.
        held_during["n"] = _held_slots_for_url(gate_root, client._base_url)
        return stream

    client.client.chat.completions.create = _create
    out = _drive(client)
    assert "hello" in "".join(out)
    assert held_during.get("n") == 2, (
        "the round loop must hold the model slot and the endpoint-cap "
        "slot while the request runs")
    assert G.EndpointGate(client._base_url).enabled
    assert _held_slots_for_url(gate_root, client._base_url) == 0


def test_slot_released_after_exception(monkeypatch, tmp_path):
    A, client, gate_root = _make_client(monkeypatch, tmp_path)

    def _create(**kwargs):
        raise RateLimited("429 too many requests")

    client.client.chat.completions.create = _create
    # RateLimited is transient -> the ladder retries to exhaustion and then
    # re-raises out of the round loop. Every attempt must have released.
    with pytest.raises(Exception):
        list(client.stream_message(
            "sys", [{"role": "user", "content": "go"}], max_tokens=100))
    assert _held_slots_for_url(gate_root, client._base_url) == 0


def test_slot_released_after_retry(monkeypatch, tmp_path):
    A, client, gate_root = _make_client(monkeypatch, tmp_path)
    calls = {"n": 0}

    def _create(**kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise APITimeoutError("upstream timed out")
        return _Stream([_Chunk(
            [_Choice(_Delta(content="second"), finish="stop")],
            usage=_Usage())])

    client.client.chat.completions.create = _create
    out = _drive(client)
    assert calls["n"] == 2
    assert "second" in "".join(out)
    assert _held_slots_for_url(gate_root, client._base_url) == 0


def test_slot_released_when_closed_at_wait_notice_yield(monkeypatch, tmp_path):
    """A consumer that closes the generator right after the round's
    waiting-for-a-slot notice yield must not leave the slot held. That yield
    sits before any stream-iteration try in the round, so it is the one path
    a naive release can miss."""
    A, client, gate_root = _make_client(monkeypatch, tmp_path)
    # Hold a slot long enough for acquire() to poll and publish a wait_notice,
    # then release it so acquire() completes and the wiring can yield that
    # notice. endpoint_gate sleeps on its own clock, so acquire really waits.
    # Per-model limits are the default, so the round waits on the MODEL key.
    gate = G.EndpointGate(client._base_url, "qwen2.5-coder:7b")
    dirpath = gate_root / gate._model_key
    dirpath.mkdir(parents=True, exist_ok=True)
    held_ev = threading.Event()

    def _hold_then_release():
        fd = os.open(str(dirpath / "slot-0"), os.O_RDWR | os.O_CREAT, 0o600)
        import fcntl
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        held_ev.set()                        # slot is now contended
        # Wait on a real timer, NOT time.sleep: _make_client no-ops the global
        # time.sleep, so the slot must be held via an Event timeout to let the
        # round's acquire() poll and publish its wait_notice.
        threading.Event().wait(0.3)
        fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)

    holder = threading.Thread(target=_hold_then_release)
    holder.start()
    assert held_ev.wait(timeout=5), "background slot holder never started"

    def _create(**kwargs):
        return _Stream([_Chunk(
            [_Choice(_Delta(content="ok"), finish="stop")], usage=_Usage())])

    client.client.chat.completions.create = _create
    gen = client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100)
    # Drive until the waiting notice is actually yielded, then close right
    # there -- the yield is inside the gate try, so closing releases the slot.
    noticed = False
    for _ in range(60):
        ev = next(gen)
        if ev.type == "notice":
            noticed = True
            break
    gen.close()
    holder.join(timeout=5)
    assert noticed, "the round did not publish a waiting-for-a-slot notice"
    assert _held_slots_for_url(gate_root, client._base_url) == 0


def test_slot_released_when_consumer_stops_early(monkeypatch, tmp_path):
    A, client, gate_root = _make_client(monkeypatch, tmp_path)
    # A stream with many chunks: the consumer takes a couple of events and
    # walks away, closing the generator mid-stream.
    stream = _Stream([
        _Chunk([_Choice(_Delta(content="a"), finish=None)], usage=_Usage())
        for _ in range(50)])

    def _create(**kwargs):
        return stream

    client.client.chat.completions.create = _create
    gen = client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100)
    for i, ev in enumerate(gen):
        if i >= 2:
            break                       # abandon mid-stream
    del gen                              # generator closed -> GeneratorExit
    assert _held_slots_for_url(gate_root, client._base_url) == 0


def test_cancel_releases_the_slot(monkeypatch, tmp_path):
    A, client, gate_root = _make_client(monkeypatch, tmp_path)

    def _create(**kwargs):
        raise RateLimited("429 too many requests")

    client.client.chat.completions.create = _create
    # The stop path: signal_stop (cuts the in-flight transport) before driving,
    # plus the engine's stop probe (the flag the loop checks -- real stopping
    # sets both), so the round error arrives "while stopped" and the loop
    # returns cleanly instead of retrying/re-raising. signal_stop REPLACES
    # self.client with a fresh SDK client (OpenAIClient.signal_stop
    # in api_client.py), so the fake is re-attached to the new one afterwards.
    client.signal_stop()
    client.should_stop = lambda: True
    client.client.chat.completions.create = _create
    events = list(client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100))
    assert any(getattr(ev, "stop_reason", "") == "stopped"
               for ev in events if ev.type == "message_delta")
    assert _held_slots_for_url(gate_root, client._base_url) == 0


def test_two_concurrent_sessions_never_exceed_the_limit(
        monkeypatch, tmp_path):
    """Two real stream_message rounds in flight at once, limit 1: the fake
    endpoint counts in-flight create() calls; the gate must keep the
    observed maximum at 1 and both sessions must still get their answer."""
    A, client1, gate_root = _make_client(monkeypatch, tmp_path)
    _, client2, _ = _make_client(monkeypatch, tmp_path, model="qwen2.5-coder:7b")

    lock = threading.Lock()
    stats = {"in_flight": 0, "max": 0}

    def _create(**kwargs):
        with lock:
            stats["in_flight"] += 1
            stats["max"] = max(stats["max"], stats["in_flight"])
        try:
            time.sleep(0.3)
            return _Stream([_Chunk(
                [_Choice(_Delta(content="answer"), finish="stop")],
                usage=_Usage())])
        finally:
            with lock:
                stats["in_flight"] -= 1

    client1.client.chat.completions.create = _create
    client2.client.chat.completions.create = _create
    results = {}

    def _run(name, cl):
        results[name] = "".join(_drive(cl))

    t1 = threading.Thread(target=_run, args=("one", client1))
    t2 = threading.Thread(target=_run, args=("two", client2))
    t1.start()
    t2.start()
    t1.join(timeout=60)
    t2.join(timeout=60)
    assert "answer" in results["one"]
    assert "answer" in results["two"]
    assert stats["max"] == 1, stats
    assert _held_slots_for_url(gate_root, client1._base_url) == 0


def test_the_slot_is_free_during_a_retry_backoff(monkeypatch, tmp_path):
    """A round that waits out a backoff before its retry holds no slot while
    it sleeps: another session may use the endpoint in that time. The
    congestion is reported once, not again by the round's finally."""
    A, client, gate_root = _make_client(monkeypatch, tmp_path)
    held_during_sleep = []
    # Only the backoff before the retry (>= 1 s, _STREAM_RETRY delay) is
    # measured: the short ticks while a round waits for its first byte run
    # with the slot held on purpose.
    def _sleep(seconds=0, *a, **k):
        if float(seconds or 0) >= 1.0:
            held_during_sleep.append(_held_slots_for_url(gate_root, client._base_url))
    monkeypatch.setattr(A.time, "sleep", _sleep)
    reports = []
    real_report = G.report_result
    monkeypatch.setattr(G, "report_result", lambda gate, **kw: (
        reports.append(kw.get("exception")), real_report(gate, **kw))[1])
    calls = {"n": 0}

    def _create(**kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise APITimeoutError("upstream timed out")
        return _Stream([_Chunk(
            [_Choice(_Delta(content="second"), finish="stop")],
            usage=_Usage())])

    client.client.chat.completions.create = _create
    out = _drive(client)
    assert "second" in "".join(out)
    assert held_during_sleep and all(n == 0 for n in held_during_sleep), held_during_sleep
    timeouts = [e for e in reports if isinstance(e, APITimeoutError)]
    assert len(timeouts) == 1, reports
    assert _held_slots_for_url(gate_root, client._base_url) == 0


# -- the wait the user reads -------------------------------------------------

def test_a_slot_wait_retry_is_named_a_local_wait():
    """A round that spent its deadline waiting for a LOCAL slot -- not the
    endpoint at all -- must never be reported as "No byte from the
    endpoint". Even at 400 s (past the zero-byte cut) it is a local wait,
    and a real endpoint stall still says "No byte from the endpoint"."""
    import delfin.agent.api_client as A
    long = 400.0                       # >= _ZERO_BYTE_CUT_S (300)
    slot = A._retry_notice(
        "SlotWaitTimeout", long, False, 2, 3, 5.0,
        slot_wait=True, base_url="https://kit.example.org")
    assert "local slot" in slot and "No byte from the endpoint" not in slot
    # A genuine endpoint stall keeps its existing message.
    stall = A._retry_notice("TimeoutError", long, False, 2, 3, 5.0)
    assert "No byte from the endpoint" in stall
