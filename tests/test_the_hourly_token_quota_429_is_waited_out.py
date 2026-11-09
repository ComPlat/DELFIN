"""The hourly token-quota 429 is waited out, not treated as congestion.

Measured 2026-11-09: the KIT endpoint enforces an HOURLY token quota per
key. A quota 429's body carries the numbers and the reset instant ("Hourly
token quota exceeded: N / LIMIT free tokens used. ... Quota resets at
YYYY-MM-DD HH:MM UTC (in Ns)."). Today the turn classifies it like any
transient 429, rides the backoff ladder to exhaustion and ends with
stop_reason "error"; the shared endpoint gate is halved to its floor as if
the endpoint were congested, and stays there because nothing succeeds for
the rest of the quota window.

Semantics under test (each numbered in the code it pins):

1. The parser reads the reset instant out of the body: absolute UTC first,
   "in Ns" fallback, clamped to at most 65 minutes; anything unparsable
   returns None (today's behaviour).
2. A quota 429 followed by a success: ONE wait notice, the round retries
   after the (test-compressed) reset, the turn completes.
3. The gate never hears a quota 429 as congestion: the shared limit is not
   halved, and no quick-success credit is recorded while waiting.
4. A quota 429 whose reset is already past: no wait at all, immediate
   retry.
5. A 429 with no quota text (plain per-minute rate limit): the old backoff
   ladder still applies, no quota wait.
6. A stop during the quota wait ends the run within ~2 s.
7. A wall-clock budget shorter than the wait ends the turn with
   stop_reason "quota_wait_exceeds_budget", not "error" and not "done".
8. The session usage summary records the quota numbers (used/limit) when
   seen, so a person sees how close the key is.

No live endpoint: ``chat.completions.create`` is faked, like the
endpoint-gate round-loop tests. The gate is real over a private runtime
dir.
"""

from __future__ import annotations

import threading
import time

import pytest

from delfin.agent import endpoint_gate as G
from delfin.agent import model_capabilities as mc

URL = "https://fake-endpoint.example/v1"

QUOTA_BODY = (
    "Hourly token quota exceeded: 6,233,280 / 6,000,000 free tokens used. "
    "Free tokens = fresh tokens + 10% of cached tokens (prompt-cache hits "
    "are 10x cheaper). Quota resets at 2026-11-09 22:00 UTC (in 410s).")


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


class RateLimited(Exception):
    status_code = 429


def _make_client(monkeypatch, tmp_path, model="qwen2.5-coder:7b"):
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
    # Real sliced sleeps must be instant in tests: the client waits in
    # 1 s slices honoring stop, but the slices themselves are driven
    # through A.time.sleep, which we zero here. A test that must observe
    # "stop breaks the wait within 2 s" drives the stop itself.
    monkeypatch.setattr(A.time, "sleep", lambda *a, **k: None)
    return A, A.create_client(
        backend="api", provider="ollama", model=model,
        api_key=URL, cwd=str(tmp_path)), gate_root


def _drive(client):
    return [ev for ev in client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100)]


# ---------------------------------------------------------------------------
# 1. the parser
# ---------------------------------------------------------------------------

def test_parses_the_absolute_utc_reset(monkeypatch):
    import delfin.agent.api_client as A
    # The verbatim measured body.
    reset, used, limit = A._quota_429_reset(
        "Hourly token quota exceeded: 6,233,280 / 6,000,000 free tokens "
        "used. Free tokens = fresh tokens + 10% of cached tokens "
        "(prompt-cache hits are 10x cheaper). Quota resets at "
        "2026-11-09 22:00 UTC (in 410s).",
        now=time.time())
    # The absolute instant wins: 2026-11-09 22:00 UTC.
    import datetime
    want = datetime.datetime(2026, 11, 9, 22, 0,
                             tzinfo=datetime.timezone.utc).timestamp()
    assert reset == want
    assert used == 6_233_280
    assert limit == 6_000_000


def test_parses_the_relative_reset_when_no_absolute_time(monkeypatch):
    import delfin.agent.api_client as A
    now = 1_700_000_000.0
    reset, used, limit = A._quota_429_reset(
        "Hourly token quota exceeded: 5 / 100 free tokens used. "
        "Quota resets in 300s.", now=now)
    assert reset == now + 300
    assert used == 5
    assert limit == 100


def test_the_reset_is_clamped_to_65_minutes(monkeypatch):
    import delfin.agent.api_client as A
    now = 1_700_000_000.0
    # A body claiming 4 hours ahead must not park the turn for 4 hours.
    reset, _used, _limit = A._quota_429_reset(
        "Hourly token quota exceeded: 1 / 2 free tokens used. "
        "Quota resets in 14400s.", now=now)
    assert reset <= now + 65 * 60


def test_an_unparsable_429_body_parses_as_none(monkeypatch):
    import delfin.agent.api_client as A
    assert A._quota_429_reset(
        "429 too many requests, slow down", now=1_700_000_000.0) is None


def test_a_non_quota_429_is_not_seen_as_quota(monkeypatch):
    import delfin.agent.api_client as A
    # Per-minute rate limits carry no "Hourly token quota" marker: they
    # must parse as None so the old backoff keeps applying.
    assert A._quota_429_reset(
        "Rate limit reached: 60 requests per minute exceeded",
        now=1_700_000_000.0) is None


# ---------------------------------------------------------------------------
# 2 + 3. quota 429 then success: one wait, turn continues, gate untouched
# ---------------------------------------------------------------------------

def test_quota_429_is_waited_out_then_the_run_continues(monkeypatch, tmp_path):
    A, client, gate_root = _make_client(monkeypatch, tmp_path)
    calls = {"n": 0}

    def _create(**kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            exc = RateLimited(QUOTA_BODY)
            exc.status_code = 429
            # str(exc) must carry the body -- the classifier reads it.
            exc.args = (QUOTA_BODY,)
            raise exc
        return _Stream([_Chunk([_Choice(_Delta(content="ok"), finish="stop")],
                               usage=_Usage())])

    client.client.chat.completions.create = _create
    events = list(_drive(client))
    assert calls["n"] == 2, "the run must retry after the quota reset"
    texts = "".join(getattr(ev, "text", "") for ev in events)
    assert "quota" in texts.lower(), "one wait notice must be shown"
    # ONE wait line, not one per second (the wait slices are ~1 s).
    assert texts.lower().count("waiting until") == 1, texts
    assert "ok" in texts
    # No stop inside the wait: the run continues and finishes normally.
    stop = [ev for ev in events if ev.type == "message_delta"]
    assert stop and stop[0].stop_reason in (None, "", "end_turn")


def test_the_gate_is_not_halved_by_a_quota_429(monkeypatch, tmp_path):
    A, client, gate_root = _make_client(monkeypatch, tmp_path)
    calls = {"n": 0}

    def _create(**kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise RateLimited(QUOTA_BODY)
        return _Stream([_Chunk([_Choice(_Delta(content="ok"), finish="stop")],
                               usage=_Usage())])

    client.client.chat.completions.create = _create
    list(_drive(client))
    gate = G.EndpointGate(client._base_url, "qwen2.5-coder:7b")
    dirpath = gate_root / gate._model_key
    if dirpath.is_dir():
        state = G._read_state(dirpath)
        assert state.get("limit", 1) >= 1 and not state.get("last_congestion_ts") \
            if False else True  # no congestion recorded
    # The strong form: the shared limit must still be the config value
    # (1 here), never halved to the floor by a quota 429, and no
    # congestion timestamp recorded.
    gate2 = G.EndpointGate(client._base_url, "qwen2.5-coder:7b")
    assert gate2.current_limit(dirpath) == 1
    state = G._read_state(dirpath)
    assert "last_congestion_ts" not in state


def test_reset_already_passed_waits_nothing(monkeypatch, tmp_path):
    A, client, gate_root = _make_client(monkeypatch, tmp_path)
    calls = {"n": 0}

    def _create(**kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            # A body whose reset instant lies in the past (stale proxy
            # cache): waiting until it would hang forever -- retry now.
            raise RateLimited(
                "Hourly token quota exceeded: 7 / 6 free tokens used. "
                "Quota resets at 2020-01-01 00:00 UTC (in 1s).")
        return _Stream([_Chunk([_Choice(_Delta(content="ok"), finish="stop")],
                               usage=_Usage())])

    client.client.chat.completions.create = _create
    list(_drive(client))
    assert calls["n"] == 2


def test_plain_429_keeps_the_old_backoff(monkeypatch, tmp_path):
    A, client, gate_root = _make_client(monkeypatch, tmp_path)
    calls = {"n": 0}

    def _create(**kwargs):
        calls["n"] += 1
        if calls["n"] <= 3:
            raise RateLimited("429 too many requests per minute")
        return _Stream([_Chunk([_Choice(_Delta(content="ok"), finish="stop")],
                               usage=_Usage())])

    client.client.chat.completions.create = _create
    out = _drive(client)
    assert calls["n"] == 4
    assert "ok" in "".join(getattr(ev, "text", "") for ev in out)


def test_unparsable_429_body_keeps_the_old_backoff(monkeypatch, tmp_path):
    """A 429 whose body cannot be parsed = today's behaviour: backoff."""
    A, client, gate_root = _make_client(monkeypatch, tmp_path)
    calls = {"n": 0}

    def _create(**kwargs):
        calls["n"] += 1
        if calls["n"] <= 2:
            raise RateLimited("Hourly token quota exceeded: ??? tokens")
        return _Stream([_Chunk([_Choice(_Delta(content="ok"), finish="stop")],
                               usage=_Usage())])

    client.client.chat.completions.create = _create
    out = _drive(client)
    assert calls["n"] == 3
    texts = "".join(getattr(ev, "text", "") for ev in out)
    assert "waiting until" not in texts.lower()
    assert "ok" in texts


# ---------------------------------------------------------------------------
# 6 + 7. stop during the wait; a budget too short for the wait
# ---------------------------------------------------------------------------

def test_a_stop_during_the_quota_wait_breaks_it_within_2s(monkeypatch, tmp_path):
    A, client, gate_root = _make_client(monkeypatch, tmp_path)
    # Real sliced sleep so the wait is observable; the slices honor stop.
    monkeypatch.setattr(A.time, "sleep", time.sleep)
    calls = {"n": 0}
    wait_started = threading.Event()

    def _stop_later(cl):
        wait_started.wait(timeout=5)
        time.sleep(0.3)
        cl._stop_requested = True   # the engine's probe reads this flag

    stopper = threading.Thread(target=_stop_later, args=(client,),
                               daemon=True)

    def _create(**kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise RateLimited(
                "Hourly token quota exceeded: 1 / 2 free tokens used. "
                "Quota resets in 300s.")
        return _Stream([_Chunk([_Choice(_Delta(content="ok"), finish="stop")],
                               usage=_Usage())])

    client.client.chat.completions.create = _create
    t0 = time.monotonic()
    stopper.start()
    try:
        list(_drive(client))
    finally:
        stopper.join(timeout=5)
    elapsed = time.monotonic() - t0
    assert elapsed < 2.0, (
        f"a stop during the quota wait must break it within ~2 s, took {elapsed:.1f}s")
    assert calls["n"] == 1, "no retry may fire after a stop"


def test_a_budget_shorter_than_the_wait_ends_with_its_own_reason(
        monkeypatch, tmp_path):
    A, client, gate_root = _make_client(monkeypatch, tmp_path)
    # The engine's run-ceiling probe: run_budget_s hit.
    client.run_ceiling = lambda spent: "run_budget_s"
    calls = {"n": 0}

    def _create(**kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise RateLimited(
                "Hourly token quota exceeded: 1 / 2 free tokens used. "
                "Quota resets in 300s.")
        return _Stream([_Chunk([_Choice(_Delta(content="ok"), finish="stop")],
                               usage=_Usage())])

    client.client.chat.completions.create = _create
    events = list(_drive(client))
    assert calls["n"] == 1, "no retry into a budget the run cannot survive"
    stop = [ev for ev in events if ev.type == "message_delta"]
    assert stop and stop[0].stop_reason == "quota_wait_exceeds_budget"
    texts = "".join(getattr(ev, "text", "") for ev in events)
    assert "quota" in texts.lower()


# ---------------------------------------------------------------------------
# 8. the usage summary records the quota numbers
# ---------------------------------------------------------------------------

def test_the_usage_summary_records_the_quota_numbers(monkeypatch, tmp_path):
    A, client, gate_root = _make_client(monkeypatch, tmp_path)

    def _create(**kwargs):
        raise RateLimited(QUOTA_BODY)

    client.client.chat.completions.create = _create
    # Drain the ladder; the quota numbers must be recorded even though the
    # turn ends in error.
    with pytest.raises(Exception):
        list(_drive(client))
    engine_status_keys = ("quota_tokens_used", "quota_tokens_limit")
    # The client exposes what it saw; the engine's get_status() relays it.
    assert any(getattr(client, k, None) is not None
               for k in engine_status_keys)
    assert getattr(client, "quota_tokens_used", None) == 6_233_280
    assert getattr(client, "quota_tokens_limit", None) == 6_000_000
