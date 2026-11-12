"""A person's interactive session never queues behind headless runs.

Seen on the Pi on 2026-11-12: in the person's phone and dashboard sessions
nearly every turn showed "waiting for a slot on …: 2/3/4 in use", once 90 s,
while headless agent runs (swarm members, ``--unattended``) shared the same
KIT key. The gate is a poll race over slot-0..N (no FIFO, no priority), so
interactive and headless contended equally.

The fix reserves the lowest slot (slot-0) for interactive requests:
``acquire(interactive=True)`` scans ONLY for slot-0; a headless request
scans slots 1..N and never takes slot-0. A person's session therefore never
waits behind headless runs, headless still progress on the rest of the cap,
and total concurrency over the endpoint stays within the configured limit.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import endpoint_gate as G


@pytest.fixture
def runtime(monkeypatch, tmp_path):
    rt = tmp_path / "rt"
    rt.mkdir()
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(rt))
    monkeypatch.delenv("DELFIN_ENDPOINT_MAX_CONCURRENCY", raising=False)
    monkeypatch.delenv("DELFIN_ENDPOINT_PER_MODEL_LIMIT", raising=False)
    monkeypatch.setattr("delfin.user_settings.load_settings", lambda: {})
    return rt


URL = "https://kit.example.org/api/v1"


def _limit(g: G.EndpointGate, n: int):
    (g._dir() / "state.json").write_text(json.dumps({"limit": n}))


def test_an_interactive_request_takes_the_reserved_low_slot(runtime):
    """A person's request never queues behind a headless holder: it gets
    slot-0 even when every other slot is taken by headless runs."""
    g = G.EndpointGate(URL, "m")
    _limit(g, 5)                     # the provider cap
    holders = [G.EndpointGate(URL, "m") for _ in range(4)]   # headless
    for h in holders:
        h.acquire()                  # headless -> slots 1..4
    interactive = G.EndpointGate(URL, "m", interactive=True)
    assert interactive.acquire(timeout_s=1.0) == 0     # the reserved slot
    assert "waiting for a slot" not in interactive.wait_notice
    interactive.release()


def test_a_headless_request_prefers_the_high_slots(runtime):
    """Headless runs take slots 1..N first, so slot-0 tends to stay free
    for a person while the total concurrency stays within the cap."""
    g = G.EndpointGate(URL, "m")
    _limit(g, 5)
    taken = [h.acquire() for h in (G.EndpointGate(URL, "m")
                                   for _ in range(4))]
    assert 0 not in taken            # no headless holder grabbed slot-0
    assert len(taken) == 4           # all four headless still progress


def test_one_reserved_slot_serves_a_person_behind_a_fixed_cap(runtime):
    """1 interactive + 4 headless callers: the interactive request is served
    at once and the headless stay within the cap -- the total held never
    exceeds the configured limit."""
    g = G.EndpointGate(URL, "m")
    _limit(g, 5)
    headless = [G.EndpointGate(URL, "m") for _ in range(4)]
    for h in headless:
        h.acquire()                  # headless take slots 1..4
    person = G.EndpointGate(URL, "m", interactive=True)
    assert person.acquire(timeout_s=1.0) == 0      # never queues behind them
    # One more headless caller cannot exceed the cap (slots 1..4 plus the
    # reserved 0 are all held).
    extra = G.EndpointGate(URL, "m")
    with pytest.raises(G.SlotWaitTimeout):
        extra.acquire(timeout_s=0.3)


def test_a_headless_request_takes_slot_0_only_as_last_resort(runtime):
    """When every preferred slot is held, a headless request still uses
    slot-0 -- the reservation must not cut headless capacity."""
    g = G.EndpointGate(URL, "m")
    _limit(g, 2)
    holders = [G.EndpointGate(URL, "m") for _ in range(2)]
    assert [h.acquire() for h in holders] == [1, 0]   # prefer 1, then 0
    for h in holders:
        h.release()


def test_per_request_interactive_override_serves_the_person(runtime):
    """The session kind can also be passed per request, so a caller that is
    not built interactive can still reserve the slot for one request."""
    g = G.EndpointGate(URL, "m")
    _limit(g, 2)
    headless = G.EndpointGate(URL, "m")
    headless.acquire()               # headless -> slot 1
    person = G.EndpointGate(URL, "m")
    assert person.acquire(timeout_s=1.0, interactive=True) == 0


def test_for_client_reads_the_session_marker(runtime, monkeypatch):
    """for_client wires the session kind: an attended session's env marker
    builds an interactive gate, a headless process (no marker) does not."""
    from delfin.agent.approval_guard import ENV_MARKER
    monkeypatch.delenv(G.INTERACTIVE_ENV, raising=False)
    assert G.for_client(URL, "m").interactive is False
    monkeypatch.setenv(G.INTERACTIVE_ENV, "1")
    assert G.for_client(URL, "m").interactive is True
    # the marker must be its own variable: it may not collide with the
    # agent-session marker used for approvals
    assert G.INTERACTIVE_ENV != ENV_MARKER
