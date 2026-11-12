"""A person's interactive session never queues behind headless runs.

Seen on the Pi on 2026-11-12: in the person's phone and dashboard sessions
nearly every turn showed "waiting for a slot on …: 2/3/4 in use", once 90 s,
while headless agent runs (swarm members, ``--unattended``) shared the same
KIT key. The gate is a poll race over slot-0..N (no FIFO, no priority), so
interactive and headless contended equally.

The fix reserves the lowest slot (slot-0) for interactive requests:
``acquire(interactive=True)`` scans ONLY for slot-0; a headless request
scans slots 1..N and never takes slot-0. A person's session therefore never
waits behind headless runs, headless still progress on N-1 slots, and total
concurrency over the endpoint stays within the configured cap.
"""

from __future__ import annotations

import json

import pytest

from conftest import child_env  # noqa: F401  (kept importable, matching the suite)
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
    slot-0 even when every headless slot is taken."""
    g = G.EndpointGate(URL, "m")
    _limit(g, 4)
    # Four headless holders fill slots 1..4; slot-0 must stay free for the
    # interactive session alone.
    holders = [G.EndpointGate(URL, "m") for _ in range(4)]
    for h in holders:
        h.acquire()                      # headless -> slots 1..4
    interactive = G.EndpointGate(URL, "m", interactive=True)
    assert interactive.acquire(timeout_s=1.0) == 0     # the reserved slot
    assert "waiting for a slot" not in interactive.wait_notice
    interactive.release()


def test_a_headless_request_never_takes_the_reserved_slot(runtime):
    """Headless runs leave slot-0 alone: they use slots 1..N, so the <=N
    cap is kept and an interactive request is always one slot away."""
    g = G.EndpointGate(URL, "m")
    _limit(g, 4)
    taken = [h.acquire() for h in (G.EndpointGate(URL, "m") for _ in range(4))]
    assert 0 not in taken               # no headless holder grabbed slot-0


def test_one_reserved_slot_serves_a_person_behind_a_fixed_cap(runtime):
    """With 1 interactive + 4 headless callers the interactive request is
    served at once while the headless stay within the 4-slot cap."""
    g = G.EndpointGate(URL, "m")
    _limit(g, 4)
    headless = [G.EndpointGate(URL, "m") for _ in range(4)]
    for h in headless:
        h.acquire()                     # headless take slots 1..4
    person = G.EndpointGate(URL, "m", interactive=True)
    assert person.acquire(timeout_s=1.0) == 0      # never queues behind them
    # One more headless caller cannot exceed the cap.
    extra = G.EndpointGate(URL, "m")
    with pytest.raises(G.SlotWaitTimeout):
        extra.acquire(timeout_s=0.3)
