"""A person's interactive session never queues behind headless runs.

Seen on the Pi on 2026-11-12: in the person's phone and dashboard sessions
nearly every turn showed "waiting for a slot on …: 2/3/4 in use", once 90 s,
while headless agent runs (swarm members, ``--unattended``) shared the same
KIT key. The gate is a poll race over slot-0..N (no FIFO, no priority), so
interactive and headless contended equally.

The fix reserves the lowest slot (slot-0) for interactive requests:
``acquire(interactive=True)`` scans slot-0 first, then any free slot; a
headless request scans slots 1..N-1 and, with more than one slot, NEVER
takes slot-0. A person's session therefore never waits behind headless
runs, two of the person's sessions no longer queue on one slot while the
rest are free, headless still progress on the rest of the cap, and total
concurrency over the endpoint stays within the configured limit.
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


def test_an_interactive_request_takes_the_reserved_lowest_slot(runtime):
    """The person's session takes slot-0 when it is free."""
    g = G.EndpointGate(URL, "m")
    _limit(g, 5)
    interactive = G.EndpointGate(URL, "m", interactive=True)
    assert interactive.acquire(timeout_s=1.0) == 0     # the reserved slot
    assert "waiting for a slot" not in interactive.wait_notice
    interactive.release()


def test_a_headless_request_takes_the_high_slots(runtime):
    """Headless runs take slots 1..N-1 first, so slot-0 tends to stay free
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


def test_a_headless_request_never_takes_slot_0_under_full_load(runtime):
    """Even when every preferred slot is held, a headless request must not
    fall back into the person's slot: under saturation -- the case the
    reservation exists for -- headless grabbing slot-0 would put the person
    behind them again. The reservation does cost headless one slot when a
    person is not competing and the cap is full."""
    g = G.EndpointGate(URL, "m")
    _limit(g, 2)
    holders = [G.EndpointGate(URL, "m") for _ in range(2)]
    # Only one slot is usable by headless when N==2, so only one headless
    # holder succeeds; the second must wait rather than take the person's
    # slot-0.
    assert holders[0].acquire(timeout_s=1.0) == 1    # headless -> slot 1
    with pytest.raises(G.SlotWaitTimeout):
        holders[1].acquire(timeout_s=0.2)
    for h in holders:
        h.release()


def test_a_single_slot_still_serves_headless(runtime):
    """With only one slot there is no second slot to reserve: headless must
    still be able to use it (a cap of 1 keeps working)."""
    g = G.EndpointGate(URL, "m")
    _limit(g, 1)
    h = G.EndpointGate(URL, "m")
    assert h.acquire(timeout_s=1.0) == 0
    h.release()


def test_two_interactive_sessions_each_take_a_free_slot(runtime):
    """Two sessions of the same person (dashboard + phone) must not queue on
    the single slot-0 while slots 1..N-1 are free: interactive scans slot-0
    first, then any free slot, so the second interactive request takes the
    next slot instead of waiting."""
    g = G.EndpointGate(URL, "m")
    _limit(g, 4)
    a = G.EndpointGate(URL, "m", interactive=True)
    b = G.EndpointGate(URL, "m", interactive=True)
    taken_a = a.acquire(timeout_s=1.0)
    taken_b = b.acquire(timeout_s=1.0)
    assert taken_a == 0                 # the first prefers slot-0
    assert taken_b != taken_a           # the second takes a different slot
    assert taken_b in (1, 2, 3)
    a.release()
    b.release()


def test_per_request_interactive_override_serves_the_person(runtime):
    """The session kind can also be passed per request, so a caller that is
    not built interactive can still reserve the slot for one request."""
    g = G.EndpointGate(URL, "m")
    _limit(g, 2)
    headless = G.EndpointGate(URL, "m")
    headless.acquire()               # headless -> slot 1
    person = G.EndpointGate(URL, "m")
    assert person.acquire(timeout_s=1.0, interactive=True) == 0


def test_for_client_reads_the_process_marker(runtime, monkeypatch):
    """for_client wires the session kind from the process marker an attended
    surface sets (endpoint_gate.set_process_interactive): a headless process
    (never set) builds a headless gate, an attended one builds an
    interactive gate. The marker is a Python module global, NOT an
    os.environ value, so it is not inherited by child processes."""
    assert G.for_client(URL, "m").interactive is False
    G.set_process_interactive(True)
    try:
        assert G.for_client(URL, "m").interactive is True
        # A caller may still pass the kind per request.
        assert G.for_client(URL, "m", interactive=False).interactive is False
    finally:
        G.set_process_interactive(False)


def test_the_marker_is_not_an_env_var_and_is_not_inherited(runtime):
    """The interactive marker must not ride os.environ: it is a per-process
    flag, so a child process that forks+execs (a delegate, a --unattended
    run, a background job) starts headless even when the attended session
    that spawned it is interactive. An env var would leak into every child."""
    import subprocess
    import sys
    G.set_process_interactive(True)
    try:
        out = subprocess.run(
            [sys.executable, "-c",
             "from delfin.agent import endpoint_gate as G; "
             "print(G.for_client(%r, 'm').interactive)" % URL],
            capture_output=True, text=True, timeout=30,
            env={k: v for k, v in __import__('os').environ.items()},
        )
    finally:
        G.set_process_interactive(False)
    # The child starts a fresh interpreter: the module global is not there,
    # so its gate is headless even though the parent's is interactive.
    assert out.stdout.strip() == "False", out.stderr
    assert "DELFIN_INTERACTIVE_SESSION" not in __import__('os').environ


def test_the_dashboard_kernel_marks_its_process_interactive():
    # The dashboard's default ("kernel") mode builds its engine in
    # tab_agent and never passes through cli._start_interactive_session:
    # without its own call the person at the page queued behind headless
    # runs exactly as before. The call follows the engine it serves.
    from pathlib import Path
    src = (Path(__file__).resolve().parents[1]
           / "delfin" / "dashboard" / "tab_agent.py").read_text()
    build = src.index("engine = AgentEngine(\n                repo_dir=repo_dir,")
    after = src[build:build + 2500]
    assert "set_process_interactive(True)" in after
