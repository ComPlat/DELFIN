"""The endpoint gate keeps N processes under one shared limit.

Measured 2026-10 in a 60-agent run: every DELFIN process on one host
called the KIT endpoint on its own; together they exceeded its capacity
and each retried independently. The gate is a host-wide semaphore over
file locks, so the tests here run REAL child processes
(``multiprocessing``, each with the suite's child-safe environment) -- a
thread test would prove nothing about sharing across processes.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
import time

import pytest

from conftest import child_env
from delfin.agent import endpoint_gate as G


@pytest.fixture
def runtime(monkeypatch, tmp_path):
    """A private runtime directory, as rundir resolves it."""
    rt = tmp_path / "rt"
    rt.mkdir()
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(rt))
    monkeypatch.delenv("DELFIN_ENDPOINT_MAX_CONCURRENCY", raising=False)
    monkeypatch.delenv("DELFIN_ENDPOINT_PER_MODEL_LIMIT", raising=False)
    monkeypatch.setattr("delfin.user_settings.load_settings", lambda: {})
    return rt


URL = "https://kit.example.org/api/v1"


# -- unit: keying, limits, state -------------------------------------------

def test_the_default_key_includes_the_model(runtime):
    """Per-model is the default: two models behind one endpoint get distinct
    adaptive keys, while the endpoint-wide cap key stays shared."""
    g = G.EndpointGate(URL, "m")
    assert g.enabled
    g2 = G.EndpointGate(URL + "/", "other")
    assert g._key != g2._key          # different MODEL: distinct limits
    assert g._cap_key == g2._cap_key  # same endpoint: one TOTAL cap
    g3 = G.EndpointGate("http://kit.example.org/api/v1", "m")
    assert g3._key != g._key          # different scheme: not the same host


def test_no_url_and_zero_disable_the_gate(runtime, monkeypatch):
    assert not G.EndpointGate("", "m").enabled
    monkeypatch.setenv("DELFIN_ENDPOINT_MAX_CONCURRENCY", "0")
    assert not G.EndpointGate(URL, "m").enabled


def test_disabled_gate_acquires_none(runtime, monkeypatch):
    monkeypatch.setenv("DELFIN_ENDPOINT_MAX_CONCURRENCY", "0")
    g = G.EndpointGate(URL, "m")
    with g.slot() as index:
        assert index is None
    g.release()                       # safe when nothing was held


def test_settings_limit_is_honoured(runtime, monkeypatch):
    monkeypatch.setattr(
        "delfin.user_settings.load_settings",
        lambda: {"agent": {"endpoint_max_concurrency": 3}})
    g = G.EndpointGate(URL, "m")
    assert g.current_limit() == 3


def test_limit_below_one_falls_back_to_default(runtime):
    g = G.EndpointGate(URL, "m")
    assert g.current_limit() == G.DEFAULT_MAX_CONCURRENCY


def test_a_lower_setting_clamps_at_once_and_raising_climbs_back(
        runtime, monkeypatch):
    """The configured value is the CEILING of the adaptive limit (the
    AIMD state may sit anywhere at or below it). Lowering the setting
    takes effect immediately; raising it does not snap a halved limit
    back up -- it climbs again on quick successes, never above the
    ceiling, and never below the halving floor."""
    g = G.EndpointGate(URL, "m")
    d = g._dir()
    (d / "state.json").write_text(json.dumps({"limit": 6}))
    monkeypatch.setenv("DELFIN_ENDPOINT_MAX_CONCURRENCY", "4")
    assert g.current_limit(d) == 4     # adaptive state 6, clamped at ceiling 4
    G.report_result(g, exception=TimeoutError(
        "no byte from the endpoint in 300s"))
    assert g.current_limit(d) == 2     # 4//2=2, the floor
    monkeypatch.setenv("DELFIN_ENDPOINT_MAX_CONCURRENCY", "16")
    for _ in range(G._QUICK_SUCCESS_RUN):
        G.report_result(g, first_byte_s=0.5)
    assert g.current_limit(d) == 3     # climbing again, ceiling now 16


# -- unit: the AIMD feedback ------------------------------------------------

def _gate_with_limit(runtime, limit):
    g = G.EndpointGate(URL, "m")
    d = g._dir()
    (d / "state.json").write_text(json.dumps({"limit": limit}))
    return g, d


def test_a_429_halves_the_shared_limit(runtime):
    g, d = _gate_with_limit(runtime, 16)
    G.report_result(g, exception=RuntimeError("Error code: 429 - rate"))
    assert g.current_limit(d) == 8


def test_a_503_halves_the_shared_limit(runtime):
    g, d = _gate_with_limit(runtime, 16)
    exc = RuntimeError("Error code: 503 - Service Unavailable")
    exc.status_code = 503
    G.report_result(g, exception=exc)
    assert g.current_limit(d) == 8


def test_a_stall_halves_the_shared_limit(runtime):
    g, d = _gate_with_limit(runtime, 16)
    G.report_result(g, exception=TimeoutError(
        "no byte from the endpoint in 300s"))
    assert g.current_limit(d) == 8


def test_a_stall_never_goes_below_the_halving_floor(runtime):
    g, d = _gate_with_limit(runtime, 1)
    G.report_result(g, exception=TimeoutError(
        "no byte from the endpoint in 300s"))
    assert g.current_limit(d) == G._HALVING_FLOOR  # floor 2, not 1


def test_a_deterministic_failure_is_no_signal(runtime):
    g, d = _gate_with_limit(runtime, 16)
    G.report_result(g, exception=RuntimeError(
        "Error code: 400 - model not found"))
    assert g.current_limit(d) == 16


def test_quick_successes_climb_back(runtime):
    g, d = _gate_with_limit(runtime, 8)
    for _ in range(G._QUICK_SUCCESS_RUN):
        G.report_result(g, first_byte_s=1.0)
    assert g.current_limit(d) == 9
    for _ in range(G._QUICK_SUCCESS_RUN):
        G.report_result(g, first_byte_s=1.0)
    assert g.current_limit(d) == 10


def test_a_slow_success_does_not_climb(runtime):
    g, d = _gate_with_limit(runtime, 8)
    for _ in range(G._QUICK_SUCCESS_RUN * 2):
        G.report_result(g, first_byte_s=60.0)
    assert g.current_limit(d) == 8


def test_a_congestion_signal_resets_the_climb(runtime):
    g, d = _gate_with_limit(runtime, 8)
    for _ in range(G._QUICK_SUCCESS_RUN - 1):
        G.report_result(g, first_byte_s=1.0)
    G.report_result(g, exception=TimeoutError(
        "no byte from the endpoint in 300s"))
    # 8 -> 4, and the near-complete climb is gone.
    assert g.current_limit(d) == 4
    assert json.loads((d / "state.json").read_text())["quick_run"] == 0


def test_the_climb_never_passes_the_configured_max(runtime, monkeypatch):
    g, d = _gate_with_limit(runtime, 15)
    monkeypatch.setenv("DELFIN_ENDPOINT_MAX_CONCURRENCY", "16")
    for _ in range(G._QUICK_SUCCESS_RUN * 10):
        G.report_result(g, first_byte_s=0.5)
    assert g.current_limit(d) == 16


def test_reporting_never_raises(runtime):
    g = G.EndpointGate(URL, "m")
    G.report_result(None, first_byte_s=1.0)      # no gate: nothing
    G.report_result(g, exception=ValueError("x"))  # non-congestion: nothing
    bad = G.EndpointGate("", "")                  # unkeyed
    G.report_result(bad, first_byte_s=1.0)


# -- unit: slots -------------------------------------------------------------

def test_a_slot_blocks_a_second_holder(runtime):
    g = G.EndpointGate(URL, "m")
    monkeypatch_env = {"limit": 1}
    (g._dir() / "state.json").write_text(json.dumps(monkeypatch_env))
    first = G.EndpointGate(URL, "m")
    assert first.acquire() == 0
    second = G.EndpointGate(URL, "m")
    with pytest.raises(G.SlotWaitTimeout):
        second.acquire(timeout_s=0.3)
    assert "waiting for a slot" in second.wait_notice
    assert "1 in use" in second.wait_notice
    first.release()
    assert second.acquire(timeout_s=1.0) == 0
    second.release()


def test_a_second_acquire_on_one_gate_is_refused(runtime):
    """One slot per gate instance: a re-entrant acquire must not
    overwrite the outer slot's descriptor (the flock would leak on
    release). A second request in the same process builds its own gate."""
    g = G.EndpointGate(URL, "m")
    (g._dir() / "state.json").write_text(json.dumps({"limit": 2}))
    assert g.acquire() == 0
    with pytest.raises(RuntimeError):
        g.acquire(timeout_s=0.3)
    g.release()


def test_release_in_a_with_block_always_happens(runtime):
    g = G.EndpointGate(URL, "m")
    (g._dir() / "state.json").write_text(json.dumps({"limit": 1}))
    with pytest.raises(ZeroDivisionError):
        with g.slot():
            assert g._held_fd is not None          # the slot is held
            1 / 0
    # The with-block's finally released the slot: acquiring works again.
    g2 = G.EndpointGate(URL, "m")
    assert g2.acquire(timeout_s=1.0) == 0
    g2.release()


def test_a_killed_holder_frees_its_slot(runtime, tmp_path):
    """The core promise: a slot held by a dead process is free, because
    the flock dies with the file descriptor."""
    import os as _os
    repo = _os.path.dirname(_os.path.dirname(
        _os.path.abspath(__file__)))
    child = textwrap.dedent(f"""
        import sys
        sys.path.insert(0, {repo!r})
        import os, time
        os.environ.setdefault("DELFIN_RUNTIME_DIR",
                              {os.environ["DELFIN_RUNTIME_DIR"]!r})
        from delfin.agent.endpoint_gate import EndpointGate
        # One slot in total, so "the child holds a slot" means
        # "no slot is free".
        g = EndpointGate({URL!r}, "m")
        (g._dir() / "state.json").write_text(
            __import__("json").dumps({{"limit": 1}}))
        i = g.acquire()
        print(i, flush=True)
        time.sleep(300)          # killed by the parent, never releases
    """)
    proc = subprocess.Popen(
        [sys.executable, "-c", child],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        env=child_env(tmp_path))
    try:
        line = proc.stdout.readline().strip()
        assert line == "0"
        g = G.EndpointGate(URL, "m")
        (g._dir() / "state.json").write_text(json.dumps({"limit": 1}))
        with pytest.raises(G.SlotWaitTimeout):
            g.acquire(timeout_s=0.5)          # held by the live child
        proc.kill()
        proc.wait(timeout=10)
        assert g.acquire(timeout_s=5.0) == 0  # the kernel freed it
        g.release()
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait(timeout=10)


# -- multi-process: the shared limit holds -----------------------------------

_WORKER = textwrap.dedent("""
    import json, os, sys, time
    url, n, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    from delfin.agent.endpoint_gate import EndpointGate
    g = EndpointGate(url, "m")
    events = []
    t0 = time.monotonic()
    idx = g.acquire()
    events.append(("in", idx, time.monotonic() - t0))
    time.sleep(0.4)
    g.release()
    events.append(("out",))
    open(out, "w").write(json.dumps(events))
""")


def _run_workers(runtime, n, url=URL):
    procs, outs = [], []
    scratch = runtime.parent / "workers"
    scratch.mkdir(exist_ok=True)
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for i in range(n):
        out = scratch / f"w{i}.json"
        outs.append(out)
        env = dict(os.environ)
        env["DELFIN_RUNTIME_DIR"] = str(runtime)
        env["PYTHONPATH"] = repo
        env.pop("DELFIN_ENDPOINT_MAX_CONCURRENCY", None)
        procs.append(subprocess.Popen(
            [sys.executable, "-c", _WORKER, url, "0", str(out)],
            env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True))
    for p in procs:
        _, err = p.communicate(timeout=120)
        assert p.returncode == 0, err
    events = []
    for out in outs:
        events.append(json.loads(out.read_text()))
    return events


def test_three_processes_never_exceed_limit_one(runtime):
    """With a limit of 1, at most one worker holds a slot at any time --
    across PROCESS boundaries, which a threads test cannot show."""
    g = G.EndpointGate(URL, "m")
    (g._dir() / "state.json").write_text(json.dumps({"limit": 1}))
    events = _run_workers(runtime, 3)
    timelines = sorted((e[2], e[2] + 0.4) for ev in events
                       for e in ev if e[0] == "in")
    overlap = 0
    for a, b in timelines:
        held = sum(1 for c, d in timelines if c <= a < d)
        overlap = max(overlap, held)
    assert overlap == 1
    assert G.EndpointGate(URL, "m").current_limit() == 1


def test_workers_wait_for_the_holder_to_release(runtime):
    """A worker that had to wait: the slot was contended and released by
    another process, and the wait notice named the count in use."""
    g = G.EndpointGate(URL, "m")
    (g._dir() / "state.json").write_text(json.dumps({"limit": 1}))
    holder = G.EndpointGate(URL, "m")
    assert holder.acquire() == 0
    scratch = runtime.parent / "workers"
    scratch.mkdir(exist_ok=True)
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env = dict(os.environ)
    env["DELFIN_RUNTIME_DIR"] = str(runtime)
    env["PYTHONPATH"] = repo
    env.pop("DELFIN_ENDPOINT_MAX_CONCURRENCY", None)
    outs = [scratch / f"w{i}.json" for i in range(2)]
    procs = [subprocess.Popen(
        [sys.executable, "-c", _WORKER, URL, "0", str(out)],
        env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True) for out in outs]
    time.sleep(1.0)                    # both workers are waiting now
    holder.release()                   # free the single slot for them
    for p in procs:
        _, err = p.communicate(timeout=60)
        assert p.returncode == 0, err
    # Each worker had to WAIT before it got the slot (the notice the
    # stream loop surfaces carries the same in-use count).
    for out in outs:
        events = json.loads(out.read_text())
        for e in events:
            if e[0] == "in":
                assert e[2] > 0.5, events   # waited for the holder


def test_credentials_in_the_url_never_reach_the_key(runtime):
    """The key names files on disk: user info, query and fragment of the
    base URL are not part of it, so a password in the URL never becomes a
    file name (and never goes into a hash)."""
    key = G._endpoint_key("https://alice:s3cr3t-pw@kit.example.org/api/v1?token=t0k#frag")
    assert "s3cr3t" not in key and "alice" not in key
    assert "t0k" not in key and "frag" not in key
    assert key == G._endpoint_key("https://kit.example.org/api/v1")
    long_url = "https://bob:pw-9x@" + "a" * 200 + ".example.org/v1"
    long_key = G._endpoint_key(long_url)
    assert "pw-9x" not in long_key and "bob" not in long_key
    assert len(long_key) <= 160
