"""The adaptive limit is per (endpoint, MODEL), not per endpoint.

Measured 2026-10: one limit per ENDPOINT let a slow model (GLM, first byte
often minutes) halve the shared limit for every model behind that endpoint,
so the fast models (DeepSeek, Mistral) were throttled by the slow one's
congestion and the limit fell to 1 for hours. This suite locks the fix:

  * per-(endpoint, model) is the DEFAULT (DELFIN_ENDPOINT_PER_MODEL_LIMIT,
    when set to a falsy value, reverts to the legacy one-limit-per-endpoint);
  * a congestion signal on one model halves ONLY that model's limit;
  * the limit recovers on its own by +1 after a quiet period without a
    congestion signal (no "quick success" needs to occur -- under load none
    ever does, which is why the limit stayed down before);
  * the halving floor is 2 per model;
  * the endpoint still has a TOTAL cap when per-model limits are on, so many
    models at their own high limits cannot together overwhelm the endpoint;
  * waiting for a local slot never counts as a congestion signal.

No live endpoint. The gate is real, over a private runtime dir; adaptive
state is written straight to state.json so the mechanism is exercised
without sleeping for the real quiet period.
"""

from __future__ import annotations

from conftest import child_env  # noqa: E402
import json
import os
import subprocess
import sys
import textwrap
import time

import pytest

from delfin.agent import endpoint_gate as G


URL = "https://kit.example.org/api/v1"


@pytest.fixture
def runtime(monkeypatch, tmp_path):
    """A private runtime dir, as rundir resolves it, private to one test so
    per-model state never leaks across tests."""
    rt = tmp_path / "rt"
    rt.mkdir()
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(rt))
    monkeypatch.delenv("DELFIN_ENDPOINT_MAX_CONCURRENCY", raising=False)
    monkeypatch.delenv("DELFIN_ENDPOINT_PER_MODEL_LIMIT", raising=False)
    monkeypatch.setattr("delfin.user_settings.load_settings", lambda: {})
    return rt


def _gate_with_limit(runtime, model, limit):
    g = G.EndpointGate(URL, model)
    d = g._dir()
    (d / "state.json").write_text(json.dumps({"limit": limit}))
    return g, d


def _stall(exc_text="no byte from the endpoint in 300s"):
    return TimeoutError(exc_text)


# -- key structure -----------------------------------------------------------

def test_per_model_is_the_default(runtime):
    """Without the env var, two models behind one endpoint get DISTINCT
    keys -> DISTINCT adaptive limits."""
    g1 = G.EndpointGate(URL, "slow")
    g2 = G.EndpointGate(URL, "fast")
    assert g1._key != g2._key
    assert g1._key.endswith("@slow")
    assert g2._key.endswith("@fast")


def test_the_env_can_revert_to_one_limit_per_endpoint(monkeypatch, runtime):
    """DELFIN_ENDPOINT_PER_MODEL_LIMIT=0 turns the per-model default off:
    one key per endpoint, so the legacy behaviour is reachable."""
    monkeypatch.setenv("DELFIN_ENDPOINT_PER_MODEL_LIMIT", "0")
    g1 = G.EndpointGate(URL, "slow")
    g2 = G.EndpointGate(URL, "fast")
    assert g1._key == g2._key
    assert "@slow" not in g1._key and "@fast" not in g2._key


# -- one model's congestion leaves the others alone --------------------------

def test_a_slow_model_halves_only_its_own_limit(runtime):
    g_slow, d_slow = _gate_with_limit(runtime, "slow", 8)
    g_fast, d_fast = _gate_with_limit(runtime, "fast", 8)
    G.report_result(g_slow, exception=_stall())
    assert g_slow.current_limit(d_slow) == 4       # slow halved
    assert g_fast.current_limit(d_fast) == 8       # fast keeps all its slots


def test_the_fast_model_keeps_its_slots(runtime):
    """Repeat halving on the slow model never touches the fast one: the
    fast model keeps its full starting limit."""
    g_slow, d_slow = _gate_with_limit(runtime, "slow", 8)
    g_fast, d_fast = _gate_with_limit(runtime, "fast", 8)
    for _ in range(6):
        G.report_result(g_slow, exception=_stall())
    assert g_slow.current_limit(d_slow) == 2       # 8 -> 4 -> 2 -> 2... (floor)
    assert g_fast.current_limit(d_fast) == 8       # still full, not halved


# -- recovery ----------------------------------------------------------------

def test_the_limit_recovers_after_the_quiet_period(runtime):
    """+1 after a quiet period without a congestion signal -- no quick
    success needed. The quiet baseline is recorded as last_congestion_ts; a
    success past it raises the limit by one."""
    g, d = _gate_with_limit(runtime, "m", 4)
    past = time.monotonic() - (G._QUIET_RECOVER_S + 60.0)
    state = {"limit": 4, "last_congestion_ts": past}
    (d / "state.json").write_text(json.dumps(state))
    G.report_result(g, first_byte_s=1.0)           # a slow-ish environment
    assert g.current_limit(d) == 5


def test_a_congestion_signal_resets_the_quiet_period(runtime):
    """A halving sets last_congestion_ts=now, so the +1 cannot fire for the
    full quiet period after it."""
    g, d = _gate_with_limit(runtime, "m", 4)
    (d / "state.json").write_text(json.dumps({"limit": 4}))
    G.report_result(g, exception=_stall())           # 4 -> 2, baseline=now
    assert g.current_limit(d) == 2
    # immediately report a success: the 120 s is NOT up, so no +1.
    G.report_result(g, first_byte_s=1.0)
    assert g.current_limit(d) == 2


# -- floor -------------------------------------------------------------------

def test_the_halving_floor_holds_at_two(runtime):
    g, d = _gate_with_limit(runtime, "m", 2)
    G.report_result(g, exception=_stall())
    assert g.current_limit(d) == 2                  # floor 2: not 1


def test_halving_three_lands_on_two_not_one(runtime):
    g, d = _gate_with_limit(runtime, "m", 3)
    G.report_result(g, exception=_stall())
    assert g.current_limit(d) == 2                  # 3//2=1 -> floor 2


# -- slot-wait timeout is not congestion -------------------------------------

def test_a_slot_wait_timeout_is_not_a_congestion_signal(runtime):
    g, d = _gate_with_limit(runtime, "m", 8)
    timeout = G.SlotWaitTimeout("no slot on https://x within 900s")
    assert not G.halve_for_exception(timeout)
    G.report_result(g, exception=timeout)
    assert g.current_limit(d) == 8                  # never halved


def test_slot_wait_timeout_from_acquire_never_halves(runtime):
    """The path the stream loop actually takes: acquire() times out AND the
    retry nothing reports a congestion -- the round's finally reports the
    exception and must not halve the limit."""
    g, d = _gate_with_limit(runtime, "m", 1)
    holder = G.EndpointGate(URL, "m")
    assert holder.acquire() == 0                    # the only slot is taken
    with pytest.raises(G.SlotWaitTimeout):
        G.EndpointGate(URL, "m").acquire(timeout_s=0.2)
    holder.release()
    assert g.current_limit(d) == 1                  # state untouched


# -- the endpoint still has a total cap --------------------------------------

_OVERLAP_WORKER = textwrap.dedent("""
    import json, os, sys, time
    url, model, cap, out = sys.argv[1:5]
    from delfin.agent.endpoint_gate import EndpointGate
    g = EndpointGate(url, model)
    events = []
    idx = g.acquire()
    t0 = time.monotonic()
    events.append(("in", idx, t0))
    time.sleep(0.5)
    g.release()
    events.append(("out",))
    open(out, "w").write(json.dumps(events))
""")


def test_per_model_limits_still_share_the_endpoints_total_cap(runtime):
    """config cap=2 over FOUR models, each requesting a slot at once: no
    matter how high each model's own limit is, the endpoint total never
    exceeds the cap. The models' adaptive state says 16 each -- the cap is
    the only thing that can bind."""
    cap = 2
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    scratch = runtime.parent / "workers2"
    scratch.mkdir(exist_ok=True)
    env = child_env(runtime.parent)
    env["DELFIN_ENDPOINT_MAX_CONCURRENCY"] = str(cap)
    env["DELFIN_RUNTIME_DIR"] = str(runtime)
    env["PYTHONPATH"] = repo
    env.pop("DELFIN_ENDPOINT_PER_MODEL_LIMIT", None)
    models = [f"m{i}" for i in range(4)]
    procs, outs = [], []
    for i, model in enumerate(models):
        out = scratch / f"w{i}.json"
        outs.append(out)
        procs.append(subprocess.Popen(
            [sys.executable, "-c", _OVERLAP_WORKER,
             URL, model, str(cap), str(out)],
            env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True))
    for p in procs:
        _, err = p.communicate(timeout=120)
        assert p.returncode == 0, err
    ins = [json.loads(o.read_text()) for o in outs]
    series = sorted((e[2], e[2] + 0.5) for ev in ins for e in ev
                    if e[0] == "in")
    overlap = 0
    for a, b in series:
        held = sum(1 for c, d in series if c <= a < d)
        overlap = max(overlap, held)
    assert overlap == cap, series                 # never above the endpoint cap


def test_a_url_shown_to_the_user_carries_no_credentials():
    from delfin.agent.endpoint_gate import display_url
    assert display_url("https://u:secret@ex.org:8443/api/v1?k=1") == "https://ex.org:8443/api/v1"
    assert "secret" not in display_url("https://u:secret@ex.org/v1")
