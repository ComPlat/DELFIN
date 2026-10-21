"""UX2 adversarial tests for the phase-3 watch/emission layer (nacht-s4).

All cases use emit=False, so they need no attention kind from the
operator's build: they are red NOW against 880bd2da if the logic is
wrong, independent of .gate/ux2_attention_kinds.patch.
"""

from delfin.agent import session_health as health_mod


def _rec(**kw):
    rec = {
        "name": "nacht-s9",
        "session_id": "abc12345",
        "paused": False,
        "pid_alive": True,
        "open_approvals": [],
        "approval_opened_at": None,
        "last_turn_end": None,
        "last_turn": {},
        "idle_s": 0,
        "api_errors_in_a_row": 0,
        "resume_planned": False,
        "open_work": False,
        "started_at": None,
        "title": "demo",
    }
    rec.update(kw)
    return rec


def _state():
    return {"seen": {}}


def test_watch_reemits_after_leaving_and_reentering():
    """The brief: dedup holds only UNTIL the session leaves the state.
    Leaving must clear the (session, verdict) entry, so re-entering emits
    again — regardless of any attention kind."""
    state = _state()
    paused = _rec(paused=True)
    working = _rec()
    again_paused = _rec(paused=True)
    health_mod.watch_events([paused], 1_000_000.0, state, emit=False)
    health_mod.watch_events([working], 1_000_001.0, state, emit=False)
    again = health_mod.watch_events(
        [again_paused], 1_000_002.0, state, emit=False)
    assert len(again) == 1, (
        "re-entering a state must emit again; the leave-path must clear "
        "the (name, verdict) dedup entry, not pop by bare name")


def test_watch_uses_settings_thresholds():
    """agent.health.stalled_min=1 must shorten the stalled threshold at
    the real call shape: watch_events passes settings, classify must see
    them. With idle 90 s and stalled_min=1 the session is stalled."""
    state = _state()
    settings = {"health": {"stalled_min": 1}}
    rec = _rec(last_turn={"stopped": True}, open_work=True, idle_s=90)
    events = health_mod.watch_events(
        [rec], 1_000_000.0, state, settings=settings, emit=False)
    verdicts = [text for _, text in events]
    assert any("stalled" in t for t in verdicts), (
        "settings thresholds must reach classify through watch_events; "
        f"got {verdicts!r}")


def test_watch_does_not_emit_before_grace_for_waiting_on_person():
    """A waiting_on_person younger than 15 min is not announced yet; the
    same record older than 15 min is. Both with emit=False."""
    state = _state()
    young = _rec(open_approvals=["a1"], approval_opened_at=999_900.0)
    events = health_mod.watch_events(
        [young], 1_000_000.0, state, emit=False)
    assert events == [], "an approval open 100 s must not be announced"
    state = _state()
    old = _rec(open_approvals=["a1"], approval_opened_at=999_999.0 - 900)
    events = health_mod.watch_events(
        [old], 1_000_000.0, state, emit=False)
    assert len(events) == 1, "an approval open 15 min must be announced"


def test_watch_stop_all_set_then_cleared_then_set_again():
    """Stop-all set -> one event; cleared (fresh session) -> one; a SECOND
    stop after the clear must announce again, not be swallowed by the
    first stop's dedup mark."""
    stop_seen = {}
    state = _state()
    e1 = health_mod.watch_events(
        [_rec()], 1_000_000.0, state, emit=False, stop_all={"at": 900.0},
        stop_all_seen=stop_seen)
    assert len(e1) == 1
    # a fresh session started after the stop: the stop is over
    e2 = health_mod.watch_events(
        [_rec(started_at=1_000.0)], 1_001.0, state, emit=False,
        stop_all={"at": 900.0}, stop_all_seen=stop_seen)
    assert any("over" in t or "allowed again" in t for _, t in e2), e2
    # the stop is lifted entirely, then set AGAIN with a NEW time
    e3 = health_mod.watch_events(
        [_rec()], 2_000_000.0, state, emit=False, stop_all=None,
        stop_all_seen=stop_seen)
    e4 = health_mod.watch_events(
        [_rec()], 2_000_001.0, state, emit=False, stop_all={"at": 1_800.0},
        stop_all_seen=stop_seen)
    assert len(e4) == 1, (
        "a second stop-all must be announced again; stop_seen must not "
        f"hold the first stop's mark (got {e4!r}, stop_seen {stop_seen!r})")
