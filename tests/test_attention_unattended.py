"""What the attention inbox owes a user who is not at the screen.

The scenario every test here is a piece of: the user starts a run at
18:00 and goes home. Overnight the agent asks one question and blocks on
it, a SLURM job dies, the permission gate refuses a command the model
then tries to walk around, the run reaches its cost ceiling and stops,
and the scheduler files a few hundred routine run notices in between. In
the morning the user opens the dashboard, answers the question, and goes
to make coffee while the agent picks the answer up.

Measured against a real inbox before these tests existed, that morning
looked like this: 893 pending records, all of them notices, oldest 15
days, with the parked question somewhere inside; the only cleanup
gesture on offer ("dismiss all") destroyed the question; the answer, if
given, never reached the agent because the session id had changed since
the question was parked, and no surface could show it any more; the
failed job was not in the inbox at all; the containment block died with
the process that recorded it; the budget stop left a silent desktop, an
empty inbox and a doctor reporting PASS; and nothing anywhere had ever
checked that a single notification transport on that machine worked.

Every test below fails on that code and passes on this one.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

import delfin.agent.attention as attention
import delfin.agent.notify as _notify_module

# Captured before any fixture mutes it — the tests that assert on the
# transport itself need the real implementation, not the mute.
_REAL_SEND_NOTIFICATION = _notify_module.send_notification


@pytest.fixture(autouse=True)
def _fake_home(tmp_path, monkeypatch):
    """Isolate ~ (inbox + settings) and mute the out-of-band transports."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", lambda: home)
    import delfin.agent.notify as notify
    monkeypatch.setattr(notify, "send_notification", lambda *a, **k: True)
    monkeypatch.setattr(
        notify, "send_remote_trigger",
        lambda *a, **k: notify.TriggerResult(sent=False, error="test"))
    return home


def _age(event_id: str, seconds: float) -> None:
    """Backdate a stored event (nothing here waits for real days)."""
    events = attention._load_events()
    for ev in events:
        if ev["id"] == event_id:
            ev["created_at"] = time.time() - seconds
    attention._write_events(events)


# ---------------------------------------------------------------------------
# 1. Notices expire; a blanket clear never destroys a blocking item
# ---------------------------------------------------------------------------

def test_blanket_clear_keeps_the_question_the_agent_is_blocked_on():
    """The one gesture for an inbox full of notices must not kill the ask.

    The agent's own timeout message tells the user the question can still
    be answered in the inbox. Clearing the notices then silently resolved
    it, with ``acknowledged=True`` so it could never be replayed either.
    """
    for i in range(3):
        attention.emit_attention("run_finished", title=f"run {i} finished")
    question = attention.emit_attention(
        "question_pending", session_id="s1", title="Which basis set?",
        options=["def2-SVP", "def2-TZVP"])

    res = attention.clear_all()
    assert res["ok"] is True
    assert res["cleared"] == 3
    assert res["kept"] == 1
    assert "blocking" in res["hint"]

    still = attention.list_pending()
    assert [ev["id"] for ev in still] == [question]
    assert still[0]["kind"] == "question_pending"


def test_blocking_item_is_clearable_by_id_or_explicit_opt_in():
    """Refusing is not the same as making it impossible."""
    a = attention.emit_attention("confirm_pending", title="confirm rm")
    b = attention.emit_attention("question_pending", title="which one?")
    assert attention.dismiss_item(a)["ok"] is True
    assert attention.clear_all("question_pending")["cleared"] == 1
    assert attention.list_pending() == []
    assert attention.drain_resolved("") == []      # neither replays

    c = attention.emit_attention("plan_pending", title="plan v3")
    assert attention.clear_all(include_blocking=True)["cleared"] == 1
    assert c not in {ev["id"] for ev in attention.list_pending()}


def test_notices_expire_and_blocking_events_never_do():
    """Nothing resolved a notice, and the prune exempted pending — so an
    inbox only ever grew. A notice is a report, not a to-do: it ages out.
    """
    finished = attention.emit_attention("run_finished", title="old finish")
    failed = attention.emit_attention("run_failed", title="old failure")
    question = attention.emit_attention("question_pending", title="old ask")
    for eid in (finished, failed, question):
        _age(eid, 15 * 86400)

    pending = {ev["id"] for ev in attention.list_pending()}
    assert finished not in pending          # 2-day TTL
    assert failed not in pending            # 7-day TTL
    assert question in pending              # never expires

    # A fresh failure is still there: expiry is about age, not kind alone.
    fresh = attention.emit_attention("run_failed", title="new failure")
    assert fresh in {ev["id"] for ev in attention.list_pending()}


def test_expired_notices_leave_the_file_at_the_next_write():
    attention.emit_attention("run_finished", title="stale")
    _age(attention._load_events()[0]["id"], 30 * 86400)
    assert len(attention._load_events()) == 1        # still on disk
    assert attention.list_pending() == []            # but not reported
    res = attention.prune_inbox()
    assert (res["ok"], res["removed"], res["kept"]) == (True, 1, 0)
    assert attention._load_events() == []


def test_pending_notices_are_capped():
    """A scheduler firing every few minutes fills a TTL window on its own,
    and every read of the inbox parses the whole file."""
    now = time.time()
    events = [{
        "id": f"att-notice-{i}", "kind": "run_finished", "session_id": "",
        "title": f"run {i}", "detail": "", "options": None, "workspace": "",
        "created_at": now - i, "status": "pending", "answer": None,
        "resolved_at": None, "acknowledged": False,
    } for i in range(attention._MAX_PENDING_NOTICES + 50)]
    attention._write_events(events)
    question = attention.emit_attention("question_pending", title="keep me")

    kept = attention._load_events()
    notices = [ev for ev in kept if ev["kind"] == "run_finished"]
    assert len(notices) == attention._MAX_PENDING_NOTICES
    assert question in {ev["id"] for ev in kept}     # cap spares the ask
    # The oldest went first.
    assert "att-notice-0" in {ev["id"] for ev in kept}
    assert f"att-notice-{attention._MAX_PENDING_NOTICES + 49}" not in {
        ev["id"] for ev in kept}


def test_clearing_a_large_inbox_rewrites_the_file_once():
    """The old clear resolved item by item, each a full read and rewrite:
    893 passes over a 352 kB file for one gesture."""
    for i in range(50):
        attention.emit_attention("run_finished", title=f"run {i}")
    writes = {"n": 0}
    real_write = attention._write_events

    def _counting(events):
        writes["n"] += 1
        real_write(events)

    attention._write_events = _counting
    try:
        assert attention.clear_all()["cleared"] == 50
    finally:
        attention._write_events = real_write
    assert writes["n"] == 1


# ---------------------------------------------------------------------------
# 2. The answer must survive a changed session id
# ---------------------------------------------------------------------------

def test_answer_reaches_a_session_that_has_a_different_id():
    """Park under one id, answer, drain under another — the whole point.

    The id is re-minted on a new cycle, cleared when the CLI backend
    restarts, and overwritten mid-turn by the backend's own id, so exact
    matching stranded the answer forever.
    """
    parked = attention.emit_attention(
        "question_pending", session_id="cycle-1-uuid",
        title="Which functional?", workspace="/ws")
    assert attention.answer_item(parked, "B3LYP")["ok"] is True

    drained = attention.drain_resolved("cycle-2-uuid", workspace="/ws")
    assert [(ev["id"], ev["answer"]) for ev in drained] == [(parked, "B3LYP")]
    assert attention.drain_resolved("cycle-2-uuid", workspace="/ws") == []


def test_answer_reaches_a_session_with_no_id_at_all():
    """The CLI backend runs with an empty id until the stream supplies one."""
    parked = attention.emit_attention(
        "question_pending", session_id="minted-uuid", title="continue?")
    attention.answer_item(parked, "yes")
    assert [ev["id"] for ev in attention.drain_resolved("")] == [parked]


def test_undelivered_answer_is_visible_to_user_and_doctor():
    """Resolved-but-unacknowledged was invisible everywhere: the listing
    filtered on ``pending``, so the user saw an empty inbox and believed
    the agent had it."""
    parked = attention.emit_attention(
        "question_pending", session_id="s1", title="Which basis set?")
    attention.answer_item(parked, "def2-TZVP")

    waiting = attention.list_undelivered()
    assert [ev["id"] for ev in waiting] == [parked]

    text = attention.render_inbox()
    assert "Answered, not yet delivered (1)" in text
    assert parked in text

    from delfin.agent import doctor
    rows = doctor._check_attention({"workspace": "", "settings": {},
                                    "fast": True})
    answers = [r for r in rows if r["check"] == "attention answers"]
    assert len(answers) == 1
    assert answers[0]["status"] == "WARN"

    # Once a session has it, both surfaces go quiet again.
    assert len(attention.drain_resolved("any-session")) == 1
    assert attention.list_undelivered() == []
    assert "not yet delivered" not in attention.render_inbox()


# ---------------------------------------------------------------------------
# 3. Retention must not delete an answer on its way to the agent
# ---------------------------------------------------------------------------

def test_retention_cap_never_drops_an_undelivered_answer():
    """The cap counted unacknowledged records and dropped the OLDEST
    first — precisely the long-parked question the user finally answered.
    """
    parked = attention.emit_attention(
        "question_pending", session_id="s1", title="answered at last")
    attention.answer_item(parked, "def2-TZVP")
    _age(parked, 86400)                       # older than every filler

    now = time.time()
    attention._write_events(attention._load_events() + [{
        "id": f"att-old-{i}", "kind": "run_finished", "session_id": "",
        "title": f"f{i}", "detail": "", "options": None, "workspace": "",
        "created_at": now, "status": "resolved", "answer": None,
        "resolved_at": now, "acknowledged": True,
    } for i in range(attention._MAX_NON_PENDING + 10)])

    attention.emit_attention("run_finished", title="one more, prune runs")

    stored = {ev["id"] for ev in attention._load_events()}
    assert parked in stored
    assert [ev["answer"] for ev in attention.drain_resolved("s1")] \
        == ["def2-TZVP"]


# ---------------------------------------------------------------------------
# 4. A failed cluster job belongs in the inbox
# ---------------------------------------------------------------------------

def test_failed_slurm_job_enters_the_inbox(monkeypatch):
    """The product's core case wrote to its own file and a popup only, so
    the inbox stayed empty and the doctor said "no pending events"."""
    from delfin.agent import job_monitor as jm
    monkeypatch.setattr(jm, "monitor_settings", lambda *a, **k: {
        "enabled": True, "interval_s": 600, "auto_diagnose": False,
        "webhook_url": "", "provider": "", "model": "", "backend": "",
    })
    finding = jm.Finding(
        job_id="4711", folder="/scratch/run", state="OUT_OF_MEMORY",
        signatures=["out-of-memory"], diagnosis_session="monitor-4711-1",
        summary="the job needed more memory than requested")
    jm.announce(finding)

    (event,) = attention.list_pending()
    assert event["kind"] == "run_failed"
    assert "4711" in event["title"]
    assert "OUT_OF_MEMORY" in event["title"] or "OUT_OF_MEMORY" in event["detail"]
    assert "out-of-memory" in event["detail"]
    assert event["workspace"] == "/scratch/run"
    assert event["delivered"] == ["desktop"]      # what actually went out


def test_announce_still_notifies_when_the_inbox_cannot_be_written(monkeypatch):
    """Durability is added on top of the popup, never in place of it."""
    from delfin.agent import job_monitor as jm
    from delfin.agent import notify

    sent: list[tuple] = []
    monkeypatch.setattr(notify, "send_notification",
                        lambda *a, **k: sent.append(a) or True)
    monkeypatch.setattr(jm, "monitor_settings", lambda *a, **k: {
        "enabled": True, "interval_s": 600, "auto_diagnose": False,
        "webhook_url": "", "provider": "", "model": "", "backend": "",
    })

    def _boom(*_a, **_k):
        raise OSError("disk full")

    monkeypatch.setattr(attention, "emit_attention", _boom)
    jm.announce(jm.Finding(job_id="9", folder="/x", state="FAILED"))
    assert len(sent) == 1


def test_announce_does_not_notify_twice(monkeypatch):
    """The inbox drives the same two transports the monitor already uses."""
    from delfin.agent import job_monitor as jm
    from delfin.agent import notify

    calls: list[str] = []
    monkeypatch.setattr(
        notify, "send_notification",
        lambda *a, **k: calls.append("desktop") or True)
    monkeypatch.setattr(
        notify, "send_remote_trigger",
        lambda payload, **k: calls.append(str(payload.get("event")))
        or notify.TriggerResult(sent=True, status_code=200))
    monkeypatch.setattr(jm, "monitor_settings", lambda *a, **k: {
        "enabled": True, "interval_s": 600, "auto_diagnose": False,
        "webhook_url": "https://example.invalid/hook",
        "provider": "", "model": "", "backend": "",
    })
    jm.announce(jm.Finding(job_id="7", folder="/x", state="TIMEOUT"))
    assert calls == ["desktop", "job_failed"]
    (event,) = attention.list_pending()
    assert sorted(event["delivered"]) == ["desktop", "webhook"]


# ---------------------------------------------------------------------------
# 5. Containment events must outlive the process that recorded them
# ---------------------------------------------------------------------------

def test_blocked_containment_event_becomes_a_durable_alert():
    """A ring buffer of 200 in the process that blocked is not a surface:
    the scheduler and the monitor record theirs inside child processes no
    dashboard is attached to."""
    from delfin.agent import security_events as se
    se.clear()
    try:
        se.record("denied_again", "bash", "rm -rf /scratch (second attempt)")
        (event,) = attention.list_pending()
        assert event["kind"] == "security_alert"
        assert "Refusal circumvented" in event["title"]
        assert "bash" in event["title"]
        assert "rm -rf /scratch" in event["detail"]
    finally:
        se.clear()


def test_security_alerts_are_deduped_per_kind_and_tool():
    from delfin.agent import security_events as se
    se.clear()
    try:
        for i in range(20):
            se.record("locked_scope_bash", "bash", f"attempt {i}")
        se.record("locked_scope_bash", "write_file", "different tool")
        se.record("secret_path", "bash", "~/.ssh/id_rsa")
        titles = [ev["title"] for ev in attention.list_pending()]
        assert len(titles) == 3
    finally:
        se.clear()


def test_panel_only_and_allowed_events_raise_nothing():
    """A refusal the USER chose, an expired approval window (the confirm
    broker parks that one itself) and an allowed action stay in the panel.
    """
    from delfin.agent import security_events as se
    se.clear()
    try:
        se.record("denied_by_user", "bash", "user clicked deny")
        se.record("approval_timeout", "bash", "window expired")
        se.record("egress", "web_fetch", "https://example.org", blocked=False)
        assert attention.list_pending() == []
        assert se.counts()["total"] == 3
    finally:
        se.clear()


def test_security_alert_failure_never_breaks_the_gate(monkeypatch):
    from delfin.agent import security_events as se
    se.clear()

    def _boom(*_a, **_k):
        raise OSError("disk full")

    monkeypatch.setattr(attention, "emit_attention", _boom)
    try:
        se.record("self_mod", "edit_file", "delfin/agent/api_client.py")
        assert se.counts()["blocked"] == 1     # the panel still has it
    finally:
        se.clear()


# ---------------------------------------------------------------------------
# 6. The budget kinds need an emitter
# ---------------------------------------------------------------------------

def _budget_engine(tmp_path, *, spent: float, budget: float = 10.0):
    """An engine object with only the budget state the notice reads."""
    from delfin.agent.engine import AgentEngine
    eng = object.__new__(AgentEngine)
    eng._budget_attention_levels = set()
    eng.session_id = "sid"
    eng.repo_dir = tmp_path
    eng.cost_usd = spent
    eng.run_budget_usd = budget
    eng.run_budget_s = 0.0
    eng._run_started_at = 0.0
    eng._unpriced_turns = 0
    eng._measured_cost_turns = 3
    eng._non_billing_turns = 0
    eng._unmeasured_budget_notice_shown = False
    eng.token_usage = {"input": 0, "output": 0}
    return eng


def test_wind_down_threshold_raises_a_budget_event(tmp_path):
    """A ceiling that only ever reached the model through a prompt block.

    Overnight the run winds down and stops; the desktop stays silent, the
    inbox stays empty and the doctor reports PASS.
    """
    from delfin.agent.engine import AgentEngine
    eng = _budget_engine(tmp_path, spent=8.5)
    block = AgentEngine._build_budget_block(eng)
    assert "wind-down" in block

    (event,) = attention.list_pending()
    assert event["kind"] == "budget_warning"
    assert "85%" in event["detail"] and "$10.00" in event["detail"]

    # Every following turn rebuilds the same block; one event is enough.
    AgentEngine._build_budget_block(eng)
    AgentEngine._build_budget_block(eng)
    assert len(attention.list_pending()) == 1


def test_exhausted_budget_raises_its_own_event(tmp_path):
    from delfin.agent.engine import AgentEngine
    eng = _budget_engine(tmp_path, spent=8.5)
    AgentEngine._build_budget_block(eng)
    eng.cost_usd = 10.5
    block = AgentEngine._build_budget_block(eng)
    assert "EXHAUSTED" in block
    titles = [ev["title"] for ev in attention.list_pending()]
    assert len(titles) == 2
    assert any("winding down" in t for t in titles)
    assert any("exhausted" in t.lower() for t in titles)


def test_budget_event_survives_a_broken_inbox(tmp_path, monkeypatch):
    from delfin.agent.engine import AgentEngine

    def _boom(*_a, **_k):
        raise OSError("disk full")

    monkeypatch.setattr(attention, "emit_attention", _boom)
    eng = _budget_engine(tmp_path, spent=9.9)
    assert "wind-down" in AgentEngine._build_budget_block(eng)


# ---------------------------------------------------------------------------
# 7. Transport health — the amplifier for everything above
# ---------------------------------------------------------------------------

def test_event_records_which_transports_delivered(_fake_home, monkeypatch):
    from delfin.agent import notify
    monkeypatch.setattr(notify, "send_notification", lambda *a, **k: False)
    eid = attention.emit_attention("question_pending", title="anyone there?")
    (stored,) = [ev for ev in attention._load_events() if ev["id"] == eid]
    assert stored["delivered"] == []

    monkeypatch.setattr(notify, "send_notification", lambda *a, **k: True)
    eid = attention.emit_attention("question_pending", title="and now?")
    (stored,) = [ev for ev in attention._load_events() if ev["id"] == eid]
    assert stored["delivered"] == ["desktop"]


def test_transport_status_reports_a_headless_login_node(_fake_home,
                                                        monkeypatch):
    from delfin.agent import notify
    monkeypatch.setattr(
        notify, "notification_backend",
        lambda: {"available": False, "backend": "",
                 "detail": "notify-send not on PATH"})
    st = attention.transport_status()
    assert st["usable"] == []
    assert "notify-send" in st["detail"]["desktop"]
    assert "not configured" in st["detail"]["webhook"]
    assert "not configured" in st["detail"]["hook"]


def test_transport_status_sees_a_configured_hook(_fake_home, monkeypatch):
    (_fake_home / ".delfin_settings.json").write_text(json.dumps({
        "agent": {"attention": {"notify_command": "/bin/echo --to phone",
                                "desktop": False}},
    }), encoding="utf-8")
    st = attention.transport_status()
    assert st["usable"] == ["hook"]
    assert "disabled" in st["detail"]["desktop"]

    (_fake_home / ".delfin_settings.json").write_text(json.dumps({
        "agent": {"attention": {
            "notify_command": "/does/not/exist --loud", "desktop": False}},
    }), encoding="utf-8")
    assert attention.transport_status()["usable"] == []


def test_notification_backend_needs_more_than_the_binary(monkeypatch):
    """``notify-send`` on a login node exits non-zero and shows nothing."""
    from delfin.agent import notify
    monkeypatch.setattr(notify, "_platform", lambda: "linux")
    monkeypatch.setattr(notify.shutil, "which",
                        lambda name: "/usr/bin/notify-send")
    for var in ("DISPLAY", "WAYLAND_DISPLAY", "DBUS_SESSION_BUS_ADDRESS"):
        monkeypatch.delenv(var, raising=False)
    probe = notify.notification_backend()
    assert probe["available"] is False
    assert "no desktop session" in probe["detail"]

    monkeypatch.setenv("DISPLAY", ":0")
    assert notify.notification_backend()["available"] is True


def test_send_notification_reports_a_refusing_notifier(monkeypatch):
    """``check=False`` turned every failure into a reported success."""
    from delfin.agent import notify

    class _Proc:
        returncode = 1

    monkeypatch.setattr(notify, "_platform", lambda: "linux")
    monkeypatch.setattr(notify.shutil, "which", lambda name: "/usr/bin/x")
    monkeypatch.setattr(notify.subprocess, "run", lambda *a, **k: _Proc())
    assert _REAL_SEND_NOTIFICATION("t", "b") is False

    _Proc.returncode = 0
    assert _REAL_SEND_NOTIFICATION("t", "b") is True


# ---------------------------------------------------------------------------
# The morning after — the whole scenario, end to end
# ---------------------------------------------------------------------------

def test_the_night_the_user_was_away(tmp_path, monkeypatch):
    """18:00 to 08:00, then the answer reaching the agent."""
    from delfin.agent import doctor, job_monitor as jm, security_events as se
    from delfin.agent.engine import AgentEngine

    monkeypatch.setattr(jm, "monitor_settings", lambda *a, **k: {
        "enabled": True, "interval_s": 600, "auto_diagnose": False,
        "webhook_url": "", "provider": "", "model": "", "backend": "",
    })
    se.clear()
    try:
        # ... the scheduler files its routine notices, night after night.
        for day in range(3):
            for i in range(20):
                eid = attention.emit_attention(
                    "run_finished", title=f"scheduled run {day}-{i} finished")
                _age(eid, (day + 3) * 86400)

        # 23:40 — the agent asks and blocks. Its own timeout message tells
        # the user the question can still be answered in the inbox.
        question = attention.emit_attention(
            "question_pending", session_id="session-of-the-evening",
            title="Which basis set for the follow-up?",
            options=["def2-SVP", "def2-TZVP"], workspace=str(tmp_path))
        # 01:12 — a cluster job dies.
        jm.announce(jm.Finding(job_id="4711", folder=str(tmp_path),
                               state="OUT_OF_MEMORY",
                               signatures=["out-of-memory"]))
        # 02:05 — the model retries something the gate already refused.
        se.record("denied_again", "bash", "curl secrets.example.invalid")
        # 03:00 — the run reaches its ceiling and stops.
        AgentEngine._build_budget_block(_budget_engine(tmp_path, spent=9.2))

        # 08:00. Nothing from three days ago is left, and everything that
        # happened tonight is here.
        text = attention.render_inbox()
        assert "1 blocking the agent" in text
        assert "scheduled run 0-0" not in text          # expired
        kinds = {ev["kind"] for ev in attention.list_pending()}
        assert kinds == {"question_pending", "run_failed", "security_alert",
                         "budget_warning"}

        # The doctor says the same thing.
        rows = {r["check"]: r for r in doctor.run_doctor(str(tmp_path))}
        assert rows["attention inbox"]["status"] == "WARN"
        assert "blocking" in rows["attention inbox"]["detail"]

        # The user clears the noise before reading — the ask survives.
        cleared = attention.clear_all()
        assert cleared["kept"] == 1
        assert [ev["id"] for ev in attention.list_pending()] == [question]

        # They answer it. The agent comes back on a NEW session id.
        assert attention.answer_item(question[:12], "def2-TZVP")["ok"] is True
        assert "Answered, not yet delivered (1)" in attention.render_inbox()
        drained = attention.drain_resolved("session-of-the-morning",
                                           workspace=str(tmp_path))
        assert [ev["answer"] for ev in drained] == ["def2-TZVP"]
        assert attention.drain_resolved("session-of-the-morning") == []
        assert attention.list_pending() == []
        assert attention.list_undelivered() == []
    finally:
        se.clear()
