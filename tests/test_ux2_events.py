"""UX2 phase-3 red tests: pause/resume emit one attention event; the
watch loop dedupes; hostile reasons and titles stay one safe line.

Every test redirects HOME (state_paths table: attention._inbox_path and
session_pause._DIR are both resolved per call) and never touches the real
~/.delfin. These events go out through attention.emit_attention only.
"""

import json
from pathlib import Path

import pytest

from delfin.agent import attention as attention_mod
from delfin.agent import session_health as health_mod
from delfin.agent import session_pause as pause_mod
from delfin.agent.state_paths import RedirectedUserState


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    with RedirectedUserState(tmp_path):
        yield tmp_path


def _inbox_events(home):
    path = home / ".delfin" / "attention_inbox.jsonl"
    if not path.exists():
        return []
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line:
            out.append(json.loads(line))
    return out


def _safe_kinds(home):
    return [e for e in _inbox_events(home)
            if e.get("kind") == "session_state"]


def _as_dicts(events):
    return [{"name": n, "title": t} for n, t in events]


# --- pause emits exactly one event, with the undo command ------------------

def test_first_pause_emits_one_event_with_undo(home):
    pause_mod.pause("nacht-s29", "anomaly brake (2 severe in 30 min)")
    events = _safe_kinds(home)
    assert len(events) == 1
    ev = events[0]
    assert "nacht-s29" in ev["title"]
    assert "anomaly brake (2 severe in 30 min)" in ev["title"]
    assert "delfin-agent resume nacht-s29" in ev["title"]
    assert ev.get("session_id") == "nacht-s29"


def test_second_pause_of_a_paused_session_emits_nothing(home):
    pause_mod.pause("nacht-s29", "brake")
    pause_mod.pause("nacht-s29", "brake again")
    assert len(_safe_kinds(home)) == 1


def test_resume_emits_one_event_only_when_a_flag_was_removed(home):
    pause_mod.pause("nacht-s29", "brake")
    pause_mod.resume("nacht-s29")
    pause_mod.resume("nacht-s29")            # unpaused: emits nothing
    pause_mod.resume("never-paused")         # unpaused: emits nothing
    events = _safe_kinds(home)
    resumed = [e for e in events if "resumed" in e["title"]]
    assert len(resumed) == 1
    assert "delfin-agent pause nacht-s29" in resumed[0]["title"]


def test_hostile_reason_stays_one_safe_line(home):
    reason = "brake\n**injected second line** <script>x</script>\x1b[31mred"
    pause_mod.pause("nacht-s29", reason)
    ev = _safe_kinds(home)[0]
    assert "\n" not in ev["title"]
    assert "\x1b" not in ev["title"]
    assert "injected second line" in ev["title"]


# --- the watch loop: one event per (session, verdict), deduped -------------

class _Clock:
    def __init__(self):
        self.t = 1_000_000.0

    def __call__(self):
        return self.t


def _rec(**kw):
    rec = {"name": "nacht-s7", "session_id": "id7", "paused": False,
           "pid_alive": True, "open_approvals": [], "approval_opened_at": None,
           "last_turn_end": None, "last_turn": {}, "idle_s": 0,
           "api_errors_in_a_row": 0, "resume_planned": False,
           "no_token_s": 0, "open_work": False, "title": "w",
           "started_at": None, "reason": None, "next": None}
    rec.update(kw)
    return rec


def _watch_state():
    return {"seen": {}}


def test_watch_emits_once_per_state_and_dedupes_across_ticks(home):
    clock = _Clock()
    record = _rec(paused=True)
    state = _watch_state()
    s1 = health_mod.watch_events([record], clock, state)
    assert len(s1) == 1
    s2 = health_mod.watch_events([record], clock, state)
    assert len(s2) == 0                       # dedup holds across two ticks


def test_watch_emits_again_after_the_session_leaves_and_reenters(home):
    clock = _Clock()
    state = _watch_state()
    rec = _rec(paused=True)
    health_mod.watch_events([rec], clock, state)
    rec2 = _rec()
    health_mod.watch_events([rec2], clock, state)   # leaves paused
    rec3 = _rec(paused=True)
    again = health_mod.watch_events([rec3], clock, state)
    assert len(again) == 1                    # re-entering emits again


def test_watch_emits_stop_all_set(home):
    clock = _Clock()
    stop = {"at": clock.t - 10, "reason": "emergency", "host": "h"}
    events = _as_dicts(health_mod.watch_events(
        [], clock, _watch_state(), stop_all=stop, stop_all_seen={}))
    assert len(events) == 1
    assert "stop-all" in events[0]["title"].lower()


def test_watch_does_not_reemit_the_same_stop_all(home):
    clock = _Clock()
    stop = {"at": clock.t - 10, "reason": "emergency", "host": "h"}
    seen = {}
    health_mod.watch_events([], clock, _watch_state(),
                            stop_all=stop, stop_all_seen=seen)
    events = health_mod.watch_events([], clock, _watch_state(),
                                     stop_all=stop, stop_all_seen=seen)
    assert len(events) == 0


def test_watch_reports_a_stop_all_clear_only_once(home):
    clock = _Clock()
    stop = {"at": clock.t - 10, "reason": "emergency", "host": "h"}
    seen = {}
    health_mod.watch_events([], clock, _watch_state(),
                            stop_all=stop, stop_all_seen=seen)
    # a fresh session started AFTER the stop: the stop is observed as over
    rec = _rec(paused=False, pid_alive=True,
               started_at=clock.t + 60)
    events = _as_dicts(health_mod.watch_events(
        [rec], clock.t + 120, _watch_state(),
        stop_all=stop, stop_all_seen=seen))
    cleared = [e for e in events if "cleared" in e["name"].lower()
               or "is over" in e["title"].lower()]
    assert len(cleared) == 1
    more = health_mod.watch_events([rec], clock.t + 240, _watch_state(),
                                   stop_all=stop, stop_all_seen=seen)
    assert not [e for e in more if "cleared" in str(e)]


# --- ops health through the CLI: read-only ----------------------------------

def _announce(key, session_id="id-s29", **kw):
    from delfin.agent import session_presence
    session_presence.announce(key, session_id=session_id,
                              title=kw.pop("title", "demo"),
                              workspace=str(kw.pop("workspace", "")),
                              **kw)


def test_ops_health_lists_a_line_and_changes_no_file(home, tmp_path):
    pause_mod.pause("nacht-s29", "brake")
    _announce("nacht-s29")
    from delfin.agent import cli
    code = cli.main(["ops", "health", "--json"])
    assert code == 0


def test_ops_health_json_marks_a_paused_session(home, capsys):
    pause_mod.pause("nacht-s29", "brake")
    _announce("nacht-s29")
    from delfin.agent import cli
    cli.main(["ops", "health", "--json"])
    out = capsys.readouterr().out
    rows = json.loads(out)
    assert isinstance(rows, list)
    names = [r.get("name") for r in rows]
    assert "nacht-s29" in names
    row = [r for r in rows if r.get("name") == "nacht-s29"][0]
    assert row["verdict"] == "paused"


def test_ops_health_text_marks_a_paused_session(home, capsys):
    pause_mod.pause("nacht-s29", "brake")
    _announce("nacht-s29")
    from delfin.agent import cli
    cli.main(["ops", "health"])
    out = capsys.readouterr().out
    assert "nacht-s29" in out
    assert "paused" in out


def test_ops_health_includes_a_paused_session_not_announced(home, capsys):
    # a session whose process died leaving the pause flag: still shown
    pause_mod.pause("nacht-s30", "brake")
    from delfin.agent import cli
    cli.main(["ops", "health", "--json"])
    rows = json.loads(capsys.readouterr().out)
    names = [r.get("name") for r in rows]
    assert "nacht-s30" in names
    row = [r for r in rows if r.get("name") == "nacht-s30"][0]
    assert row["verdict"] == "paused"


def test_ops_health_wait_subcommand_emits_the_events(home, capsys, monkeypatch):
    from delfin.agent import cli
    clock = _Clock()
    rec = _rec(paused=True)
    rec["name"] = "nacht-s29"
    rec["session_id"] = "id29"
    import delfin.agent.session_health as sh
    monkeypatch.setattr(sh, "pause_watch", lambda seconds: None)
    ticks = [0]

    def _collect(*a, **k):
        ticks[0] += 1
        if ticks[0] >= 2:
            raise KeyboardInterrupt()
        return [dict(rec, started_at=None)]

    monkeypatch.setattr(sh, "collect", _collect)
    code = cli.main(["ops", "health", "--wait", "--interval", "0"])
    assert code == 0
    events = _safe_kinds(home)
    assert len(events) == 1
    assert "nacht-s29" in events[0]["title"]
