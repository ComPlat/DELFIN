"""The progress watchdog (hg1): a stall's class, its age, and the report.

Input: ``delfin.agent.progress_watchdog.classify`` over activity
snapshots (stand-ins with phase/phase_s, on a scripted clock — no real
sleeps), and ``ProgressWatchdog.tick`` driven per tick.
Asserted: each waiting phase classifies to its own stall class with the
duration and the safe action per class (approval-wait keeps waiting;
model-wait reports and leaves the retry to the existing deadline;
gate-queue outranks model-wait when the gate's wait notice is up;
tool-running never proposes a kill); idle is not a stall; a tick emits
one stream-json event and one status line, and stays silent while
nothing changes; a per-class limit marks limit_crossed; the settings
resolver parses ``agent.watchdog`` defensively.
"""

from __future__ import annotations

import pytest

from delfin.agent import progress_watchdog as pw


class _Snap:
    def __init__(self, phase, phase_s=0.0):
        self.phase = phase
        self.phase_s = phase_s


class _Clock:
    def __init__(self):
        self.t = 1000.0

    def __call__(self):
        return self.t


# ---------------------------------------------------------------------------
# classify: the phase -> class mapping
# ---------------------------------------------------------------------------

def test_an_approval_wait_classifies_and_keeps_waiting():
    row = pw.classify(_Snap("person", phase_s=780.0))
    assert row == {"stall": True, "class": "approval-wait",
                   "duration_s": 780.0, "action": pw.KEEP_WAITING}


def test_a_model_wait_classifies_with_the_retry_action():
    row = pw.classify(_Snap("model_wait", phase_s=120.0))
    assert row["class"] == "model-wait"
    assert row["stall"] and row["duration_s"] == 120.0
    assert row["action"] == pw.REPORT_MODEL_RETRY


def test_a_gate_queue_outranks_the_plain_model_wait():
    row = pw.classify(_Snap("model_wait", phase_s=300.0),
                      gate_notice="waiting for a slot on https://x: 2 in use")
    assert row["class"] == "gate-queue"
    assert row["stall"] and row["action"] == pw.KEEP_WAITING


def test_an_empty_gate_notice_is_no_gate_queue():
    for notice in ("", "   "):
        row = pw.classify(_Snap("model_wait", 5.0), gate_notice=notice)
        assert row["class"] == "model-wait"


def test_a_running_tool_is_a_stall_that_is_never_killed():
    row = pw.classify(_Snap("tool", phase_s=900.0))
    assert row == {"stall": True, "class": "tool-running",
                   "duration_s": 900.0, "action": pw.KEEP_WAITING}
    assert "kill" not in str(row["action"])


@pytest.mark.parametrize("phase", ["idle", "done", "", "streaming"])
def test_non_waiting_phases_are_not_stalls(phase):
    row = pw.classify(_Snap(phase, phase_s=50.0))
    assert row == {"stall": False, "class": "idle",
                   "duration_s": 50.0, "action": pw.NONE}


def test_a_broken_snapshot_classifies_as_no_stall():
    assert pw.classify(None)["stall"] is False
    assert pw.classify(object())["stall"] is False


# ---------------------------------------------------------------------------
# status line + stream-json event
# ---------------------------------------------------------------------------

def test_the_status_line_names_the_class_and_the_clock():
    line = pw.status_line({"stall": True, "class": "approval-wait",
                           "duration_s": 780.0, "action": pw.KEEP_WAITING})
    assert "approval-wait" in line and "13 min" in line


def test_the_event_shape_is_replacement_not_append():
    ev = pw.event({"stall": True, "class": "gate-queue",
                   "duration_s": 90.0, "action": pw.KEEP_WAITING})
    assert ev["type"] == "stall"
    assert ev["class"] == "gate-queue" and ev["duration_s"] == 90.0
    assert ev["action"] == pw.KEEP_WAITING


def test_no_status_line_when_not_a_stall():
    assert pw.status_line({"stall": False, "class": "idle",
                           "duration_s": 0.0, "action": pw.NONE}) == ""


# ---------------------------------------------------------------------------
# ProgressWatchdog.tick: emit-on-entry, silence between reports, the limit
# ---------------------------------------------------------------------------

def test_the_first_tick_emits_and_quiet_ticks_do_not():
    clock = _Clock()
    wd = pw.ProgressWatchdog(now=clock, limits={"approval-wait": 0.0},
                             tick_s=30.0)
    ev, line = wd.tick(_Snap("person", 5.0))
    assert ev and ev["type"] == "stall" and line
    clock.t += 5
    assert wd.tick(_Snap("person", 10.0)) == (None, "")
    clock.t += 30
    ev2, line2 = wd.tick(_Snap("person", 40.0))
    assert ev2 and ev2["duration_s"] == 40.0 and line2


def test_a_class_change_emits_even_before_the_tick_interval():
    clock = _Clock()
    wd = pw.ProgressWatchdog(now=clock, tick_s=60.0)
    wd.tick(_Snap("model_wait", 10.0))
    clock.t += 2
    ev, line = wd.tick(_Snap("person", 3.0))
    assert ev["class"] == "approval-wait" and line


def test_a_duration_never_loses_its_maximum():
    clock = _Clock()
    wd = pw.ProgressWatchdog(now=clock, tick_s=1.0)
    wd.tick(_Snap("model_wait", 100.0))
    clock.t += 2
    # The phase clock reset (a new request) must not shorten the reported
    # stall of the class that is still open.
    ev, _ = wd.tick(_Snap("model_wait", 1.0))
    assert ev["duration_s"] == 100.0


def test_a_limit_marks_limit_crossed_and_zero_means_never():
    clock = _Clock()
    wd = pw.ProgressWatchdog(now=clock,
                             limits={"tool-running": 600.0}, tick_s=1.0)
    ev, _ = wd.tick(_Snap("tool", 500.0))
    assert ev and wd.limit_crossed is False
    clock.t += 2
    wd.tick(_Snap("tool", 900.0))
    assert wd.limit_crossed is True


def test_leaving_the_stall_resets_the_state():
    clock = _Clock()
    wd = pw.ProgressWatchdog(now=clock, tick_s=1.0)
    wd.tick(_Snap("tool", 900.0))
    assert wd.tick(_Snap("idle", 0.0)) == (None, "")
    clock.t += 1
    ev, _ = wd.tick(_Snap("tool", 5.0))
    assert ev["duration_s"] == 5.0, "the old stall must not leak in"


def test_the_limit_crossed_never_reports_a_kill_action():
    wd = pw.ProgressWatchdog(limits={"approval-wait": 60.0}, tick_s=1.0)
    ev, line = wd.tick(_Snap("person", 120.0))
    assert wd.limit_crossed is True
    assert ev["action"] == pw.KEEP_WAITING
    assert "kill" not in line.lower()


# ---------------------------------------------------------------------------
# settings parsing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cfg, key, expected", [
    (None, None, pw.DEFAULT_LIMITS_S),
    ({}, None, pw.DEFAULT_LIMITS_S),
    ({"watchdog": {"model-wait_s": 120}}, "model-wait", 120.0),
    ({"watchdog": {"approval-wait_s": "90"}}, "approval-wait", 90.0),
    ({"watchdog": {"tool-running_s": 0}}, "tool-running", 0.0),
    ({"watchdog": {"tool-running_s": -5}}, "tool-running", 0.0),
    ({"watchdog": {"tool-running_s": "soon"}}, "tool-running", 0.0),
    ({"watchdog": {"tool-running_s": True}}, "tool-running", 0.0),
])
def test_resolve_limits_parses_defensively(cfg, key, expected):
    limits = pw.resolve_limits(cfg)
    if key is None:
        assert limits == pw.DEFAULT_LIMITS_S
    else:
        assert limits[key] == expected


def test_resolve_tick_s_defaults_and_bounds():
    assert pw.resolve_tick_s(None) == 30.0
    assert pw.resolve_tick_s({"watchdog": {"tick_s": 5}}) == 5.0
    assert pw.resolve_tick_s({"watchdog": {"tick_s": "x"}}) == 30.0
    assert pw.resolve_tick_s({"watchdog": {"tick_s": -1}}) == 30.0
