"""Resuming reset two session counters to zero, and one of them was a budget.

``_run_started_at`` is the start of the wall-clock run budget --
``_run_budget_status`` measures ``time.time() - _run_started_at`` against
``run_budget_s`` and, at 100%, tells the agent to stop starting work. It
was in neither the export nor the declaration, so it was re-stamped to
"now" by the constructor on every resume. A run that was out of time got
its whole time budget back by being resumed, and could do that as often
as it liked: resume was a way to launder a time budget on exactly the
unattended runs the budget exists for.

``_last_outcome_cost`` is the baseline the next outcome's delta is
measured from: ``cost_delta = cost_usd - _last_outcome_cost``. Zero after
a resume means the first outcome record books the WHOLE session's spend
as that one turn -- on top of the records already written for the turns
that spent it, which the activity view then sums.

Both are exported now. The clock is carried as ELAPSED rather than as a
timestamp, and the restore sets the start to now minus that, so the same
clock keeps running across a process boundary.

The third thing a headless resume did not do: re-bind the task session
id. The dashboard patched ``kit_permissions.task_session_id`` from
outside after calling ``restore_state``; the headless path had no such
line, so a resumed unattended run filtered the task list on the id the
constructor had minted and saw zero open tasks -- for a session whose
whole point was to continue them.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent.engine import AgentEngine


class _Perms:
    def __init__(self):
        self.task_session_id = "a-throwaway-id"


def _engine(perms=None) -> AgentEngine:
    eng = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(eng, spec.attr, spec.reset())
    eng.mode = "solo"
    eng.route = []
    eng.messages = []
    eng.session_id = ""
    eng.client = None
    eng._perms_for_test = perms
    return eng


@pytest.fixture(autouse=True)
def _perms_property(monkeypatch):
    """kit_permissions is a read-only property on the real engine."""
    monkeypatch.setattr(
        AgentEngine, "kit_permissions",
        property(lambda self: getattr(self, "_perms_for_test", None)),
        raising=False)


# ---------------------------------------------------------------------------
# The wall clock
# ---------------------------------------------------------------------------

def test_the_run_clock_survives_a_resume():
    src = _engine()
    src._run_started_at = time.time() - 3600.0    # an hour into the run
    resumed = _engine()
    resumed.restore_state(src.export_state())
    elapsed = time.time() - resumed._run_started_at
    assert 3595.0 < elapsed < 3610.0, (
        f"the resumed run thinks it has been going {elapsed:.0f}s -- a "
        f"resume that resets the clock hands back the whole time budget")


def test_a_resumed_run_that_is_out_of_time_stays_out_of_time():
    """The consequence, at the level the budget is read."""
    src = _engine()
    src._run_started_at = time.time() - 7200.0
    resumed = _engine()
    resumed.restore_state(src.export_state())
    resumed.run_budget_usd = 0.0
    resumed.run_budget_s = 3600.0
    _frac, exhausted = resumed._run_budget_status()
    assert exhausted, "the run budget came back unspent after a resume"


def test_the_export_carries_elapsed_seconds_not_a_wall_clock_stamp():
    """A timestamp from another machine (or another timezone's clock skew)
    is not a quantity this can trust; the elapsed time is."""
    src = _engine()
    src._run_started_at = time.time() - 120.0
    state = src.export_state()
    assert 118.0 <= state["run_elapsed_s"] <= 125.0


def test_a_fresh_cycle_starts_the_clock_again():
    eng = _engine()
    eng._run_started_at = time.time() - 5000.0
    eng._reset_session_fields()
    assert time.time() - eng._run_started_at < 5.0


# ---------------------------------------------------------------------------
# The cost baseline
# ---------------------------------------------------------------------------

def test_the_cost_baseline_survives_a_resume():
    src = _engine()
    src.cost_usd = 4.25
    src._last_outcome_cost = 4.25       # every cent already booked
    resumed = _engine()
    resumed.restore_state(src.export_state())
    assert resumed._last_outcome_cost == 4.25, (
        "the next outcome would book the whole session's spend as one "
        "turn, on top of the records that already booked it")


def test_the_first_outcome_after_a_resume_books_only_the_new_spend():
    src = _engine()
    src.cost_usd = 4.25
    src._last_outcome_cost = 4.25
    resumed = _engine()
    resumed.restore_state(src.export_state())
    resumed.cost_usd = 4.30             # one more turn, five cents
    delta = max(0.0, resumed.cost_usd - resumed._last_outcome_cost)
    assert delta == pytest.approx(0.05)


# ---------------------------------------------------------------------------
# The task list the resumed session is supposed to continue
# ---------------------------------------------------------------------------

def test_a_restore_rebinds_the_task_session_id():
    perms = _Perms()
    eng = _engine(perms)
    eng.restore_state({"mode": "solo", "engine_messages": [],
                       "session_id": "the-real-session"})
    assert perms.task_session_id == "the-real-session", (
        "the resumed run filters its task list on the id the constructor "
        "minted, so every task it is meant to continue is invisible")


def test_a_restore_without_permissions_does_not_raise():
    eng = _engine(None)
    eng.restore_state({"mode": "solo", "engine_messages": [],
                       "session_id": "s1"})
    assert eng.session_id == "s1"
