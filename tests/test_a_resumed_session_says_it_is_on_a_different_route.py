"""``route`` was exported by every save and read by no restore.

``export_state`` writes it. ``restore_state`` re-derives the route from
the mode via ``_load_mode`` and never looks at the saved value. Harmless
while the pack holds still — and a field that looks carried and is not is
the shape this codebase keeps finding bugs in, so it gets one answer or
the other.

Reinstating the saved route would be wrong: a role that the pack no
longer defines cannot run, and resuming into it would fail later and
further from the cause. The current pack decides what runs.

What was missing is the sentence. The exported value is the only record
of the route the session ACTUALLY ran, so a pack edited between two
sessions is exactly the case where it is worth something — and a restore
that quietly continues on a different chain of roles is a difference the
person resuming should be told about, not one they discover from the
output.
"""

from __future__ import annotations

import pytest

from delfin.agent.engine import AgentEngine


def _blank_engine():
    """An engine with the declared fields reset and nothing else — the
    same construction the other resume tests use."""
    eng = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(eng, spec.attr, spec.reset())
    eng.mode, eng.route, eng.messages = "solo", [], []
    eng.session_id, eng.client = "", None
    return eng


def _restore(saved_route, current_route):
    eng = _blank_engine()
    eng.route = list(current_route)
    return eng.restore_state({
        "mode": "solo", "route": list(saved_route),
        "engine_messages": [], "session_id": "s1",
    })


# ---------------------------------------------------------------------------
# The difference is named
# ---------------------------------------------------------------------------

def test_a_route_that_changed_under_the_session_is_reported():
    report = _restore(["a", "b", "c"], ["a", "c"])
    assert not report.complete
    joined = " ".join(report.failed)
    assert "route" in joined
    assert "a→b→c" in joined          # what the session ran
    assert "a→c" in joined            # what it will run now


def test_a_pack_that_lost_the_mode_entirely_is_reported():
    report = _restore(["a", "b"], [])
    assert "nothing" in " ".join(report.failed)


def test_the_summary_carries_it():
    report = _restore(["a", "b"], ["a"])
    assert "route" in report.summary()


# ---------------------------------------------------------------------------
# An unchanged route is not a complaint
# ---------------------------------------------------------------------------

def test_an_unchanged_route_restores_quietly():
    report = _restore(["a", "b"], ["a", "b"])
    assert report.failed == ()
    assert "route" in report.restored


def test_a_session_saved_without_a_route_says_nothing_about_one():
    """Older session files have no such key, and inventing a complaint
    from an absent field would be worse than the silence it replaces."""
    eng = _blank_engine()
    eng.route = ["a"]
    report = eng.restore_state({"mode": "solo", "engine_messages": [],
                                "session_id": "s1"})
    assert not any("route" in f for f in report.failed)
    assert "route" not in report.restored


def test_the_pack_still_decides_what_runs():
    """Named, not reinstated: resuming into roles the pack no longer
    defines would fail later and further from the cause."""
    eng = _blank_engine()
    eng.route = ["a", "c"]
    eng.restore_state({"mode": "solo", "route": ["a", "b", "c"],
                       "engine_messages": [], "session_id": "s1"})
    assert eng.route == ["a", "c"]


# ---------------------------------------------------------------------------
# It survives the round trip it is exported for
# ---------------------------------------------------------------------------

def test_export_still_carries_the_route():
    eng = _blank_engine()
    eng.route = ["session_manager", "builder_agent"]
    assert eng.export_state()["route"] == ["session_manager", "builder_agent"]


@pytest.mark.parametrize("saved", [None, [], ()])
def test_an_empty_saved_route_is_not_a_drift(saved):
    eng = _blank_engine()
    eng.route = ["a"]
    report = eng.restore_state({"mode": "solo", "route": saved,
                                "engine_messages": [], "session_id": "s"})
    assert not any("route" in f for f in report.failed)
