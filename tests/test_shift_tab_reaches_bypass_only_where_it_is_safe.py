"""shift+tab's ladder: default -> acceptEdits -> plan -> bypassPermissions
-> default, in the REPL and in the window (the worker's session.setMode).

Bypass is on the ladder only where the managed policy allows it and a
filesystem layer holds each command; elsewhere the key steps over it. The
first step onto it in a session switches nothing: it is offered, and only
a deliberate ``y`` (REPL) or a setMode marked ``confirmed`` (the window's
y) switches. A line that is not y, or no answer, leaves the posture where
it was. After a yes the key passes through bypass without asking again.
"""

from __future__ import annotations

import io

import pytest

from delfin.agent import posture_cycle as PC
from delfin.agent import repl, repl_render as rr

PLAIN = rr.Theme(enabled=False)


class _Perms:
    def __init__(self, mode="plan"):
        self.mode = mode
        self.unisolated_fallback = False


class _Engine:
    def __init__(self, mode="plan"):
        self.kit_permissions = _Perms(mode)
        self.token_usage = {"input": 0, "output": 0}
        self.client = None
        self.switches: list[str] = []

    def set_kit_permission_mode(self, mode, announce=True):
        self.kit_permissions.mode = mode
        self.switches.append(mode)
        return True


class _Tty(io.StringIO):
    def isatty(self):
        return True


@pytest.fixture
def caged(monkeypatch):
    """A host whose layer holds every command, and no policy: the cage
    bypass relies on, the same on every machine."""
    from delfin.agent import isolation_posture as IP
    from delfin.agent import managed_policy as MP
    state = {"isolated": True, "policy": None}
    monkeypatch.setattr(IP, "isolation_status", lambda perms=None: IP.IsolationStatus(
        state["isolated"], layer="bwrap" if state["isolated"] else "",
        reason="" if state["isolated"] else "no layer works on this host"))
    monkeypatch.setattr(MP, "mode_refusal", lambda mode: state["policy"] if mode == PC.BYPASS else None)
    monkeypatch.setattr(PC, "_DELIBERATE_S", 0.0)
    return state


def _agent(engine):
    err = _Tty()
    agent = repl.TerminalAgent(engine, out=io.StringIO(), err=err)
    agent.transcript.theme = PLAIN
    return agent, err


# -- the ladder ---------------------------------------------------------------

def test_the_order_with_the_cage(caged):
    assert PC.ladder(_Engine()) == ["default", "acceptEdits", "plan", "bypassPermissions"]


def test_a_running_turn_steps_over_plan(caged):
    # main's rule (mode_cycle): while a turn runs the ladder passes over
    # plan, so a tool call landing mid-cycle is not refused.
    eng = _Engine()
    eng.turn_in_flight = True
    assert PC.step(eng, "acceptEdits", confirmed=True) == ("bypassPermissions", False)
    assert PC.step(eng, "plan", confirmed=True)[0] == "bypassPermissions"
    eng.turn_in_flight = False
    assert PC.step(eng, "acceptEdits", confirmed=True) == ("plan", False)


def test_bypass_is_stepped_over_without_isolation(caged):
    caged["isolated"] = False
    engine = _Engine()
    assert PC.ladder(engine) == ["default", "acceptEdits", "plan"]
    assert "isolation is not active" in PC.bypass_refusal(engine)


def test_bypass_is_stepped_over_when_the_policy_disallows_it(caged):
    caged["policy"] = "bypassPermissions is not available: the administrator's managed policy disallows it"
    assert PC.BYPASS not in PC.ladder(_Engine())


def test_bypass_is_stepped_over_in_the_unisolated_fallback(caged):
    engine = _Engine()
    engine.kit_permissions.unisolated_fallback = True
    assert PC.BYPASS not in PC.ladder(engine)


# -- the REPL -----------------------------------------------------------------

def test_the_repl_walks_default_accept_plan(caged):
    caged["isolated"] = False
    engine = _Engine("default")
    agent, _err = _agent(engine)
    seen = []
    for _ in range(4):
        agent._cycle_mode()
        seen.append(engine.kit_permissions.mode)
    assert seen == ["acceptEdits", "plan", "default", "acceptEdits"]


def test_the_first_step_onto_bypass_asks_and_switches_nothing(caged):
    engine = _Engine("plan")
    agent, err = _agent(engine)
    agent._cycle_mode()
    assert engine.kit_permissions.mode == "plan" and engine.switches == []
    assert "y and enter" in err.getvalue()
    assert agent._answer_bypass_offer("y") is True
    assert engine.kit_permissions.mode == "bypassPermissions"


def test_an_unanswered_offer_stays_on_the_previous_mode(caged):
    engine = _Engine("plan")
    agent, err = _agent(engine)
    agent._cycle_mode()
    # Another line is no answer: it is handled as what it is.
    assert agent._answer_bypass_offer("look at the logs") is False
    assert engine.kit_permissions.mode == "plan"
    assert "still plan" in err.getvalue()
    # shift+tab again, the offer open: on past bypass, never onto it.
    agent._cycle_mode()
    agent._cycle_mode()
    assert engine.kit_permissions.mode == "default"
    assert "bypassPermissions" not in engine.switches


def test_a_y_sooner_than_the_pause_is_no_answer(caged, monkeypatch):
    monkeypatch.setattr(PC, "_DELIBERATE_S", 60.0)
    engine = _Engine("plan")
    agent, _err = _agent(engine)
    agent._cycle_mode()
    assert agent._answer_bypass_offer("y") is True       # taken, so not sent as a prompt
    assert engine.kit_permissions.mode == "plan"


def test_after_a_yes_the_ladder_passes_through_bypass(caged):
    engine = _Engine("plan")
    agent, _err = _agent(engine)
    agent._cycle_mode()
    agent._answer_bypass_offer("y")
    agent._cycle_mode()
    assert engine.kit_permissions.mode == "default"
    for _ in range(3):
        agent._cycle_mode()
    assert engine.kit_permissions.mode == "bypassPermissions"
    assert agent._bypass_offer is None


def test_the_cage_gone_between_offer_and_yes_refuses(caged):
    engine = _Engine("plan")
    agent, err = _agent(engine)
    agent._cycle_mode()
    caged["isolated"] = False
    agent._answer_bypass_offer("y")
    assert engine.kit_permissions.mode == "plan"
    assert "isolation is not active" in err.getvalue()


def test_entering_bypass_is_recorded(caged, monkeypatch):
    from delfin.agent import security_events
    seen = []
    monkeypatch.setattr(security_events, "record", lambda kind, tool, detail, blocked=True: seen.append((kind, detail)))
    engine = _Engine("plan")
    assert PC.enter_bypass(engine, via="test") is None
    assert seen and seen[0][0] == "permission_mode" and "bypassPermissions" in seen[0][1]


# -- the window (the worker's session.setMode) ---------------------------------

class _Worker:
    """The worker's setMode and the facts it announces, on a stub engine."""

    def __init__(self, engine):
        from delfin.agent.service import worker as W
        self.w = W.Worker.__new__(W.Worker)
        self.w.engine = engine
        self.w.agent = type("A", (), {"post": lambda self, fn: None})()
        self.w._tasks, self.w._tasks_key = [], ()
        self.w.workspace = "."
        self.w._bypass_confirmed = False

    def set(self, mode, **extra):
        return self.w._m_session_setMode({"mode": mode, **extra})


def test_the_window_switches_to_bypass_only_when_confirmed(caged):
    from delfin.agent.service.worker import _Refusal
    engine = _Engine("plan")
    w = _Worker(engine)
    with pytest.raises(_Refusal):
        w.set("bypassPermissions")
    assert engine.kit_permissions.mode == "plan"
    assert w.set("bypassPermissions", confirmed=True)["approval"] == "bypassPermissions"
    w.set("default")
    # Confirmed once in this session: the ladder passes through.
    assert w.set("bypassPermissions")["approval"] == "bypassPermissions"


def test_the_window_cannot_reach_bypass_the_policy_or_the_host_rules_out(caged):
    from delfin.agent.service.worker import _Refusal
    caged["isolated"] = False
    engine = _Engine("plan")
    w = _Worker(engine)
    with pytest.raises(_Refusal) as refused:
        w.set("bypassPermissions", confirmed=True)
    assert "isolation is not active" in str(refused.value)
    assert engine.kit_permissions.mode == "plan"


def test_the_window_is_told_whether_the_ladder_may_stop_on_bypass(caged):
    engine = _Engine("plan")
    w = _Worker(engine)
    facts = w.w._session_facts()
    assert facts["bypassBlocked"] == "" and facts["bypassConfirmed"] is False
    caged["isolated"] = False
    assert "isolation" in w.w._session_facts()["bypassBlocked"]


def test_the_session_info_a_window_opens_with_carries_both(caged):
    """session.create's answer and a replay header: info(), not only the
    updates (a window that is told nothing steps over bypass)."""
    w = _Worker(_Engine("plan")).w
    w.emitter = type("E", (), {"session_id": "s1"})()
    w.title, w.banner, w.notices = "", [], []
    info = w.info()
    assert info["bypassBlocked"] == "" and info["bypassConfirmed"] is False
