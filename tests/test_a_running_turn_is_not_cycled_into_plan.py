"""Cycling modes during a running turn does not land in plan; a restored
pending plan is announced once.

Measured on the shared bug archive, 2026-09-16 .. 10-14: 15 calls refused
by plan mode in 9 sessions, and 10 "A pending plan was restored from the
saved session" notices in 2 sessions, all of them run in all_free. In one
session the user cycled bypass -> plan -> default -> diff -> acceptEdits
-> bypass with Shift+Tab while a turn ran; each step took effect at once,
a push landed while the session stood in plan and was refused, and the
model stopped with "I'm back in plan mode".

  idle                    the cycle is unchanged, plan included
  a turn running          the cycle passes over plan; an explicit
                          /mode plan or the dropdown still reaches it
  restored pending plan   announced on the first restore of that plan,
                          not on every later one
"""

from __future__ import annotations

import io
import re
from pathlib import Path

import pytest

from delfin.agent import repl, repl_render as rr


ORDER = ("plan", "default", "diff_approval", "acceptEdits",
         "bypassPermissions")


def test_idle_the_cycle_is_unchanged():
    from delfin.agent.mode_cycle import next_in_cycle
    assert next_in_cycle(ORDER, "bypassPermissions",
                         turn_running=False) == "plan"
    assert next_in_cycle(ORDER, "plan", turn_running=False) == "default"


def test_during_a_turn_the_cycle_passes_over_plan():
    from delfin.agent.mode_cycle import next_in_cycle
    assert next_in_cycle(ORDER, "bypassPermissions",
                         turn_running=True) == "default"
    # Walking the whole ring mid-turn never lands in plan.
    cur, seen = "default", []
    for _ in range(8):
        cur = next_in_cycle(ORDER, cur, turn_running=True)
        seen.append(cur)
    assert "plan" not in seen


def test_an_unknown_mode_starts_the_ring_safely():
    from delfin.agent.mode_cycle import next_in_cycle
    assert next_in_cycle(ORDER, "weird", turn_running=True) == "default"


def test_a_ring_of_plan_alone_is_left_alone():
    from delfin.agent.mode_cycle import next_in_cycle
    assert next_in_cycle(("plan",), "plan", turn_running=True) == "plan"


class _Perms:
    def __init__(self, mode):
        self.mode = mode


class _Engine:
    def __init__(self, mode, running):
        self.kit_permissions = _Perms(mode)
        self.token_usage = {"input": 0, "output": 0}
        self.client = None
        self.turn_in_flight = running

    def set_kit_permission_mode(self, mode):
        self.kit_permissions.mode = mode


class _Tty(io.StringIO):
    def isatty(self):
        return True


def _agent(engine):
    agent = repl.TerminalAgent(engine, out=io.StringIO(), err=_Tty())
    agent.transcript.theme = rr.Theme(enabled=False)
    return agent


def test_the_terminal_cycle_skips_plan_while_a_turn_runs():
    engine = _Engine("acceptEdits", running=True)
    agent = _agent(engine)
    agent._cycle_mode()
    assert engine.kit_permissions.mode == "default"


def test_the_terminal_cycle_reaches_plan_when_idle():
    engine = _Engine("acceptEdits", running=False)
    agent = _agent(engine)
    agent._cycle_mode()
    assert engine.kit_permissions.mode == "plan"


def test_the_engine_says_whether_a_turn_runs(tmp_path):
    from unittest.mock import MagicMock, patch
    from delfin.agent import engine as E
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                            model="kit.qwen3.5-397b-A17b", mode="solo")
    assert eng.turn_in_flight is False
    eng._turn_in_flight = True
    assert eng.turn_in_flight is True


def _tab_source() -> str:
    import delfin.dashboard.tab_agent as T
    src = Path(T.__file__).read_text(encoding="utf-8")
    return "\n".join(l for l in src.splitlines()
                     if not l.lstrip().startswith("#"))


def test_both_dashboard_cycles_use_the_turn_aware_ring():
    """The chip button and Shift+Tab (/perm-cycle) both step through
    next_in_cycle with the engine's turn state."""
    code = _tab_source()
    chip = code[code.index("def _cycle_kit_mode"):]
    chip = chip[:chip.index("\n    def ", 10)]
    assert "next_in_cycle(" in chip and "turn_running=" in chip
    perm = code[code.index('if cmd == "/perm-cycle":'):]
    perm = perm[:perm.index("return True")]
    assert "next_in_cycle(" in perm and "turn_running=" in perm


# -- a restored pending plan --------------------------------------------------

@pytest.fixture
def store(tmp_path, monkeypatch):
    from delfin.agent import session_store as S
    monkeypatch.setattr(S, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(S, "_ensure_dir",
                        lambda: (tmp_path / "sessions").mkdir(
                            parents=True, exist_ok=True)
                        or tmp_path / "sessions")
    monkeypatch.setattr(S, "acquire_session_lock", lambda sid: None)
    return S


def test_a_pending_plan_is_announced_once(store):
    S = store
    S.save_session("s1", pending_plan_body="1. do x\n2. do y")
    assert S.pending_plan_needs_notice(S.load_session("s1")) is True
    assert S.mark_pending_plan_noticed("s1") is True
    assert S.pending_plan_needs_notice(S.load_session("s1")) is False
    # A later save of the same plan keeps the mark.
    S.save_session("s1", pending_plan_body="1. do x\n2. do y")
    assert S.pending_plan_needs_notice(S.load_session("s1")) is False


def test_a_new_plan_is_announced_again(store):
    S = store
    S.save_session("s1", pending_plan_body="1. do x")
    S.mark_pending_plan_noticed("s1")
    S.save_session("s1", pending_plan_body="1. do z instead")
    assert S.pending_plan_needs_notice(S.load_session("s1")) is True


def test_no_plan_needs_no_notice(store):
    S = store
    S.save_session("s1", pending_plan_body="")
    assert S.pending_plan_needs_notice(S.load_session("s1")) is False


def test_the_restore_announces_through_the_mark():
    code = _tab_source()
    at = code.index("A pending plan was restored")
    window = code[max(0, at - 1500):at]
    assert "pending_plan_needs_notice(" in window
    assert "mark_pending_plan_noticed(" in code[at - 1500:at + 800]
    assert re.search(r"pending_plan_needs_notice\(", code)
