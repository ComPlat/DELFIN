"""Every turn-level improvement reaches the dashboard, not only the terminal.

Waves 12 and 13 were driven from the terminal, so their features were wired
into repl.py only: the announced-work follow-up, pause/resume, and the fresh
restart over the token budget. Most users work in the dashboard. This file
checks that both surfaces call the same hooks, and that the shared fresh-start
helper behaves the same for either caller.
"""
from __future__ import annotations

import inspect

from delfin.agent import repl
from delfin.dashboard import tab_agent

# hook -> what it is; each must appear in BOTH surfaces' source.
_HOOKS = {
    "pending_turn_continuation": "announced-work follow-up (wave 12, R1)",
    "pause_key": "delfin-agent pause reaches the session (wave 12, R1)",
    "wake_blocked": "no turn starts on its own while paused (wave 12, R1)",
}


def test_both_surfaces_wire_every_turn_hook():
    terminal = inspect.getsource(repl)
    dashboard = inspect.getsource(tab_agent)
    missing = {h: why for h, why in _HOOKS.items()
               if h not in terminal or h not in dashboard}
    assert not missing, f"wired on one surface only: {missing}"


def test_the_dashboard_uses_the_shared_fresh_restart():
    assert "prepare_fresh_turn" in inspect.getsource(tab_agent)


class _Engine:
    def __init__(self, n_messages, chars):
        self.messages = [{"role": "user", "content": "x" * chars}
                         for _ in range(n_messages)]
        self.session_id = ""
        self.token_usage = {"input": 0, "output": 0}
        self.context_budget = 1_000
        self.start_fresh = False


def test_under_budget_the_prompt_is_sent_unchanged():
    eng = _Engine(1, 40)
    assert repl.prepare_fresh_turn(eng, "next step") == "next step"
    assert not eng.start_fresh


def test_over_budget_the_turn_starts_fresh_with_the_prompt_kept(monkeypatch):
    archived = []
    monkeypatch.setattr(repl, "_archive_cut_history",
                        lambda engine, history: archived.append(len(history)))
    eng = _Engine(40, 4_000)                      # far over 1_000 tokens
    out = repl.prepare_fresh_turn(eng, "next step")
    assert out.endswith("next step") and out != "next step"
    assert eng.start_fresh
    assert archived == [40]


def test_a_broken_engine_degrades_to_the_plain_prompt():
    assert repl.prepare_fresh_turn(object(), "p") == "p"
