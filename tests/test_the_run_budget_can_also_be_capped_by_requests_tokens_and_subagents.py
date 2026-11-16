"""The optional run-budget dimensions: requests, tokens, sub-agents.

Each cap is OFF by default and counted from state that already exists:
requests from the per-session request metrics, tokens from the
session-cumulative token counter PLUS every delegated run's tokens
(counted towards the parent's cap), sub-agents from the delegated-spend
ledger's session count. ``_run_budget_detail`` names the dimension that
is furthest spent, exactly as for the USD / wall-clock / turn caps.
"""

from __future__ import annotations

import pytest


class _StandIn:
    """Only the attributes the budget code reads."""

    def __init__(self, **kw):
        self.cost_usd = 0.0
        self.run_budget_usd = 0.0
        self.run_budget_s = 0.0
        self.run_budget_turns = 0
        self._run_started_at = 0.0
        self._unpriced_turns = 0
        self._measured_cost_turns = 0
        self._non_billing_turns = 0
        self._unmeasured_budget_notice_shown = False
        self._budget_attention_levels = set()
        self.token_usage = {"input": 0, "output": 0}
        self.session_id = "sid"
        self.repo_dir = ""
        for k, v in kw.items():
            setattr(self, k, v)


def _ledger(engine):
    from delfin.agent.budget_ledger import BudgetLedger
    return BudgetLedger(engine)


def _with_delegates(eng, count=0, t_in=0, t_out=0, usd=0.0):
    from delfin.agent.agent_metrics import DelegateLedger, DelegateSpend
    led = DelegateLedger()
    led.session = DelegateSpend(count=count, input_tokens=t_in,
                                output_tokens=t_out, cost_usd=usd)
    eng._delegate_spend = led
    return eng


def test_the_request_cap_counts_the_sessions_model_requests(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    from delfin.agent import turn_metrics
    entries = [{"purpose": "", "input_tokens": 100, "output_tokens": 10},
               {"purpose": "", "input_tokens": 200, "output_tokens": 20},
               {"purpose": "classifier", "input_tokens": 1,
                "output_tokens": 1}]
    monkeypatch.setattr(turn_metrics, "read_requests",
                        lambda session, **kw: entries)
    led = _ledger(_StandIn(run_budget_requests=2, session_id="s1"))
    frac, exhausted, dim = led._run_budget_detail()
    assert (frac, dim) == (pytest.approx(1.0), "agent.run_budget_requests")
    assert exhausted is True


def test_the_request_cap_holds_an_empty_session_at_zero(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    from delfin.agent import turn_metrics
    monkeypatch.setattr(turn_metrics, "read_requests", lambda s, **kw: [])
    led = _ledger(_StandIn(run_budget_requests=2))
    # A dimension at 0 spent names no dim, the same convention the turn
    # cap follows (share 0.0 is never "> frac 0.0").
    assert led._run_budget_detail() == (0.0, False, "")


def test_the_request_cap_reads_the_attribute_before_the_setting(monkeypatch):
    """Same precedence as the other three dimensions."""
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"run_budget_requests": 40}})
    led = _ledger(_StandIn())
    assert led._run_budget_requests() == 40
    led._engine.run_budget_requests = 7
    assert led._run_budget_requests() == 7


def test_the_token_cap_counts_parent_and_delegated_tokens(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    eng = _with_delegates(
        _StandIn(run_budget_tokens=1000,
                 token_usage={"input": 300, "output": 200}),
        count=1, t_in=400, t_out=100)
    frac, exhausted, dim = _ledger(eng)._run_budget_detail()
    assert (frac, dim) == (pytest.approx(1.0), "agent.run_budget_tokens")
    assert exhausted is True


def test_the_token_cap_reads_the_attribute_before_the_setting(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"run_budget_tokens": 5000}})
    led = _ledger(_StandIn())
    assert led._run_budget_tokens() == 5000
    led._engine.run_budget_tokens = 900
    assert led._run_budget_tokens() == 900


def test_the_token_cap_tolerates_a_missing_counter(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    eng = _StandIn(run_budget_tokens=1000)
    del eng.token_usage                     # engines built bare have none
    assert _ledger(eng)._run_budget_tokens_spent() == 0


def test_the_subagent_cap_counts_the_delegates_towards_the_parent(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    eng = _with_delegates(_StandIn(run_budget_subagents=3), count=3,
                          t_in=100, t_out=100)
    frac, exhausted, dim = _ledger(eng)._run_budget_detail()
    assert (frac, exhausted, dim) == (pytest.approx(1.0), True,
                                      "agent.run_budget_subagents")


def test_the_subagent_cap_holds_without_delegates(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    eng = _StandIn(run_budget_subagents=3)
    assert _ledger(eng)._run_budget_detail() == (0.0, False, "")


def test_the_worst_dimension_is_reported_across_all_six(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    from delfin.agent import turn_metrics
    monkeypatch.setattr(turn_metrics, "read_requests", lambda s, **kw: [])
    eng = _with_delegates(
        _StandIn(run_budget_usd=10.0, run_budget_turns=10,
                 run_budget_requests=50, run_budget_tokens=100_000,
                 run_budget_subagents=10, _measured_cost_turns=2,
                 token_usage={"input": 20_000, "output": 10_000}),
            count=2, t_in=5_000, t_out=2_000)
    # usd 0/10, turns 0.2, requests 0, tokens 0.37, subagents 0.2 -> tokens
    frac, exhausted, dim = _ledger(eng)._run_budget_detail()
    assert dim == "agent.run_budget_tokens"
    assert frac == pytest.approx((20_000 + 10_000 + 5_000 + 2_000) / 100_000)
    assert exhausted is False


def test_a_parent_with_two_children_stops_at_the_subagent_cap(monkeypatch):
    """The task's scenario: a fake parent with two delegates.

    The parent has 2 of its 2-delegate cap used by two children; the
    detail tuple says the run is stopped on the sub-agent cap, and the
    budget block instructs the agent to wind down with no new work.
    """
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    eng = _with_delegates(_StandIn(run_budget_subagents=2,
                                   token_usage={"input": 100, "output": 50}),
                          count=2, t_in=900, t_out=300)
    led = _ledger(eng)
    frac, exhausted, dim = led._run_budget_detail()
    assert (frac, exhausted, dim) == (pytest.approx(1.0), True,
                                      "agent.run_budget_subagents")
    assert led._run_budget_status() == (pytest.approx(1.0), True)
    eng._emit_budget_attention = lambda *a, **k: None
    block = led._build_budget_block()
    assert block.startswith("# Run budget EXHAUSTED")
    assert "No new work." in block
    # The stop reason the turn gate records is the setting's short name,
    # the same convention the other dimensions follow.
    assert dim.rsplit(".", 1)[-1] == "run_budget_subagents"
