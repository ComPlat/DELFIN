"""With no new budget setting present, _run_budget_detail is unchanged.

Characterization pin for the run-budget extension: the three OPTIONAL
dimensions (requests, tokens, sub-agents) default OFF, and with the new
settings absent the detail tuple is exactly what it was before they
existed -- same fraction, same exhausted flag, same dimension string,
for every state the existing pins cover plus the delegated-spend states.
"""

from __future__ import annotations


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


def test_no_new_setting_detail_tuple_is_todays(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    led = _ledger(_StandIn())
    assert led._run_budget_detail() == (0.0, False, "")
    assert led._run_budget_status() == (0.0, False)


def test_no_new_setting_usd_state_is_todays(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    led = _ledger(_StandIn(run_budget_usd=10.0, _measured_cost_turns=5,
                           cost_usd=5.0))
    assert led._run_budget_detail() == (0.5, False, "agent.run_budget_usd")


def test_no_new_setting_turn_state_is_todays(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    led = _ledger(_StandIn(run_budget_turns=10, _measured_cost_turns=10))
    assert led._run_budget_detail() == (1.0, True, "agent.run_budget_turns")


def test_no_new_setting_with_delegates_present_is_todays(monkeypatch):
    """A session that delegated must not read as spending a cap that is
    not configured."""
    from delfin.agent.agent_metrics import DelegateLedger, DelegateSpend
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    eng = _StandIn()
    led = DelegateLedger()
    led.session = DelegateSpend(count=2, input_tokens=5000, output_tokens=900)
    eng._delegate_spend = led
    assert _ledger(eng)._run_budget_detail() == (0.0, False, "")


def test_no_new_setting_turn_budget_unaffected_by_token_count(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    led = _ledger(_StandIn(run_budget_turns=4, _measured_cost_turns=3,
                           token_usage={"input": 90_000, "output": 9_000}))
    assert led._run_budget_detail() == (0.75, False, "agent.run_budget_turns")


def test_no_new_setting_settings_keys_read_only_the_old_three(monkeypatch):
    """The settings file may hold unrelated agent keys; nothing new fires
    unless one of the three new keys is present."""
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"run_budget_turns": 2},
                                 "other": 1})
    led = _ledger(_StandIn(_measured_cost_turns=2))
    assert led._run_budget_detail() == (1.0, True, "agent.run_budget_turns")


def test_the_wind_down_block_is_unchanged_without_new_settings(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    eng = _StandIn(run_budget_usd=10.0, cost_usd=8.5, _measured_cost_turns=3)
    eng._emit_budget_attention = lambda *a, **k: None
    block = _ledger(eng)._build_budget_block()
    assert block == (
        "# Run budget wind-down (auto-injected)\n"
        "- 85% of this run's budget is spent. Start wrapping up NOW: complete "
        "or checkpoint the current step, commit/save state, mark tasks, and "
        "summarise open work. If a follow-up run makes sense, note exactly "
        "where to resume. Do not start new large work items.")
