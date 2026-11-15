"""The budget surface is the same after the move into BudgetLedger.

Characterization tests for the extraction of the run-budget methods from
``delfin/agent/engine.py`` into ``delfin/agent/budget_ledger.py``. They
pin TODAY'S behaviour on three fixed engine states (no budget; a USD
budget at 50%; the turn budget reached), on the detail tuples, on the
>110 % early return, and on the export/restore round trip of the budget
fields -- so the move shows up here if it changes anything.

They run against the ENGINE (the only surface a caller may use), plus
one direct test of the new module against a stand-in engine object, so
budget_ledger.py has to work on its own, without AgentEngine.
"""

from __future__ import annotations

import time

import pytest


# ---------------------------------------------------------------------------
# A minimal stand-in engine: only the attributes the budget code reads.
# ---------------------------------------------------------------------------

class _StandIn:
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


# ---------------------------------------------------------------------------
# 1. The new module works on its own, on a stand-in engine
# ---------------------------------------------------------------------------

def test_the_ledger_reads_the_settings_fallback_through_the_engine():
    eng = _StandIn()
    assert _ledger(eng)._run_budget() == (0.0, 0.0)
    eng.run_budget_usd = 10.0
    assert _ledger(eng)._run_budget() == (10.0, 0.0)


def test_the_ledger_reports_the_worst_dimension():
    eng = _StandIn(run_budget_usd=10.0, run_budget_turns=10,
                   _measured_cost_turns=9)
    assert _ledger(eng)._run_budget_detail() == (0.9, False,
                                                 "agent.run_budget_turns")
    eng.cost_usd = 9.5
    assert _ledger(eng)._run_budget_detail() == (0.95, False,
                                                 "agent.run_budget_usd")


def test_the_ledger_counts_every_turn_not_just_the_priced_ones():
    eng = _StandIn(_unpriced_turns=2, _non_billing_turns=1,
                   _measured_cost_turns=3)
    assert _ledger(eng)._turns_taken() == 6


def test_the_ledger_does_not_copy_engine_state():
    """The ledger reads through the reference: a later change on the
    engine is visible in the next read, and a write through the ledger's
    bookkeeping lands on the engine."""
    eng = _StandIn()
    led = _ledger(eng)
    eng.cost_usd = 4.0
    assert led._run_budget_detail()[0] == 0.0      # no budget configured
    eng.run_budget_usd = 8.0
    assert led._run_budget_detail()[0] == pytest.approx(0.5)
    eng._unpriced_turns = 1
    assert led._usd_budget_enforced() is False


# ---------------------------------------------------------------------------
# 2. The engine surface: same method names, same answers
# ---------------------------------------------------------------------------

def _engine_via_new(**kw) -> object:
    """An AgentEngine without __init__ (which opens a backend), holding
    exactly the state the budget methods read."""
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine.__new__(AgentEngine)
    eng.cost_usd = 0.0
    eng.run_budget_usd = 0.0
    eng.run_budget_s = 0.0
    eng.run_budget_turns = 0
    eng._run_started_at = 0.0
    eng._unpriced_turns = 0
    eng._measured_cost_turns = 0
    eng._non_billing_turns = 0
    eng._unmeasured_budget_notice_shown = False
    eng._budget_attention_levels = set()
    eng.token_usage = {"input": 0, "output": 0}
    eng.session_id = "sid"
    eng.repo_dir = ""
    for k, v in kw.items():
        setattr(eng, k, v)
    return eng


def test_engine_budget_detail_tuples_are_unchanged():
    eng = _engine_via_new(run_budget_usd=10.0, _measured_cost_turns=5,
                          cost_usd=5.0)
    assert eng._run_budget_detail() == (0.5, False, "agent.run_budget_usd")
    assert eng._run_budget_status() == (0.5, False)

    eng2 = _engine_via_new(run_budget_turns=10, _measured_cost_turns=10)
    assert eng2._run_budget_detail() == (1.0, True, "agent.run_budget_turns")

    eng3 = _engine_via_new()
    assert eng3._run_budget_detail() == (0.0, False, "")
    assert eng3._run_budget_status() == (0.0, False)


def test_engine_usd_budget_enforced_is_unchanged():
    eng = _engine_via_new(run_budget_usd=10.0)
    assert eng._usd_budget_enforced() is True
    eng._unpriced_turns = 1
    assert eng._usd_budget_enforced() is False
    assert _engine_via_new()._usd_budget_enforced() is False  # no budget


def test_engine_cost_hard_cap_reads_the_setting(monkeypatch):
    eng = _engine_via_new()
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"cost_hard_limit_usd": 1.5}})
    assert eng._cost_hard_cap() == 1.5
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})
    assert eng._cost_hard_cap() == 50.0


def test_engine_run_budget_turns_precedence_is_unchanged(monkeypatch):
    eng = _engine_via_new()
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"run_budget_turns": 6}})
    assert eng._run_budget_turns() == 6
    eng.run_budget_turns = 3
    assert eng._run_budget_turns() == 3          # the attribute wins
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"run_budget_turns": -5}})
    eng.run_budget_turns = 0
    assert eng._run_budget_turns() == 0          # negative is off


# ---------------------------------------------------------------------------
# 3. The rendered budget block for three fixed states
# ---------------------------------------------------------------------------

def test_the_budget_block_is_silent_without_a_budget():
    eng = _engine_via_new(cost_usd=999.0)
    assert eng._build_budget_block() == ""


def test_the_budget_block_at_50_percent_is_silent():
    eng = _engine_via_new(run_budget_usd=10.0, cost_usd=5.0,
                          _measured_cost_turns=3)
    assert eng._build_budget_block() == ""


WIND_DOWN = (
    "# Run budget wind-down (auto-injected)\n"
    "- 85% of this run's budget is spent. Start wrapping up NOW: complete "
    "or checkpoint the current step, commit/save state, mark tasks, and "
    "summarise open work. If a follow-up run makes sense, note exactly "
    "where to resume. Do not start new large work items."
)

EXHAUSTED = (
    "# Run budget EXHAUSTED\n"
    "- 105% of this run's budget is spent. Finish your CURRENT sentence "
    "of work only: save state, mark task statuses honestly, and summarise "
    "what remains. No new work."
)


def test_the_budget_block_wind_down_text_is_unchanged():
    eng = _engine_via_new(run_budget_usd=10.0, cost_usd=8.5,
                          _measured_cost_turns=3)
    assert eng._build_budget_block() == WIND_DOWN


def test_the_budget_block_exhausted_text_is_unchanged():
    eng = _engine_via_new(run_budget_usd=10.0, cost_usd=10.5,
                          _measured_cost_turns=3)
    assert eng._build_budget_block() == EXHAUSTED


def test_the_budget_block_names_every_ceiling_in_force():
    eng = _engine_via_new(run_budget_usd=10.0, run_budget_s=600.0,
                          run_budget_turns=20, cost_usd=10.5,
                          _measured_cost_turns=21)
    block = eng._build_budget_block()
    assert block.startswith("# Run budget EXHAUSTED")
    assert "$10.00" not in block                 # the ceiling list goes to
    assert block == EXHAUSTED                    # the attention detail, not


def test_the_budget_attention_goes_to_the_inbox_once_per_level():
    seen = []
    eng = _engine_via_new(run_budget_usd=10.0, cost_usd=8.5,
                          _measured_cost_turns=3)
    from unittest.mock import patch
    with patch("delfin.agent.attention.emit_attention",
               side_effect=lambda kind, **k: seen.append((kind, k))):
        eng._build_budget_block()
        eng._build_budget_block()
    assert len(seen) == 1
    kind, kw = seen[0]
    assert kind == "budget_warning"
    assert kw["title"] == "Run budget at 80% — the agent is winding down"
    assert "85% of the run budget ($10.00)" in kw["detail"]


# ---------------------------------------------------------------------------
# 4. The unmeasured-USD notice
# ---------------------------------------------------------------------------

def test_the_unmeasured_notice_text_is_unchanged():
    eng = _engine_via_new(run_budget_usd=10.0, _unpriced_turns=1,
                          _measured_cost_turns=0,
                          token_usage={"input": 1000, "output": 500})
    block = eng._build_budget_block()
    assert block.startswith("# Run budget NOT enforceable in USD "
                            "(auto-injected, once)")
    assert "- 1 of 1 turn(s) so far ran on a model with no published " \
           "rate" in block
    assert "- Measured instead: 1 turn(s), 1,000 input + 500 output " \
           "tokens." in block
    assert "For a ceiling that can be enforced, set " \
           "agent.run_budget_turns or agent.run_budget_s, or pick a " \
           "model id with a published rate." in block


def test_the_unmeasured_notice_says_what_is_in_force():
    eng = _engine_via_new(run_budget_usd=10.0, run_budget_turns=20,
                          _unpriced_turns=1)
    block = eng._build_budget_block()
    assert "- The 20 turn(s) ceiling IS measured and still ends this " \
           "run." in block


def test_the_unmeasured_notice_fires_once():
    eng = _engine_via_new(run_budget_usd=10.0, _unpriced_turns=1)
    assert "NOT enforceable" in eng._build_budget_block()
    assert eng._build_budget_block() == ""


# ---------------------------------------------------------------------------
# 5. The early return at the cap
# ---------------------------------------------------------------------------

def test_the_110_percent_gate_refuses_the_turn_with_the_right_dial():
    """stream_response's budget gate: past 110% no turn starts, and the
    refusal names the setting that stopped it."""
    from unittest.mock import MagicMock, patch
    from delfin.agent.engine import AgentEngine
    client = MagicMock()
    client.stream_message = MagicMock(
        side_effect=lambda *a, **k: iter(()))
    calls = client.stream_message
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine.__new__(AgentEngine)
    # The gate reads the session fields it needs; only the budget ones
    # and the turn bookkeeping matter for this assertion.
    for attr in ("cost_usd", "messages"):
        if not hasattr(eng, attr):
            setattr(eng, attr, 0.0 if attr == "cost_usd" else [])
    eng.run_budget_usd = 10.0
    eng.run_budget_s = 0.0
    eng.run_budget_turns = 4
    eng.cost_usd = 0.0
    eng._unpriced_turns = 0
    eng._measured_cost_turns = 5                 # 5/4 = 125%
    eng._non_billing_turns = 0
    # The gate needs a little more engine scaffolding; assert at the
    # level the gate itself uses instead of driving a whole turn:
    frac, exhausted, dim = eng._run_budget_detail()
    assert (frac, exhausted, dim) == (1.25, True, "agent.run_budget_turns")
    assert frac >= 1.1
    assert calls.call_count == 0                 # nothing ran to get here


# ---------------------------------------------------------------------------
# 6. export_state / restore_state round-trips the budget fields
# ---------------------------------------------------------------------------

def test_the_budget_fields_round_trip_through_the_session_state():
    from delfin.agent.engine import AgentEngine

    def _blank():
        eng = AgentEngine.__new__(AgentEngine)
        for spec in AgentEngine._SESSION_FIELDS:
            setattr(eng, spec.attr, spec.reset())
        return eng

    src = _blank()
    src.mode, src.route, src.messages = "solo", [], []
    src.session_id, src.client = "sid", None
    src.cost_usd = 2.5
    src._measured_cost_turns = 4
    src._unpriced_turns = 7
    src._non_billing_turns = 1
    src._run_started_at = time.time() - 600.0
    src._delegate_spend = None                   # a broken ledger must
    # survive the export the same way it always did

    resumed = _blank()
    resumed.mode, resumed.route, resumed.messages = "solo", [], []
    resumed.session_id, resumed.client = "", None
    resumed.restore_state(src.export_state())

    assert resumed.cost_usd == pytest.approx(2.5)
    assert resumed._measured_cost_turns == 4
    assert resumed._unpriced_turns == 7
    assert resumed._non_billing_turns == 1
    assert 595.0 < time.time() - resumed._run_started_at < 605.0
    # The gates agree the ceiling is still spent after the resume.
    resumed.run_budget_usd = 10.0
    frac, exhausted = resumed._run_budget_status()
    assert frac == pytest.approx(0.25)
    assert exhausted is False
    resumed.run_budget_usd = 0.0
    resumed.run_budget_turns = 10
    frac, exhausted = resumed._run_budget_status()
    assert exhausted is True                     # 12 turns / 10


# ---------------------------------------------------------------------------
# 7. The delegated-spend surface moves with the ledger
# ---------------------------------------------------------------------------

def test_the_delegated_spend_ledger_is_reached_through_the_engine(monkeypatch):
    import threading

    from delfin.agent import agent_metrics as am
    eng = _engine_via_new()
    eng._lock = threading.Lock()
    monkeypatch.setattr(type(eng), "kit_permissions", property(
        lambda self: getattr(self, "_kit_perms_for_test", None)))
    eng._kit_perms_for_test = None
    eng._delegate_spend = am.DelegateLedger()
    assert eng.delegate_spend() is eng._delegate_spend
    # A broken ledger is replaced, exactly as before.
    eng._delegate_spend = "nonsense"
    assert isinstance(eng.delegate_spend(), am.DelegateLedger)


def test_a_delegate_is_booked_to_the_turn_that_waited_for_it():
    import threading

    from delfin.agent import agent_metrics as am
    eng = _engine_via_new(provider="openai", backend="api")
    eng._lock = threading.Lock()
    eng._delegate_spend = am.DelegateLedger()
    eng._turn_in_flight = True
    eng._turn_id = 7
    eng._book_delegate_spend({"input_tokens": 100, "output_tokens": 10,
                              "model": "gpt-4.1"}, started_turn=7)
    spend = eng.delegate_spend()
    assert spend.turn.count == 1 and spend.session.count == 1
    assert spend.background.count == 0
    eng._turn_in_flight = False
    eng._book_delegate_spend({"input_tokens": 50, "output_tokens": 5,
                              "model": "gpt-4.1"}, started_turn=7)
    assert eng.delegate_spend().background.count == 1
    assert eng.delegate_spend().session.count == 2


def test_the_meter_wraps_the_runner_and_unwraps_itself(monkeypatch):
    import threading

    from delfin.agent import agent_metrics as am

    calls = []

    def runner(*args, **kw):
        calls.append(kw)
        return {"input_tokens": 10, "output_tokens": 1,
                "model": "gpt-4.1"}

    class _Perms:
        subagent_runner = runner

    eng = _engine_via_new(provider="openai", backend="api")
    eng._lock = threading.Lock()
    eng._delegate_spend = am.DelegateLedger()
    eng._turn_in_flight = True
    eng._turn_id = 1
    monkeypatch.setattr(type(eng), "kit_permissions", property(
        lambda self: getattr(self, "_kit_perms_for_test", None)))
    eng._kit_perms_for_test = _Perms()
    assert eng._meter_delegate_spend() is True
    eng.kit_permissions.subagent_runner(subagent_type="explore")
    assert eng.delegate_spend().turn.count == 1
    # Re-installing does not stack counters.
    assert eng._meter_delegate_spend() is True
    eng.kit_permissions.subagent_runner(subagent_type="explore")
    assert eng.delegate_spend().turn.count == 2
    assert eng.delegate_spend().session.count == 2


def test_a_backend_without_a_runner_has_no_meter(monkeypatch):
    eng = _engine_via_new()
    monkeypatch.setattr(type(eng), "kit_permissions", property(
        lambda self: getattr(self, "_kit_perms_for_test", None)))
    eng._kit_perms_for_test = None
    assert eng._meter_delegate_spend() is False


def test_the_price_state_lands_on_the_engine():
    from delfin.agent import pricing

    class _Client:
        model = "kit.qwen3.5-397b-A17b"

    eng = _engine_via_new(provider="kit", backend="api", client=_Client())
    eng._note_turn_price_state(0.0)
    assert eng._non_billing_turns == 1
    assert eng._last_turn_price_state == "non_billing"
    assert pricing.resolve("kit.qwen3.5-397b-A17b", "kit").state \
        == pricing.NON_BILLING


def test_an_observed_spend_counts_as_measured_even_unpriced():
    class _Client:
        model = "a-model-nobody-priced"

    eng = _engine_via_new(provider="openai", backend="api", client=_Client())
    eng._note_turn_price_state(0.42)
    assert eng._measured_cost_turns == 1
    assert eng._unpriced_turns == 0


def test_the_engine_budget_methods_keep_their_signatures():
    """The hand-over script (scratch/apply_engine.py) rewrites these
    methods into delegates; the def lines, parameter lists, return
    annotations and decorators must survive it exactly.  Pinned here at
    the ast level so a regenerated signature fails before the diff is
    reviewed."""
    import ast
    import delfin.agent.engine
    from pathlib import Path

    engine_path = Path(delfin.agent.engine.__file__)
    src = engine_path.read_text(encoding="utf-8")
    tree = ast.parse(src)
    cls = next(n for n in ast.walk(tree)
               if isinstance(n, ast.ClassDef) and n.name == "AgentEngine")
    methods = {n.name: n for n in cls.body
               if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))}
    names = ["_run_budget", "_run_ceiling_hit", "_run_budget_turns",
             "_turns_taken", "_note_turn_price_state",
             "_usd_budget_enforced", "_unmeasured_budget_block",
             "_run_budget_status", "_run_budget_detail",
             "_emit_budget_attention", "_build_budget_block",
             "_cost_hard_cap", "delegate_spend", "_meter_delegate_spend",
             "_book_delegate_spend"]
    for name in names:
        node = methods.get(name)
        assert node is not None, name
        # self is the first positional parameter of every method
        assert node.args.args and node.args.args[0].arg == "self", name
        # return annotations survive where they exist today
        expected = {
            "_run_budget": "tuple[float, float]",
            "_run_budget_turns": "int",
            "_run_budget_status": "tuple[float, bool]",
            "_run_budget_detail": "tuple[float, bool, str]",
            "_cost_hard_cap": "float",
            "delegate_spend": "'_agent_metrics.DelegateLedger'",
            "_meter_delegate_spend": "bool",
            "_book_delegate_spend": "None",
            "_note_turn_price_state": "None",
            "_usd_budget_enforced": "bool",
        }
        got = ast.unparse(node.returns) if node.returns else ""
        if name in expected:
            assert got == expected[name], (name, got)
        # keyword-only markers survive
        if name == "_book_delegate_spend":
            assert [x.arg for x in node.args.kwonlyargs] == ["started_turn"]
