"""The ceiling that holds when the price cannot be measured.

``test_a_budget_gate_cannot_stop_a_cost_it_never_measured`` records why
every USD gate reads a bare 0 on a model with no published rate. It ends
by offering "the ceilings that are always measured -- turns and tokens".
Three field sessions then showed that the offer alone bounds nothing,
because the surface they ran on never reads it and the dimension it named
was off by default:

======================  ==============  ========  =============
model                   input tokens    cost_usd  outcome
======================  ==============  ========  =============
kit.glm-5.3-flash           27,333,309       0.0  finished
kit.glm-5.3-flash          124,117,533       0.0  finished
Blablador GLM-5.3-Flash    363,450,615       0.0  session corrupt
======================  ==============  ========  =============

515 M input tokens, every one of them recorded as costing nothing. The
broken run made 1,696 tool calls and compacted 34 times; its last turn
alone carried 1,031,096 input tokens over 10.7 minutes, which is also why
the wall clock the CLI warning used to offer is not a substitute: 10.7
minutes is a cheap time budget and a very expensive context.

What these tests pin:

* the token dimension is computed where the USD fraction is, so the
  >110 % refusal, the 80 % wind-down and the mid-turn check reach it;
* its default would have refused the broken session and left the two
  that worked alone;
* a provider that returns no input count makes the dimension report
  UNMEASURABLE and name itself -- never a silent 0, which is the defect
  the whole file exists for;
* every figure appears with its denominator.
"""

from __future__ import annotations

import textwrap
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import budget_ledger as BL
from delfin.agent.api_client import StreamEvent


# ---------------------------------------------------------------------------
# Stand-ins
# ---------------------------------------------------------------------------

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
        self._measured_token_turns = 0
        self._untokened_turns = 0
        self._untokened_provider = ""
        self._unmeasured_budget_notice_shown = False
        self._unmeasured_token_notice_shown = False
        self._budget_attention_levels = set()
        self.token_usage = {"input": 0, "output": 0}
        self.provider = "kit"
        self.session_id = "sid"
        self.repo_dir = ""
        for k, v in kw.items():
            setattr(self, k, v)

    def _emit_budget_attention(self, *a, **k):
        pass

    def _unmeasured_budget_block(self):
        return BL.BudgetLedger(self)._unmeasured_budget_block()


def _ledger(engine):
    return BL.BudgetLedger(engine)


@pytest.fixture(autouse=True)
def _no_host_settings(monkeypatch):
    """Never read the host's own agent configuration: a ceiling the
    machine happens to configure would make these pins machine-specific.
    """
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {}})


@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# solo mode")
    (lite_dir / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


def _engine(agent_tree, *, model="opus", provider="claude",
            input_tokens=1_000, output_tokens=500):
    """An engine whose every turn reports the given token counts and a
    cost of 0.0 -- what the client hands over for an unpriced model."""
    from delfin.agent.engine import AgentEngine

    def stream(system, messages, max_tokens=4096, session_id="",
               thinking_budget=0, **kw):
        yield StreamEvent(type="text_delta", text="ok")
        yield StreamEvent(type="message_delta", input_tokens=input_tokens,
                          output_tokens=output_tokens, cost_usd=0.0)

    client = MagicMock()
    client.model = model
    client.stream_message = MagicMock(side_effect=stream)
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine(repo_dir=agent_tree, backend="api",
                          provider=provider, model=model,
                          mode="quick", pack_dir=agent_tree)
    eng.client = client
    return eng


# ---------------------------------------------------------------------------
# 1. The default, against the sessions it was derived from
# ---------------------------------------------------------------------------

#: The three field figures, as the reports' Context tables record them.
WORKED_SMALL = 27_333_309
WORKED_LARGE = 124_117_533
BROKE = 363_450_615


def test_the_default_refuses_the_session_that_broke():
    """>110 % of the ceiling is where no further turn starts."""
    led = _ledger(_StandIn(token_usage={"input": BROKE, "output": 505_054}))
    frac, exhausted, dim = led._run_budget_detail()
    assert dim == "agent.run_budget_tokens"
    assert exhausted is True
    assert frac >= 1.1, "the 363 M session must reach the hard refusal"


def test_the_default_leaves_the_two_sessions_that_worked_alone():
    """Not refused, not even wound down: 124 M is 62 % of 200 M."""
    for total in (WORKED_SMALL, WORKED_LARGE):
        led = _ledger(_StandIn(token_usage={"input": total, "output": 0}))
        frac, exhausted, _dim = led._run_budget_detail()
        assert exhausted is False, total
        assert frac < 0.8, (total, frac)
        assert led._build_budget_block() == "", total


def test_the_broken_session_is_wound_down_before_it_is_refused():
    """The 80 % band is what converts it into a saved run, not a dead
    one: 160 M of the 1,696 tool calls in, the agent is told to stop
    starting work."""
    led = _ledger(_StandIn(token_usage={"input": 170_000_000, "output": 0}))
    block = led._build_budget_block()
    assert block.startswith("# Run budget wind-down")
    assert "170.0 M of 200.0 M input+output tokens (85 %)" in block


def test_the_default_is_stated_where_it_is_derived():
    """A number in a gate with no basis beside it is a guess.

    An int carries no ``__doc__``, so the basis lives in the comment
    block above the assignment -- and this reads that block, because a
    requirement to state the basis that nothing checks is a requirement
    the next edit drops.
    """
    import pathlib
    assert BL.DEFAULT_RUN_BUDGET_TOKENS == 200_000_000
    source = pathlib.Path(BL.__file__).read_text(encoding="utf-8")
    basis = source[:source.index("DEFAULT_RUN_BUDGET_TOKENS =")]
    for figure in ("27.3 M", "124.1 M", "363.5 M", "2,328"):
        assert figure in basis, figure


# ---------------------------------------------------------------------------
# 2. It is a dial, and it can be lifted
# ---------------------------------------------------------------------------

def test_the_ceiling_is_on_without_any_setting():
    assert _ledger(_StandIn())._run_budget_tokens() == (
        BL.DEFAULT_RUN_BUDGET_TOKENS)


def test_an_explicit_zero_setting_lifts_it(monkeypatch):
    """0 in the settings is an answer, not an absent one."""
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"run_budget_tokens": 0}})
    led = _ledger(_StandIn())
    assert led._run_budget_tokens() == 0
    assert led._run_budget_detail() == (0.0, False, "")


def test_a_negative_attribute_lifts_it_for_this_run_only():
    eng = _StandIn(run_budget_tokens=-1)
    assert _ledger(eng)._run_budget_tokens() == 0
    eng.run_budget_tokens = 0
    assert _ledger(eng)._run_budget_tokens() == BL.DEFAULT_RUN_BUDGET_TOKENS


def test_the_attribute_still_beats_the_setting(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"run_budget_tokens": 5_000}})
    eng = _StandIn()
    assert _ledger(eng)._run_budget_tokens() == 5_000
    eng.run_budget_tokens = 900
    assert _ledger(eng)._run_budget_tokens() == 900


def test_a_run_stopped_on_tokens_is_resumable_by_raising_the_same_dial():
    """A refusal that sends its reader to the wrong flag is not a
    refusal they can act on."""
    from delfin.agent.run_limits import resume_hint
    hint = resume_hint(_StandIn(), "agent.run_budget_tokens")
    assert "--max-run-tokens" in hint


# ---------------------------------------------------------------------------
# 3. Mirrored, not parallel: the same three behaviours, same places
# ---------------------------------------------------------------------------

def test_the_fraction_is_computed_where_the_usd_fraction_is():
    """One function answers for every dimension, so a change to the USD
    term is a change to the token term."""
    import inspect
    src = inspect.getsource(BL.BudgetLedger._run_budget_detail)
    assert "agent.run_budget_usd" in src
    assert "agent.run_budget_tokens" in src


def test_the_110_percent_gate_refuses_a_turn_on_the_token_dimension(
        agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.token_usage["input"] = BROKE
    out = eng.stream_response("carry on")
    assert "Run budget exhausted (>110%)" in out
    assert "agent.run_budget_tokens" in out
    assert eng.last_turn_stop_reason == "run_budget_tokens"
    # ...and the client was never asked to produce another turn.
    assert eng.client.stream_message.call_count == 0


def test_the_wind_down_block_fires_on_the_token_dimension(agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.token_usage["input"] = 180_000_000
    block = eng._build_budget_block()
    assert "# Run budget wind-down" in block
    assert "180.0 M of 200.0 M input+output tokens (90 %)" in block


def test_a_running_turn_is_stopped_between_rounds_not_only_before_it(
        agent_tree):
    """The session that broke was ONE turn of 1,696 tool calls.

    A ceiling consulted only before a turn would have watched every one
    of them go past, which is why the wall clock and the spend are
    already checked between rounds. The token term is checked there too.
    """
    eng = _engine(agent_tree, model="opus")
    eng.token_usage["input"] = 199_000_000
    # Nothing yet: the finished turns are under the ceiling.
    assert eng._run_ceiling_hit(0.0, 0) == ""
    # The round that crosses it ends the turn, with the dial as the
    # reason, and the work of the rounds before it is kept.
    assert eng._run_ceiling_hit(0.0, 2_000_000) == "run_budget_tokens"


def test_the_client_passes_the_running_turns_tokens_to_the_probe():
    """A probe that is only ever handed 0 cannot stop a running turn."""
    from delfin.agent import api_client as A
    client = A.OpenAIClient.__new__(A.OpenAIClient)
    seen = []

    def probe(spent, tokens=0):
        seen.append((spent, tokens))
        return "run_budget_tokens" if tokens > 10 else ""

    client.run_ceiling = probe
    assert client._run_ceiling_reached(0.0, 11) == "run_budget_tokens"
    assert seen[-1] == (0.0, 11)


def test_an_older_one_argument_probe_still_works():
    """An engine or test double built against the previous signature
    must not have its turns ended by a TypeError."""
    from delfin.agent import api_client as A
    client = A.OpenAIClient.__new__(A.OpenAIClient)
    client.run_ceiling = lambda spent: "run_budget_usd" if spent > 1 else ""
    assert client._run_ceiling_reached(2.0, 999) == "run_budget_usd"
    assert client._run_ceiling_reached(0.0, 999) == ""


def test_the_tool_loop_reads_the_token_ceiling_inside_itself():
    """A cap nobody reads inside the loop is the defect being fixed, one
    dimension across."""
    import pathlib
    from delfin.agent import api_client as A
    src = pathlib.Path(A.__file__).read_text(encoding="utf-8")
    code = "\n".join(line for line in src.splitlines()
                     if not line.lstrip().startswith("#"))
    start = code.index("for _round in range(_MAX_TOOL_ROUNDS + 1):")
    window = code[start:start + 4000]
    assert "_turn_token_cap" in window
    # ...and in the same window as the cost cap it mirrors, so the two
    # cannot drift apart into one place that checks and one that does not.
    assert "_turn_cost_cap" in window


# ---------------------------------------------------------------------------
# 4. The per-turn breaker
# ---------------------------------------------------------------------------

def test_the_loop_can_ask_for_the_per_turn_token_cap():
    from delfin.agent import api_client as A
    client = A.OpenAIClient.__new__(A.OpenAIClient)
    assert client._turn_token_cap() == 0, "no probe installed means no cap"
    client.turn_token_cap = lambda: 1_234
    assert client._turn_token_cap() == 1_234


def test_a_broken_token_cap_probe_does_not_end_the_turn():
    from delfin.agent import api_client as A
    client = A.OpenAIClient.__new__(A.OpenAIClient)

    def boom():
        raise RuntimeError("engine gone")

    client.turn_token_cap = boom
    assert client._turn_token_cap() == 0


def test_the_engine_installs_the_per_turn_token_cap(agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.stream_response("hi")
    probe = getattr(eng.client, "turn_token_cap", None)
    assert callable(probe)
    assert probe() == BL.DEFAULT_TURN_TOKEN_CAP


def test_the_per_turn_cap_is_read_once_per_turn(agent_tree):
    """The same contract as the cost cap: snapshot at turn start, so a
    settings change mid-turn cannot move the ceiling under the turn."""
    eng = _engine(agent_tree, model="opus")
    eng.stream_response("hi")
    assert eng._token_cap_value == BL.DEFAULT_TURN_TOKEN_CAP
    assert eng._token_cap_hit is False
    assert eng._turn_start_input_tokens == 0


def test_the_per_turn_cap_can_be_lifted_or_lowered(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"token_hard_limit_in": 0}})
    assert _ledger(_StandIn())._token_hard_cap() == 0
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"token_hard_limit_in": 1_000}})
    assert _ledger(_StandIn())._token_hard_cap() == 1_000


def test_the_per_turn_cap_default_is_above_every_turn_that_worked():
    """19.4 M was the largest single turn in the three sessions, in the
    one that broke. A per-turn cap low enough to catch a 1 M turn would
    cut a quarter of the turns recorded on this host, so the run total is
    what separates a broken session from a working one and this cap is
    the runaway backstop only."""
    assert BL.DEFAULT_TURN_TOKEN_CAP == 50_000_000
    assert BL.DEFAULT_TURN_TOKEN_CAP > 19_416_805


def test_a_pre_send_figure_exists_but_is_an_estimate():
    """Requirement measured, not assumed.

    The assembled request CAN be measured before it goes out --
    ``estimate_schema_tokens`` is chars // 4 over the message list, which
    the context-length 400 handler already reports. It is a house
    estimate, so the breaker refuses on the counts the provider returned
    and therefore stops the NEXT round, not the one that crossed the cap.
    """
    from delfin.agent.tool_catalog import estimate_schema_tokens
    msgs = [{"role": "user", "content": "x" * 4_000}]
    assert estimate_schema_tokens(msgs[0]) == pytest.approx(1_000, rel=0.05)
    # ...and it is an estimate: the figure is derived from characters,
    # not from a tokeniser, so it may not be the number billed.
    assert estimate_schema_tokens({"role": "u", "content": "日本語" * 100}) > 0


# ---------------------------------------------------------------------------
# 5. Unavailable must never read as zero
# ---------------------------------------------------------------------------

def test_a_productive_turn_with_no_input_count_is_unmeasured_not_free():
    led = _ledger(_StandIn(provider="blablador"))
    led._note_turn_token_state(0, produced=True)
    assert led._engine._untokened_turns == 1
    assert led._engine._measured_token_turns == 0
    assert led._token_budget_enforced() is False


def test_a_turn_that_reported_its_size_is_measured():
    led = _ledger(_StandIn())
    led._note_turn_token_state(12_345, produced=True)
    assert led._engine._measured_token_turns == 1
    assert led._engine._untokened_turns == 0
    assert led._token_budget_enforced() is True


def test_a_turn_that_produced_nothing_claims_nothing_either_way():
    """21 of the three sessions' turns are 400s and empty stops. Counting
    them as unmeasured would raise the notice on every run that ever saw
    a failed request, which is the false statement in the other
    direction."""
    led = _ledger(_StandIn())
    led._note_turn_token_state(0, produced=False)
    assert led._engine._untokened_turns == 0
    assert led._engine._measured_token_turns == 0
    assert led._token_budget_enforced() is True


def test_the_notice_names_the_provider_that_left_the_count_out():
    eng = _StandIn(provider="blablador",
                   token_usage={"input": 1_000, "output": 10})
    led = _ledger(eng)
    led._note_turn_token_state(0, produced=True)
    block = led._unmeasured_token_block()
    assert "UNMEASURABLE" in block
    assert "blablador" in block
    assert "Their size is unknown, NOT zero." in block
    # ...and with its denominator, so the reader knows what the partial
    # figure is a part of.
    assert "1.0 k of 200.0 M (0 %)" in block


def test_the_unmeasurable_notice_is_said_once():
    eng = _StandIn(provider="kit")
    led = _ledger(eng)
    led._note_turn_token_state(0, produced=True)
    assert led._unmeasured_token_block() != ""
    assert led._unmeasured_token_block() == ""


def test_the_block_says_it_without_a_notice_being_due():
    eng = _StandIn()
    assert _ledger(eng)._unmeasured_token_block() == ""


def test_the_dimension_list_marks_the_partly_unmeasured_ceiling():
    eng = _StandIn(provider="kit", token_usage={"input": 1_000, "output": 0})
    led = _ledger(eng)
    led._note_turn_token_state(0, produced=True)
    line = next(ln for ln in led._budget_dimension_lines() if "tokens" in ln)
    assert line.endswith("— PARTLY UNMEASURED")


def test_the_engine_records_the_token_state_per_turn(agent_tree):
    """The engine's own turn end, not a hand-called helper."""
    eng = _engine(agent_tree, model="opus")
    eng.stream_response("hi")
    assert eng._measured_token_turns == 1
    assert eng._untokened_turns == 0
    assert eng._last_turn_token_state == "measured"


def test_a_turn_whose_provider_reports_no_tokens_is_caught_by_the_engine(
        agent_tree):
    eng = _engine(agent_tree, model="opus", input_tokens=0, output_tokens=0)
    eng.stream_response("hi")
    assert eng._untokened_turns == 1
    assert eng._token_budget_enforced() is False
    assert "UNMEASURABLE" in eng._build_budget_block()


# ---------------------------------------------------------------------------
# 5b. The zero that wore a measurement's clothes
#
# ``_note_turn_price_state`` counts a reported spend as measured whatever
# the rate card says, and that rule is right: a subscription's $0.00 is a
# real observation. It is also how the broken session escaped every USD
# gate -- it ran ``backend: cli``, so it counted ZERO unpriced turns,
# reported its ceiling as enforced, and would have shown
# ``$0.00 of $10.00 (0 %)`` through 363 M tokens. The two zeros are
# identical by value, so the dimension that IS measured decides.
# ---------------------------------------------------------------------------

def test_a_cli_backend_zero_after_a_material_run_is_not_a_measurement():
    """The field session's own shape: measured turns, cost 0.00, 363 M
    tokens. The pinned rule stands; the ceiling is not called enforced."""
    eng = _StandIn(run_budget_usd=10.0, cost_usd=0.0, backend="cli",
                   _measured_cost_turns=114,
                   token_usage={"input": BROKE, "output": 505_054})
    led = _ledger(eng)
    assert led._usd_budget_enforced() is False
    assert led._usd_zero_is_unexplained() is True


def test_a_small_zero_cost_run_is_left_alone():
    """Below the floor the zero is still allowed to mean thrift: a run
    that has carried 1,500 tokens tells nothing either way, which is why
    the spec file's cli-backend pin keeps passing."""
    eng = _StandIn(run_budget_usd=10.0, cost_usd=0.0, backend="cli",
                   _measured_cost_turns=1,
                   token_usage={"input": 1_000, "output": 500})
    led = _ledger(eng)
    assert led._usd_budget_enforced() is True
    assert led._usd_zero_is_unexplained() is False
    assert led._unmeasured_budget_block() == ""


def test_a_non_zero_report_is_measured_however_large_the_run():
    """The pinned rule, unchanged: an observation outranks the rate card
    when the backend actually observed something."""
    eng = _StandIn(run_budget_usd=10.0, cost_usd=0.01, backend="cli",
                   _measured_cost_turns=114,
                   token_usage={"input": BROKE, "output": 0})
    assert _ledger(eng)._usd_budget_enforced() is True


def test_a_quota_funded_run_is_never_cross_checked():
    """A provider that bills nothing produces a zero that is a fact.
    Checked BEFORE the token count, so a 363 M-token KIT run stays
    silent -- calling it unmeasured is the false statement in the other
    direction."""
    eng = _StandIn(run_budget_usd=10.0, cost_usd=0.0, provider="kit",
                   _non_billing_turns=114,
                   token_usage={"input": BROKE, "output": 0})
    led = _ledger(eng)
    assert led._usd_budget_enforced() is True
    assert led._usd_zero_is_unexplained() is False
    assert led._unmeasured_budget_block() == ""


def test_the_floor_is_derived_from_the_rate_card_and_the_sessions():
    """2 M tokens is $0.20 at the cheapest published input rate and
    $30.00 at the most expensive, so an exact 0.00 past it is not a
    priced run being cheap."""
    from delfin.agent import pricing
    assert BL.DEFAULT_ZERO_COST_TOKEN_FLOOR == 2_000_000
    rates = [r[0] for r in list(pricing.ANTHROPIC_RATES.values())
             + list(pricing.OPENAI_RATES.values()) if r[0] > 0]
    cheapest = min(rates)
    assert cheapest == pytest.approx(0.10)
    floor_cost = BL.DEFAULT_ZERO_COST_TOKEN_FLOOR / 1_000_000 * cheapest
    assert floor_cost >= 0.05, "the floor must clear the $0.005 rounding"
    # ...and the basis is written where the number is.
    import pathlib
    source = pathlib.Path(BL.__file__).read_text(encoding="utf-8")
    basis = source[:source.index("DEFAULT_ZERO_COST_TOKEN_FLOOR =")]
    for figure in ("$0.10/MTok", "363.5 M", "2,331", "21,211"):
        assert figure in basis, figure


def test_the_notice_names_the_backend_and_not_the_rate_card():
    """Two reasons a ceiling cannot fire, one sentence each: a notice
    that said "no published rate" here would send its reader to the model
    picker when the fix is the backend's own reporting."""
    eng = _StandIn(run_budget_usd=10.0, cost_usd=0.0, backend="cli",
                   _measured_cost_turns=114,
                   token_usage={"input": BROKE, "output": 0})
    block = _ledger(eng)._unmeasured_budget_block()
    assert block.startswith("# Run budget NOT enforceable in USD")
    assert "The cli backend reported exactly $0.00" in block
    assert "no published rate" not in block
    assert "ABSENCE of a measurement, not a cheap run" in block
    # ...with the figure that proves it, and the ceiling that holds.
    assert "363.5 M tokens" in block
    assert "$36.35" in block            # 363.5 M at the cheapest rate
    assert "agent.run_budget_tokens" in block


def test_the_zero_cost_notice_is_said_once():
    eng = _StandIn(run_budget_usd=10.0, cost_usd=0.0, backend="cli",
                   _measured_cost_turns=114,
                   token_usage={"input": BROKE, "output": 0})
    led = _ledger(eng)
    assert led._unmeasured_budget_block() != ""
    assert led._unmeasured_budget_block() == ""


def test_the_unpriced_reason_still_wins_when_both_apply():
    """A run with unpriced turns gets the rate-card sentence: that is the
    first thing to fix, and the block is said once."""
    eng = _StandIn(run_budget_usd=10.0, cost_usd=0.0, backend="cli",
                   _unpriced_turns=3, _measured_cost_turns=114,
                   token_usage={"input": BROKE, "output": 0})
    block = _ledger(eng)._unmeasured_budget_block()
    assert "no published rate" in block
    assert "reported exactly $0.00" not in block


def test_the_dimension_list_marks_the_unexplained_zero_too():
    """The USD line already carried — UNMEASURED for an unpriced run; the
    cross-checked zero reaches the same suffix, through the same
    predicate, so the two cannot drift."""
    eng = _StandIn(run_budget_usd=10.0, cost_usd=0.0, backend="cli",
                   _measured_cost_turns=114,
                   token_usage={"input": BROKE, "output": 0})
    line = next(ln for ln in _ledger(eng)._budget_dimension_lines()
                if "$" in ln)
    assert line == "$0.00 of $10.00 (0 %) — UNMEASURED"


# ---------------------------------------------------------------------------
# 6. Every figure with its denominator
# ---------------------------------------------------------------------------

def test_every_dimension_in_force_is_reported_as_spent_of_ceiling():
    eng = _StandIn(run_budget_usd=10.0, cost_usd=2.5, run_budget_turns=50,
                   _measured_cost_turns=10,
                   token_usage={"input": 41_000_000, "output": 0})
    lines = _ledger(eng)._budget_dimension_lines()
    assert "$2.50 of $10.00 (25 %)" in lines
    assert "10 of 50 turns (20 %)" in lines
    assert "41.0 M of 200.0 M input+output tokens (20 %)" in lines


def test_a_usd_ceiling_that_cannot_be_measured_says_so_in_the_list():
    eng = _StandIn(run_budget_usd=10.0, _unpriced_turns=1)
    line = next(ln for ln in _ledger(eng)._budget_dimension_lines()
                if "$" in ln)
    assert line == "$0.00 of $10.00 (0 %) — UNMEASURED"


def test_a_token_figure_never_appears_as_a_bare_run_of_digits():
    assert BL.format_tokens(41_000_000) == "41.0 M"
    assert BL.format_tokens(1_500) == "1.5 k"
    assert BL.format_tokens(900) == "900"
    assert BL.format_tokens(None) == "0"


def test_a_lifted_ceiling_is_a_dash_not_a_division_by_zero():
    assert BL.BudgetLedger._fraction_pct(5, 0) == "—"
    assert BL.BudgetLedger._fraction_pct(5, None) == "—"


# ---------------------------------------------------------------------------
# 7. The counters travel with the total they qualify
# ---------------------------------------------------------------------------

def test_a_resume_does_not_turn_an_unmeasured_run_into_a_measured_one():
    """The same laundering the price counters were carried for: a resume
    that dropped these would restore a token total covering only the
    turns whose size was reported and hand it to the gates as if it
    covered the run."""
    from delfin.agent.engine import AgentEngine
    src = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(src, spec.attr, spec.reset())
    src._measured_token_turns = 4
    src._untokened_turns = 7
    src._untokened_provider = "blablador"

    resumed = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(resumed, spec.attr, spec.reset())
    resumed.mode, resumed.route, resumed.messages = "solo", [], []
    resumed.session_id, resumed.client = "", None
    resumed.restore_state(src.export_state())

    assert resumed._untokened_turns == 7
    assert resumed._measured_token_turns == 4
    assert resumed._untokened_provider == "blablador"
    assert resumed._token_budget_enforced() is False
