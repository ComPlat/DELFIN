"""A turn's cost did not include the work it handed to somebody else.

The parent model's tokens went into ``engine.cost_usd``. A delegated run's
tokens went into ``~/.delfin/subagent_telemetry.jsonl``, a file that
carries no parent, no session and no turn field -- so nothing could join
the two, and nothing tried. Every figure a user reads for a turn (the live
terminal line, ``/cost``, ``/usage``, ``get_status``) was the cost of the
part of the turn the parent typed itself. A turn that handed the work to
five sub-agents showed a few cents.

What is asserted here is the REPORTED NUMBERS, never a flag: a turn is
given delegates with known token counts on a model with a known published
rate, and the total, the split and the un-double-counted-ness of both are
read back off the same surfaces a user reads.

Nothing here may touch the real home directory, a real backend or a real
sub-agent: those are properties of the machine, not of the mechanism, and
a test that asserts them passes here and fails on CI.
"""

from __future__ import annotations

import pytest

from delfin.agent import agent_metrics as am
from delfin.agent import repl_commands as rc
from delfin.agent.engine import AgentEngine


# ``gpt-4.1`` has a published rate in pricing.OPENAI_RATES, so the USD
# figures below are arithmetic, not a guess. Provider "openai" is passed
# explicitly everywhere so the numbers do not depend on inference.
PRICED_MODEL = "gpt-4.1"
PRICED_PROVIDER = "openai"


def _rate() -> tuple[float, float]:
    from delfin.agent import pricing

    price = pricing.resolve(PRICED_MODEL, PRICED_PROVIDER)
    assert price.state == pricing.PRICED, (
        "this test prices its delegates from the one price table; a model "
        "that lost its rate makes every figure below meaningless")
    return price.input_per_mtok, price.output_per_mtok


def _usd(input_tokens: int, output_tokens: int) -> float:
    rin, rout = _rate()
    return (input_tokens * rin + output_tokens * rout) / 1_000_000


def _payload(input_tokens: int, output_tokens: int,
             model: str = PRICED_MODEL) -> dict:
    """One delegate's return value, in the shape the runner hands back."""
    return {
        "subagent_type": "explore",
        "description": "find the callers",
        "model": model,
        "result": "found them",
        "input_tokens": input_tokens,
        "output_tokens": output_tokens,
    }


class _Perms:
    """The permissions object the client attaches the runner to."""

    def __init__(self, runner=None):
        self.subagent_runner = runner


class _Client:
    """The least a client can be and still hold a permissions object."""

    def __init__(self, perms):
        self._permissions = perms
        self.model = PRICED_MODEL


def _engine(runner=None, *, provider: str = PRICED_PROVIDER) -> AgentEngine:
    """An engine with the delegate ledger and the turn identity, and
    nothing else -- ``__init__`` opens a backend."""
    eng = AgentEngine.__new__(AgentEngine)
    import threading

    eng._lock = threading.Lock()
    eng._delegate_spend = am.DelegateLedger()
    eng._turn_in_flight = False
    eng._turn_serial = 0
    eng._turn_id = 0
    eng.provider = provider
    eng.cost_usd = 0.0
    eng.token_usage = {"input": 0, "output": 0, "cached": 0}
    eng.mode = "solo"
    eng.backend = "api"
    eng.route = ["solo_agent"]
    eng.current_role_index = 0
    eng.session_id = "sid"
    eng.client = _Client(_Perms(runner))
    return eng


def _open_turn(eng: AgentEngine) -> None:
    """What stream_response does behind the turn gate."""
    eng._turn_in_flight = True
    eng._turn_serial += 1
    eng._turn_id = eng._turn_serial
    eng.delegate_spend().turn = am.DelegateSpend()
    eng._meter_delegate_spend()


def _close_turn(eng: AgentEngine) -> None:
    eng._turn_in_flight = False


# ---------------------------------------------------------------------------
# The join: which turn does a delegate belong to?
# ---------------------------------------------------------------------------

def test_a_foreground_delegate_is_charged_to_the_turn_that_waited_for_it():
    eng = _engine(lambda **kw: _payload(40_000, 3_000))
    _open_turn(eng)
    eng.kit_permissions.subagent_runner(
        subagent_type="explore", description="d", prompt="p")
    _close_turn(eng)

    spend = eng.delegate_spend().turn
    assert spend.count == 1
    assert spend.input_tokens == 40_000
    assert spend.output_tokens == 3_000
    assert spend.cost_usd == pytest.approx(_usd(40_000, 3_000))


def test_an_earlier_turns_delegate_is_not_charged_to_this_one():
    """The failure the whole join exists to prevent: a per-turn figure
    that reports whatever the log happened to contain."""
    eng = _engine(lambda **kw: _payload(40_000, 3_000))
    _open_turn(eng)
    eng.kit_permissions.subagent_runner(
        subagent_type="explore", description="d", prompt="p")
    _close_turn(eng)
    first = eng.delegate_spend().turn.cost_usd
    assert first > 0

    _open_turn(eng)          # a second turn that delegates nothing
    _close_turn(eng)
    assert eng.delegate_spend().turn.count == 0
    assert eng.delegate_spend().turn.cost_usd == 0.0
    # ...and the session still remembers the first one.
    assert eng.delegate_spend().session.cost_usd == pytest.approx(first)


def test_a_delegate_that_outlives_its_turn_is_charged_to_no_turn():
    """A backgrounded delegate returns after the turn that started it.

    Charging it to the turn running by then would charge a turn that
    delegated nothing; charging it to the turn that started it would
    change a number already printed. It is session spend, and says so.
    """
    ledger = {}

    def _slow_runner(**kw):
        # The turn ends while the delegate is still running, exactly as a
        # background thread does.
        _close_turn(ledger["engine"])
        return _payload(120_000, 9_000)

    eng = _engine(_slow_runner)
    ledger["engine"] = eng
    _open_turn(eng)
    eng.kit_permissions.subagent_runner(
        subagent_type="explore", description="d", prompt="p")

    assert eng.delegate_spend().turn.count == 0, (
        "a delegate that finished after the turn cannot be that turn's")
    assert eng.delegate_spend().background.count == 1
    assert eng.delegate_spend().background.cost_usd == pytest.approx(
        _usd(120_000, 9_000))
    assert eng.delegate_spend().session.cost_usd == pytest.approx(
        _usd(120_000, 9_000))


def test_a_delegate_that_returns_during_a_later_turn_is_not_that_turns():
    """The sharper form: the background thread comes back mid-turn-two."""
    state = {}

    def _runner(**kw):
        if kw.get("description") == "background":
            eng = state["engine"]
            _close_turn(eng)          # turn one ends...
            _open_turn(eng)           # ...and turn two starts, while it runs
            return _payload(120_000, 9_000)
        return _payload(1_000, 100)

    eng = _engine(_runner)
    state["engine"] = eng
    _open_turn(eng)
    eng.kit_permissions.subagent_runner(
        subagent_type="explore", description="background", prompt="p")

    assert eng.delegate_spend().turn.count == 0, (
        "turn two delegated nothing and must not inherit turn one's bill")
    assert eng.delegate_spend().background.count == 1


# ---------------------------------------------------------------------------
# The split, and the total
# ---------------------------------------------------------------------------

def test_the_status_reports_the_direct_spend_the_delegated_spend_and_the_total():
    eng = _engine(lambda **kw: _payload(40_000, 3_000))
    eng.cost_usd = 0.1701
    _open_turn(eng)
    eng.kit_permissions.subagent_runner(
        subagent_type="explore", description="d", prompt="p")
    _close_turn(eng)

    st = eng.get_status()
    delegated = _usd(40_000, 3_000)
    assert st["cost_usd"] == pytest.approx(0.1701), (
        "the direct figure must keep meaning exactly what it meant")
    assert st["delegated_cost_usd"] == pytest.approx(delegated)
    assert st["delegate_count"] == 1
    assert st["delegated_input_tokens"] == 40_000
    assert st["delegated_output_tokens"] == 3_000
    assert st["total_cost_usd"] == pytest.approx(0.1701 + delegated)


def test_reading_the_numbers_twice_does_not_bill_the_delegate_twice():
    """The defect a log-scanning join would have: read it again, pay again."""
    eng = _engine(lambda **kw: _payload(40_000, 3_000))
    eng.cost_usd = 0.1701
    _open_turn(eng)
    eng.kit_permissions.subagent_runner(
        subagent_type="explore", description="d", prompt="p")
    _close_turn(eng)

    first = eng.get_status()
    second = eng.get_status()
    third = eng.get_status()
    for key in ("delegated_cost_usd", "delegate_count",
                "delegated_input_tokens", "total_cost_usd"):
        assert second[key] == first[key], key
        assert third[key] == first[key], key
    assert eng.delegate_spend().turn.count == 1


def test_re_installing_the_meter_does_not_stack_two_counters():
    """The client re-binds the runner on every set_permissions, so the
    meter is re-installed every turn. Nested meters would double every
    delegate's bill."""
    eng = _engine(lambda **kw: _payload(40_000, 3_000))
    _open_turn(eng)
    for _ in range(5):
        eng._meter_delegate_spend()
    eng.kit_permissions.subagent_runner(
        subagent_type="explore", description="d", prompt="p")
    _close_turn(eng)

    assert eng.delegate_spend().turn.count == 1
    assert eng.delegate_spend().turn.cost_usd == pytest.approx(
        _usd(40_000, 3_000))


def test_a_fresh_runner_from_the_client_is_metered_again():
    """set_permissions replaces the runner outright; the next turn has to
    wrap the new one or the whole turn goes uncounted."""
    eng = _engine(lambda **kw: _payload(1_000, 100))
    _open_turn(eng)
    _close_turn(eng)
    # ...as _attach_subagent_runner does on the next set_permissions.
    eng.kit_permissions.subagent_runner = lambda **kw: _payload(40_000, 3_000)
    _open_turn(eng)
    eng.kit_permissions.subagent_runner(
        subagent_type="explore", description="d", prompt="p")
    _close_turn(eng)

    assert eng.delegate_spend().turn.count == 1
    assert eng.delegate_spend().turn.cost_usd == pytest.approx(
        _usd(40_000, 3_000))


def test_five_delegates_in_one_turn_add_up():
    eng = _engine(lambda **kw: _payload(40_000, 3_000))
    _open_turn(eng)
    for _ in range(5):
        eng.kit_permissions.subagent_runner(
            subagent_type="explore", description="d", prompt="p")
    _close_turn(eng)

    spend = eng.delegate_spend().turn
    assert spend.count == 5
    assert spend.input_tokens == 200_000
    assert spend.cost_usd == pytest.approx(5 * _usd(40_000, 3_000))


# ---------------------------------------------------------------------------
# An unpriced delegate is not a free one
# ---------------------------------------------------------------------------

def test_a_delegate_with_no_published_rate_is_counted_apart_not_as_zero():
    eng = _engine(lambda **kw: _payload(40_000, 3_000,
                                        model="a-model-nobody-priced"))
    eng.provider = "openai"
    _open_turn(eng)
    eng.kit_permissions.subagent_runner(
        subagent_type="explore", description="d", prompt="p")
    _close_turn(eng)

    spend = eng.delegate_spend().turn
    assert spend.count == 1
    assert spend.unpriced == 1
    assert spend.cost_usd == 0.0, "no rate means no figure, not a zero one"
    assert spend.input_tokens == 40_000, "the tokens are still measured"


def test_a_runner_that_raises_does_not_charge_the_turn():
    def _boom(**kw):
        raise RuntimeError("the delegate died")

    eng = _engine(_boom)
    _open_turn(eng)
    with pytest.raises(RuntimeError):
        eng.kit_permissions.subagent_runner(
            subagent_type="explore", description="d", prompt="p")
    _close_turn(eng)
    assert eng.delegate_spend().turn.count == 0


def test_a_malformed_payload_is_not_booked_as_a_free_run():
    eng = _engine(lambda **kw: "not a dict")
    _open_turn(eng)
    eng.kit_permissions.subagent_runner(
        subagent_type="explore", description="d", prompt="p")
    _close_turn(eng)
    assert eng.delegate_spend().turn.count == 0


def test_a_backend_with_no_runner_simply_has_no_meter():
    eng = _engine(None)
    assert eng._meter_delegate_spend() is False
    eng.client = None
    assert eng._meter_delegate_spend() is False


# ---------------------------------------------------------------------------
# What the turn line prints
# ---------------------------------------------------------------------------

def test_the_turn_line_names_the_split_and_the_total():
    delegated = _usd(40_000, 3_000)
    spend = am.DelegateSpend(count=3, input_tokens=40_000,
                             output_tokens=3_000, cost_usd=delegated)
    line = am.format_turn_delegation(0.1701, spend)
    assert "3 sub-agents" in line
    assert "$0.1701 direct" in line
    assert f"${delegated:.4f} delegated" in line
    assert f"${0.1701 + delegated:.4f}" in line
    assert "40,000" in line and "3,000" in line


def test_a_turn_that_delegated_nothing_prints_nothing():
    """A "delegated $0.00" line on every turn trains a reader to skip the
    line the real number appears on."""
    assert am.format_turn_delegation(0.1701, am.DelegateSpend()) == ""
    assert am.format_turn_delegation(0.1701, None) == ""


def test_the_turn_line_quotes_tokens_when_there_is_no_usd_to_quote():
    """A KIT/local run bills a quota, not dollars. Printing $0.0000 there
    reads as free rather than as unbilled."""
    spend = am.DelegateSpend(count=1, input_tokens=40_000,
                             output_tokens=3_000, cost_usd=0.0, unpriced=1)
    line = am.format_turn_delegation(0.0, spend)
    assert "40,000" in line
    assert "$" not in line.split("·")[0]
    assert "no published rate" in line


def test_the_turn_line_says_one_sub_agent_in_the_singular():
    spend = am.DelegateSpend(count=1, input_tokens=10, output_tokens=1)
    assert "1 sub-agent:" in am.format_turn_delegation(0.0, spend)


# ---------------------------------------------------------------------------
# What /cost and /usage print
# ---------------------------------------------------------------------------

class _StatusEngine:
    """An engine reduced to the one question these commands ask."""

    def __init__(self, **status):
        base = {"provider": "openai", "role": "solo_agent",
                "input_tokens": 1_200, "output_tokens": 340,
                "cached_tokens": 0, "cost_usd": 0.1701}
        base.update(status)
        self._status = base
        self.messages: list = []
        self.client = type("C", (), {"model": PRICED_MODEL})()
        self.session_id = "sid"
        self.kit_permissions = None

    def get_status(self):
        return dict(self._status)


def _delegating_status(**over) -> dict:
    delegated = _usd(120_000, 9_000)
    out = {
        "delegated_cost_usd": delegated,
        "delegated_input_tokens": 120_000,
        "delegated_output_tokens": 9_000,
        "delegate_count": 5,
        "delegates_unpriced": 0,
        "background_delegated_cost_usd": 0.0,
        "total_cost_usd": 0.1701 + delegated,
    }
    out.update(over)
    return out


def _run(name: str, engine, tmp_path) -> str:
    return rc.BUILTINS[name].handler(
        rc.ReplContext(engine=engine, workspace=tmp_path), "").output


def test_cost_reports_the_direct_the_delegated_and_the_sum(tmp_path):
    delegated = _usd(120_000, 9_000)
    out = _run("/cost", _StatusEngine(**_delegating_status()), tmp_path)
    assert "$0.1701" in out
    assert f"${delegated:.4f}" in out
    assert "5 sub-agent(s)" in out
    assert f"= ${0.1701 + delegated:.4f}" in out


def test_cost_stays_one_line_when_nothing_was_delegated(tmp_path):
    out = _run("/cost", _StatusEngine(), tmp_path)
    assert out.count("\n") == 0
    assert "delegated" not in out


def test_cost_names_the_background_share_as_belonging_to_no_turn(tmp_path):
    out = _run("/cost", _StatusEngine(**_delegating_status(
        background_delegated_cost_usd=_usd(60_000, 4_000))), tmp_path)
    assert "background" in out
    assert "charged to no turn" in out


def test_cost_says_when_a_delegate_could_not_be_priced(tmp_path):
    out = _run("/cost", _StatusEngine(**_delegating_status(
        delegates_unpriced=2)), tmp_path)
    assert "no published rate" in out
    assert "NOT in the figure" in out


def test_usage_reports_the_delegated_tokens_and_the_total(tmp_path):
    delegated = _usd(120_000, 9_000)
    out = _run("/usage", _StatusEngine(**_delegating_status()), tmp_path)
    assert "delegated" in out
    assert "120,000" in out and "9,000" in out
    assert f"${delegated:.4f} delegated" in out
    assert f"= ${0.1701 + delegated:.4f}" in out


def test_usage_says_nothing_about_delegation_when_there_was_none(tmp_path):
    out = _run("/usage", _StatusEngine(), tmp_path)
    assert "delegated" not in out


def test_a_status_without_the_new_keys_still_answers(tmp_path):
    """An engine that predates these fields must not break the command."""
    engine = _StatusEngine()
    engine._status.pop("delegate_count", None)
    assert "$0.1701" in _run("/cost", engine, tmp_path)
    assert "rate" in _run("/usage", engine, tmp_path)


# ---------------------------------------------------------------------------
# The record the log keeps
# ---------------------------------------------------------------------------

@pytest.fixture()
def log_path(tmp_path, monkeypatch):
    """Never the real ~/.delfin: the log is a property of the machine."""
    p = tmp_path / "agent_metrics.jsonl"
    monkeypatch.setattr(am, "_LOG_PATH", p)
    return p


def test_a_recorded_turn_carries_both_halves_and_a_total(log_path):
    am.record_turn(am.TurnMetrics(
        model="A", cost_usd=0.1701, delegate_count=5,
        delegated_input_tokens=120_000, delegated_output_tokens=9_000,
        delegated_cost_usd=0.43, delegates_unpriced=1,
    ), path=log_path)
    rec = am.read_turns(path=log_path)[0]
    assert rec["cost_usd"] == pytest.approx(0.1701)
    assert rec["delegated_cost_usd"] == pytest.approx(0.43)
    assert rec["delegate_count"] == 5
    assert rec["delegates_unpriced"] == 1
    assert am.TurnMetrics(cost_usd=0.1701,
                          delegated_cost_usd=0.43).total_cost_usd == (
        pytest.approx(0.6001))


def test_the_aggregate_keeps_the_two_halves_apart(log_path):
    """Direct spend falling while the total rises is what "it delegated
    the work" looks like, and one column cannot show it."""
    am.record_turn(am.TurnMetrics(model="A", cost_usd=0.10,
                                  delegated_cost_usd=0.40,
                                  delegate_count=2), path=log_path)
    am.record_turn(am.TurnMetrics(model="A", cost_usd=0.10), path=log_path)
    agg = am.aggregate_by_model(path=log_path)["A"]
    assert agg["total_cost_usd"] == pytest.approx(0.20)
    assert agg["total_delegated_cost_usd"] == pytest.approx(0.40)
    assert agg["total_cost_usd_incl_delegated"] == pytest.approx(0.60)
    assert agg["avg_delegated_cost_usd"] == pytest.approx(0.20)
    assert agg["delegating_turns"] == 1


def test_a_record_written_before_these_columns_existed_still_aggregates(
        log_path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text('{"model": "A", "cost_usd": 0.02}\n',
                        encoding="utf-8")
    agg = am.aggregate_by_model(path=log_path)["A"]
    assert agg["total_delegated_cost_usd"] == 0.0
    assert agg["total_cost_usd_incl_delegated"] == pytest.approx(0.02)


# ---------------------------------------------------------------------------
# The ledger survives a resume without carrying a turn across it
# ---------------------------------------------------------------------------

def test_the_session_total_survives_a_resume_and_the_turn_bucket_does_not():
    ledger = am.DelegateLedger()
    ledger.turn.merge(am.DelegateSpend(count=1, input_tokens=10,
                                       output_tokens=1, cost_usd=0.01))
    ledger.session.merge(am.DelegateSpend(count=4, input_tokens=400,
                                          output_tokens=40, cost_usd=0.44))
    ledger.background.merge(am.DelegateSpend(count=2, input_tokens=200,
                                             output_tokens=20, cost_usd=0.22))
    back = am.DelegateLedger.from_dict(ledger.as_dict())
    assert back.session.cost_usd == pytest.approx(0.44)
    assert back.background.cost_usd == pytest.approx(0.22)
    assert back.turn.count == 0, (
        "a resumed session has no turn in flight, so a carried per-turn "
        "bucket would charge its first turn for another process's work")


def test_a_session_saved_before_the_ledger_existed_loads_empty():
    fresh = am.DelegateLedger.from_dict(None)
    assert fresh.session.count == 0 and fresh.session.cost_usd == 0.0
    assert am.DelegateLedger.from_dict({"session": "nonsense"}).session.count == 0


def test_the_ledger_is_declared_session_state():
    """Export, restore and reset all read one declaration; a field that is
    not in it comes back missing after every resume."""
    declared = {spec.attr for spec in AgentEngine._SESSION_FIELDS}
    assert "_delegate_spend" in declared


# ---------------------------------------------------------------------------
# The one writer of the record
# ---------------------------------------------------------------------------

def test_the_only_writer_of_turnmetrics_fills_the_delegation_fields():
    """`TurnMetrics` has exactly one call site.

    Adding five columns and not passing them there is the shape this
    whole area keeps producing: the mechanism exists, one hop drops it,
    and the record carries the fields with no figure in them — a turn
    that delegated five agents reading as though it spent only what it
    spent itself.
    """
    import inspect

    from delfin.dashboard import tab_agent

    src = inspect.getsource(tab_agent)
    assert src.count("record_turn(TurnMetrics(") == 1, (
        "a second writer needs the same fields, and this test to be widened")
    assert "_delegation_fields(engine)" in src


def test_the_dashboard_handover_reads_the_turn_bucket_not_the_session():
    """Session spend accumulates; a per-turn record needs the per-turn
    bucket, or every row after the first over-reports."""
    import inspect

    from delfin.dashboard import tab_agent

    src = inspect.getsource(tab_agent)
    start = src.index("def _delegation_fields")
    body = src[start:start + 900]
    assert "delegate_spend().turn" in body
    assert ".session" not in body


def test_a_missing_ledger_costs_the_row_nothing():
    """Best-effort, like the rest of that block: a metrics record is worth
    less than the turn it describes."""
    import inspect

    from delfin.dashboard import tab_agent

    src = inspect.getsource(tab_agent)
    start = src.index("def _delegation_fields")
    body = src[start:start + 900]
    assert "except Exception" in body
    assert "return {}" in body
