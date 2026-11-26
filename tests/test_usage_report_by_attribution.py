"""The `usage` report: what a molecule/workflow/task/agent cost.

Two properties are what these controls exist for, and each is a defect that
was measured elsewhere before it was written:

* **A `usage` question must never become a billed chat turn.** On the tree
  before this verb, `usage` is not in ``_subcommand_names()``, so
  ``_route_argv`` prepends ``chat`` and the invocation is answered by the
  model. That is why the router test below asserts on the parser and the
  routing function rather than only on the report: a regression that
  un-registers the verb does not raise, it silently spends money.
* **An empty table is not the same answer as "nobody attributed these".**
  The store today carries no molecule/workflow/task/agent column at all, so
  a report that grouped by a missing field and printed nothing would look
  like healthy output over a store it cannot read. It must say so.

The whole file drives ``delfin.agent.cli.main([...])`` -- parser, routing,
dispatch, exit code, stdout -- because a test that calls
``report_costs.collect_usage`` directly measures that the function exists,
not that a person can reach it.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import agent_metrics
from delfin.agent import cli
from delfin.agent import report_costs


def _write_store(tmp_path, records):
    p = tmp_path / "agent_metrics.jsonl"
    p.write_text(
        "".join(json.dumps(r) + "\n" for r in records), encoding="utf-8")
    return p


def _row(ts, *, cost=0.0, price_state="unknown", in_tok=0, out_tok=0,
         tool_calls=0, duration_s=0.0, model="m1", **attribution):
    r = {
        "model": model, "ts": ts, "cost_usd": cost,
        "price_state": price_state,
        "input_tokens": in_tok, "output_tokens": out_tok,
        "tool_calls": tool_calls, "duration_s": duration_s,
        "delegated_cost_usd": 0.0,
    }
    r.update(attribution)
    return r


@pytest.fixture(autouse=True)
def _never_start_a_chat(monkeypatch):
    """Any route into chat fails loudly instead of billing a turn.

    ``build_parser`` reads ``cmd_chat`` from the module at call time, so
    replacing it here covers the real dispatch path. Without this, the
    base tree would answer ``usage`` with a model call -- which is exactly
    the defect under test, and not something a test should trigger.

    It also clears the process-level attribution. ``cli.main`` adopts
    ``DELFIN_SWARM_RUN``/``DELFIN_AGENT_LABEL`` into a module global that
    outlives the call (it has to: it must survive the thread hop), so a
    suite started from inside a swarm member would stamp that member's ids
    on every row written here and these tests would report the harness
    instead of the code. Clearing the environment is NOT enough -- the
    adopted global is a separate thing from the variables it came from.
    """
    from delfin.agent import attribution

    attribution.reset_process()
    try:
        def _refuse(args):  # pragma: no cover - only reached on the defect
            raise AssertionError(
                "'usage' was routed into a chat session instead of the "
                "usage verb -- the token became a billed prompt")
        monkeypatch.setattr(cli, "cmd_chat", _refuse)
        yield
    finally:
        attribution.reset_process()


def test_usage_is_a_verb_and_is_not_routed_into_a_chat_prompt():
    """RED on the base tree: `usage` is swallowed by the chat router."""
    parser = cli.build_parser()
    known = cli._subcommand_names(parser)

    assert "usage" in known, (
        "the usage verb is not registered, so `delfin-agent usage --by "
        "molecule` becomes a chat prompt")
    # The router leaves a real subcommand exactly as typed.
    assert cli._route_argv(["usage", "--by", "molecule"], known) == [
        "usage", "--by", "molecule"]


def test_usage_verb_reaches_the_report_through_cli_main(
        tmp_path, capsys, monkeypatch):
    """The product entry point, not the helper: cli.main builds the report."""
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 30, cost=0.25, price_state="measured",
             in_tok=1000, out_tok=100, tool_calls=3, duration_s=12.0,
             molecule_id="mol-A", workflow_id="wf-1", role="solo"),
        _row(now - 20, cost=0.10, price_state="measured",
             in_tok=500, out_tok=50, tool_calls=1, duration_s=8.0,
             molecule_id="mol-A", workflow_id="wf-1", role="solo"),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    rc = cli.main(["usage", "--by", "molecule", "--json"])
    out = capsys.readouterr().out

    assert rc == 0
    data = json.loads(out)
    assert data["by"] == "molecule"
    assert data["n_rows"] == 2
    group = data["groups"]["mol-A"]
    assert group["cost_total_usd"] == pytest.approx(0.35)
    assert group["input_tokens"] == 1500
    assert group["output_tokens"] == 150
    assert group["tool_calls"] == 4
    assert group["wall_time_s"] == pytest.approx(20.0)
    assert group["cost_meaning"] == "measured"


def test_a_store_with_no_attribution_columns_reports_that_not_an_empty_table(
        tmp_path, capsys, monkeypatch):
    """Every row unattributed: said plainly, and the rows are still counted.

    This is the shape of the real store today. The failure this guards
    against is a report that groups by a missing key, finds nothing, and
    prints an empty table that reads like a clean bill of health.
    """
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 40, cost=0.20, price_state="measured", in_tok=10),
        _row(now - 30, cost=0.30, price_state="measured", in_tok=20),
        _row(now - 20, cost=0.10, price_state="measured", in_tok=30),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    rc = cli.main(["usage", "--by", "molecule"])
    out = capsys.readouterr().out

    assert rc == 0
    # Not an empty table: the turns are counted and reported as unattributed.
    assert "No agent turn records" not in out
    assert "(unattributed)" in out
    assert "all 3 turns (100%) are unattributed" in out
    # And it names what the writers must record for this to work.
    assert "molecule_id" in out


def test_unattributed_share_is_fractional_when_only_some_rows_carry_it(
        tmp_path, capsys, monkeypatch):
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 40, cost=0.1, price_state="measured", molecule_id="mol-A"),
        _row(now - 30, cost=0.1, price_state="measured", molecule_id="mol-A"),
        _row(now - 20, cost=0.1, price_state="measured"),
        _row(now - 10, cost=0.1, price_state="measured"),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    rc = cli.main(["usage", "--by", "molecule", "--json"])
    data = json.loads(capsys.readouterr().out)

    assert rc == 0
    assert data["n_rows"] == 4
    assert data["n_unattributed"] == 2
    assert data["unattributed_share"] == pytest.approx(0.5)
    assert data["columns_present"] == ["molecule_id"]


def test_a_zero_cost_over_unpriced_turns_says_unpriced_not_free(
        tmp_path, capsys, monkeypatch):
    """A 0.00 with nothing behind it must not read as a measured zero."""
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 30, cost=0.0, price_state="unknown", molecule_id="mol-U"),
        _row(now - 20, cost=0.0, price_state="unknown", molecule_id="mol-U"),
        _row(now - 10, cost=0.0, price_state="unknown", molecule_id="mol-U"),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    rc = cli.main(["usage", "--by", "molecule"])
    out = capsys.readouterr().out

    assert rc == 0
    assert "unpriced" in out
    # The bare "0.00" never appears without its meaning attached.
    assert "$0.00 unpriced" in out


def test_a_non_billing_group_is_not_reported_as_unpriced(
        tmp_path, capsys, monkeypatch):
    """The other real zero: the provider charges no USD at all."""
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 30, cost=0.0, price_state="non_billing", molecule_id="m"),
        _row(now - 20, cost=0.0, price_state="non_billing", molecule_id="m"),
        _row(now - 10, cost=0.0, price_state="non_billing", molecule_id="m"),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    rc = cli.main(["usage", "--by", "molecule"])
    out = capsys.readouterr().out

    assert rc == 0
    assert "non-billing" in out
    assert "unpriced" not in out


def test_a_group_below_the_threshold_is_shown_and_marked(
        tmp_path, capsys, monkeypatch):
    """Small groups are shown, not hidden -- but marked."""
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 90, cost=0.1, price_state="measured", molecule_id="big"),
        _row(now - 80, cost=0.1, price_state="measured", molecule_id="big"),
        _row(now - 70, cost=0.1, price_state="measured", molecule_id="big"),
        _row(now - 60, cost=0.1, price_state="measured", molecule_id="tiny"),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    rc = cli.main(["usage", "--by", "molecule", "--json"])
    data = json.loads(capsys.readouterr().out)

    assert rc == 0
    assert data["min_rows"] == report_costs.MIN_GROUP_ROWS
    assert data["groups"]["big"]["below_threshold"] is False
    assert data["groups"]["tiny"]["below_threshold"] is True
    # Shown, not dropped.
    assert data["groups"]["tiny"]["n_rows"] == 1

    text_rc = cli.main(["usage", "--by", "molecule"])
    text = capsys.readouterr().out
    assert text_rc == 0
    assert "tiny" in text
    assert "fewer than" in text and "not comparable" in text


def test_absent_cached_and_reasoning_tokens_are_not_reported_as_zero(
        tmp_path, capsys, monkeypatch):
    """Absent is "nobody recorded it", not "measured zero".

    Neither field is on the store's records today. A sum that folded
    absence into 0 would invent a measurement, which is the same defect as
    an unpriced turn reported as free.
    """
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 30, cost=0.1, price_state="measured", molecule_id="mol-A"),
        _row(now - 20, cost=0.1, price_state="measured", molecule_id="mol-A"),
        _row(now - 10, cost=0.1, price_state="measured", molecule_id="mol-A"),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    rc = cli.main(["usage", "--by", "molecule", "--json"])
    data = json.loads(capsys.readouterr().out)

    assert rc == 0
    assert data["groups"]["mol-A"]["cached_tokens"] is None
    assert data["groups"]["mol-A"]["reasoning_tokens"] is None

    cli.main(["usage", "--by", "molecule"])
    text = capsys.readouterr().out
    assert "not recorded by any row, so not shown as 0" in text
    assert "cached_tokens" in text and "reasoning_tokens" in text


def test_cached_tokens_are_summed_when_the_store_does_carry_them(
        tmp_path, capsys, monkeypatch):
    """The control for the test above: present values ARE summed."""
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 30, cost=0.1, price_state="measured", molecule_id="mol-A",
             cached_tokens=400, reasoning_tokens=70),
        _row(now - 20, cost=0.1, price_state="measured", molecule_id="mol-A",
             cached_tokens=600, reasoning_tokens=30),
        _row(now - 10, cost=0.1, price_state="measured", molecule_id="mol-A",
             cached_tokens=0, reasoning_tokens=0),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    rc = cli.main(["usage", "--by", "molecule", "--json"])
    group = json.loads(capsys.readouterr().out)["groups"]["mol-A"]

    assert rc == 0
    assert group["cached_tokens"] == 1000
    assert group["reasoning_tokens"] == 100


def test_output_carries_no_identifier_beyond_the_grouped_id(
        tmp_path, capsys, monkeypatch):
    """A usage report must not leak session or profile identifiers."""
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 30, cost=0.1, price_state="measured",
             molecule_id="mol-A", session_id="SESSION-SECRET-123",
             profile_name="PROFILE-SECRET-9"),
        _row(now - 20, cost=0.1, price_state="measured",
             molecule_id="mol-A", session_id="SESSION-SECRET-123"),
        _row(now - 10, cost=0.1, price_state="measured",
             molecule_id="mol-A", session_id="SESSION-SECRET-123"),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    for argv in (["usage", "--by", "molecule"],
                 ["usage", "--by", "molecule", "--json"]):
        assert cli.main(argv) == 0
        out = capsys.readouterr().out
        assert "mol-A" in out
        assert "SESSION-SECRET-123" not in out
        assert "PROFILE-SECRET-9" not in out


def test_by_agent_groups_on_the_delegate_then_the_role(
        tmp_path, capsys, monkeypatch):
    """"Which agent spent the tokens": the delegate, else the role -- never model.

    Falling back to ``model`` would report every row as attributed, which
    is the opposite of what the share-of-unattributed column is for.

    The rows are WRITTEN, not hand-built. An earlier version of this test
    composed dicts carrying ``subagent=`` and ``role=``, which are the
    reader's names, not the writer's -- so it was green while ``--by agent``
    grouped 2 of 2 real rows as unattributed. The writer emits
    ``subagent_name``/``agent_role`` (``attribution.COLUMNS``), and only a
    row that went through it can show the two halves agree.
    """
    from delfin.agent import attribution

    store = tmp_path / "agent_metrics.jsonl"
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    cases = [
        ({"role": "subagent:explore", "subagent_name": "explore"}, 0.1),
        ({"role": "solo"}, 0.1),
        ({"role": "solo", "subagent_name": "explore"}, 0.1),
        ({}, 0.1),
    ]
    for attrib, cost in cases:
        with attribution.attribution_scope(**attrib):
            agent_metrics.record_turn(agent_metrics.TurnMetrics(
                model="m1", cost_usd=cost, price_state="measured"))

    rc = cli.main(["usage", "--by", "agent", "--json"])
    data = json.loads(capsys.readouterr().out)

    assert rc == 0
    assert set(data["groups"]) == {"explore", "solo", "(unattributed)"}
    assert data["groups"]["explore"]["n_rows"] == 2
    assert data["groups"]["solo"]["n_rows"] == 1
    assert data["n_unattributed"] == 1
    # "m1" is the model every row names; it must not appear as an agent.
    assert "m1" not in data["groups"]
    # The columns the writer put there are the ones the report found: if the
    # two halves drift again, this is what says so. Only the two ids these
    # rows were written under are present, and both are the writer's names --
    # ``agent_role``/``subagent_name``, not the reader's ``role``/``subagent``.
    assert data["columns_present"] == ["agent_role", "subagent_name"], (
        "the report read columns no writer emits: by=agent would be dead "
        "against real rows even though this fixture looks fine")


def test_since_limits_the_window_and_a_typo_is_refused(
        tmp_path, capsys, monkeypatch):
    """The window is the flag's, and it is narrower than the default.

    An earlier version of this test only asserted ``"recent" in groups``
    after ``--since 2020-01-01`` on a fixture holding a 30-day-old and a
    60-second-old row. Both rows fall inside a six-year window AND inside
    the default seven-day window, so the assertion held with ``--since``
    wired to nothing at all -- measured: replacing ``_parse_since`` with a
    constant left it green. The assertions below separate the two windows:
    the 30-day row is OUT under the default and IN under a 2020 cut.
    """
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 30 * 86_400, cost=5.0, price_state="measured",
             molecule_id="old"),
        _row(now - 60, cost=0.2, price_state="measured",
             molecule_id="recent"),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    # No flag: the default window (the last 7 days) excludes the 30-day row.
    assert cli.main(["usage", "--by", "molecule", "--json"]) == 0
    default = json.loads(capsys.readouterr().out)
    assert "recent" in default["groups"]
    assert "old" not in default["groups"], (
        "the default window did not exclude a 30-day-old row, so the "
        "--since assertion below cannot tell the two windows apart")
    assert default["n_rows"] == 1

    # A 2020 cut reaches back far enough to include it -- the flag widened
    # the window, which only a working --since can do.
    rc = cli.main(["usage", "--by", "molecule", "--since", "2020-01-01",
                   "--json"])
    data = json.loads(capsys.readouterr().out)
    assert rc == 0
    assert "recent" in data["groups"]
    assert "old" in data["groups"], (
        "--since 2020-01-01 did not widen the window, so --since is not "
        "reaching the report")
    assert data["n_rows"] == 2

    # A value that is not a date is refused, not read as "now".
    bad_rc = cli.main(["usage", "--by", "molecule", "--since", "yesterday-ish"])
    captured = capsys.readouterr()
    assert bad_rc == 2
    assert "is not a date" in captured.err
    assert captured.out == ""


def test_an_empty_store_says_there_are_no_records(
        tmp_path, capsys, monkeypatch):
    """Empty is a different statement from "no attribution"."""
    store = tmp_path / "agent_metrics.jsonl"
    store.write_text("", encoding="utf-8")
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    rc = cli.main(["usage", "--by", "molecule"])
    out = capsys.readouterr().out

    assert rc == 0
    assert "No agent turn records" in out
    assert "unattributed" not in out


def test_an_unknown_grouping_is_rejected_by_the_parser():
    """--by is a closed set, so a mistyped grouping never reaches a store.

    Both directions are asserted: rejecting everything would also make the
    second half pass, which is why the valid grouping is parsed too.
    """
    parser = cli.build_parser()
    valid = parser.parse_args(["usage", "--by", "molecule"])
    assert valid.by == "molecule"
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["usage", "--by", "molekuel"])
    assert exc.value.code == 2


def test_a_group_with_unpriced_delegates_is_not_called_measured(
        tmp_path, capsys, monkeypatch):
    """RED before the fix: unpriced delegate spend was silently dropped.

    ``agent_metrics`` records ``delegates_unpriced`` next to
    ``delegated_cost_usd`` for exactly one reason (agent_metrics.py:252):
    a delegate whose model has no published rate contributes 0.0 to the
    measured sum, and counting it as free is the same defect as reading an
    unpriced turn's 0.00 as free. A report that adds the priced delegates
    and stays silent about the rest certifies a group as "measured" while
    part of what it spent is unknown.

    The witness is one row: its own turn is priced, its delegates are
    priced 0.00 over five unpriced delegates. The group must NOT read as a
    clean "measured" $0.01.
    """
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 30, cost=0.01, price_state="measured",
             molecule_id="mol-D", delegates_unpriced=5),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    rc = cli.main(["usage", "--by", "molecule", "--json"])
    data = json.loads(capsys.readouterr().out)

    assert rc == 0
    group = data["groups"]["mol-D"]
    # The number itself is honest -- the priced part is still $0.01 ...
    assert group["cost_total_usd"] == pytest.approx(0.01)
    # ... but what it MEANS is not "measured": five delegates went unpriced.
    assert group["delegates_unpriced"] == 5
    assert group["cost_meaning"] == "mixed", (
        "a group whose delegates include 5 unpriced ones reported "
        f"{group['cost_meaning']!r}; its cost figure is partly unknown")

    rc_text = cli.main(["usage", "--by", "molecule"])
    out = capsys.readouterr().out
    assert rc_text == 0
    assert "mixed" in out, "the table hid that part of the spend is unpriced"


def test_a_group_whose_delegates_are_all_priced_stays_measured(
        tmp_path, capsys, monkeypatch):
    """The control case: the fix must not relabel every group "mixed".

    Same shape as the test above with ``delegates_unpriced=0`` -- this one
    passes BEFORE the fix and after it, so a fix that marks everything
    "mixed" is caught here rather than passing the red test by accident.

    It deliberately asserts nothing about the ``delegates_unpriced`` key:
    that field is new, and depending on it here would make the control red
    on the base too, which is the one thing a control must not be.
    """
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 30, cost=0.02, price_state="measured",
             molecule_id="mol-E", delegates_unpriced=0),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    rc = cli.main(["usage", "--by", "molecule", "--json"])
    data = json.loads(capsys.readouterr().out)

    assert rc == 0
    group = data["groups"]["mol-E"]
    assert group["cost_meaning"] == "measured"


def test_a_window_that_excludes_everything_says_so_instead_of_nothing_recorded(
        tmp_path, capsys, monkeypatch):
    """An empty store and an empty window are different statements.

    Both used to print "No agent turn records in the recent window.", so a
    reader holding a full store and a too-narrow ``--since`` was told that
    nothing had ever been recorded. The distinguishing input is a store
    with rows, all older than the cut.
    """
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 400 * 86_400, cost=1.0, price_state="measured",
             molecule_id="ancient"),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    rc = cli.main(["usage", "--by", "molecule"])
    out = capsys.readouterr().out

    assert rc == 0
    assert "No agent turn records in the window" in out
    assert "1 older records exist" in out
    # It must NOT claim the store was empty.
    assert "No agent turn records in the recent window." not in out

    # The control: a genuinely empty store keeps the original sentence.
    empty = tmp_path / "empty.jsonl"
    empty.write_text("", encoding="utf-8")
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", empty)
    assert cli.main(["usage", "--by", "molecule"]) == 0
    assert "No agent turn records in the recent window." in (
        capsys.readouterr().out)


def test_min_rows_from_the_cli_reaches_the_mark(tmp_path, capsys, monkeypatch):
    """--min-rows is wired to the report, not pinned to the constant.

    Nothing asserted this before, so the handler could pass
    ``report_costs.MIN_GROUP_ROWS`` instead of ``args.min_rows`` and every
    test would still have passed. The witness is a group of 1 row, which is
    marked at the default (3) and unmarked at ``--min-rows 1`` -- the two
    runs differ only in the flag.

    ``--min-rows 0`` is asserted on the JSON only: at 0 no group can be
    marked, so the table renders identically whether the value reached it
    or was coerced back to 3. An assertion on that text could not fail.
    """
    now = time.time()
    store = _write_store(tmp_path, [
        _row(now - 60, cost=0.1, price_state="measured", molecule_id="solo"),
    ])
    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)

    assert cli.main(["usage", "--by", "molecule", "--json"]) == 0
    default = json.loads(capsys.readouterr().out)
    assert default["min_rows"] == report_costs.MIN_GROUP_ROWS
    assert default["groups"]["solo"]["below_threshold"] is True

    assert cli.main(["usage", "--by", "molecule", "--min-rows", "1",
                     "--json"]) == 0
    widened = json.loads(capsys.readouterr().out)
    assert widened["min_rows"] == 1
    assert widened["groups"]["solo"]["below_threshold"] is False

    # And the table follows the same flag: the footnote disappears.
    assert cli.main(["usage", "--by", "molecule", "--min-rows", "1"]) == 0
    assert "not comparable" not in capsys.readouterr().out
    assert cli.main(["usage", "--by", "molecule"]) == 0
    assert "not comparable" in capsys.readouterr().out


def test_a_store_written_by_the_writer_comes_back_with_all_its_numbers(tmp_path,
                                                                      capsys,
                                                                      monkeypatch):
    """The report reads the fields the WRITER writes, not a shape of my choosing.

    Written after an integrator (me) guessed the schema wrong twice while
    checking the report by hand: ``ts`` as an ISO string instead of the epoch
    float ``record_turn`` stores, then ``wall_s``/``cost_meaning`` instead of
    ``duration_s``/``price_state``. Both times the report was right and the
    fixture was wrong, and both times the symptom looked like a report bug --
    "no records", "0.0 seconds", "unstated".

    The other tests build rows by hand, so they inherit whatever names the
    author believed. This one writes through ``record_turn``, the only writer
    of the store, and asserts the report reads back the values that went in.
    Any future rename on the writer side breaks it here instead of in a
    person's terminal.

    It deliberately groups by ``agent`` and does not set ``molecule_id``:
    this test is about the numeric columns, not attribution, so it uses only
    the fields it needs and reads the group back by whatever key the map
    holds. (The attribution columns ARE on ``TurnMetrics`` now --
    ``attribution.COLUMNS`` -- and
    ``test_by_agent_groups_on_the_delegate_then_the_role`` covers them.)
    """
    metrics = agent_metrics.TurnMetrics(
        ts=time.time() - 60,
        session_id="s-1",
        model="m-1",
        input_tokens=1000,
        output_tokens=200,
        cost_usd=0.01,
        price_state="measured",
        tool_calls=4,
        duration_s=12.5,
    )
    store = tmp_path / "agent_metrics.jsonl"
    agent_metrics.record_turn(metrics, path=store)

    monkeypatch.setattr(agent_metrics, "_LOG_PATH", store)
    assert cli.main(["usage", "--by", "agent", "--json"]) == 0
    data = json.loads(capsys.readouterr().out)

    assert data["n_rows"] == 1, (
        "the row the writer just wrote was not read back — if this is 0 the "
        "reader and the writer disagree about a field name")
    assert len(data["groups"]) == 1
    group = next(iter(data["groups"].values()))
    assert group["n_rows"] == 1
    assert group["input_tokens"] == 1000
    assert group["output_tokens"] == 200
    assert group["tool_calls"] == 4
    assert group["wall_time_s"] == pytest.approx(12.5), (
        "wall time is read from a different key than record_turn writes")
    assert group["cost_meaning"] == "measured", (
        "the writer's price_state did not reach the report's cost meaning")
