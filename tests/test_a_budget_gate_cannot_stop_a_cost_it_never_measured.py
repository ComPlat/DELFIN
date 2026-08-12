"""The run budget passed every check on a run whose cost nobody measured.

``pricing.resolve`` answers with one of three states -- PRICED,
NON_BILLING, UNKNOWN -- and the clients collapse the third into a number
before the engine ever sees it::

    return pricing.cost_usd(self.model, "claude", i, o) or 0.0

So a turn on a model with no published rate adds ``0.0`` to
``engine.cost_usd``. Every USD gate reads that bare float:

* ``_run_budget_status`` divides it by the configured ceiling, so the
  fraction stays at 0% for a run of any size;
* the >110% refusal at the top of ``stream_response`` therefore never
  fires;
* the 80% wind-down block stays silent;
* the per-turn cost circuit-breaker compares the same delta against its
  hard cap and never trips.

The tier aliases the dropdown actually ships -- ``opus``, ``sonnet``,
``haiku`` -- are exactly the ids with no rate on record (they have meant
both 15/75 and 5/25 per MTok, so a rate there would be a guess). The
default Anthropic run is therefore the one where the ceiling is inert,
and it is inert *silently*: a $10.00 budget and a 0% reading look
identical to a cheap run.

A genuinely free provider must NOT be swept up in this. A KIT-Toolbox or
locally served turn costs no money at all; its zero is a measurement, and
calling it "unmeasured" would be a second false statement in the other
direction.

These tests pin the distinction: the gates now carry the price STATE
alongside the number, say once that no monetary ceiling is in force, and
offer the ceilings that are always measured -- turns and tokens.
"""

from __future__ import annotations

import textwrap
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent.api_client import StreamEvent


@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# solo mode")
    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)
    return tmp_path


def _answering_client(model: str):
    """A client that answers with one billed message_delta."""
    def stream(system, messages, max_tokens=4096, session_id="",
               thinking_budget=0, **kw):
        yield StreamEvent(type="text_delta", text="ok")
        # What the client reports for an UNKNOWN-priced model: the
        # tokens are real, the cost could not be computed.
        yield StreamEvent(type="message_delta", input_tokens=1_000,
                          output_tokens=500, cost_usd=0.0)
    client = MagicMock()
    client.model = model
    client.stream_message = MagicMock(side_effect=stream)
    return client


def _engine(agent_tree, *, model: str, provider: str = "claude"):
    from delfin.agent.engine import AgentEngine
    client = _answering_client(model)
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine(repo_dir=agent_tree, backend="api",
                          provider=provider, model=model,
                          mode="quick", pack_dir=agent_tree)
    eng.client = client
    return eng


# ---------------------------------------------------------------------------
# Where the 0.0 comes from
# ---------------------------------------------------------------------------

def test_the_client_hands_the_engine_a_zero_for_an_unpriced_model():
    """The gates never see the state: the client resolves it away.

    ``_estimate_cost`` ends in ``or 0.0``, so UNKNOWN and "this turn was
    free" reach ``engine.cost_usd`` as the same number.
    """
    from types import SimpleNamespace

    from delfin.agent import pricing
    from delfin.agent.api_client import APIClient

    assert pricing.cost_usd("opus", "claude", 1_000_000, 1_000_000) is None
    collapsed = APIClient._estimate_cost(
        SimpleNamespace(model="opus"), 1_000_000, 1_000_000)
    assert collapsed == 0.0


# ---------------------------------------------------------------------------
# The state travels with the number
# ---------------------------------------------------------------------------

def test_a_turn_on_an_unpriced_model_leaves_the_usd_ceiling_unenforced(
        agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.stream_response("hi")
    assert eng.cost_usd == 0.0                  # nothing measurable arrived
    assert eng._unpriced_turns == 1
    assert eng._usd_budget_enforced() is False


def test_a_priced_model_keeps_an_enforceable_ceiling(agent_tree):
    eng = _engine(agent_tree, model="claude-opus-4-20250514")
    eng.run_budget_usd = 10.0
    eng.stream_response("hi")
    assert eng._unpriced_turns == 0
    assert eng._usd_budget_enforced() is True
    assert eng._build_budget_block() == ""


def test_a_quota_funded_run_is_a_measured_zero_not_an_unmeasured_one(
        agent_tree):
    """KIT bills a token quota; the zero is a fact and must stay silent."""
    eng = _engine(agent_tree, model="kit.qwen3.5-397b-A17b", provider="kit")
    eng.run_budget_usd = 10.0
    eng.stream_response("hi")
    assert eng._non_billing_turns == 1
    assert eng._unpriced_turns == 0
    assert eng._usd_budget_enforced() is True
    assert eng._build_budget_block() == ""


def test_a_locally_served_run_says_nothing_either(agent_tree):
    eng = _engine(agent_tree, model="qwen3-coder:32b", provider="ollama")
    eng.run_budget_usd = 10.0
    eng.stream_response("hi")
    assert eng._usd_budget_enforced() is True
    assert eng._build_budget_block() == ""


def test_a_backend_that_reports_its_own_spend_is_measured_whatever_the_table_says(
        agent_tree):
    """An observation outranks the rate card.

    The CLI backend passes the provider's ``total_cost_usd`` straight
    through, so an id with no rate on record still produced a measured
    turn — and warning that nothing was measured would be the false
    statement here.
    """
    eng = _engine(agent_tree, model="opus")
    eng.backend = "cli"
    eng.run_budget_usd = 10.0
    eng.stream_response("hi")
    assert eng._unpriced_turns == 0
    assert eng._usd_budget_enforced() is True
    assert eng._build_budget_block() == ""


def test_an_observed_spend_on_an_unpriced_id_counts_as_measured(agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng._note_turn_price_state(0.42)            # the backend reported spend
    assert eng._unpriced_turns == 0
    assert eng._measured_cost_turns == 1


# ---------------------------------------------------------------------------
# It says so once, and says what it CAN enforce
# ---------------------------------------------------------------------------

def test_the_run_states_plainly_that_no_monetary_ceiling_is_in_force(
        agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.stream_response("hi")
    block = eng._build_budget_block()
    assert "NOT enforceable" in block
    assert "no published rate" in block


def test_the_notice_offers_the_ceilings_that_were_measured(agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.stream_response("hi")
    block = eng._build_budget_block()
    assert "1 turn" in block                    # turns are always countable
    assert "1,000 input" in block               # so are the token totals
    assert "500 output" in block


def test_the_run_says_it_once_and_does_not_nag_every_turn(agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.stream_response("hi")
    assert "NOT enforceable" in eng._build_budget_block()
    eng.stream_response("again")
    assert eng._build_budget_block() == ""


def test_a_run_without_a_usd_budget_has_nothing_to_warn_about(agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 0.0
    eng.run_budget_s = 0.0
    eng.stream_response("hi")
    assert eng._unpriced_turns == 1
    assert eng._build_budget_block() == ""


def test_a_wall_clock_ceiling_is_named_as_the_one_still_in_force(agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.run_budget_s = 600.0
    eng.stream_response("hi")
    assert "600s" in eng._build_budget_block()


def test_the_wind_down_block_survives_alongside_the_notice(agent_tree):
    """The time budget still winds the run down while USD cannot."""
    import time as _t
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.run_budget_s = 100.0
    eng.stream_response("hi")
    eng._run_started_at = _t.time() - 90.0      # 90% of the time budget
    block = eng._build_budget_block()
    assert "NOT enforceable" in block
    assert "wind-down" in block


def test_a_resume_does_not_turn_an_unmeasured_run_into_a_measured_one():
    """The counts travel with the total they qualify.

    ``cost_usd`` is carried across a resume. Dropping the states behind it
    would restore a spend figure that covers only the priced turns and
    hand it to the gates as if it covered the run -- the same laundering
    the run clock was fixed for, one level up.
    """
    from delfin.agent.engine import AgentEngine
    src = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(src, spec.attr, spec.reset())
    src.cost_usd = 2.0
    src._measured_cost_turns = 4
    src._unpriced_turns = 7
    src._non_billing_turns = 1

    resumed = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(resumed, spec.attr, spec.reset())
    resumed.mode, resumed.route, resumed.messages = "solo", [], []
    resumed.session_id, resumed.client = "", None
    resumed.restore_state(src.export_state())

    assert resumed._unpriced_turns == 7
    resumed.run_budget_usd = 10.0
    assert resumed._usd_budget_enforced() is False


def test_the_notice_never_claims_a_ceiling_the_run_does_not_have(agent_tree):
    """With no time budget it must point at a ceiling, not invent one."""
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.run_budget_s = 0.0
    eng.stream_response("hi")
    block = eng._build_budget_block()
    assert "run_budget_s" in block
    assert "EXHAUSTED" not in block
