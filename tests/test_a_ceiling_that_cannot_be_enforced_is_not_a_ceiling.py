"""Saying that the ceiling is not in force is not the same as having one.

The first half of this was already shipped: a run on a model with no
published rate now states ONCE, in the system prompt, that its USD
ceiling cannot be enforced, and reports the turn and token counts that
were measured instead. Correct, and still not a ceiling.

What was measured before touching it, on the shipped code:

* ``_run_budget_status`` takes the worst fraction over the enabled
  dimensions, and there are two: USD and wall clock. On an unpriced run
  the USD term is ``0.0 / ceiling`` for a run of any size, so with only
  ``run_budget_usd`` configured the fraction never leaves 0% -- the 80%
  wind-down never speaks and the >110% refusal in ``stream_response``
  never fires. The run is unbounded.
* The notice that says so goes into the SYSTEM PROMPT. An unattended run
  -- which is the entire reason the run budget exists -- shows a silent
  desktop, an empty attention inbox and a green doctor while the one
  ceiling its owner configured does nothing.
* The docstring of the older test file already promised "the ceilings
  that are always measured -- turns and tokens". The notice printed the
  counts; nothing anywhere could stop on them. The rule was a sentence
  and the mechanism permitted its opposite.

So: a turn ceiling that is enforced through the same status function as
the other two, and the unenforceable-USD fact raised once where an
absent operator can see it.

Off by default. A run that configures no turn ceiling behaves exactly as
it did, and a priced run is not told anything new -- both are pinned
below, because a ceiling nobody asked for would be the regression.
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
    def stream(system, messages, max_tokens=4096, session_id="",
               thinking_budget=0, **kw):
        yield StreamEvent(type="text_delta", text="ok")
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
# The gap, stated as a measurement
# ---------------------------------------------------------------------------

def test_a_usd_only_ceiling_on_an_unpriced_run_bounds_nothing(agent_tree):
    """The state before the turn ceiling existed, kept as the reference.

    Ten turns, a $10 ceiling, and a fraction that has not moved. This is
    not a bug in the arithmetic -- 0/10 IS 0% -- it is a gate reading a
    number that was never measured.
    """
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.run_budget_s = 0.0
    for _ in range(10):
        eng.stream_response("hi")
    frac, exhausted = eng._run_budget_status()
    assert eng._unpriced_turns == 10
    assert frac == 0.0
    assert exhausted is False


# ---------------------------------------------------------------------------
# A ceiling that is always measurable
# ---------------------------------------------------------------------------

def test_the_turn_ceiling_is_spent_by_turns_that_no_price_could_measure(
        agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.run_budget_turns = 10
    for _ in range(8):
        eng.stream_response("hi")
    frac, exhausted = eng._run_budget_status()
    assert frac == pytest.approx(0.8)
    assert exhausted is False


def test_the_wind_down_speaks_on_a_run_the_usd_gate_cannot_see(agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.run_budget_turns = 10
    for _ in range(8):
        eng.stream_response("hi")
    block = eng._build_budget_block()
    assert "wind-down" in block
    assert "80%" in block


def test_the_run_stops_instead_of_continuing_forever(agent_tree):
    """Past 110% no further turn starts, and the refusal says which
    ceiling stopped it -- naming only the USD one would send the user to
    a setting that had nothing to do with it."""
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.run_budget_turns = 4
    for _ in range(5):                       # 5/4 = 125%
        eng.stream_response("hi")
    out = eng.stream_response("one more")
    assert "budget" in out.lower()
    assert "run_budget_turns" in out
    # ... and it really did not run.
    assert eng.client.stream_message.call_count == 5


def test_a_counted_turn_is_any_turn_that_ran_priced_or_not(agent_tree):
    """The count must not be the unpriced count. A run that switches to a
    priced model half way would otherwise stop spending its ceiling."""
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_turns = 10
    eng.stream_response("hi")
    eng.client.model = "claude-opus-4-20250514"
    eng.model = "claude-opus-4-20250514"
    eng.stream_response("hi")
    assert eng._unpriced_turns == 1
    frac, _ = eng._run_budget_status()
    assert frac == pytest.approx(0.2)


# ---------------------------------------------------------------------------
# Nobody who did not ask for it is affected
# ---------------------------------------------------------------------------

def test_no_turn_ceiling_configured_changes_nothing(agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.run_budget_turns = 0
    for _ in range(50):
        eng.stream_response("hi")
    frac, exhausted = eng._run_budget_status()
    assert frac == 0.0
    assert exhausted is False


def test_a_priced_run_is_not_told_about_a_ceiling_it_does_not_need(
        agent_tree):
    eng = _engine(agent_tree, model="claude-opus-4-20250514")
    eng.run_budget_usd = 10.0
    eng.stream_response("hi")
    assert eng._build_budget_block() == ""


def test_the_worst_dimension_wins_and_the_others_do_not_hold_it_back(
        agent_tree):
    """Turns are added to the max, not averaged into it: a run at 90% of
    one ceiling and 5% of another is at 90%."""
    import time as _t
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_s = 1000.0
    eng.run_budget_turns = 10
    eng._run_started_at = _t.time() - 50.0       # 5% of the time budget
    for _ in range(9):                           # 90% of the turn budget
        eng.stream_response("hi")
    frac, _ = eng._run_budget_status()
    assert frac == pytest.approx(0.9, abs=0.02)


# ---------------------------------------------------------------------------
# The notice points at the ceiling that now exists
# ---------------------------------------------------------------------------

def test_the_notice_names_the_ceiling_that_can_be_enforced(agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.run_budget_s = 0.0
    eng.stream_response("hi")
    block = eng._build_budget_block()
    assert "run_budget_turns" in block


def test_a_run_that_already_set_a_turn_ceiling_is_told_it_is_in_force(
        agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    eng.run_budget_turns = 20
    eng.stream_response("hi")
    block = eng._build_budget_block()
    assert "20 turn" in block


# ---------------------------------------------------------------------------
# An absent operator finds out
# ---------------------------------------------------------------------------

def test_an_unenforceable_ceiling_reaches_the_attention_inbox(agent_tree):
    """The prompt notice is read by the model and nobody else. The run
    budget exists for runs nobody is watching."""
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    seen: list[tuple] = []
    with patch("delfin.agent.attention.emit_attention",
               side_effect=lambda *a, **k: seen.append((a, k))):
        eng.stream_response("hi")
        eng._build_budget_block()
    assert len(seen) == 1
    kinds = [a[0] for a, _ in seen]
    assert kinds == ["budget_warning"]
    detail = seen[0][1].get("detail", "") + seen[0][1].get("title", "")
    assert "not in force" in detail or "not enforce" in detail


def test_it_reaches_the_inbox_once_however_long_the_run_goes_on(agent_tree):
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_usd = 10.0
    seen: list = []
    with patch("delfin.agent.attention.emit_attention",
               side_effect=lambda *a, **k: seen.append(k)):
        for _ in range(5):
            eng.stream_response("hi")
            eng._build_budget_block()
    assert len(seen) == 1


def test_a_measured_run_raises_nothing(agent_tree):
    eng = _engine(agent_tree, model="claude-opus-4-20250514")
    eng.run_budget_usd = 10.0
    seen: list = []
    with patch("delfin.agent.attention.emit_attention",
               side_effect=lambda *a, **k: seen.append(k)):
        eng.stream_response("hi")
        eng._build_budget_block()
    assert seen == []


def test_a_quota_funded_run_raises_nothing_either(agent_tree):
    eng = _engine(agent_tree, model="kit.qwen3.5-397b-A17b", provider="kit")
    eng.run_budget_usd = 10.0
    seen: list = []
    with patch("delfin.agent.attention.emit_attention",
               side_effect=lambda *a, **k: seen.append(k)):
        eng.stream_response("hi")
        eng._build_budget_block()
    assert seen == []


# ---------------------------------------------------------------------------
# The setting is read from where every other budget setting is read
# ---------------------------------------------------------------------------

def test_the_turn_ceiling_can_be_configured_in_settings(agent_tree,
                                                        monkeypatch):
    eng = _engine(agent_tree, model="opus")
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"run_budget_turns": 6}})
    assert eng._run_budget_turns() == 6


def test_an_explicit_attribute_outranks_the_setting(agent_tree, monkeypatch):
    """Same precedence as the other two: the scheduler daemon sets these
    per entry and must not be overridden by a global default."""
    eng = _engine(agent_tree, model="opus")
    eng.run_budget_turns = 3
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"run_budget_turns": 99}})
    assert eng._run_budget_turns() == 3


def test_a_broken_settings_file_does_not_invent_a_ceiling(agent_tree,
                                                          monkeypatch):
    eng = _engine(agent_tree, model="opus")

    def _boom():
        raise OSError("settings unreadable")

    monkeypatch.setattr("delfin.user_settings.load_settings", _boom)
    assert eng._run_budget_turns() == 0


def test_a_negative_ceiling_is_off_not_immediately_exhausted(agent_tree,
                                                             monkeypatch):
    eng = _engine(agent_tree, model="opus")
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: {"agent": {"run_budget_turns": -5}})
    assert eng._run_budget_turns() == 0


# ---------------------------------------------------------------------------
# The ceiling survives what the spend total survives
# ---------------------------------------------------------------------------

def test_a_resume_does_not_hand_the_ceiling_a_fresh_turn_count():
    """``cost_usd`` is carried across a resume and so are the turn counts
    it is qualified by; a restored run that reset them to zero would give
    itself the whole ceiling back."""
    from delfin.agent.engine import AgentEngine
    src = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(src, spec.attr, spec.reset())
    src._unpriced_turns = 9

    resumed = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(resumed, spec.attr, spec.reset())
    resumed.mode, resumed.route, resumed.messages = "solo", [], []
    resumed.session_id, resumed.client = "", None
    resumed.restore_state(src.export_state())
    resumed.run_budget_turns = 10

    frac, _ = resumed._run_budget_status()
    assert frac == pytest.approx(0.9)


# ---------------------------------------------------------------------------
# The other half of the same task: a model id the backend cannot accept
# ---------------------------------------------------------------------------

def test_a_tier_alias_is_refused_before_the_first_request():
    """``opus`` is not an id the Messages API knows.

    The dropdown's fallback list for this provider ships exactly the three
    tier aliases, and they are valid for the CLI backend -- but this
    client puts the string straight into the request body, so the API
    backend answered every such run with a vendor 404 after the session
    was already built. Refusing at construction turns an opaque failure
    on the first turn into one sentence at the point where the contract
    is known.

    Resolving the alias to a dated id instead would be a guess about
    WHICH model the user gets -- the same guess the pricing table refuses
    to make about the rate.
    """
    from delfin.agent.api_client import APIClient
    for alias in ("opus", "sonnet", "haiku"):
        with pytest.raises(ValueError) as caught:
            APIClient(api_key="test-key", model=alias)
        assert alias in str(caught.value)
        assert "claude-" in str(caught.value)      # names a usable id


def test_a_pinned_id_is_accepted():
    from delfin.agent.api_client import APIClient
    client = APIClient(api_key="test-key", model="claude-sonnet-4-20250514")
    assert client.model == "claude-sonnet-4-20250514"


def test_the_default_model_of_this_client_is_one_it_accepts():
    from delfin.agent.api_client import APIClient
    client = APIClient(api_key="test-key")
    assert client.model == APIClient.DEFAULT_MODEL


def test_switching_to_an_alias_mid_run_is_refused_too():
    """A model switch is the same contract one turn later."""
    from delfin.agent.api_client import APIClient
    client = APIClient(api_key="test-key", model="claude-sonnet-4-20250514")
    with pytest.raises(ValueError):
        client.switch_model("opus")
    assert client.model == "claude-sonnet-4-20250514"


def test_the_cli_backend_still_takes_the_aliases_it_was_built_for():
    """The aliases are correct there; refusing them everywhere would
    remove a working path to fix a broken one."""
    from delfin.agent.api_client import CLIClient
    assert CLIClient.DEFAULT_MODEL == "sonnet"
