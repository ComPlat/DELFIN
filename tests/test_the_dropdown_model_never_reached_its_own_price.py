"""The model the user picked was never the model that got priced.

``APIClient._PRICING`` was keyed by dated API ids -- ``claude-opus-4-20250514``
and friends -- but the string that arrives at the client is the dropdown
value, and the dropdown ships ``("Opus", "opus")``, ``("Sonnet", "sonnet")``,
``("Haiku", "haiku")`` (``tab_agent`` ``_PROVIDER_MODELS_FALLBACK["claude"]``).
``dict.get("opus")`` is None for all three, so **every** Anthropic call fell
through to the ``(3.0, 15.0)`` fallback. Measured against the shipped code:

    APIClient._estimate_cost(1_000_000, 1_000_000)
        "opus"                      -> $18.00
        "sonnet"                    -> $18.00
        "haiku"                     -> $18.00
        "claude-opus-4-20250514"    -> $90.00   (the table's own opus row)

Three different models, one number; and the one row that could have priced
an opus turn was unreachable from the UI. The same number then drove the
per-turn circuit breaker and the run-budget refusal, so the guard on the
most expensive tier ran five times too loose.

``/usage`` made it visible on one screen: its own local table quoted
``$15.0/MTok`` for ``opus`` while the ``Total:`` line beneath came from an
engine total computed at (3.0, 15.0) -- the block contradicted itself by a
factor of five, and never noticed.

Now: one table, no fallback. A pinned id gets its published rate. A tier
alias -- which has meant both (15.0, 75.0) and (5.0, 25.0) in this line, so
no single rate is true of it -- resolves to "unknown" and is printed as
unmeasured, never as a rate and never as zero.
"""

from __future__ import annotations

from delfin.agent import pricing
from delfin.agent.api_client import APIClient
from delfin.dashboard import tab_agent


_MTOK = 1_000_000


def _estimate(model: str) -> float:
    client = APIClient.__new__(APIClient)
    client.model = model
    return client._estimate_cost(_MTOK, _MTOK)


# ---------------------------------------------------------------------------
# The alias no longer collects somebody else's price
# ---------------------------------------------------------------------------


def test_a_tier_alias_is_not_charged_at_sonnets_rate():
    """Each of the three shipped aliases used to bill $18.00/MTok-pair."""
    for alias in ("opus", "sonnet", "haiku"):
        assert pricing.price_for(alias, "claude") is None, alias
        assert _estimate(alias) == 0.0, alias


def test_the_alias_says_why_it_has_no_price():
    """Unknown is a decision on record, not an id nobody classified."""
    price = pricing.resolve("opus", "claude")
    assert price.state == pricing.UNKNOWN
    assert price.declared is True
    assert "alias" in price.reason


def test_a_pinned_model_id_is_charged_its_own_published_rate():
    """The rows that were unreachable from the dropdown now resolve."""
    assert pricing.price_for("claude-opus-4-20250514", "claude") == (15.0, 75.0)
    assert _estimate("claude-opus-4-20250514") == 90.0
    assert pricing.price_for("claude-opus-5", "claude") == (5.0, 25.0)
    assert pricing.price_for("claude-sonnet-5", "claude") == (3.0, 15.0)


def test_the_opus_and_sonnet_rows_are_no_longer_the_same_number():
    """The defect in one line: two tiers, one price."""
    opus = _estimate("claude-opus-4-20250514")
    sonnet = _estimate("claude-sonnet-4-20250514")
    assert opus != sonnet
    assert opus == 5 * sonnet


def test_an_unlisted_anthropic_id_is_unknown_rather_than_sonnet():
    """There is no fallback rate left to inherit."""
    assert pricing.price_for("claude-something-unreleased", "claude") is None
    assert _estimate("claude-something-unreleased") == 0.0


# ---------------------------------------------------------------------------
# The dashboard no longer quotes a rate it cannot stand behind
# ---------------------------------------------------------------------------


def test_the_status_row_does_not_invent_a_dollar_figure_for_an_alias():
    """It used to print "~$0.018" for 1k/1k on any non-KIT API run."""
    out = tab_agent._estimate_cost_str(
        "api", 1_000, 1_000, provider="claude", model="opus",
    )
    assert "$" not in out
    assert "unmeasured" in out
    assert "2,000 tokens" in out


def test_the_status_row_still_prices_a_model_it_knows():
    out = tab_agent._estimate_cost_str(
        "api", _MTOK, _MTOK, provider="claude", model="claude-opus-4-20250514",
    )
    assert out == "~$90.000"


def test_the_status_row_keeps_saying_subscription_for_the_cli_backend():
    out = tab_agent._estimate_cost_str(
        "cli", 1_000, 1_000, provider="claude", model="opus",
    )
    assert out == "included in subscription"


def test_the_status_bar_asks_for_the_model_not_just_the_provider():
    """Within one provider the rate is per model, so the row needs it.

    ``_render_status`` used to be handed only the provider and picked one
    of two hardcoded rates from it.
    """
    priced = tab_agent._render_status(
        "solo", "api", "solo_agent", 0, 1, _MTOK, _MTOK, 0.0,
        provider="claude", model="claude-opus-4-20250514",
    )
    unknown = tab_agent._render_status(
        "solo", "api", "solo_agent", 0, 1, _MTOK, _MTOK, 0.0,
        provider="claude", model="opus",
    )
    assert "$90.000" in priced
    assert "unmeasured" in unknown
    assert "$" not in unknown.split('class="tokens-info"')[1]
