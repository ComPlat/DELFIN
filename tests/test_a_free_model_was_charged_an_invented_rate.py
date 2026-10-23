"""Models that bill nothing were charged, and the branch that knew better
could never run.

``OpenAIClient._estimate_cost`` stripped a ``kit.``/``azure.`` prefix and
looked the remainder up in the GPT rate card. ``kit.qwen3.5-397b-A17b``
became ``qwen3.5-397b-A17b``, which is in no table, so it took the
``(2.0, 8.0)`` fallback. Measured against the shipped code, per 1M input +
1M output tokens:

    kit.qwen3.5-397b-A17b   -> $10.00
    kit.gpt-oss-120b        -> $10.00
    kit.minimax-m2.7-229b   -> $10.00
    azure.gpt-5.1           -> $10.00   (GPT rates, on a KIT quota call)
    qwen3-coder:32b         -> $10.00   (a model running on the user's box)
    llama3.3:70b            -> $10.00

None of that money exists. KIT-Toolbox runs against a token quota; an
Ollama/vLLM endpoint runs on the user's own hardware.

The dashboard already had the honest answer -- ``_estimate_cost_str``'s KIT
branch prints "KIT quota: N tokens" rather than faking a price -- but the
status row only reaches it when the reported cost is ``<= 0``, and the
fallback above guaranteed the cost was never 0. The one place in the repo
that refused to invent a price was unreachable *because* of the invention
it was there to avoid.

It also fed a hard stop: ``benchmark_runner`` sets
``engine.run_budget_usd = max(1.0, task.max_cost_usd * 4)``, and the engine
refuses further turns once ``self.cost_usd`` crosses it. A long local-model
benchmark run could be cut off by a budget of money nobody spent -- at
$10 per 1M+1M tokens, a $1 floor is reached after ~100k+100k tokens.
"""

from __future__ import annotations

import pytest

from delfin.agent import pricing
from delfin.agent.api_client import OpenAIClient
from delfin.dashboard import tab_agent


_MTOK = 1_000_000

_KIT_MODELS = [
    "kit.qwen3.5-397b-A17b",
    "kit.gpt-oss-120b",
    "kit.gemma4-31b-it",
    "kit.minimax-m2.7-229b",
    "kit.mistral-small-4-119b-a8b",
    "azure.gpt-5.1",
    "azure.gpt-4.1-mini",
    "azure.o4-mini",
]

_LOCAL_MODELS = [
    "qwen3-coder:32b",
    "qwen2.5-coder:7b",
    "llama3.3:70b",
    "deepseek-coder-v2:lite",
    "mistral-nemo:12b",
]


def _estimate(model: str, provider: str) -> float:
    client = OpenAIClient.__new__(OpenAIClient)
    client.model = model
    client._provider = provider
    return client._estimate_cost(_MTOK, _MTOK)


# ---------------------------------------------------------------------------
# Nothing is charged for a call that costs nothing
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("model", _KIT_MODELS)
def test_a_kit_model_is_not_charged_gpt_rates(model):
    assert _estimate(model, "kit") == 0.0
    assert pricing.resolve(model, "kit").state == pricing.NON_BILLING


@pytest.mark.parametrize("model", _LOCAL_MODELS)
def test_a_locally_served_model_is_not_charged_at_all(model):
    assert _estimate(model, "ollama") == 0.0
    assert pricing.resolve(model, "ollama").state == pricing.NON_BILLING


def test_the_provider_outranks_the_model_name():
    """A GPT deployment behind the KIT gateway bills KIT's quota.

    Prefix-stripping made ``azure.gpt-5.1`` look like an OpenAI call and
    charged the user OpenAI's rate card for it.
    """
    assert _estimate("azure.gpt-5.1", "kit") == 0.0
    assert _estimate("gpt-5.1", "openai") == 10.0


def test_a_zero_from_a_free_provider_is_a_measured_zero():
    """Free and unknown must not collapse into the same 0.0."""
    free = pricing.resolve("kit.qwen3.5-397b-A17b", "kit")
    unknown = pricing.resolve("some-model-nobody-listed", "openai")
    assert free.measurable is True and free.rates == (0.0, 0.0)
    assert unknown.measurable is False and unknown.rates is None


def test_an_unlisted_openai_id_is_unknown_rather_than_two_dollars():
    assert pricing.price_for("gpt-9-turbo", "openai") is None
    assert _estimate("gpt-9-turbo", "openai") == 0.0


def test_a_listed_openai_id_keeps_its_rate():
    assert pricing.price_for("gpt-4.1", "openai") == (2.0, 8.0)
    assert pricing.price_for("gpt-5-nano", "openai") == (0.10, 0.40)


# ---------------------------------------------------------------------------
# The branch that refuses to fake a price can finally run
# ---------------------------------------------------------------------------


def test_the_quota_line_is_reachable_now_that_cost_is_really_zero():
    """It was guarded on cost <= 0, which the invented rate prevented."""
    client = OpenAIClient.__new__(OpenAIClient)
    client.model = "kit.qwen3.5-397b-A17b"
    client._provider = "kit"
    cost = client._estimate_cost(120_000, 30_000)
    assert cost == 0.0

    shown = (
        tab_agent._fmt_cost(cost) if cost > 0
        else tab_agent._estimate_cost_str(
            "api", 120_000, 30_000, provider="kit", model=client.model,
        )
    )
    assert shown == "KIT quota: 150,000 tokens"


def test_a_local_run_reports_tokens_not_dollars():
    shown = tab_agent._estimate_cost_str(
        "api", 4_000, 1_000, provider="ollama", model="qwen3-coder:32b",
    )
    assert shown == "Local quota: 5,000 tokens"
    assert "$" not in shown


def test_the_run_budget_is_no_longer_burned_by_imaginary_spend():
    """A $1 run budget used to be exhausted by ~200k free tokens."""
    client = OpenAIClient.__new__(OpenAIClient)
    client.model = "kit.qwen3.5-397b-A17b"
    client._provider = "kit"
    spent = sum(client._estimate_cost(100_000, 100_000) for _ in range(20))
    assert spent == 0.0
