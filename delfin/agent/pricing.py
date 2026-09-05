"""One source of truth for what a token costs, and for admitting when we
do not know.

Four independent price tables used to live in this repo -- two in the API
clients, one in the dashboard's ``/usage`` block, one in the dashboard
status row -- and each carried its own fallback rate for anything it did
not recognise. Three different guesses, none of them a published price:

* The Anthropic client's table was keyed by dated ids
  (``claude-opus-4-20250514``). The string that actually arrives is the
  dropdown value -- ``"opus"``, ``"sonnet"``, ``"haiku"`` -- so no key was
  ever reached and every call took the ``(3.0, 15.0)`` fallback. Measured:
  1M in + 1M out reported $18.00 for all three, while the table's own
  ``claude-opus-4-20250514`` row says $90.00.
* The OpenAI-compatible client fell back to ``(2.0, 8.0)``, so
  ``kit.qwen3.5-397b-A17b`` and ``llama3.3:70b`` -- neither of which bills
  a cent -- were charged $10.00 per 1M in + 1M out of money nobody spent.
  The dashboard has a branch that refuses to fake a price for those
  providers, but it only runs when the reported cost is 0, which that
  fallback guaranteed it never was.
* The dashboard's own table priced ``opus`` at 15/75 while the total on
  the same screen was computed at (3.0, 15.0) -- one block contradicting
  itself by a factor of five.

There is no fallback rate here. :func:`resolve` answers with one of three
states and never with a number it cannot source:

``PRICED``
    A published per-million-token rate.
``NON_BILLING``
    The provider charges no USD for this call (KIT-Toolbox quota, a local
    Ollama/vLLM/LM-Studio endpoint). Zero here is a *measured* zero and
    may be treated as such.
``UNKNOWN``
    No rate on record. This is not zero and it is not a guess. Consumers
    must render it as unmeasured -- never as ``$0.00``, never as a rate --
    and must not let it stand in for a budget check that was never made.

Adding a model to a shipped dropdown without adding it here (or declaring
it unpriced in :data:`DECLARED_UNPRICED`) fails the coverage test, so a
new model cannot quietly inherit somebody's guess again.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "PRICED",
    "NON_BILLING",
    "UNKNOWN",
    "Price",
    "resolve",
    "price_for",
    "cost_usd",
    "is_non_billing",
    "infer_provider",
    "NON_BILLING_PROVIDERS",
    "DECLARED_UNPRICED",
    "ANTHROPIC_RATES",
    "OPENAI_RATES",
]


PRICED = "priced"
NON_BILLING = "non_billing"
UNKNOWN = "unknown"


# ---------------------------------------------------------------------------
# Published rates -- USD per million tokens, (input, output)
# ---------------------------------------------------------------------------

# Anthropic API model ids. Read from the vendor's model/pricing reference
# on 2026-08-12. Only exact API ids appear: an id is a billing contract, a
# tier name is not (see DECLARED_UNPRICED).
ANTHROPIC_RATES: dict[str, tuple[float, float]] = {
    "claude-fable-5": (10.0, 50.0),
    "claude-opus-5": (5.0, 25.0),
    "claude-opus-4-8": (5.0, 25.0),
    "claude-opus-4-7": (5.0, 25.0),
    "claude-opus-4-6": (5.0, 25.0),
    "claude-sonnet-5": (3.0, 15.0),
    "claude-sonnet-4-6": (3.0, 15.0),
    "claude-haiku-4-5": (1.0, 5.0),
    # Dated ids. The first two are carried over unchanged from the table
    # this module replaces; the third is corrected -- that table said
    # (0.80, 4.0) for the 4-5 haiku line, the published rate is (1.0, 5.0),
    # so the old entry under-reported real spend by 20%.
    "claude-opus-4-20250514": (15.0, 75.0),
    "claude-sonnet-4-20250514": (3.0, 15.0),
    "claude-haiku-4-5-20251001": (1.0, 5.0),
}

# OpenAI / Azure-OpenAI model ids, carried over unchanged from the table
# this module replaces. Nothing here is re-derived and nothing is added:
# the change is that an id NOT in this table is now unknown instead of
# silently taking (2.0, 8.0).
OPENAI_RATES: dict[str, tuple[float, float]] = {
    # GPT-5 family
    "gpt-5.4": (2.0, 8.0),
    "gpt-5.4-mini": (0.40, 1.60),
    "gpt-5.3-codex": (2.0, 8.0),
    "gpt-5.2-codex": (2.0, 8.0),
    "gpt-5.2": (2.0, 8.0),
    "gpt-5.1": (2.0, 8.0),
    "gpt-5.1-codex-max": (2.0, 8.0),
    "gpt-5.1-codex-mini": (0.40, 1.60),
    "gpt-5": (2.0, 8.0),
    "gpt-5-mini": (0.40, 1.60),
    "gpt-5-nano": (0.10, 0.40),
    # GPT-4 family
    "gpt-4.1": (2.0, 8.0),
    "gpt-4.1-mini": (0.40, 1.60),
    "gpt-4.1-nano": (0.10, 0.40),
    # o-series reasoning
    "o4-mini": (1.10, 4.40),
    "o3": (2.0, 8.0),
}


# ---------------------------------------------------------------------------
# Providers that bill no USD per call
# ---------------------------------------------------------------------------

# KIT-Toolbox is provided by KIT against a token quota; Ollama / vLLM /
# LM-Studio / llama.cpp run on the user's own hardware. For all of these
# the per-token USD rate is genuinely zero -- which is a fact, not a
# fallback -- and the honest display is tokens spent, not dollars.
NON_BILLING_PROVIDERS: frozenset[str] = frozenset({
    "kit", "ollama", "local", "vllm", "lmstudio", "llamacpp", "llama.cpp",
})

# Server-side aliases the KIT gateway routes to a default model. They
# carry no ``kit.`` prefix, so name them explicitly.
_NON_BILLING_MODEL_IDS: frozenset[str] = frozenset({
    "standard-extern", "standard-external", "standard-local",
})

# ``kit.`` is KIT's own serving prefix. ``azure.`` looks like a vendor
# prefix but in this app it appears ONLY in the KIT provider's model list
# (the KIT gateway fronts those deployments), so an ``azure.`` id is a KIT
# call and costs the user nothing per token.
_NON_BILLING_PREFIXES: tuple[str, ...] = ("kit.", "azure.")


# ---------------------------------------------------------------------------
# Ids that are deliberately left unpriced
# ---------------------------------------------------------------------------

_TIER_ALIAS_REASON = (
    "tier alias, not a billable model id -- the same alias has meant both "
    "(15.0, 75.0) and (5.0, 25.0) per MTok, so any rate here would be a "
    "guess; pick a pinned model id to get a priced run"
)

#: Ids that ship in a dropdown but deliberately have no rate, mapped to the
#: reason. Being in here is what separates "we decided this cannot be
#: priced" from "somebody forgot" -- the coverage test accepts the first
#: and rejects the second.
DECLARED_UNPRICED: dict[str, str] = {
    "opus": _TIER_ALIAS_REASON,
    "sonnet": _TIER_ALIAS_REASON,
    "haiku": _TIER_ALIAS_REASON,
}


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Price:
    """What one (model, provider) pair costs, or why we cannot say."""

    state: str                          # PRICED | NON_BILLING | UNKNOWN
    input_per_mtok: float = 0.0
    output_per_mtok: float = 0.0
    provider: str = ""
    reason: str = ""
    #: True when UNKNOWN was a decision recorded in DECLARED_UNPRICED
    #: rather than an id nobody has classified yet.
    declared: bool = False

    @property
    def measurable(self) -> bool:
        """True when a USD figure for this model means something.

        ``NON_BILLING`` is measurable: its answer is a real zero.
        ``UNKNOWN`` is not -- a zero there is the absence of a number.
        """
        return self.state in (PRICED, NON_BILLING)

    @property
    def rates(self) -> tuple[float, float] | None:
        """``(input, output)`` per MTok, or ``None`` when unknown."""
        if self.state == UNKNOWN:
            return None
        return (self.input_per_mtok, self.output_per_mtok)

    def cost(self, input_tokens: int, output_tokens: int) -> float | None:
        """USD for this many tokens, or ``None`` when unknown."""
        if self.state == UNKNOWN:
            return None
        return (int(input_tokens) * self.input_per_mtok
                + int(output_tokens) * self.output_per_mtok) / 1_000_000


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------


def infer_provider(model: str) -> str:
    """Best-effort provider for callers that only carry a model id.

    Returns ``""`` when the id does not identify its provider; that is not
    an error, it just means the model tables decide on their own.
    """
    mid = (model or "").strip().lower()
    if not mid:
        return ""
    if mid.startswith(_NON_BILLING_PREFIXES) or mid in _NON_BILLING_MODEL_IDS:
        return "kit"
    # Ollama tag convention is ``name:tag`` (``qwen3-coder:32b``). A colon
    # in a bare id is a locally served model; a path-style id is not.
    if ":" in mid and "/" not in mid:
        return "ollama"
    if mid.startswith("claude") or mid in DECLARED_UNPRICED:
        return "claude"
    if mid.startswith(("gpt-", "o1", "o3", "o4")):
        return "openai"
    return ""


def is_non_billing(model: str = "", provider: str = "") -> bool:
    """True when this call costs no USD at all."""
    return resolve(model, provider).state == NON_BILLING


def resolve(model: str = "", provider: str = "") -> Price:
    """Resolve one (model, provider) pair to a :class:`Price`.

    ``provider`` wins over the model id: a GPT deployment served through
    the KIT gateway bills KIT's quota, not OpenAI's rate card, so the
    caller's provider is asked first and only then the model tables.
    """
    mid = (model or "").strip()
    key = mid.lower()
    prov = (provider or "").strip().lower() or infer_provider(mid)

    if prov in NON_BILLING_PROVIDERS:
        return Price(
            state=NON_BILLING, provider=prov,
            reason=("provider bills a token quota, not USD" if prov == "kit"
                    else "model runs on local hardware, no per-token USD cost"),
        )

    if key in DECLARED_UNPRICED:
        return Price(state=UNKNOWN, provider=prov,
                     reason=DECLARED_UNPRICED[key], declared=True)

    rates = ANTHROPIC_RATES.get(key) or OPENAI_RATES.get(key)
    if rates is not None:
        return Price(state=PRICED, input_per_mtok=rates[0],
                     output_per_mtok=rates[1], provider=prov)

    return Price(
        state=UNKNOWN, provider=prov,
        reason=(f"no published rate on record for {mid!r}" if mid
                else "no model named, so no rate to apply"),
    )


def price_for(model: str = "", provider: str = "") -> tuple[float, float] | None:
    """``(input, output)`` USD per million tokens, or ``None``.

    ``None`` means exactly one thing: **we do not know this model's
    price**. It never means free. A provider that bills no USD answers
    ``(0.0, 0.0)`` -- a rate we can stand behind -- and callers that need
    to tell the two apart ask :func:`resolve` for the state.
    """
    return resolve(model, provider).rates


def cost_usd(model: str, provider: str,
             input_tokens: int, output_tokens: int) -> float | None:
    """USD for a turn, or ``None`` when the model has no known rate."""
    return resolve(model, provider).cost(input_tokens, output_tokens)
