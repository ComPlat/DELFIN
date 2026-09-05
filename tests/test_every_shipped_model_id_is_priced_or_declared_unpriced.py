"""A model could be added to a dropdown and silently inherit a guess.

Prices lived in four independent tables with three different fallbacks, so
"this model is not in the table" and "this model costs $2/$8" were the same
event. Nothing anywhere compared the shipped model lists against the price
tables, and the shipped lists had drifted a long way from them:

    of the 5 KIT models, 10 Azure models and 10 Ollama tags the dropdown
    ships, 0 had a price entry -- all 25 took the invented (2.0, 8.0)
    OpenAI fallback,
    and of the 3 Anthropic entries, 0 had one either -- all 3 took the
    invented (3.0, 15.0) fallback.

Only the 11 plain OpenAI ids actually resolved. 25 of 39 selectable models
were being billed at a number nobody published.

This test closes that door. Every id the shipped lists produce must land in
exactly one of three states, and "unknown" only counts when it was written
down as a decision in ``pricing.DECLARED_UNPRICED``:

    priced        a published rate is on record
    non_billing   the provider charges no USD for this call
    unknown       declared, with the reason, in DECLARED_UNPRICED

An id that resolves to an *undeclared* unknown fails here. Adding a model
to a dropdown therefore forces a deliberate answer about what it costs --
which may be "we don't know", as long as somebody says so.

Scope note: this covers the named model-list assignments the UI ships plus
the client defaults and the static capability registry -- the places that
produce a concrete id. Prefix tables (``model_capabilities._STATIC_PREFIX``)
are family patterns, not ids, and free-typed ids are unknowable by
construction; both correctly resolve to unknown and are shown as unmeasured.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from delfin.agent import model_capabilities, pricing
from delfin.agent.api_client import APIClient, CodexCLIClient, OpenAIClient

_REPO = Path(__file__).resolve().parents[1]
_TAB_AGENT = _REPO / "delfin" / "dashboard" / "tab_agent.py"
_TAB_SETTINGS = _REPO / "delfin" / "dashboard" / "tab_settings.py"


def _code_only(path: Path) -> str:
    """Source with comments, strings and whitespace removed.

    A price table has to be recognised in the code, not in the prose that
    describes the one we removed.
    """
    import io
    import tokenize
    kept = []
    with open(path, "rb") as fh:
        for tok in tokenize.tokenize(io.BytesIO(fh.read()).readline):
            if tok.type in (tokenize.COMMENT, tokenize.STRING,
                            tokenize.NL, tokenize.NEWLINE,
                            tokenize.INDENT, tokenize.DEDENT):
                continue
            kept.append(tok.string)
    return "".join(kept)


def _literal_assignment(path: Path, name: str):
    """Read one literal assignment out of a module, nested or not.

    The dropdown lists are locals inside ``build_agent_tab``, so they
    cannot be imported -- but they are plain literals, so they can be read
    from the source. Parsing keeps this test honest: it sees whatever the
    UI actually ships, not a copy maintained alongside it.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == name:
                    return ast.literal_eval(node.value)
    raise AssertionError(f"{name} not found in {path.name}")


def _shipped_ids() -> dict[str, str]:
    """Every concrete model id the app offers -> where it came from."""
    found: dict[str, str] = {}

    def add(model, where):
        mid = (model or "").strip()
        if mid:
            found.setdefault(mid, where)

    lists = _literal_assignment(_TAB_AGENT, "_PROVIDER_MODELS_FALLBACK")
    for provider, entries in lists.items():
        for _label, mid in entries:
            add(mid, f"tab_agent dropdown[{provider}]")

    for name in ("_PROVIDER_DEFAULTS", "_PROVIDER_CHEAP"):
        for provider, mid in _literal_assignment(_TAB_AGENT, name).items():
            add(mid, f"tab_agent {name}[{provider}]")

    jobmon = _literal_assignment(_TAB_SETTINGS, "_JOBMON_MODEL_SUGGESTIONS")
    for provider, entries in jobmon.items():
        for mid in entries:
            add(mid, f"tab_settings jobmon[{provider or 'default'}]")

    for mid in model_capabilities._STATIC:
        add(mid, "model_capabilities._STATIC")

    add(model_capabilities.KIT_BEST_MODEL, "model_capabilities.KIT_BEST_MODEL")
    add(APIClient.DEFAULT_MODEL, "APIClient.DEFAULT_MODEL")
    add(OpenAIClient.DEFAULT_MODEL, "OpenAIClient.DEFAULT_MODEL")
    add(CodexCLIClient.DEFAULT_MODEL, "CodexCLIClient.DEFAULT_MODEL")
    return found


_SHIPPED = _shipped_ids()


def test_the_shipped_lists_were_actually_found():
    """Guard the guard: an empty sweep would pass everything."""
    assert len(_SHIPPED) >= 35, sorted(_SHIPPED)
    for expected in ("opus", "kit.qwen3.5-397b-A17b", "gpt-4.1",
                     "qwen3-coder:32b", "azure.gpt-5.1"):
        assert expected in _SHIPPED


@pytest.mark.parametrize("model", sorted(_SHIPPED))
def test_a_shipped_model_is_priced_free_or_declared_unpriced(model):
    price = pricing.resolve(model)
    if price.state == pricing.UNKNOWN:
        assert price.declared, (
            f"{model!r} ({_SHIPPED[model]}) has no price, is not a "
            f"non-billing provider, and is not declared in "
            f"pricing.DECLARED_UNPRICED. Add a sourced rate, or declare "
            f"it unpriced with the reason -- do not leave it to a guess."
        )
        assert price.reason


def test_a_priced_model_has_a_rate_above_zero():
    """A published rate of exactly zero would be a table typo, not a gift."""
    for model in sorted(_SHIPPED):
        price = pricing.resolve(model)
        if price.state == pricing.PRICED:
            assert price.input_per_mtok > 0, model
            assert price.output_per_mtok > 0, model


def test_a_new_model_is_unknown_until_somebody_prices_it():
    """The default for an unrecognised id is 'we do not know'."""
    price = pricing.resolve("brand-new-frontier-model", "openai")
    assert price.state == pricing.UNKNOWN
    assert price.declared is False
    assert price.rates is None
    assert price.cost(1_000_000, 1_000_000) is None


def test_there_is_one_price_table_and_not_four():
    """The root cause: four tables meant four chances to disagree.

    Two lived in ``api_client`` (``_PRICING``, plus the Codex client
    aliasing one of them), one in the ``/usage`` chat command, one in the
    dashboard status row -- and on one screen ``/usage`` quoted a rate
    from its own table beside a Total computed from another. Rates belong
    in exactly one module; everything else asks.
    """
    for path in (
        _REPO / "delfin" / "agent" / "api_client.py",
        _REPO / "delfin" / "agent" / "benchmark.py",
        _TAB_AGENT,
        _TAB_SETTINGS,
    ):
        code = _code_only(path)
        assert "_PRICING" not in code, f"{path.name} grew a price table again"
        for guess in ("(3.0,15.0)", "(2.0,8.0)", "(15.0,75.0)",
                      "*3.0+", "*2.0+", "*15.0+"):
            assert guess not in code, (
                f"{path.name} contains {guess!r} -- rates and fallback "
                f"rates belong in delfin/agent/pricing.py"
            )


def test_the_usage_command_reads_the_shared_table():
    """``/usage`` prints per-MTok rates, so it must not keep its own."""
    text = _TAB_AGENT.read_text(encoding="utf-8")
    usage = text[text.index('if cmd == "/usage":'):]
    usage = usage[:usage.index('if cmd == "/status":')]
    assert "_pricing.resolve(" in usage
    assert "/MTok" in usage          # it still shows a rate when it has one
    assert "unmeasured" in usage     # and says so when it has none


def test_none_never_means_free():
    """The contract the consumers rely on, stated once."""
    free = pricing.resolve("kit.gpt-oss-120b", "kit")
    unknown = pricing.resolve("brand-new-frontier-model", "openai")
    assert free.rates == (0.0, 0.0) and free.cost(10, 10) == 0.0
    assert unknown.rates is None and unknown.cost(10, 10) is None
    assert free.measurable and not unknown.measurable
