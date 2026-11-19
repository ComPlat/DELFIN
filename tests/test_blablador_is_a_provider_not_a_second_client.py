"""Helmholtz Blablador is reached through the client DELFIN already has.

Blablador (Forschungszentrum Jülich) serves an OpenAI-compatible API, so
adding it is an endpoint and a key, not a second client. That matters
beyond tidiness: the OpenAI-compatible client is where the tool gate, the
workspace sandbox, the transient-error backoff and the secret handling
live. A separate implementation would be a second place for each of them
to be forgotten -- the same argument the shell escape was rebuilt around.

What this file pins is therefore mostly *absence of a copy*: the same
class, the same sandbox, the same retry ladder, with only the base URL
and the credential name differing.

Its model ids are server-side aliases (alias-fast, alias-code,
alias-large, alias-huge) that point at changing concrete models. Nothing
ships a concrete name; the list is read from /models.
"""

from __future__ import annotations

import inspect

import pytest

from delfin.agent import api_client, backend_parity, pricing
from delfin.agent.model_capabilities import _BLABLADOR_BASE_URL

KEY_NAME = "BLABLADOR_API_KEY"
STAND_IN = "stand-in-value-never-a-real-key-04b1"


def test_the_endpoint_is_the_juelich_one_and_is_not_repeated_anywhere():
    assert _BLABLADOR_BASE_URL == "https://api.blablador.fz-juelich.de/v1"

    src = inspect.getsource(api_client)
    assert src.count("blablador.fz-juelich.de") == 0, (
        "the endpoint belongs in model_capabilities and is imported from "
        "there; a second copy is a second thing to change")


def test_the_client_is_the_openai_compatible_one_with_its_own_endpoint(
        tmp_path):
    client = api_client.create_client(
        backend="api", provider="blablador", model="alias-code",
        api_key=STAND_IN, cwd=str(tmp_path))

    assert type(client).__name__ == "OpenAIClient", (
        "an OpenAI-compatible service must not get a client of its own")
    assert getattr(client, "_provider", "") == "blablador"
    assert _BLABLADOR_BASE_URL in str(getattr(client, "base_url", "")
                                      or getattr(client, "_base_url", ""))


def test_it_is_sandboxed_like_the_other_chat_apis(tmp_path):
    """Built without a permissions policy, such a client offered
    read_file on any absolute path. The sandbox is not optional."""
    client = api_client.create_client(
        backend="api", provider="blablador", model="alias-fast",
        api_key=STAND_IN, cwd=str(tmp_path))

    assert getattr(client, "_permissions", None) is not None


def test_the_key_is_read_from_the_store_not_only_the_environment(
        tmp_path, monkeypatch):
    """A key typed into the Settings tab has to work without a restart,
    the same as every other provider's."""
    monkeypatch.delenv(KEY_NAME, raising=False)
    seen = {}

    def _fake_load(name):
        seen["name"] = name
        return STAND_IN

    monkeypatch.setattr("delfin.agent.credentials.load_credential",
                        _fake_load)
    api_client.create_client(backend="api", provider="blablador",
                            model="alias-fast", cwd=str(tmp_path))

    assert seen.get("name") == KEY_NAME


def test_the_service_bills_nothing_and_says_so_rather_than_guessing():
    """A model that resolves to an undeclared unknown is billed at an
    invented rate. Blablador is a free research service, so zero here is
    measured, not assumed."""
    for alias in ("alias-fast", "alias-code", "alias-large", "alias-huge"):
        assert pricing.infer_provider(alias) == "blablador", alias
        assert pricing.is_non_billing(alias), alias

    assert "blablador" in pricing.NON_BILLING_PROVIDERS


def test_an_unrelated_alias_is_not_claimed():
    """``alias-`` is generic enough to belong to somebody else. The four
    names are listed, not matched on a prefix."""
    assert pricing.infer_provider("alias-of-some-other-service") != "blablador"


def test_the_parity_matrix_knows_what_this_backend_can_do():
    """The capability report is what a person is shown when a backend
    cannot do something. An unknown provider silently falls back to the
    CLI row, which would describe the wrong thing entirely."""
    assert backend_parity.canonical_backend(provider="blablador") == \
        backend_parity.BACKEND_OPENAI

    row = backend_parity.parity_matrix()[backend_parity.BACKEND_OPENAI]
    assert row["bash"]["support"] == "conditional"


def test_it_inherits_the_retry_ladder_rather_than_bringing_its_own():
    """A shared research service will answer 429. The OpenAI-compatible
    client already backs off with jitter and bounded retries; this test
    is here so a future change does not give Blablador its own."""
    src = inspect.getsource(api_client)

    assert "_TRANSIENT_API_STATUS" in src
    assert "429" in src
    start = src.index("_TRANSIENT_API_STATUS")
    assert "blablador" not in src[start:start + 2000].lower(), (
        "the backoff must stay provider-agnostic")


@pytest.mark.parametrize("surface", ["tab_agent", "tab_settings"])
def test_a_person_can_pick_it_and_store_its_key(surface):
    """Two halves of one thing: the provider is offered where sessions
    are started, and its key can be stored where keys are stored. Either
    alone is a dead end -- a provider nobody can authenticate, or a key
    for a provider nobody can choose."""
    import importlib

    module = importlib.import_module(f"delfin.dashboard.{surface}")
    src = inspect.getsource(module.create_tab)

    assert KEY_NAME in src, f"{surface} does not mention {KEY_NAME}"
    if surface == "tab_agent":
        assert '"blablador"' in src
        assert "Helmholtz Blablador" in src, "name it as the service is named"
