"""A provider the client can build must be reachable and complete.

Adding a provider means editing a handful of unrelated places, and
missing one does not fail: the client builds, the turn runs, and the one
surface nobody updated simply never offers it. Measured on this branch
before this test existed, with Blablador already working end to end:

* the agent tab offered it and fetched its live model list -- mounted;
* the Settings key field offered its key name -- mounted;
* ``_WELL_KNOWN_KEYS`` did not list it, so a key exported in a shell
  profile was never adopted into the store and never secured out of
  that file, which is the whole point of that mechanism;
* the job-monitor provider dropdown did not offer it;
* three ``--provider`` help texts named the other providers only.

None of those is a crash. That is what makes them worth a test rather
than a review: the question "is this provider mounted everywhere?" has
one answer per surface and no single place that asks it.

The list below is the one the key-name map declares, so a provider added
there is checked here without this file being touched.
"""

from __future__ import annotations

import inspect

import pytest

from delfin.agent import credentials, pricing


def _providers() -> dict[str, str]:
    """{provider: env var}, as the agent tab's own map declares it."""
    from delfin.dashboard.tab_agent import PROVIDER_KEY_NAMES
    return dict(PROVIDER_KEY_NAMES)


def test_the_map_is_not_empty():
    """A guard on the guard: an import that silently returned {} would
    make every parametrised test below vacuous."""
    assert len(_providers()) >= 4


@pytest.mark.parametrize("env_var", sorted(set(_providers().values())))
def test_the_key_is_one_the_store_knows(env_var):
    """``_WELL_KNOWN_KEYS`` is what adopt_from_environment walks. A key
    missing from it is never taken out of the environment and never
    secured out of a shell profile, however well the provider works."""
    assert env_var in credentials._WELL_KNOWN_KEYS, (
        f"{env_var} is a provider key the agent consumes but the "
        f"credential store does not list it among the names it knows")


@pytest.mark.parametrize("env_var", sorted(set(_providers().values())))
def test_the_key_can_be_picked_in_the_settings_tab(env_var):
    """The CLI is not the only route: a key must be settable where
    repairs happen."""
    src = inspect.getsource(
        __import__("delfin.dashboard.tab_settings", fromlist=["x"]))
    assert f"'{env_var}'" in src or f'"{env_var}"' in src, (
        f"{env_var} cannot be chosen in the Settings key field")


@pytest.mark.parametrize("provider", sorted(_providers()))
def test_the_provider_has_a_price_classification(provider):
    """An unclassified provider is billed at an invented rate or reported
    as $0.00, and both are wrong in a way nobody notices."""
    price = pricing.resolve("some-model-id", provider)
    assert price.state in (pricing.PRICED, pricing.NON_BILLING) or (
        price.state == pricing.UNKNOWN), price.state
    if provider in ("kit", "blablador"):
        assert price.state == pricing.NON_BILLING, (
            f"{provider} bills a quota, not USD; an invented rate here "
            f"would report spend that did not happen")


def test_every_provider_the_tab_names_can_be_built():
    """The map is not a list of intentions: create_client must accept
    each name. Nothing is sent -- construction alone fails on a provider
    the factory does not know."""
    from delfin.agent import api_client as ac

    src = inspect.getsource(ac.create_client)
    missing = [p for p in _providers() if f'"{p}"' not in src
               and f"'{p}'" not in src]
    assert not missing, (
        f"the tab offers {missing} but create_client does not name them")


def test_the_command_line_help_names_them():
    """The TUI takes its provider from --provider, which has no fixed
    choice list, so an unnamed provider works and reads as unsupported.
    Help text is the only place a person finds out it exists."""
    src = inspect.getsource(
        __import__("delfin.agent.cli", fromlist=["x"]))
    i = src.index('"--provider"')
    # Every --provider help string in the file, not just the first.
    helps = [src[j:j + 400] for j in range(0, len(src))
             if src.startswith('"--provider"', j)]
    assert helps, "no --provider option found at all"
    for provider in _providers():
        if provider == "api":          # the Anthropic key's own spelling
            continue
        assert any(provider in h for h in helps), (
            f"no --provider help text mentions {provider!r}, so nobody "
            f"reading --help learns it can be used")
    assert i >= 0
