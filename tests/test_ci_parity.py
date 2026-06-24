"""CI-parity guard: the local test environment must be key-free, like CI.

CI runs fully mocked with no secrets (.github/workflows/ci.yml). The autouse
``_ci_parity_no_secrets`` fixture in conftest strips API keys from BOTH the
environment and the ``~/.delfin/credentials.json`` store, so a test cannot
quietly depend on an ambient key — which would pass on a developer box that
has one but fail in CI ("lokal grün, CI rot"). These tests pin that contract.
"""

from __future__ import annotations

import os

import pytest

from delfin.agent import credentials

# Kept in sync with conftest._SECRET_KEY_VARS (duplicated to avoid importing the
# conftest module, which pytest loads specially rather than as a package).
_SECRET_KEY_VARS = ("KIT_TOOLBOX_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY")


@pytest.mark.parametrize("var", _SECRET_KEY_VARS)
def test_secret_key_absent_from_environment(var):
    assert os.environ.get(var) in (None, "")


@pytest.mark.parametrize("var", _SECRET_KEY_VARS)
def test_secret_key_unreachable_from_credential_store(var):
    # With the env stripped AND the credential store isolated, resolution must
    # yield nothing — exactly what a client sees in CI.
    assert credentials.load_credential(var) == ""


def test_building_a_kit_client_without_a_key_fails_like_in_ci():
    """A KIT client needs a key; in the key-free test env it must fail the same
    way it would in CI — so a test that builds one without providing a key is
    caught locally, not after the push."""
    from delfin.agent.api_client import create_client

    with pytest.raises(Exception):
        create_client(backend="api", provider="kit", api_key="",
                      model="kit.qwen3.5-397b-A17b")
