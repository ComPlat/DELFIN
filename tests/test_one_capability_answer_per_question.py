"""One capability answer per question — pinned on the base commit.

This module pins, for every model id the repository names, the answer the
capability layer gives *today* for the questions a caller can ask:

  - "is this a weak model?"   -> ``model_capabilities._is_weak``
  - "how many auto tool rounds?" -> ``ModelProfile.max_tool_rounds``
  - "what output-token cap does this role carry?" ->
    ``engine._ROLE_MAX_TOKENS`` / ``engine._DEFAULT_MAX_TOKENS``

and the one URL every module that talks to the KIT Toolbox endpoint resolves:

  - ``model_capabilities._KIT_BASE_URL``

The values are literals recorded from the base (a later move of the weak-model
heuristic must keep every cell identical; a difference is a behaviour change
and has to be argued separately). The pin also cross-checks the weak-model
question against ``PromptLoader._is_weak_model`` — the two public answers must
agree for every id the repository names, otherwise there are two answers to
one question and the merge is rejected.

A source-text check makes sure no *unprotected* module spells the KIT endpoint
literal any more: the constant lives only in model_capabilities (its owner).
"""

from __future__ import annotations

import pytest

from delfin.agent import engine
from delfin.agent import model_capabilities as mc
from delfin.agent.model_profiles import list_profiles
from delfin.agent.prompt_loader import PromptLoader

# The KIT Toolbox endpoint the capability layer owns, in one place.
_EXPECTED_KIT_BASE_URL = "https://ki-toolbox.scc.kit.edu/api/v1"
_EXPECTED_KIT_MODELS_URL = _EXPECTED_KIT_BASE_URL + "/models"

# Every registered model profile: (is_weak, max_tool_rounds, compact, core_only).
# Recorded from the shipped registry as it stands. A change to any cell is a
# behaviour change, not a refactor.
_PINNED_PROFILES = {
    "azure.gpt-5":           (False, 100, False, False),
    "azure.gpt-5-mini":      (False, 100, False, False),
    "azure.gpt-5-nano":      (False, 100, False, False),
    "azure.gpt-5.1":         (False, 100, False, False),
    "azure.gpt-5.4":         (False, 100, False, False),
    "haiku":                 (False, 100, False, False),
    "kit.deepseek-v4-flash": (False, 100, False, False),
    "kit.gemma4-31b-it":     (False,  10, False, False),
    "kit.glm-5.3":           (False, 100, False, False),
    "kit.gpt-oss-120b":      (False,  15, False, True),
    "kit.qwen3.5-397b-A17b": (False, 100, False, False),
    "opus":                  (False, 100, False, False),
    "sonnet":                (False, 100, False, False),
}

# Non-registered model ids the repository names (static table, prefix table,
# existing tests). Only the weak-model answer is pinned for these.
_PINNED_WEAK_ONLY = {
    "claude-opus-5":             False,
    "codellama:7b":              True,
    "deepseek-r1:32b":           False,
    "deepseek-r1:671b":          False,
    "deepseek-r1:7b":            True,
    "deepseek-v3":               False,
    "gemma2:2b":                 True,
    "gemma3:2b":                 True,
    "gemma4:7b":                 True,
    "gemma:27b":                 False,
    "gpt-4.1":                   False,
    "gpt-5":                     False,
    "gpt-5.4":                   False,
    "kit.glm-5.3-flash":         False,
    "kit.mistral-small-4-119b-a8b": False,
    "llama3.1:8b":               True,
    "llama3.2:3b":               True,
    "llama3.2:8b":               True,
    "llama3.3:70b":              False,
    "llava:7b":                  True,
    "mistral-7b":                True,
    "mistral-small:22b":         False,
    "o3":                        False,
    "o4-mini":                   False,
    "phi-3.5":                   True,
    "qwen2.5-coder:14b":         True,
    "qwen2.5-coder:32b":         False,
    "qwen2.5-coder:7b":          True,
    "qwen3-coder:32b":           False,
    "qwen3-vl:4b":               True,
    "qwen3:0.6b":                True,
    "qwen3:14b":                 True,
    "qwq:32b":                   False,
    "standard-extern":           False,
    "standard-local":            False,
    "totally-unknown-model":     False,
}

# Per role: the output-token cap a request carries.
_PINNED_ROLE_MAX_TOKENS = {
    "builder_agent": 32768,
    "chief_agent": 8192,
    "critic_agent": 16384,
    "dashboard_agent": 16384,
    "research_agent": 8192,
    "reviewer_agent": 16384,
    "runtime_agent": 8192,
    "session_manager": 16384,
    "solo_agent": 32768,
    "test_agent": 8192,
}
_PINNED_DEFAULT_MAX_TOKENS = 8192

# Modules that talk to the KIT *endpoint* and must all resolve one URL string.
# api_client.py and tab_agent.py are protected core files (missing here on the
# branch; the handover script makes them agree and this is asserted on the
# integration tree where the script has run).
_KIT_ENDPOINT_MODULES = (
    "delfin/doctor.py",
    "delfin/agent/model_capabilities.py",
)

_ALL_IDS = sorted(set(_PINNED_PROFILES) | set(_PINNED_WEAK_ONLY))


def _profile_dict():
    return dict(list_profiles())


# ---------------------------------------------------------------------------
# The weak-model question: one answer, both public entries agree
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("model", _ALL_IDS)
def test_is_weak_answers_agree_across_the_two_entries(model):
    """model_capabilities._is_weak and PromptLoader._is_weak_model must not
    be two independent answers: they are one heuristic, exposed two ways."""
    expected_prof = _PINNED_PROFILES.get(model)
    expected_weak = (
        expected_prof[0] if expected_prof is not None
        else _PINNED_WEAK_ONLY[model]
    )
    assert mc._is_weak(model) is expected_weak
    assert PromptLoader()._is_weak_model(model) is expected_weak


def test_every_name_we_pin_is_still_known_to_the_layer():
    """Every pinned id still reaches the same profile / answer path: the
    registry has not dropped a profile, and no pin names a model the layer
    cannot answer for."""
    names = {name for name, _ in list_profiles()}
    missing = sorted(set(_PINNED_PROFILES) - names)
    assert not missing, f"profiles dropped from the registry: {missing}"


# ---------------------------------------------------------------------------
# The tool-rounds question (per registered profile)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("model", sorted(_PINNED_PROFILES))
def test_max_tool_rounds_is_pinned_per_profile(model):
    weak, rounds, compact, core_only = _PINNED_PROFILES[model]
    profile = dict(list_profiles())[model]
    assert mc._is_weak(model) is weak
    assert profile.compact_prompt is compact
    assert profile.core_tools_only is core_only
    assert profile.max_tool_rounds == rounds


# ---------------------------------------------------------------------------
# The output-cap question (per role)
# ---------------------------------------------------------------------------


def test_the_per_role_output_token_caps_are_pinned():
    assert engine._ROLE_MAX_TOKENS == _PINNED_ROLE_MAX_TOKENS


def test_the_default_output_token_cap_is_pinned():
    assert engine._DEFAULT_MAX_TOKENS == _PINNED_DEFAULT_MAX_TOKENS


# ---------------------------------------------------------------------------
# The KIT endpoint URL question
# ---------------------------------------------------------------------------


def test_the_kit_endpoint_url_is_one_pinned_value():
    assert mc._KIT_BASE_URL == _EXPECTED_KIT_BASE_URL


def test_every_kit_endpoint_module_resolves_the_same_url():
    """Every module that talks to the KIT endpoint builds the URL from the one
    constant, so a change to the endpoint travels to all of them at once."""
    from delfin import doctor

    assert mc._KIT_BASE_URL + "/models" == _EXPECTED_KIT_MODELS_URL
    assert doctor._KIT_MODELS_URL == _EXPECTED_KIT_MODELS_URL


def test_no_unprotected_module_spells_the_kit_endpoint_literal():
    """The endpoint literal ('.../api/v1' or '.../api/v1/models') may live only
    in model_capabilities, which owns it. Protected core files are rewritten by
    the handover script and are asserted on the integration tree."""
    from pathlib import Path

    root = Path(__file__).resolve().parent.parent / "delfin"
    # model_capabilities owns the constant (it spells the literal by design).
    owner = root / "agent" / "model_capabilities.py"
    # Protected files the handover script (not this test tree) rewrites.
    protected = {
        root / "agent" / "api_client.py",
        root / "dashboard" / "tab_agent.py",
    }
    literal = "https://ki-toolbox.scc.kit.edu/api/v1"
    offenders = []
    for path in sorted(root.rglob("*.py")):
        if path == owner or path in protected:
            continue
        text = path.read_text(encoding="utf-8")
        if literal in text:
            offenders.append(str(path.relative_to(root)))
    assert not offenders, f"endpoint literal still hardcoded in: {offenders}"
