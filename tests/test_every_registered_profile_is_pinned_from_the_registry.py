"""The registered profile layer, pinned from the registry itself.

Unlike the hand-written pin tables in
``test_the_model_capability_layer_decisions_are_pinned.py`` and
``test_one_capability_answer_per_question.py``, this file derives its
cases from ``model_profiles.list_profiles()``: a profile added to the
registry is covered by these assertions automatically, with no edit
here.

Pinned, for every registered profile:
  - the weak/strong decision comes from the single heuristic
    (``model_capabilities._is_weak``) and is a plain bool,
  - the profile carries a positive integer ``max_tool_rounds``,
  - ``get_profile`` returns exactly the registered profile.

Pinned, for the output-token cap:
  - ``ModelProfile`` carries NO ``max_output_tokens`` field: the cap a
    request actually carries is decided per role by
    ``engine._ROLE_MAX_TOKENS`` / ``engine._DEFAULT_MAX_TOKENS``
    (``ModelProfile`` never carried it; the write-only field of the
    same name that ``ModelCapabilities`` once had was removed — grep
    evidence in the PR). These pins keep it that way.
"""

from __future__ import annotations

from delfin.agent import engine
from delfin.agent import model_capabilities as mc
from delfin.agent.model_profiles import (
    ModelProfile,
    get_profile,
    list_profiles,
)


def test_the_registry_is_not_empty():
    assert list_profiles(), "profile registry is empty"


def test_get_profile_returns_the_registered_profile_for_every_name():
    registered = dict(list_profiles())
    for name, profile in registered.items():
        assert get_profile(name) is profile, name


def test_every_registered_profile_gets_one_weak_answer():
    """The point of the single heuristic: every caller of the question
    gets the same answer for every registered name.

    Asserting only that the answer is a ``bool`` would pass whatever the
    heuristic said, so it would not notice a second copy drifting. The
    prompt loader's delegate is the other caller that decides prompt
    size, so the two are compared here for the whole registry instead of
    for the hand-written list of ids the older files carry.
    """
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader()
    for name in dict(list_profiles()):
        decision = mc._is_weak(name)
        assert isinstance(decision, bool), name
        assert loader._is_weak_model(name) == decision, name


def test_a_profiles_tool_round_cap_is_the_one_the_registry_holds():
    """``get_profile`` is what the engine asks, so what it hands back has
    to be the registered number, not the dataclass default."""
    for name, profile in list_profiles():
        assert get_profile(name).max_tool_rounds == profile.max_tool_rounds, name


def test_every_registered_profile_has_a_positive_tool_round_cap():
    for name, profile in list_profiles():
        assert isinstance(profile.max_tool_rounds, int), name
        assert profile.max_tool_rounds > 0, name


def test_no_profile_field_is_named_max_output_tokens():
    """The applied output cap is engine's per-role table, not a profile
    field. A field of this name must not come back on ModelProfile."""
    fields = set(ModelProfile.__dataclass_fields__)
    assert "max_output_tokens" not in fields


def test_the_applied_output_cap_lives_in_the_engine():
    assert isinstance(engine._ROLE_MAX_TOKENS, dict) and engine._ROLE_MAX_TOKENS
    for role, cap in engine._ROLE_MAX_TOKENS.items():
        assert isinstance(cap, int) and cap > 0, role
    assert isinstance(engine._DEFAULT_MAX_TOKENS, int)
    assert engine._DEFAULT_MAX_TOKENS > 0
