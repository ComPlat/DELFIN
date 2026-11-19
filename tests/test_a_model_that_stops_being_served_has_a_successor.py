"""A turn does not end because one model went away -- on any provider.

The recovery existed and was written for one provider. `_kit_fallback_model`
opened with

    if str(provider or "").lower() != "kit":
        return ""

so an unavailable model on any other provider had no successor and the
turn ended. The gap was invisible as long as only KIT was used.

It matters for Blablador in particular: those models are served from a
shared research cluster (JUWELS Booster) whose nodes are drained and
restarted, so "this model is not being served right now" is a normal
event there rather than an incident. An unattended run that meets it at
night would simply stop.

The line is not KIT versus the rest, and not free versus paid as a
matter of taste: switching the model of a PAID provider spends the
person's money on a model they did not choose. That is a decision for
them, not a repair the client may make, so openai and claude keep
answering "" and this file pins that too.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import model_capabilities as mc
from delfin.agent import pricing


# --- the gap that was closed ---------------------------------------------


def test_blablador_has_a_successor(monkeypatch):
    """The case that used to answer "" and end the turn."""
    monkeypatch.setattr(mc, "_fetch_openai_models", lambda *a, **k: None)

    alt = A._kit_fallback_model("alias-large", "blablador")

    assert alt, "an unavailable Blablador model had no successor"
    assert alt != "alias-large", "and the successor is not the model just left"


def test_the_successor_is_one_the_endpoint_lists(monkeypatch):
    """What is served outranks the ranking: a ranked model that is not
    listed is the very situation being recovered from."""
    monkeypatch.setattr(mc, "_fetch_openai_models", lambda *a, **k: [
        {"id": "alias-code"}, {"id": "alias-large"}])

    alt = A._kit_fallback_model("alias-large", "blablador")

    assert alt == "alias-code", (
        "the first ranked model the endpoint actually lists wins")


def test_a_model_already_left_is_not_handed_back(monkeypatch):
    """Without this a turn bounces A->B->A until the switch budget runs
    out, having made no progress."""
    monkeypatch.setattr(mc, "_fetch_openai_models", lambda *a, **k: None)
    ranked = mc.MODEL_PREFERENCE_BY_PROVIDER["blablador"]

    alt = A._kit_fallback_model(ranked[1], "blablador",
                                already_left={ranked[0]})

    assert alt not in (ranked[0], ranked[1])


def test_an_unreadable_listing_still_yields_a_successor(monkeypatch):
    """A roster that cannot be read is not a reason to end the turn. The
    ranking is still better than stopping, with the current model out."""
    def _boom(*_a, **_k):
        raise OSError("the gateway is not answering")
    monkeypatch.setattr(mc, "_fetch_openai_models", _boom)

    assert A._kit_fallback_model("alias-huge", "blablador") == ""


def test_an_empty_listing_falls_back_to_the_ranking(monkeypatch):
    monkeypatch.setattr(mc, "_fetch_openai_models", lambda *a, **k: [])
    ranked = mc.MODEL_PREFERENCE_BY_PROVIDER["blablador"]

    assert A._kit_fallback_model(ranked[0], "blablador") == ranked[1]


# --- the line, and why it is where it is ----------------------------------


@pytest.mark.parametrize("provider", ["openai", "claude", "anthropic", "", "x"])
def test_a_paid_or_unknown_provider_gets_no_successor(provider):
    """Switching a paid provider's model would spend money on a model the
    person did not pick."""
    assert A._kit_fallback_model("gpt-5.4", provider) == ""


def test_every_provider_with_a_roster_bills_no_usd():
    """The map is the mechanism, so the rule lives there rather than in a
    comment: a paid provider must not acquire a successor list by being
    added to it."""
    for provider in mc.MODEL_PREFERENCE_BY_PROVIDER:
        price = pricing.resolve("", provider)
        assert price.state == pricing.NON_BILLING, (
            f"{provider} has a fallback roster but is not a non-billing "
            f"provider; switching its model spends the person's money")


def test_kit_is_unchanged():
    """The provider this was written for keeps its behaviour."""
    assert mc.MODEL_PREFERENCE_BY_PROVIDER["kit"] is mc._KIT_PREFERENCE


# --- the roster itself ----------------------------------------------------


def test_the_blablador_roster_holds_only_measured_models():
    """Ranked from a measurement on 2026-11-19, not from the names. The
    two candidates that scored 2/4 on the code-reasoning items are absent
    on purpose: continuing a turn on a model that misreads the code is
    not a recovery, and `alias-kimi-k3-1m` is absent for its 452 s.
    """
    ranked = mc.MODEL_PREFERENCE_BY_PROVIDER["blablador"]

    assert "alias-glm-5.3-flash" not in ranked
    assert "alias-qwen3.8-flash-next" not in ranked
    assert "alias-kimi-k3-1m" not in ranked
    assert ranked[0] == "deepseek-v4.1-flash-at-hzdr", (
        "the one that scored 4/4 at 6.5 s with the highest request limit "
        "leads; a fallback is taken mid-turn, so latency is part of it")


def test_no_roster_offers_an_embedding_or_image_model():
    """A successor must be able to answer a chat turn and drive tools."""
    bad = ("embedding", "rerank", "whisper", "image", "flux", "ocr",
           "video", "voice")
    for provider, ranked in mc.MODEL_PREFERENCE_BY_PROVIDER.items():
        for model in ranked:
            assert not any(b in model.lower() for b in bad), (
                f"{provider} would continue a turn on {model}")


# --- one model must be enough --------------------------------------------
#
# A standing requirement: DELFIN runs on a single model -- a local server,
# one key, one roster entry -- and gains from several without requiring
# them. Several mechanisms here PREFER another model: the classifier
# confirms an allow with a second family, sub-agents route read-only work
# to a cheap tier, and a turn whose model went away looks for a successor.
# Each must degrade to "the one model decides" rather than to a refusal.
#
# The failure this guards against is not a crash. It is an unattended run
# on one model where every allow needs a second opinion that cannot
# arrive, so each action is put to a person who is not there. Measured
# 2026-11-19: with one served model the confirmation chain is empty and
# the primary decides alone, which is the behaviour these tests pin.


@pytest.mark.parametrize("served", [
    {"kit.glm-5.3-flash"},          # exactly one, and it IS the primary
    {"kit.deepseek-v4-flash"},      # exactly one, a different one
])
def test_one_served_model_decides_alone(served):
    """An allow needs two families only where two exist. With one, the
    primary decides -- otherwise nothing could ever be allowed."""
    primary = mc.classifier_choice(
        served=served, session_model=next(iter(served)),
        cfg_model="kit.deepseek-v4-flash", cfg_explicit=False)
    chain = mc.classifier_confirmation(
        served=served, primary=primary, cfg_model="kit.deepseek-v4-flash")

    assert primary, "there is still a classifier"
    assert tuple(chain) == (), (
        "a second family was asked for although only one model is served; "
        "every allow would then wait for a verdict that cannot come")


@pytest.mark.parametrize("served", [None, set()])
def test_an_unreadable_roster_does_not_demand_a_second_family(served):
    """A gateway that will not list its models is not a reason to stop
    allowing actions."""
    primary = mc.classifier_choice(served=served, session_model="",
                                   cfg_model="kit.deepseek-v4-flash",
                                   cfg_explicit=False)
    assert primary
    assert tuple(mc.classifier_confirmation(
        served=served, primary=primary,
        cfg_model="kit.deepseek-v4-flash")) == ()


def test_a_single_model_roster_yields_no_successor_and_does_not_raise(
        monkeypatch):
    """With one model there is nothing to switch to, and that must be a
    quiet "" rather than an exception on the turn's error path."""
    monkeypatch.setattr(mc, "MODEL_PREFERENCE_BY_PROVIDER",
                        {"blablador": ("alias-code",)})
    monkeypatch.setattr(mc, "_fetch_openai_models",
                        lambda *a, **k: [{"id": "alias-code"}])

    assert A._kit_fallback_model("alias-code", "blablador") == ""


def test_the_classifier_choice_never_answers_empty():
    """Whatever the roster says, there is a classifier. An empty answer
    would disable the watcher silently, which is worse than a weak one."""
    for served in (None, set(), {"kit.glm-5.3-flash"},
                   {"something-unranked"}):
        assert mc.classifier_choice(
            served=served, session_model="", cfg_model="kit.deepseek-v4-flash",
            cfg_explicit=False), f"no classifier for served={served!r}"
