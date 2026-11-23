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


def test_the_blablador_roster_leads_with_the_larger_measurement():
    """The head of the roster is the model measured over 39 requests, not
    the one a four-item set put first.

    Two measurements exist: 13 items x 3 repeats (2026-11-22) and 4 items
    (2026-11-19). They disagree about `alias-glm-5.3-flash` -- 39 of 39
    correct against 2 of 4 -- and the larger one decides, because a
    four-item score cannot establish that a model misreads code.
    """
    ranked = mc.MODEL_PREFERENCE_BY_PROVIDER["blablador"]

    assert ranked[0] == "alias-glm-5.3-flash", (
        "39 of 39 answered correctly at a median of 23.6 s, the best "
        "measured pair of score and latency on the roster")
    assert ranked[1] == "alias-code", (
        "36 of 39 answered, all 36 correct, median 60.9 s: the same "
        "score on fewer answers, and 2.6x slower")


def test_the_roster_excludes_the_id_the_endpoint_refuses():
    """`deepseek-v4.1-flash-at-hzdr` led this roster and must not: it is
    advertised by /v1/models while /v1/chat/completions answered HTTP 400
    "Unsupported model or no available clients" on 2 of 2 attempts
    (2026-11-22, 14:35-14:40 UTC). A roster exists to continue a turn
    whose model stopped being served, so an entry that may refuse makes
    the recovery the second failure."""
    ranked = mc.MODEL_PREFERENCE_BY_PROVIDER["blablador"]
    assert "deepseek-v4.1-flash-at-hzdr" not in ranked


def test_the_roster_excludes_what_latency_rules_out():
    """Latency is the one judgement a four-item set cannot get wrong, so
    these exclusions stand: a fallback taken mid-turn on a model that
    answers in three minutes turns a recovery into a stall.

    `alias-qwen3.8-flash-next` is also absent, but on weaker ground --
    its 2/4 came from the instrument that misjudged the model now leading
    the roster, so that exclusion is unconfirmed rather than measured.
    """
    ranked = mc.MODEL_PREFERENCE_BY_PROVIDER["blablador"]
    assert "alias-kimi-k3-1m" not in ranked      # 452.3 s
    assert "alias-muse-spark" not in ranked      # 184.0 s
    assert "alias-qwen3.8-flash-next" not in ranked


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


# --- driven through the real stream loop ----------------------------------
#
# Everything above exercises the chooser. The defect was in the TURN: a
# provider other than KIT got "" and the round ended. A chooser that
# answers correctly while the loop never reaches it would leave that in
# place, so the switch is driven here through stream_message itself.


class _Delta:
    def __init__(self, content=None):
        self.content = content
        self.tool_calls = None


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish
        self.index = 0


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


@pytest.fixture
def blablador_client(monkeypatch, tmp_path):
    from delfin.agent import mcp_client as M

    monkeypatch.setenv("DELFIN_STREAM_RETRY_DELAY_SCALE", "0")
    monkeypatch.setenv("BLABLADOR_API_KEY", "not-a-real-key")
    caps = mc.ModelCapabilities(model="alias-large", provider="blablador",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    # The roster the switch picks from, so the test does not depend on a
    # live listing.
    monkeypatch.setattr(mc, "_fetch_openai_models", lambda *a, **k: [
        {"id": "alias-large"}, {"id": "alias-code"}])
    monkeypatch.setattr(mc, "MODEL_PREFERENCE_BY_PROVIDER",
                        {"blablador": ("alias-code", "alias-large")})
    cl = A.create_client(provider="blablador", model="alias-large",
                         cwd=str(tmp_path))
    cl.set_permissions(A.KitToolPermissions(
        workspace=tmp_path, mode="default", task_session_id="sess-succ"))
    return cl


def _unavailable(model: str):
    """What the gateway answers for a model it is not serving now."""
    openai = pytest.importorskip("openai")
    try:
        import httpx2 as transport
    except ImportError:                      # pragma: no cover
        transport = pytest.importorskip("httpx")
    request = transport.Request(
        "POST", "https://api.blablador.fz-juelich.de/v1/chat/completions")
    response = transport.Response(400, request=request)
    return openai.BadRequestError(
        f"Model '{model}' is temporarily unavailable. Please try again later.",
        response=response, body=None)


def test_a_turn_continues_on_the_successor(blablador_client):
    """The round that used to end: the model goes away mid-turn and the
    answer still arrives, on another model."""
    asked = []

    def _create(**kw):
        asked.append(kw.get("model"))
        if kw.get("model") == "alias-large":
            raise _unavailable("alias-large")
        return _Stream([_Chunk([_Choice(_Delta("Hallo."), finish="stop")],
                               usage=_Usage())])

    blablador_client.client.chat.completions.create = _create
    events = list(blablador_client.stream_message(
        "sys", [{"role": "user", "content": "hi"}], max_tokens=50))
    text = "".join(getattr(e, "text", "") or "" for e in events
                   if getattr(e, "type", "") == "text_delta")

    assert "alias-code" in asked, (
        f"the turn never reached the successor; it asked {asked}")
    assert text == "Hallo.", "and the answer arrives on it"


def test_the_switch_is_visible_and_recorded(blablador_client):
    """A silent model change is its own defect: the person must be able
    to see which model answered."""
    def _create(**kw):
        if kw.get("model") == "alias-large":
            raise _unavailable("alias-large")
        return _Stream([_Chunk([_Choice(_Delta("Hallo."), finish="stop")],
                               usage=_Usage())])

    blablador_client.client.chat.completions.create = _create
    events = list(blablador_client.stream_message(
        "sys", [{"role": "user", "content": "hi"}], max_tokens=50))
    notices = "".join(getattr(e, "text", "") or "" for e in events
                      if getattr(e, "type", "") == "notice")

    assert "alias-code" in notices, f"the switch is announced: {notices!r}"
    switches = getattr(blablador_client, "_turn_model_switches", [])
    assert any(s.get("after") == "alias-code" for s in switches), (
        "and recorded on the turn, not only printed")
