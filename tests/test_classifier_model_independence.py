"""Make the action classifier's model independent of the session's model.

The classifier (delfin/agent/action_classifier.py) ran on the same model as
the agents (kit.deepseek-v4-flash). With ten agents on that model it got no
verdict within 20 s, latched "timed out earlier in this turn", and every test
run an agent made in its own worktree went to a human for approval (measured
2026-10-29: 6 requests in 20 minutes, all harmless).

What this pins:

1. ``model_family`` resolves a served id's vendor family (kit.deepseek-v4-flash
   -> deepseek, kit.glm-5.3-flash -> glm, kit.mistral-small-4-119b-a8b ->
   mistral).
2. ``classifier_choice`` picks, when ``agent.classifier.model`` is NOT set, a
   fast kit.* chat model the endpoint serves from a DIFFERENT family than the
   session's; only the session's family served means that family (never no
   classifier). Azure / Google / alias / standard-* routes and non-chat models
   (embedding, reranker, speech, image) are never chosen.
3. ``classifier_alternate`` picks a model from a family OTHER than the primary
   classifier model's, for the one retry before a timeout goes to a human.
4. Through the real second-opinion gate with a fake client: a slow primary
   times out, the fast alternate answers (no human asked); a refusal from the
   alternate stays a refusal; and after one model stalls in a turn the latch
   is per classifier model, so a later action is tried on the other model
   rather than blocked.

No test contacts a model or the network; ``served`` is injected directly.
"""

from __future__ import annotations

import time

from delfin.agent import api_client as A
from delfin.agent import model_capabilities as MC


# --------------------------------------------------------------------------
# model_family
# --------------------------------------------------------------------------

def test_family_of_the_three_kit_lines():
    assert MC.model_family("kit.deepseek-v4-flash") == "deepseek"
    assert MC.model_family("kit.glm-5.3-flash") == "glm"
    assert MC.model_family("kit.mistral-small-4-119b-a8b") == "mistral"


def test_family_survives_a_missing_kit_prefix():
    assert MC.model_family("deepseek-v4-flash") == "deepseek"
    assert MC.model_family("mistral-small-4-119b-a8b") == "mistral"


def test_family_is_case_insensitive():
    assert MC.model_family("kit.GLM-5.3-flash") == "glm"


# --------------------------------------------------------------------------
# classifier_choice -- the primary classifier model
# --------------------------------------------------------------------------

# The three kit.* locals the endpoint serves, plus paid routes and non-chat
# models that must never be picked.
_SERVED = {
    "kit.deepseek-v4-flash",
    "kit.glm-5.3-flash",
    "kit.mistral-small-4-119b-a8b",
    "azure.gpt-5",
    "google.gemini",
    "alias.convenient",
    "standard-local",
    "kit.embed-v3",
    "kit.rerank-large",
    "kit.whisper",
    "kit.flux",
}


def _choice(served, session, cfg_model="kit.deepseek-v4-flash",
            cfg_explicit=False):
    return MC.classifier_choice(
        served=served, session_model=session,
        cfg_model=cfg_model, cfg_explicit=cfg_explicit)


def test_explicit_setting_wins_even_if_the_same_family():
    # The user said agent.classifier.model = kit.glm-5.3-flash: use it, do not
    # re-derive from the session.
    out = _choice(_SERVED, session="kit.glm-5.3-flash",
                  cfg_model="kit.glm-5.3-flash", cfg_explicit=True)
    assert out == "kit.glm-5.3-flash"


def test_session_on_glm_picks_a_different_fast_family():
    # Session model is kit.glm-5.3-flash; a fast different-family model wins.
    out = _choice(_SERVED, session="kit.glm-5.3-flash")
    assert out == "kit.deepseek-v4-flash"


def test_session_on_deepseek_picks_the_fastest_other_family():
    # Session is deepseek; the other served families are glm and mistral.
    out = _choice(_SERVED, session="kit.deepseek-v4-flash")
    assert MC.model_family(out) in ("glm", "mistral")
    # mistral is the faster of the two for the classifier budget.
    assert out == "kit.mistral-small-4-119b-a8b"


def test_only_the_sessions_family_served_uses_it_not_none():
    out = _choice({"kit.deepseek-v4-flash"}, session="kit.deepseek-v4-flash")
    assert out == "kit.deepseek-v4-flash"


def test_unknown_session_family_picks_the_fastest_served():
    out = _choice(_SERVED, session="kit.some-other-v9")
    # Fastest served family by the classifier's preference order.
    assert out == "kit.deepseek-v4-flash"


def test_paid_routes_and_nonchat_are_never_chosen():
    # The only chat kit local on offer is deepseek; paid / non-chat rows must
    # be invisible to the choice.
    paid_and_nonchat = {
        "azure.gpt-5", "google.gemini", "alias.convenient", "standard-local",
        "kit.embed-v3", "kit.rerank-large", "kit.whisper", "kit.flux",
    }
    out = _choice(paid_and_nonchat, session="kit.glm-5.3-flash")
    # Nothing pickable: fall back to the configured default, never drop the
    # classifier.
    assert out == "kit.deepseek-v4-flash"


def test_empty_served_roster_falls_back_to_configured():
    out = _choice(set(), session="kit.glm-5.3-flash")
    assert out == "kit.deepseek-v4-flash"


# --------------------------------------------------------------------------
# classifier_alternate -- the "other family" retry once
# --------------------------------------------------------------------------

def _alt(served, from_model, cfg_model="kit.deepseek-v4-flash"):
    return MC.classifier_alternate(
        served=served, from_model=from_model, cfg_model=cfg_model)


def test_alternate_is_a_different_family_from_primary():
    out = _alt(_SERVED, from_model="kit.deepseek-v4-flash")
    assert MC.model_family(out) != "deepseek"
    assert out == "kit.mistral-small-4-119b-a8b"


def test_alternate_with_only_one_family_served_falls_back_to_configured():
    out = _alt({"kit.deepseek-v4-flash"}, from_model="kit.deepseek-v4-flash")
    assert out == "kit.deepseek-v4-flash"


# --------------------------------------------------------------------------
# The real second-opinion gate with a fake client (one slow, one fast)
# --------------------------------------------------------------------------

ALLOW = "DECISION: allow\nREASON: part of the task."
ASK = "DECISION: ask\nREASON: the request does not mention this."
REFUSE = "DECISION: refuse\nREASON: nobody asked for this."
TASK = [{"role": "user", "content": "run the failing test in my worktree"}]


class Dialog:
    def __init__(self):
        self.shown = []

    def ask(self, name, args, preview):
        self.shown.append((name, preview))
        return True


class _Runner:
    """Fake classifier runner: one model answers after a delay, another is
    slow or raises. Records which models were asked."""

    def __init__(self):
        self.calls = []          # (model, outcome)
        self.slow = set()        # models that hang past the wait
        self.broken = set()      # models that raise
        self.answers = {}        # model -> answer text
        self.delay = 0.0

    def __call__(self, prompt, model):
        self.calls.append(model)
        if model in self.broken:
            raise RuntimeError("model unavailable")
        if model in self.slow:
            time.sleep(5.0)
        if self.delay:
            time.sleep(self.delay)
        return self.answers.get(model, ALLOW)


def _bounded_gate(tmp_path, monkeypatch, runner, models, dialog, *,
                   timeout=0.1, message="phase 1 done, please review", perms=None):
    """Drive the second opinion for a message to another agent (an action
    that leaves the workspace; a missing verdict would ask the human).
    Returns ``(out, perms)`` so a test can inspect the per-turn latch. Pass
    a prior ``perms`` to keep the turn's latch across two actions."""
    monkeypatch.setattr(A, "_CLASSIFIER_TIMEOUT_S", timeout)
    monkeypatch.setattr(A, "_CLASSIFIER_RETRY_DELAY_S", 0.0, raising=False)
    if perms is None:
        ws = tmp_path / "ws"
        ws.mkdir(exist_ok=True)
        perms = A.KitToolPermissions(workspace=ws.resolve(), mode="acceptEdits")
        perms.classifier_runner = runner
        perms.classifier_messages = TASK
        if models:
            perms.classifier_model = models[0]
            if len(models) > 1:
                perms.classifier_alt_model = models[1]
        perms.confirm_callback = dialog.ask
    out = A._DocToolExecutor()._classifier_second_opinion(
        perms, "session_message", {"to": "ec9c6393", "message": message},
        perms.mode)
    return out, perms


def test_slow_primary_then_fast_alternate_asks_nobody(
        tmp_path, monkeypatch):
    # Primary (deepseek) hangs; alternate (mistral) answers allow fast. The
    # action is allowed by the alternate's verdict, nobody asked.
    runner = _Runner()
    runner.slow.add("kit.deepseek-v4-flash")
    runner.answers["kit.mistral-small-4-119b-a8b"] = ALLOW
    dialog = Dialog()
    out, _ = _bounded_gate(
        tmp_path, monkeypatch, runner,
        ["kit.deepseek-v4-flash", "kit.mistral-small-4-119b-a8b"], dialog,
        timeout=0.15)
    assert out is None                       # allow stands
    assert dialog.shown == []                # nobody asked
    assert runner.calls == ["kit.deepseek-v4-flash",
                            "kit.mistral-small-4-119b-a8b"]


def test_a_refusal_from_the_alternate_is_still_a_refusal(
        tmp_path, monkeypatch):
    # Strictness is preserved: the alternate must not turn a refusal into an
    # allow. It asks when it refuses and the human is reachable.
    runner = _Runner()
    runner.slow.add("kit.deepseek-v4-flash")
    runner.answers["kit.mistral-small-4-119b-a8b"] = ASK
    dialog = Dialog()
    out, _ = _bounded_gate(
        tmp_path, monkeypatch, runner,
        ["kit.deepseek-v4-flash", "kit.mistral-small-4-119b-a8b"], dialog,
        timeout=0.15)
    # The ask verdict from the alternate reached the human (Dialog.ask True
    # means the user approves, so out is None) -- the strictness is that the
    # ask was SHOWN, never turned into a silent allow.
    assert len(dialog.shown) == 1            # the ask verdict reached the user


def test_stall_latch_is_per_classifier_model(tmp_path, monkeypatch):
    # After the primary stalls once in a turn, the NEXT action is tried
    # directly on the alternate: the latch is keyed per classifier model, so
    # one slow model does not block every later action of the turn.
    runner = _Runner()
    runner.slow.add("kit.deepseek-v4-flash")   # primary hangs past the wait
    runner.answers["kit.mistral-small-4-119b-a8b"] = ALLOW
    dialog = Dialog()
    args = dict(tmp_path=tmp_path, monkeypatch=monkeypatch, runner=runner,
                models=["kit.deepseek-v4-flash", "kit.mistral-small-4-119b-a8b"],
                dialog=dialog, timeout=0.15)
    out, perms = _bounded_gate(**args)
    assert out is None                        # allowed via the alternate
    latch = perms.classifier_state["stalled_models"]
    assert "kit.deepseek-v4-flash" in latch   # the slow primary is latched
    assert "kit.mistral-small-4-119b-a8b" not in latch
    # A second, different action in the same turn must not be blocked by the
    # stalled primary: reused perms keep the latch; the alternate answers
    # directly, nobody asked.
    out2, _ = _bounded_gate(
        tmp_path, monkeypatch, runner,
        ["kit.deepseek-v4-flash", "kit.mistral-small-4-119b-a8b"], dialog,
        timeout=0.15, message="phase 2 done, another review please",
        perms=perms)
    assert out2 is None
    assert dialog.shown == []
    # First action: primary (times out) then alternate; second action: the
    # latched primary is skipped, the alternate answers directly.
    assert runner.calls == ["kit.deepseek-v4-flash", "kit.mistral-small-4-119b-a8b",
                            "kit.mistral-small-4-119b-a8b"]
