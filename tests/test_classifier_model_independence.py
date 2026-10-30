"""The action classifier's model is independent of the session's model.

The classifier (delfin/agent/action_classifier.py) ran on the same model as
the agents (kit.deepseek-v4-flash). With ten agents on that model it got no
verdict within 20 s, latched "timed out earlier in this turn", and every test
run an agent made in its own worktree went to a human for approval (measured
2026-10-29: 6 requests in 20 minutes, all harmless).

What this pins:

1. ``model_family`` resolves a served id's vendor family (kit.deepseek-v4-flash
   -> deepseek, kit.glm-5.3-flash -> glm, kit.mistral-small-4-119b-a8b ->
   mistral).
2. ``classifier_choice`` picks, when ``agent.classifier.model`` is NOT set, the
   first served model in the single bench table ``_CLASSIFIER_MODELS`` (the
   primary). Azure / Google / alias / standard-* routes and non-chat models
   (embedding, reranker, speech, image) are never chosen.
3. ``classifier_confirmation`` returns the ALLOW-confirmation chain: served
   bench models of a family OTHER than the primary's, in bench order. Empty
   exactly when only the primary's family is served (it decides alone).
4. Through the real second-opinion gate with a fake client, "allow needs two
   families": the primary's refuse/ask stands at once; an ALLOW is confirmed
   by a model of another family. An allow whose confirming model refuses
   becomes a refuse; an allow whose confirming model times out asks (never
   allowed on one family's vote); a primary refuse never calls the confirm.
   A per-model stall latch keeps one slow model from blocking every later
   action, and a timed-out thread cannot pollute a later family's verdict.

No test contacts a model or the network; ``served`` is injected directly.
"""

from __future__ import annotations

import time

from delfin.agent import action_classifier as _ac
from delfin.agent import api_client as A
from delfin.agent import model_capabilities as MC


# --------------------------------------------------------------------------
# one table: the bench lineup is fixed and ordered
# --------------------------------------------------------------------------

def test_the_bench_lineup_is_one_ordered_table():
    # The ranking lives in ONE table so it can be re-set from a future
    # orchestrator measurement without touching code. Order is the bench:
    # primary first, then the ALLOW-confirmation chain.
    assert MC._CLASSIFIER_MODELS == (
        "kit.mistral-small-4-119b-a8b",
        "kit.deepseek-v4-flash",
        "kit.glm-5.3-flash",
    )
    # And the three families differ, so a confirm is always another family.
    fams = [MC.model_family(m) for m in MC._CLASSIFIER_MODELS]
    assert len(set(fams)) == len(fams)


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
# classifier_choice -- the PRIMARY classifier model
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


def test_primary_is_the_first_served_bench_model():
    # When all three locals are served, the primary is mistral (first in the
    # bench table), not deepseek.
    got = MC.classifier_choice(
        served=_SERVED, session_model="kit.deepseek-v4-flash",
        cfg_model=MC._CLASSIFIER_MODELS[1], cfg_explicit=False)
    assert got == "kit.mistral-small-4-119b-a8b"


def test_primary_skips_to_the_next_served_bench_model():
    # mistral not served -> deepseek (next in the table).
    got = MC.classifier_choice(
        served={"kit.deepseek-v4-flash", "kit.glm-5.3-flash"},
        session_model="kit.deepseek-v4-flash",
        cfg_model="kit.glm-5.3-flash", cfg_explicit=False)
    assert got == "kit.deepseek-v4-flash"


def test_explicit_classifier_model_is_respected_as_is():
    got = MC.classifier_choice(
        served=_SERVED, session_model="kit.deepseek-v4-flash",
        cfg_model="kit.other", cfg_explicit=True)
    assert got == "kit.other"


def test_paid_routes_and_nonchat_are_never_chosen():
    # Only the kit.* chat locals are candidates; a roster that lists only
    # paid / non-chat entries must fall back, never pick one of them.
    bad = {"azure.gpt-5", "google.gemini", "alias.convenient",
           "standard-local", "kit.embed-v3", "kit.rerank-large",
           "kit.whisper", "kit.flux"}
    got = MC.classifier_choice(
        served=bad, session_model="kit.deepseek-v4-flash",
        cfg_model="kit.fallback", cfg_explicit=False)
    assert got == "kit.fallback"


def test_empty_served_roster_falls_back_to_configured():
    got = MC.classifier_choice(
        served=set(), session_model="kit.deepseek-v4-flash",
        cfg_model="kit.fallback", cfg_explicit=False)
    assert got == "kit.fallback"


def test_unknown_served_family_decides_alone():
    # A roster of a family no bench covers: the endpoint decides alone
    # (never no classifier).
    got = MC.classifier_choice(
        served={"kit.frankenstein-1"}, session_model="kit.deepseek-v4-flash",
        cfg_model="kit.fallback", cfg_explicit=False)
    assert got == "kit.frankenstein-1"


# --------------------------------------------------------------------------
# classifier_confirmation -- the ALLOW second opinion
# --------------------------------------------------------------------------

def test_confirmation_is_the_other_families_in_bench_order():
    # Primary mistral -> deepseek then glm (deepseek is faster in the bench,
    # so it confirms first).
    chain = MC.classifier_confirmation(
        served={"kit.mistral-small-4-119b-a8b", "kit.deepseek-v4-flash",
                "kit.glm-5.3-flash"},
        primary="kit.mistral-small-4-119b-a8b", cfg_model="kit.glm-5.3-flash")
    assert chain == ("kit.deepseek-v4-flash", "kit.glm-5.3-flash")


def test_confirmation_excludes_the_primary_family():
    chain = MC.classifier_confirmation(
        served={"kit.mistral-small-4-119b-a8b", "kit.deepseek-v4-flash",
                "kit.glm-5.3-flash"},
        primary="kit.deepseek-v4-flash", cfg_model="kit.glm-5.3-flash")
    assert "kit.deepseek-v4-flash" not in chain
    assert chain == ("kit.mistral-small-4-119b-a8b", "kit.glm-5.3-flash")


def test_confirmation_empty_when_only_one_family_served():
    # Only the primary's family can be served: it decides alone (today's
    # behaviour) -- an empty confirmation chain. A single served family is a
    # legitimate one-vote decision; an invented second model would only 404.
    chain = MC.classifier_confirmation(
        served={"kit.mistral-small-4-119b-a8b"},
        primary="kit.mistral-small-4-119b-a8b", cfg_model="kit.deepseek-v4-flash")
    assert chain == ()


def test_confirmation_filters_foreign_and_nonchat():
    chain = MC.classifier_confirmation(
        served={"kit.mistral-small-4-119b-a8b", "azure.gpt-5"},
        primary="kit.mistral-small-4-119b-a8b", cfg_model="kit.deepseek-v4-flash")
    # azure.gpt-5 is never a confirm.
    assert "azure.gpt-5" not in chain


# --------------------------------------------------------------------------
# The real second-opinion gate with fake clients (one slow, one fast)
# --------------------------------------------------------------------------

ALLOW = "DECISION: allow\nREASON: part of the task."
ASK = "DECISION: ask\nREASON: the request does not mention this."
REFUSE = "DECISION: refuse\nREASON: nobody asked for this."
TASK = [{"role": "user", "content": "run the failing test in my worktree"}]

PRIMARY = "kit.mistral-small-4-119b-a8b"
DEEPSEEK = "kit.deepseek-v4-flash"
GLM = "kit.glm-5.3-flash"


class Dialog:
    def __init__(self):
        self.shown = []

    def ask(self, name, args, preview):
        self.shown.append((name, preview))
        return True


class _Runner:
    """Fake classifier runner: each model answers from a table; a model can
    be ``slow`` (hangs past the wait) or ``broken`` (raises). Records which
    models were asked."""

    def __init__(self):
        self.calls = []
        self.slow = set()        # models that hang past the wait
        self.broken = set()      # models that raise
        self.answers = {}        # model -> answer text
        self.delay = 0.0

    def __call__(self, prompt, model):
        self.calls.append(model)
        if model in self.broken:
            raise RuntimeError("model unavailable")
        if model in self.slow:
            # Long enough that a 0.1-0.3 s wait always times out first, and
            # long enough to still be alive when a later model is asked (the
            # per-model box binding is what keeps its late write out).
            time.sleep(5.0)
        if self.delay:
            time.sleep(self.delay)
        return self.answers.get(model, ALLOW)


def _bounded_gate(tmp_path, monkeypatch, runner, models, dialog, *,
                    timeout=0.12, message="phase 1 done, please review",
                    perms=None, stalled_models=None):
    """Drive the second opinion for a message to another agent (an action
    that leaves the workspace; a missing verdict would ask the human).

    ``models[0]`` is the primary, ``models[1:]`` the ALLOW-confirmation chain.
    Returns ``(out, perms)`` so a test can inspect the per-turn latch. Pass a
    prior ``perms`` to keep the turn's latch across two actions."""
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
            perms.classifier_confirm_models = tuple(models[1:])
        if stalled_models:
            # Pre-seed the per-turn latch AND a recent attempt marker, so the
            # confirmation-with-a-latched-confirm path is reached on the FIRST
            # action deterministically (otherwise the probe window notion of
            # "recent" would ignore the latch on a fresh action).
            # The gate starts a fresh per-turn state (latch cleared) when the
            # turn key changes, so the seeded latch carries this turn's key.
            users = A._classifier_user_texts(TASK)
            perms.classifier_state.update(
                turn=(len(users), users[-1][:200] if users else ""),
                allowed=set(), refused=set(), refusals=0, calls=[],
                stalled_models=set(stalled_models))
            perms.classifier_health["last_attempt_at"] = (
                A._classifier_clock())
        perms.confirm_callback = dialog.ask
    out = A._DocToolExecutor()._classifier_second_opinion(
        perms, "session_message", {"to": "ec9c6393", "message": message},
        perms.mode)
    return out, perms


def test_primary_allow_is_confirmed_by_a_second_family(tmp_path, monkeypatch):
    # Primary allows; the confirming family also allows -> the action runs,
    # nobody asked.
    runner = _Runner()
    runner.answers[PRIMARY] = ALLOW
    runner.answers[DEEPSEEK] = ALLOW
    dialog = Dialog()
    out, _ = _bounded_gate(
        tmp_path, monkeypatch, runner,
        [PRIMARY, DEEPSEEK], dialog)
    assert out is None           # allowed after two families
    assert dialog.shown == []
    assert runner.calls == [PRIMARY, DEEPSEEK]


def test_primary_allow_then_confirm_refuses_is_a_refusal(tmp_path, monkeypatch):
    # "nothing refused today is allowed": an allow that a second family
    # refuses becomes a refusal -- never allowed on one vote.
    runner = _Runner()
    runner.answers[PRIMARY] = ALLOW
    runner.answers[DEEPSEEK] = REFUSE
    dialog = Dialog()
    out, _ = _bounded_gate(
        tmp_path, monkeypatch, runner,
        [PRIMARY, DEEPSEEK], dialog)
    assert out is not None             # refused
    assert "refused it" in str(out)    # the classifier refused, nothing allowed
    assert dialog.shown == []          # never silently allowed; held, not asked
    assert runner.calls == [PRIMARY, DEEPSEEK]


def test_primary_allow_confirm_times_out_asks(tmp_path, monkeypatch):
    # The confirming family times out: the allow is NOT granted; it asks and
    # the human is shown the reason (never allowed on one family's vote).
    runner = _Runner()
    runner.answers[PRIMARY] = ALLOW
    runner.slow.add(DEEPSEEK)
    dialog = Dialog()
    out, _ = _bounded_gate(
        tmp_path, monkeypatch, runner,
        [PRIMARY, DEEPSEEK], dialog)
    assert len(dialog.shown) == 1      # the ask reached the human
    assert runner.calls == [PRIMARY, DEEPSEEK]


def test_primary_refuse_does_not_call_the_confirm(tmp_path, monkeypatch):
    # The primary's own refuse stands at once; the confirming model is never
    # asked.
    runner = _Runner()
    runner.answers[PRIMARY] = REFUSE
    runner.answers[DEEPSEEK] = ALLOW
    dialog = Dialog()
    out, _ = _bounded_gate(
        tmp_path, monkeypatch, runner,
        [PRIMARY, DEEPSEEK], dialog)
    assert out is not None             # refused by the primary
    assert runner.calls == [PRIMARY]   # the confirm was not called


def test_confirm_tries_the_next_family_and_uses_its_verdict(
        tmp_path, monkeypatch):
    # primary allows; the first confirm (deepseek) hangs past the wait; the
    # next confirm (glm) refuses. The verdict is glm's refuse -- the timed-out
    # thread for deepseek must not get the chance to override it with its own
    # late allow. This is exactly the per-model box binding (_work binds its
    # own box and model) under test.
    runner = _Runner()
    runner.answers[PRIMARY] = ALLOW
    runner.slow.add(DEEPSEEK)          # hangs, then would write allow late
    runner.answers[GLM] = REFUSE
    dialog = Dialog()
    out, _ = _bounded_gate(
        tmp_path, monkeypatch, runner,
        [PRIMARY, DEEPSEEK, GLM], dialog)
    # glm's refuse wins, not deepseek's hypothetical late allow.
    assert out is not None
    assert "not run" in out
    assert runner.calls == [PRIMARY, DEEPSEEK, GLM]


def test_per_model_stall_latch_lets_a_later_action_use_another_model(
        tmp_path, monkeypatch):
    # After the primary times out once in a turn, the latch is keyed per
    # classifier model: the slow PRIMARY is latched but the (fast, other
    # family) confirm is NOT, so a later action is still judged rather than
    # blocked turn-wide. The estimator of "one slow model does not latch all"
    # is that only the slow model lands in stalled_models.
    runner = _Runner()
    runner.slow.add(PRIMARY)
    runner.answers[PRIMARY] = ALLOW       # would allow, but hangs first
    runner.answers[DEEPSEEK] = ALLOW
    dialog = Dialog()
    out, perms = _bounded_gate(
        tmp_path, monkeypatch, runner,
        [PRIMARY, DEEPSEEK], dialog)
    # The primary hung -> the turn asks (the confirm is not reached for an
    # absent primary verdict) -- the strict, safe outcome; nobody was silently
    # allowed.
    latch = perms.classifier_state["stalled_models"]
    assert PRIMARY in latch               # the slow primary is latched
    assert DEEPSEEK not in latch          # the other family is not


def test_primary_allow_with_the_only_confirm_latched_asks(tmp_path, monkeypatch):
    # pt-7 security rule: a second family was served (_confirm non-empty from
    # the start) but its only confirm is now latched (stalled_models). The
    # primary's allow must NOT stand on one vote -- it asks.
    runner = _Runner()
    runner.answers[PRIMARY] = ALLOW       # primary allows, fast
    runner.slow.add(DEEPSEEK)             # the only confirm, latches on ask
    dialog = Dialog()
    # Action 1: primary allows; confirm (deepseek) is slow -> times out and is
    # latched -> "not confirmed" -> ask (the dialog is shown; Dialog.ask True
    # means the user approves, so out is None -- the evidence of the ask is
    # that the dialog was SHOWN, never that the allow silently stood).
    # One action: the only confirm (deepseek) is pre-latched (stalled_models
    # already holds it), so _confirms is empty but the ORIGINAL _confirm was
    # non-empty (a second family is served). The primary's allow must NOT
    # stand on one vote -- it asks.
    out, perms = _bounded_gate(
        tmp_path, monkeypatch, runner,
        [PRIMARY, DEEPSEEK], dialog, stalled_models={DEEPSEEK})
    assert len(dialog.shown) == 1         # the lone allow never stood; it asked
    assert runner.calls == [PRIMARY]      # deepseek was never re-asked


def test_process_kill_command_is_asked_without_the_model():
    # A deterministic pkill/killall screen: the classifier must refuse/ask on
    # a process kill without touching the model (a killing the model would
    # otherwise let through, bench pkill -f python). The runner would ALLOW,
    # but the model must never be called.
    called = []
    runner = lambda p, m: called.append(m) or ALLOW
    verdict = _ac.judge(
        {"tool": "bash", "command": "pkill -f python"},
        {"task": "run the failing test"}, model="m", runner=runner,
        default_on_error="ask")
    assert verdict.decision == "ask"
    assert called == []                   # no model round trip for a kill


def test_attach_classifier_runner_stows_primary_and_confirm_from_the_roster(
        tmp_path, monkeypatch):
    # pt-6: the client stows classifier_model (primary) and
    # classifier_confirm_models (the allow-confirm chain) from ONE served
    # roster, and end-to-end a primary-allow + confirm-refuse becomes a refuse.
    served = {"kit.mistral-small-4-119b-a8b", "kit.deepseek-v4-flash",
              "kit.glm-5.3-flash", "azure.gpt-5", "kit.embed-v3"}
    monkeypatch.setattr(MC, "kit_models_served", lambda *a, **k: served)
    monkeypatch.setattr(A, "_classifier_settings",
                        lambda: {"model": PRIMARY, "explicit_model": False,
                                 "rules": [], "broad": False})

    class _Stub:
        _base_url = "kit://proxy"
        _api_key = "k"
        model = "kit.deepseek-v4-flash"   # the session's own model
        _classifier_complete = lambda self, prompt, model: ""  # noqa: E731

    ws = tmp_path / "ws"
    ws.mkdir()
    perms = A.KitToolPermissions(workspace=ws.resolve(), mode="acceptEdits")
    A.OpenAIClient._attach_classifier_runner(_Stub(), perms)
    assert perms.classifier_model == PRIMARY          # mistral, not the session
    assert perms.classifier_confirm_models == (DEEPSEEK, GLM)

    # End-to-end: the stowed models are used. Primary allows, confirm refuses.
    runner = _Runner()
    runner.answers[PRIMARY] = ALLOW
    runner.answers[DEEPSEEK] = REFUSE
    perms.classifier_runner = runner
    perms.classifier_messages = TASK
    dialog = Dialog()
    perms.confirm_callback = dialog.ask
    out = A._DocToolExecutor()._classifier_second_opinion(
        perms, "session_message", {"to": "ec9c6393", "message": "review"},
        perms.mode)
    assert out is not None               # refused, not allowed on one vote
    assert dialog.shown == []
    assert runner.calls == [PRIMARY, DEEPSEEK]
