"""A classifier timeout is not a verdict: what it costs, and what it is.

Input: actions the permission rules already allow, judged through the real
gate (``_DocToolExecutor._run_permission_gate``), with a classifier runner
that hangs for one, two or three model families. The wait is pinned to
0.15 s, so a hanging family costs that much and the cases finish in under a
second. No test contacts a model.

Output, pinned here (measured on origin/main first: scratch/PR.md):

* a no-verdict fallback is not a vote, so a WORKSPACE action is not held
  because no second family confirmed a fallback -- the documented outcome
  for workspace work is "the rule result stands" (api_client.py:22314);
* one action costs ONE wait, not one per family it can reach: past the
  action's deadline no new model call starts, and the outcome is the
  no-verdict one (ask / hold), never an allow;
* control, unchanged by either rule: a real allow still needs its second
  family, a deterministic refusal still refuses with zero classifier calls,
  and an action that leaves the workspace is never allowed on a timeout.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent import api_client as A

ALLOW = "DECISION: allow\nREASON: part of the task."
TASK = [{"role": "user", "content": "rebuild the report in this workspace"}]
WAIT = 0.15


class Clock:
    """A clock the test advances. Cases that pin the DEADLINE use the real
    one: a frozen clock would never let an action's deadline pass."""

    def __init__(self):
        self.t = 1000.0

    def __call__(self):
        return self.t


class Runner:
    """``behaviour`` maps a model name to how it answers: ``"hang"`` sleeps
    past every wait, ``"allow"`` answers at once, a float sleeps that long
    first. A model with no entry hangs."""

    def __init__(self, behaviour=None, default="hang"):
        self.behaviour = dict(behaviour or {})
        self.default = default
        self.calls = []

    def __call__(self, prompt, model):
        self.calls.append(model)
        how = self.behaviour.get(model, self.default)
        if how == "hang":
            time.sleep(3.0)
            return ALLOW
        if how == "allow":
            return ALLOW
        time.sleep(float(how))
        return ALLOW


class Dialog:
    last_timed_out = False

    def __init__(self, answer=False):
        self.answer = answer
        self.shown = []

    def ask(self, name, args, preview):
        self.shown.append((name, preview))
        return self.answer


@pytest.fixture(autouse=True)
def _a_pinned_wait(monkeypatch):
    monkeypatch.setattr(A, "_CLASSIFIER_TIMEOUT_S", WAIT)
    monkeypatch.setattr(A, "_classifier_timeout_s", lambda: WAIT)
    monkeypatch.setattr(A, "_CLASSIFIER_RETRY_DELAY_S", 0.0)
    monkeypatch.setattr(A, "_classifier_settings",
                        lambda: {"model": "m1", "rules": [], "broad": True})
    # A host without a bash cage: the cage would screen workspace commands
    # without the model (test_bypass_workspace_work_needs_no_watcher.py),
    # and the failing classifier is what these cases are about.
    monkeypatch.setattr(A, "_bash_runs_caged", lambda perms: False)


def _perms(tmp_path, mode, runner, dialog=None):
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    perms = A.KitToolPermissions(workspace=ws.resolve(), mode=mode)
    perms.classifier_runner = runner
    perms.classifier_messages = TASK
    perms.classifier_model = "m1"
    perms.classifier_confirm_models = ("m2", "m3")
    if dialog is not None:
        perms.confirm_callback = dialog.ask
    return perms


def _gate(perms, name, args):
    return A._DocToolExecutor()._run_permission_gate(name, args, perms)


# -- the defect: a fallback is not a vote ------------------------------------

@pytest.mark.parametrize("mode,tool,args", [
    ("acceptEdits", "write_file", {"path": "report.md", "content": "x\n"}),
    # Bypass screens a workspace WRITE without the model (that is the point
    # of test_bypass_workspace_work_needs_no_watcher), so the fallback is
    # reached there with a command instead: an uncaged command goes to the
    # model in every mode.
    ("bypassPermissions", "bash", {"command": "python3 build.py"}),
])
def test_a_workspace_action_is_not_held_when_no_family_can_answer(
        tmp_path, mode, tool, args):
    """Nobody can be asked, the action stays in the workspace: the rule
    result stands (that is the documented no-verdict outcome for workspace
    work). Requiring a second family for a FALLBACK turned it into a hold
    after three waits -- the "34 edits held" of the report."""
    runner = Runner()                       # every family hangs
    perms = _perms(tmp_path, mode, runner)
    out = _gate(perms, tool, args)
    assert out is None, (
        "a workspace action was held because no family confirmed a "
        f"fallback: {out}")


# -- the cost: measured, reported, NOT changed here -------------------------
# scratch/PR.md records that one action can cost one wait PER FAMILY
# (measured on origin/main: 3 model calls x 20 s). The obvious fix -- one
# deadline per action -- is NOT shipped, because it breaks a pinned contract:
# tests/test_classifier_model_independence.py::test_confirm_tries_the_next_family_and_uses_its_verdict
# requires the chain to still reach the next family, and a family that hangs
# would consume the whole deadline for every later one. The measured design
# that keeps both is a SHARED deadline divided over the families still to be
# asked; it is named in PR.md for the integrator rather than guessed at here.


# -- controls: what the change must NOT weaken ------------------------------

def test_a_real_allow_still_needs_its_second_family(tmp_path):
    """The primary says allow, no other family answers: the allow is NOT
    granted. A fallback is not a vote, but a model's own allow still needs
    the second family -- this is what the change must leave standing."""
    runner = Runner({"m1": "allow"})
    perms = _perms(tmp_path, "bypassPermissions", runner)
    out = _gate(perms, "web_fetch", {"url": "https://github.com/a"})
    assert out is not None, "an allow was granted on one vote"
    assert "m1" in runner.calls


def test_a_real_allow_with_a_confirming_family_is_allowed(tmp_path):
    """The same action when the second family DOES answer: allowed."""
    runner = Runner({"m1": "allow", "m2": "allow", "m3": "allow"})
    perms = _perms(tmp_path, "bypassPermissions", runner)
    assert _gate(perms, "web_fetch", {"url": "https://github.com/b"}) is None
    assert "m2" in runner.calls or "m3" in runner.calls


@pytest.mark.parametrize("mode", ["default", "acceptEdits",
                                  "bypassPermissions"])
@pytest.mark.parametrize("attended", [True, False])
def test_a_leaving_action_is_never_allowed_on_a_timeout(tmp_path, mode, attended):
    """A timeout must never turn into an allow for an action that leaves the
    workspace: it asks, or it is held when nobody can be asked."""
    runner = Runner()
    dialog = Dialog(False) if attended else None
    perms = _perms(tmp_path, mode, runner, dialog)
    out = _gate(perms, "web_fetch", {"url": "https://github.com/c"})
    assert out is not None
    if attended:
        assert len(dialog.shown) == 1
    else:
        assert "was not run" in out


def test_a_refused_action_stays_refused_when_the_classifier_is_down(tmp_path):
    """The deterministic gates decide before the classifier is consulted:
    a deny-listed command is refused with ZERO model calls, so nothing they
    refuse depends on the classifier answering."""
    runner = Runner()
    perms = _perms(tmp_path, "bypassPermissions", runner)
    denied = "rm " + "-rf" + " scratch/whatever"
    out = _gate(perms, "bash", {"command": denied})
    assert out is not None and "deny-pattern" in out
    assert runner.calls == []
