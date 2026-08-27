"""A plan already awaiting approval must not be waited for again.

The approval callback blocks for the whole approval window -- 180 s in the
dashboard. When it expires, ``exit_plan_mode`` answers "still awaiting
approval, do NOT resubmit". That sentence is all that stood between the
agent and a second full wait, and a sentence is not a mechanism.

It failed twice, on the record, in the same code:

  2026-06-25  two exit_plan_mode calls x 10 min = a 21-minute hang. The
              fix was to stop calling a timeout a "rejection" and to ask
              the model to stop. The window was also cut to 3 minutes.
  2026-08-27  same session, same shape: exit_plan_mode 180 297 ms, then
              exit_plan_mode 180 450 ms. Six minutes of a user's session
              spent waiting a second time on a click that had already not
              come once. Four bug reports in ninety minutes, one asking
              plainly whether the agent had hung.

Two calls agreeing to 150 ms is not a person deciding twice; it is a
timeout, twice. So the second submission is now answered from recorded
state and never reaches the callback at all.
"""

from __future__ import annotations

import json

from delfin.agent.api_client import KitToolPermissions, _doc_executor


class _Approver:
    """Stands in for the dashboard's blocking wait, and counts the waits."""

    def __init__(self, *, approves=False, timed_out=True):
        self.calls = 0
        self.approves = approves
        self.timed_out = timed_out

    def approve_plan(self, plan: str) -> dict:
        self.calls += 1
        if self.approves:
            return {"approved": True, "new_mode": "default"}
        return {"approved": False, "timed_out": self.timed_out}


def _perms(tmp_path, approver=None):
    p = KitToolPermissions(mode="plan", workspace=str(tmp_path))
    if approver is not None:
        p.plan_approval_callback = approver.approve_plan
    return p


def _submit(perms, plan="# Plan\n\n1. Look\n2. Report\n"):
    return json.loads(_doc_executor.execute(
        "exit_plan_mode", {"plan": plan}, perms))


# ---------------------------------------------------------------------------
# The mechanism
# ---------------------------------------------------------------------------

def test_the_second_submission_never_reaches_the_approval_wait(tmp_path):
    """The whole finding in one assertion: the callback is called ONCE."""
    approver = _Approver()
    perms = _perms(tmp_path, approver)

    first = _submit(perms)
    assert first["status"] == "awaiting_approval"
    assert approver.calls == 1

    second = _submit(perms)
    assert second["status"] == "awaiting_approval"
    assert second.get("resubmitted") is True
    assert approver.calls == 1, "the second call must not wait again"


def test_a_third_and_fourth_submission_do_not_wait_either(tmp_path):
    """The recorded session made two. Nothing caps it at two."""
    approver = _Approver()
    perms = _perms(tmp_path, approver)
    for _ in range(4):
        _submit(perms)
    assert approver.calls == 1


def test_a_different_plan_also_does_not_start_a_second_wait(tmp_path):
    """Rewriting the plan is the obvious way around a guard keyed on the
    text. The guard is keyed on "something is pending", not on the words."""
    approver = _Approver()
    perms = _perms(tmp_path, approver)
    _submit(perms, "# Plan A\n\n1. one\n")
    out = _submit(perms, "# Plan B — completely different\n\n1. two\n")
    assert approver.calls == 1
    assert out.get("resubmitted") is True


def test_the_reply_says_the_earlier_plan_is_the_one_on_screen(tmp_path):
    """A silent no-op would leave the model believing its new plan is what
    the user is looking at."""
    approver = _Approver()
    perms = _perms(tmp_path, approver)
    _submit(perms, "# Plan A\n\n1. one\n")
    out = _submit(perms, "# Plan B\n\n1. two\n")
    assert "did not" in out["message"]
    assert "still the one" in out["message"]


def test_an_identical_resubmission_is_named_as_identical(tmp_path):
    approver = _Approver()
    perms = _perms(tmp_path, approver)
    plan = "# Plan\n\n1. one\n"
    _submit(perms, plan)
    out = _submit(perms, plan)
    assert "identical" in out["message"]


# ---------------------------------------------------------------------------
# What must still get through
# ---------------------------------------------------------------------------

def test_the_first_submission_is_untouched(tmp_path):
    """The guard may not cost the user the plan they are meant to see."""
    approver = _Approver(approves=True)
    perms = _perms(tmp_path, approver)
    out = _submit(perms)
    assert out["status"] == "approved"
    assert approver.calls == 1
    assert perms.mode == "default"


def test_approval_clears_the_block(tmp_path):
    """Otherwise a session gets one plan and never another."""
    approver = _Approver(approves=True)
    perms = _perms(tmp_path, approver)
    _submit(perms)
    assert perms.plan_awaiting_approval == ""


def test_an_explicit_rejection_clears_the_block(tmp_path):
    """A rejection is an answer. The next plan is a fresh proposal and has
    to reach the user -- a guard that survives a rejection would silence
    every plan after the first one the user turned down."""
    approver = _Approver(approves=False, timed_out=False)
    perms = _perms(tmp_path, approver)
    out = _submit(perms)
    assert out["status"] == "rejected"
    assert perms.plan_awaiting_approval == ""

    approver.approves = True
    again = _submit(perms)
    assert again["status"] == "approved"
    assert approver.calls == 2, "a fresh proposal must reach the wait"


def test_leaving_plan_mode_clears_the_block(tmp_path):
    """The user approved in the UI, which flips the mode. Whatever was
    pending is settled, and a later return to plan mode must start clean."""
    approver = _Approver()
    perms = _perms(tmp_path, approver)
    _submit(perms)
    assert perms.plan_awaiting_approval != ""

    perms.mode = "acceptEdits"
    _submit(perms)                      # refused: wrong mode
    assert perms.plan_awaiting_approval == ""

    perms.mode = "plan"
    approver.approves = True
    out = _submit(perms)
    assert out["status"] == "approved"


def test_a_session_with_no_approval_channel_is_unchanged(tmp_path):
    """Headless has no callback at all: it must still refuse to
    self-approve, and it must not acquire a block it can never clear."""
    perms = _perms(tmp_path, None)
    out = _submit(perms)
    assert out["status"] == "awaiting_approval"
    assert perms.mode == "plan", "plan mode is never self-exitable"


# ---------------------------------------------------------------------------
# A delegate is not its parent
# ---------------------------------------------------------------------------

def test_a_delegate_does_not_inherit_the_parents_pending_plan(tmp_path):
    """dataclasses.replace copies field values, so without this the child's
    own first submission would be answered from a plan the user never saw
    -- the same defect the read_tracker copy next to it exists to prevent."""
    from delfin.agent.subagents import _derive_perms

    approver = _Approver()
    parent = _perms(tmp_path, approver)
    _submit(parent)
    assert parent.plan_awaiting_approval != ""

    child = _derive_perms(parent, "plan")
    assert child is not None
    assert getattr(child, "plan_awaiting_approval", "") == ""


# ---------------------------------------------------------------------------
# The instruction still matches what the mechanism does
# ---------------------------------------------------------------------------

def test_the_timeout_reply_no_longer_promises_only_a_rule(tmp_path):
    """It used to say "do NOT resubmit" and stop there. It now says what
    happens if you do, which is the part that is true whatever the model
    decides."""
    approver = _Approver()
    perms = _perms(tmp_path, approver)
    out = _submit(perms)
    assert "returns immediately" in out["message"]


def test_the_dashboard_wait_is_still_bounded():
    """The 180 s window is what makes a second wait expensive enough to be
    worth preventing. If it ever became unbounded, this guard would be
    hiding a freeze rather than preventing a repeat."""
    import inspect

    from delfin.dashboard import tab_agent

    src = inspect.getsource(tab_agent)
    assert "ev.wait(timeout=180)" in src
