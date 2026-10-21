"""A sub-agent's exit_plan_mode never reaches the user's plan dialog.

Input: a parent KitToolPermissions with a plan_approval_callback (the
dashboard's blocking "Plan awaiting approval" wait), a child derived from
it with ``_derive_perms`` (read-only presets run in mode "plan"), and an
``exit_plan_mode`` call made with the child's permissions.

Output: the child's call returns at once, without calling the callback;
the child stays in plan mode and is told to hand its plan back as its
report. The parent's own call still reaches the callback.

Semantics: the plan dialog belongs to the session the user is talking to.
``dataclasses.replace`` copied the parent's callback into every child, so
a read-only explore child that called exit_plan_mode put its plan on the
parent's screen, blocked for the dashboard's 180 s window, and the plan
stayed "pending" in the saved session and was announced again at every
restore. Measured in one all_free session: two such dialogs, four restore
notices, two children whose 300 s wall budget ran out waiting.

Rejected alternative: clearing the callback in ``_derive_perms`` only.
The executor is the one place every child's call passes through, and a
callback attached to a child later (re-wiring after an engine rebuild)
would reopen the route; the depth check holds regardless of wiring.
"""

from __future__ import annotations

import json

from delfin.agent.api_client import KitToolPermissions, _doc_executor
from delfin.agent.subagents import _derive_perms


class _Approver:
    def __init__(self):
        self.calls = 0

    def approve_plan(self, plan: str) -> dict:
        self.calls += 1
        return {"approved": True, "new_mode": "default"}


def _parent(tmp_path, mode="bypassPermissions"):
    approver = _Approver()
    p = KitToolPermissions(mode=mode, workspace=str(tmp_path))
    p.plan_approval_callback = approver.approve_plan
    return p, approver


def _submit(perms, plan="# Plan\n\n1. Read\n2. Report\n"):
    return json.loads(_doc_executor.execute(
        "exit_plan_mode", {"plan": plan}, perms))


def test_a_read_only_child_does_not_reach_the_parents_dialog(tmp_path):
    parent, approver = _parent(tmp_path)
    child = _derive_perms(parent, "plan")
    assert child.mode == "plan"
    out = _submit(child)
    assert approver.calls == 0
    assert child.mode == "plan"
    assert child.plan_awaiting_approval == ""
    assert "error" not in out


def test_the_child_is_told_to_return_its_plan_as_its_report(tmp_path):
    parent, _ = _parent(tmp_path)
    child = _derive_perms(parent, "plan")
    out = _submit(child)
    text = json.dumps(out).lower()
    assert "sub-agent" in text
    assert "report" in text


def test_a_child_of_a_plan_mode_parent_cannot_be_approved_out_of_it(tmp_path):
    parent, approver = _parent(tmp_path, mode="plan")
    child = _derive_perms(parent, "default")
    assert child.mode == "plan"
    _submit(child)
    assert approver.calls == 0
    assert child.mode == "plan"


def test_the_parents_own_plan_still_reaches_the_dialog(tmp_path):
    parent, approver = _parent(tmp_path, mode="plan")
    out = _submit(parent)
    assert approver.calls == 1
    assert parent.mode == "default"
    assert "error" not in out
