"""A subagent must never be more permissive than its parent.

Closes a plan-mode escape: a parent pinned to ``plan`` (the human-approval
gate) could call ``subagent(subagent_type="general-purpose")`` whose preset
carries ``mode="default"`` and obtain a child that writes/bashes on its behalf.
The child mode is now clamped to the more restrictive of parent/preset.
"""

from __future__ import annotations

from delfin.agent.subagents import _clamp_child_mode, _derive_perms
from delfin.agent.api_client import KitToolPermissions


def test_clamp_plan_parent_never_widens():
    assert _clamp_child_mode("plan", "default") == "plan"
    assert _clamp_child_mode("plan", "acceptEdits") == "plan"
    assert _clamp_child_mode("plan", "bypassPermissions") == "plan"


def test_clamp_preserves_a_stricter_preset():
    # A read-only explore preset (plan) stays plan even from a default parent.
    assert _clamp_child_mode("default", "plan") == "plan"
    assert _clamp_child_mode("acceptEdits", "plan") == "plan"


def test_clamp_same_mode_is_identity():
    assert _clamp_child_mode("default", "default") == "default"
    assert _clamp_child_mode("plan", "plan") == "plan"


def test_clamp_child_never_more_permissive_than_parent():
    # preset default is not widened past a default parent
    assert _clamp_child_mode("bypassPermissions", "default") == "default"
    # a bypass preset is clamped down to a default parent
    assert _clamp_child_mode("default", "bypassPermissions") == "default"


def test_derive_perms_plan_parent_clamps_child(tmp_path):
    parent = KitToolPermissions(workspace=tmp_path, mode="plan")
    child = _derive_perms(parent, "default")
    assert child.mode == "plan", "plan parent must not spawn a writable child"


def test_derive_perms_default_parent_keeps_preset(tmp_path):
    parent = KitToolPermissions(workspace=tmp_path, mode="default")
    assert _derive_perms(parent, "default").mode == "default"
    # explore preset (plan) preserved from a default parent
    assert _derive_perms(parent, "plan").mode == "plan"


def test_derive_perms_inherits_deny_and_sandbox(tmp_path):
    # The clamp only touches mode; deny-list/workspace still inherit.
    parent = KitToolPermissions(workspace=tmp_path, mode="plan")
    child = _derive_perms(parent, "default")
    assert child.workspace == parent.workspace
    assert child.bash_deny_patterns == parent.bash_deny_patterns
