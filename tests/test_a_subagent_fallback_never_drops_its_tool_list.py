"""A sub-agent's permissions are never wider than assigned, on any path.

Input: ``subagents._derive_perms`` with a parent permissions object whose
type cannot take the fields the clone wants to set (a stand-in that
predates them), and with a narrowing step that raises. Output: the
child's permissions, or ``SubagentPermissionsError``.

Semantics: the role's tool list is a restriction. A fallback that cannot
carry it does not return a child without it: the derivation refuses
(``SubagentPermissionsError``) and ``run_subagent`` reports the refusal
as the result's error instead of running a child. A narrowing step that
raises leaves the child able to call nothing, not everything. Returning
the parent's own permissions object (no copy, the parent's mode) is
never a fallback either.
"""
from __future__ import annotations

import dataclasses
import types

import pytest

from delfin.agent import subagents as SA


@dataclasses.dataclass
class _OldPerms:
    """A permissions type without session_allowed_tools / subagent_depth."""
    mode: str = "default"
    workspace: str = ""


def test_a_type_that_cannot_carry_the_tool_list_is_refused(tmp_path):
    parent = _OldPerms(workspace=str(tmp_path))
    with pytest.raises(SA.SubagentPermissionsError) as exc:
        SA._derive_perms(parent, "plan", preset_tools=("read_file", "grep_file"))
    assert "tool" in str(exc.value).lower()


def test_a_narrowing_that_raises_leaves_the_child_with_nothing(tmp_path, monkeypatch):
    from delfin.agent import api_client as A
    parent = A.KitToolPermissions(workspace=str(tmp_path), mode="default")

    def boom(*_a, **_k):
        raise RuntimeError("narrowing failed")

    monkeypatch.setattr(SA, "_narrow_allowed_tools", boom)
    child = SA._derive_perms(parent, "default", preset_tools=("read_file",))
    assert child.session_allowed_tools == frozenset({SA._PRESET_ALLOWS_NOTHING})


def test_the_parent_object_is_never_handed_out(tmp_path):
    parent = types.SimpleNamespace(mode="default", workspace=str(tmp_path),
                                   subagent_depth=0)
    with pytest.raises(SA.SubagentPermissionsError):
        SA._derive_perms(parent, "plan")


def test_a_type_without_the_fields_still_works_when_nothing_is_narrowed(tmp_path):
    parent = _OldPerms(workspace=str(tmp_path))
    child = SA._derive_perms(parent, "plan")
    assert child is not parent and child.mode == "plan"


def test_run_subagent_reports_the_refusal_instead_of_running(tmp_path):
    parent = _OldPerms(workspace=str(tmp_path))
    ran = []

    class _Client:
        def __getattr__(self, name):
            ran.append(name)
            raise AttributeError(name)

    res = SA.run_subagent(subagent_type="explore", description="d",
                          prompt="look at the code", parent_client=_Client(),
                          parent_perms=parent)
    assert res.error and "permission" in res.error.lower()
    assert res.final_text == ""
