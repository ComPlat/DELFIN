"""The reactivity department may spawn the experts its own text names.

Input:  the SHIPPED preset files in delfin/agent/pack/agents/, read by the
        registry the product builds at start, driven through the real
        executor path (``KitToolExecutor.execute`` on the subagent tool
        with a real ``KitToolPermissions``).
Output: whether a depth-1 delegation to a named expert passes the
        whitelist guard, whether one to an unnamed preset is refused, and
        whether the configuration agrees with the text.
Semantics: ``may_delegate_to`` is the whitelist ``_execute_subagent``
        enforces for depth > 1. Measured 2026-11-21: the department's
        description said it "routes a reactivity question to the right
        expert -- the ts-finder for a new transition state or barrier,
        the irc-checker to verify a TS the user already has", its prompt
        gave the two routing rules, and its own result contract asks for
        "which preset you delegated to" -- while ``may_delegate_to`` was
        ``[]``, so it could spawn nobody. Documentation promising a route
        the configuration forbids; the specialist framework was built and
        its department left empty.

Nothing here runs a model: the guard must fire before a runner is needed,
which is what the missing-runner error in the pass case proves.

Run: ~/micromamba/envs/delfin-ci/bin/python -m pytest \
     tests/test_the_department_can_reach_the_experts_it_names.py -q
"""

from __future__ import annotations

import dataclasses
import json
import re

import pytest

from delfin.agent import api_client as A
from delfin.agent import subagents as SA

DEPARTMENT = "reactivity-department"
NAMED_EXPERTS = ("ts-finder", "irc-checker")


def _perms(tmp_path, depth: int, preset_name: str):
    perms = A.KitToolPermissions(workspace=str(tmp_path), mode="default")
    perms = dataclasses.replace(
        perms,
        subagent_depth=depth,
        session_allowed_tools=frozenset({"subagent", "read_file"}),
    )
    # What _derive_perms carries on the child: the spawner's preset name,
    # which the guard reads with getattr.
    perms._preset_name = preset_name
    return perms


def _call(perms, sa_type: str) -> dict:
    executor = A._DocToolExecutor()
    return json.loads(executor.execute(
        "subagent",
        {"subagent_type": sa_type, "description": "d",
         "prompt": "please do the nested task, at least this long."},
        perms))


def _hook_applied() -> bool:
    """The call-side guard is a protected hook the operator applies."""
    import inspect
    return "may_delegate_to" in inspect.getsource(
        A._DocToolExecutor._execute_subagent)


@pytest.fixture(autouse=True)
def department_on(monkeypatch):
    """The department's switch must be on, or its presets are not in the
    registry at all (enabled_setting: agent.departments.reactivity)."""
    if DEPARTMENT not in SA.SUBAGENT_PRESETS:
        pytest.skip("the reactivity department is switched off in this "
                    "installation's settings")


def test_the_whitelist_names_exactly_the_experts_the_text_names():
    """The configuration against the prose, which is where it broke: the
    description and the routing rules name two experts, so the whitelist
    has to carry those two -- not more (a department may not reach an
    expert it never describes) and not fewer (an empty list made the
    documented route impossible)."""
    preset = SA.SUBAGENT_PRESETS[DEPARTMENT]
    assert tuple(preset.may_delegate_to) == NAMED_EXPERTS, (
        preset.may_delegate_to)

    text = f"{preset.description}\n{preset.system_prompt}"
    for expert in NAMED_EXPERTS:
        assert expert in text, (
            f"{expert} is in the whitelist but the department's own text "
            "does not say when to route to it")
    # And the other direction: no preset is named as a route in the text
    # without being in the whitelist.
    named = {n for n in SA.SUBAGENT_PRESETS
             if re.search(rf"\b{re.escape(n)}\b", text)} - {DEPARTMENT}
    assert named <= set(NAMED_EXPERTS), sorted(named)


def test_the_experts_themselves_delegate_to_nobody():
    """The depth-2 rule: a department reaches its experts, an expert
    reaches nobody. Without this the whitelist would be a chain."""
    for expert in NAMED_EXPERTS:
        assert SA.SUBAGENT_PRESETS[expert].may_delegate_to == (), expert


@pytest.mark.skipif(not _hook_applied(),
                    reason="SA0 whitelist hook not applied yet "
                           "(.gate/px copies await the operator build)")
@pytest.mark.parametrize("expert", NAMED_EXPERTS)
def test_the_department_may_spawn_a_named_expert(tmp_path, expert):
    """Driven through the real executor: the guard lets it through, and
    the error that remains is the missing runner -- which is the proof
    that the refusal did not come from the whitelist."""
    data = _call(_perms(tmp_path, 1, DEPARTMENT), expert)
    assert "whitelist" not in json.dumps(data).lower(), data
    assert "error" in data, data


@pytest.mark.skipif(not _hook_applied(),
                    reason="SA0 whitelist hook not applied yet "
                           "(.gate/px copies await the operator build)")
def test_the_department_may_not_spawn_an_expert_it_does_not_name(tmp_path):
    """The other direction, so the fix is not "let everything through":
    a preset outside the whitelist is still refused, by name."""
    outsider = "chemistry-reviewer"
    assert outsider in SA.SUBAGENT_PRESETS, "the control preset is gone"
    data = _call(_perms(tmp_path, 1, DEPARTMENT), outsider)
    blob = json.dumps(data).lower()
    assert "may_delegate_to" in blob and "whitelist" in blob, data
