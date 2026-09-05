"""Tests for plan-mode wiring: addendum injection + dropdown surfacing.

Plan-mode is solo_agent role + permission_profile="plan" (read-only)
+ a markdown addendum that tells the model to use exit_plan_mode for
approval. The addendum file lives at
``delfin/agent/pack/shared/plan_mode_addendum.md`` and the prompt
loader picks it up when ``mode_id == "plan"`` OR the active
``permission_mode == "plan"`` — plan is a permission profile now, so
setting Perms = Plan (in the Code mode) gets the full plan experience.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.prompt_loader import PromptLoader


def test_plan_mode_addendum_file_exists():
    p = Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack" / "shared" / "plan_mode_addendum.md"
    assert p.is_file(), f"missing addendum: {p}"
    body = p.read_text(encoding="utf-8")
    assert "Plan Mode" in body
    assert "exit_plan_mode" in body


def test_solo_prompt_includes_plan_addendum_when_mode_is_plan():
    loader = PromptLoader()
    prompt = loader.build_system_prompt(
        role_id="solo_agent",
        mode_id="plan",
        mode_description="plan mode test",
        task_text="figure out how to add feature X",
    )
    assert "Plan Mode" in prompt
    assert "exit_plan_mode" in prompt


def test_solo_prompt_includes_plan_addendum_when_permission_is_plan():
    # Plan is a permission profile now: Code mode + Perms=plan must still get
    # the full plan addendum (read-only-first → exit_plan_mode).
    loader = PromptLoader()
    prompt = loader.build_system_prompt(
        role_id="solo_agent",
        mode_id="solo",            # Code mode, NOT the legacy "plan" mode
        permission_mode="plan",    # but the active permission profile is plan
        task_text="figure out how to add feature X",
    )
    assert "Plan Mode" in prompt
    assert "exit_plan_mode" in prompt


def test_solo_prompt_skips_plan_addendum_in_other_modes():
    loader = PromptLoader()
    for mode in ("solo", "dashboard", "quick"):
        prompt = loader.build_system_prompt(
            role_id="solo_agent",
            mode_id=mode,
            mode_description=f"{mode} mode test",
            task_text="anything",
        )
        # No plan mode_id and no plan permission → addendum must not appear.
        assert "# Plan Mode\n" not in prompt, mode


def test_plan_mode_addendum_documents_exit_plan_mode_handoff():
    p = Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack" / "shared" / "plan_mode_addendum.md"
    body = p.read_text(encoding="utf-8")
    # Key contract elements must be present
    assert "exit_plan_mode" in body
    assert "approve" in body.lower()
    assert "acceptEdits" in body


# ---------------------------------------------------------------------------
# Plan approval decision (button vs /plan approve|reject must not drift)
# ---------------------------------------------------------------------------

from delfin.dashboard.tab_agent import _finalize_plan_decision  # noqa: E402


def _fresh_state(plan_body="step 1\nstep 2"):
    return {
        "_plan_approval_result": {"approved": False, "new_mode": "default"},
        "_pending_plan_body": plan_body,
    }


def _fake_save(plan_body, repo_root):
    # Deterministic path without touching disk.
    return f"{repo_root}/.claude/plans/saved.md"


def test_approve_flips_mode_and_persists():
    st = _fresh_state()
    msgs = _finalize_plan_decision(
        st, approved=True, repo_root="/repo", save_plan_fn=_fake_save,
    )
    res = st["_plan_approval_result"]
    assert res["approved"] is True
    assert res["new_mode"] == "acceptEdits"
    assert res["plan_path"] == "/repo/.claude/plans/saved.md"
    assert "_pending_plan_body" not in st          # body was consumed
    assert any("Plan saved" in m for m in msgs)


def test_reject_resets_mode_and_does_not_persist():
    saved = []
    st = _fresh_state()
    msgs = _finalize_plan_decision(
        st, approved=False, repo_root="/repo",
        save_plan_fn=lambda *a, **k: saved.append(1),
    )
    res = st["_plan_approval_result"]
    assert res["approved"] is False
    assert res["new_mode"] == "default"
    assert "plan_path" not in res
    assert saved == []                             # save_plan never called
    assert msgs == []


def test_button_and_slash_paths_are_identical():
    # Both call sites pass the same args to the helper — prove the helper
    # is deterministic so the button == /plan approve.
    st_btn = _fresh_state()
    st_cmd = _fresh_state()
    m_btn = _finalize_plan_decision(
        st_btn, approved=True, repo_root="/repo", save_plan_fn=_fake_save,
    )
    m_cmd = _finalize_plan_decision(
        st_cmd, approved=True, repo_root="/repo", save_plan_fn=_fake_save,
    )
    assert st_btn["_plan_approval_result"] == st_cmd["_plan_approval_result"]
    assert m_btn == m_cmd


def test_no_pending_result_is_noop():
    assert _finalize_plan_decision({}, approved=True, repo_root="/r") == []


def test_save_failure_surfaces_message_but_keeps_approval():
    def _boom(plan_body, repo_root):
        raise RuntimeError("disk full")

    st = _fresh_state()
    msgs = _finalize_plan_decision(
        st, approved=True, repo_root="/repo", save_plan_fn=_boom,
    )
    res = st["_plan_approval_result"]
    assert res["approved"] is True                 # approval still stands
    assert res["new_mode"] == "acceptEdits"
    assert "plan_path" not in res
    assert any("Plan save failed" in m for m in msgs)


def test_empty_plan_body_skips_save():
    saved = []
    st = _fresh_state(plan_body="   ")
    msgs = _finalize_plan_decision(
        st, approved=True, repo_root="/repo",
        save_plan_fn=lambda *a, **k: saved.append(1),
    )
    assert saved == []                             # nothing to persist
    assert st["_plan_approval_result"]["approved"] is True
    assert msgs == []


# ---------------------------------------------------------------------------
# Plan -> durable state bridge (approval persists the plan + scaffolds tasks)
# ---------------------------------------------------------------------------

import json


def _approving_perms(tmp_path):
    from delfin.agent.api_client import KitToolPermissions
    perms = KitToolPermissions(workspace=tmp_path, mode="plan")
    perms.plan_approval_callback = lambda plan: {
        "approved": True, "new_mode": "acceptEdits"}
    perms.task_session_id = "plan-bridge-test"
    return perms


_PLAN = """# Build the spectrum exporter

## Steps
1. Create exporter module with CSV writer
2. Wire exporter into the dashboard download button
3. Add regression tests for the exporter
"""


def test_approved_plan_is_persisted_and_scaffolds_tasks(tmp_path, monkeypatch):
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    from delfin.agent import api_client as A
    ws = tmp_path / "ws"; ws.mkdir()
    perms = _approving_perms(ws)
    out = json.loads(A._doc_executor._execute_exit_plan_mode(
        {"plan": _PLAN}, perms))
    assert out["status"] == "approved"
    assert out["tasks_created"] == 3
    assert "ALREADY in your task list" in out["instruction"]
    # Plan persisted to the workspace plans store.
    from delfin.agent.memory_store import list_plans
    plans = list_plans(ws)
    assert plans and "spectrum exporter" in str(plans).lower()
    # Tasks scaffolded from the numbered steps.
    from delfin.agent.agent_tasks import get_store
    tasks = get_store(ws).list(session_id="plan-bridge-test")
    subjects = [t["subject"] for t in tasks]
    assert any("exporter module" in s for s in subjects)
    assert len(tasks) == 3


def test_plan_without_steps_creates_no_tasks(tmp_path, monkeypatch):
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    from delfin.agent import api_client as A
    ws = tmp_path / "ws"; ws.mkdir()
    perms = _approving_perms(ws)
    out = json.loads(A._doc_executor._execute_exit_plan_mode(
        {"plan": "Just prose, no numbered steps in this plan at all."},
        perms))
    assert out["status"] == "approved"
    assert out["tasks_created"] == 0
    from delfin.agent.agent_tasks import get_store
    assert get_store(ws).list(session_id="plan-bridge-test") == []


def test_plan_steps_parser_handles_checkboxes_and_cap():
    from delfin.agent.api_client import _DocToolExecutor
    plan = "\n".join([f"- [ ] step number {i} does a thing"
                      for i in range(20)])
    steps = _DocToolExecutor._plan_steps(plan)
    assert len(steps) == 12                      # capped
    assert steps[0].startswith("step number 0")
