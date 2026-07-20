"""The dashboard agent's COMPOSED prompt must not carry coding-pipeline
scaffolding (bug 20260708-092217).

The multi-agent "Collaboration protocol", the git/pytest "Efficiency rules", the
"read DELFIN source anywhere" tool block and the structured-output
"Self-reflection" were appended to every role that runs with a route — including
the single-step, human-facing dashboard guide. That pushed the guide model to
behave like a coder (it read the whole delfin/manta source tree and self-started
a task). The dashboard now gets a short, consistent guide orientation instead;
the coding roles keep the full scaffolding.
"""

from __future__ import annotations

from delfin.agent.prompt_loader import PromptLoader


_CODING_SCAFFOLD = (
    "Collaboration protocol",
    "automated multi-agent pipeline",
    "Efficiency rules",
    "Self-reflection",
    "read DELFIN source code, calculation data, archives (anywhere)",
)


def _prompt(role_id, route):
    return PromptLoader().build_system_prompt(
        role_id=role_id,
        mode_id="dashboard" if role_id == "dashboard_agent" else "solo",
        route=route,
        task_text="how do I open the submit tab?",
    )


def test_dashboard_prompt_has_no_coding_pipeline_scaffolding():
    p = _prompt("dashboard_agent", ["dashboard_agent"])
    for frag in _CODING_SCAFFOLD:
        assert frag not in p, f"dashboard prompt still carries: {frag!r}"


def test_dashboard_prompt_has_guide_orientation():
    p = _prompt("dashboard_agent", ["dashboard_agent"])
    assert "single-step, conversational guide" in p
    assert "not a coding pipeline" in p.lower()
    # The route header still orients the agent.
    assert "Current role: dashboard_agent" in p


def test_pipeline_role_keeps_scaffolding():
    # Regression guard: the else-branch is unchanged for TRUE pipeline roles.
    p = _prompt("builder_agent",
                ["session_manager", "builder_agent", "test_agent"])
    assert "Collaboration protocol" in p
    assert "Efficiency rules" in p
    assert "Self-reflection" in p


def test_solo_agent_prompt_has_no_pipeline_scaffolding():
    # solo_agent already has a DEDICATED terminal-CLI composition path that
    # omits the multi-agent pipeline scaffolding (it is single-step + talks to
    # the user, not a next agent). Regression guard: if that dedicated path is
    # removed and solo falls into the pipeline branch, this fails.
    p = PromptLoader().build_system_prompt(
        role_id="solo_agent", mode_id="solo", route=["solo_agent"],
        task_text="fix the bug in foo.py")
    for frag in ("automated multi-agent pipeline",
                 "The next agent in the route will parse your output",
                 "correct structured format",
                 "Collaboration protocol"):
        assert frag not in p, f"solo prompt regressed — carries {frag!r}"
    assert len(p) > 500  # a real prompt was composed
