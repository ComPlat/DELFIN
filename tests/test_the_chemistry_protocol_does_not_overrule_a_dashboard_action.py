"""A dashboard action is not a chemistry question.

Benchmark run of 2026-09-25 (kit.glm-5.3): "stell den orca functional auf
b3lyp ein" was classified as chemistry (keywords orca, functional), the
injected protocol said "Use search_docs BEFORE answering", and the
dashboard agent searched twice and set nothing -- against its own rule 1,
"dashboard action first". The protocol now says, in the dashboard, that
it is for questions and that an action comes first.
"""

from __future__ import annotations

from delfin.agent.prompt_loader import PromptLoader


def test_the_dashboard_protocol_puts_the_action_first():
    prompt = PromptLoader().build_system_prompt(
        role_id="dashboard_agent", mode_id="dashboard",
        task_text="stell den orca functional auf b3lyp ein")
    start = prompt.index("--- Chemistry Protocol ---")
    block = prompt[start:start + 600]
    assert "ACTION line first" in block
    assert block.index("ACTION line first") < block.index("search_docs")


def test_the_solo_protocol_is_unchanged():
    prompt = PromptLoader().build_system_prompt(
        role_id="solo_agent", mode_id="solo",
        task_text="which orca functional should I use for a Fe complex")
    if "--- Chemistry Protocol ---" in prompt:
        start = prompt.index("--- Chemistry Protocol ---")
        assert "ACTION line first" not in prompt[start:start + 600]
