"""Solo agent orchestration CONTRACT — keep it working like Claude Code.

The solo_agent role prompt already codifies the Claude-Code working method
(understand → plan → read → implement → verify → report; live task list; probe
before assuming; verify after modifying; parallel tool calls; change tactic on
repeated failure; ask when genuinely uncertain). This test PINS that method as
a CI-enforced contract so a future prompt trim can't silently erode the
systematics — the point is that the agent orchestrates tasks the SAME way over
time, not that it drifts.

It intentionally asserts on the shipped prompt text (like the token-budget
guard), because the methodology lives in the prompt and that is the single
source of truth the model actually receives.
"""

from __future__ import annotations

from pathlib import Path

import pytest

_SOLO = (
    Path(__file__).resolve().parents[1]
    / "delfin" / "agent" / "pack" / "agents" / "solo_agent.md"
)


def _text() -> str:
    return _SOLO.read_text(encoding="utf-8")


# The canonical Claude-Code working loop, expressed as required phrases. Each
# maps to a systematic the agent must keep. Update this list CONSCIOUSLY when
# the method itself changes — never to silence a failure from a careless trim.
_WORKING_LOOP = [
    "Understand first",          # read the request before acting
    "Plan before acting",        # brief approach for non-trivial work
    "Verify your work",          # don't declare done without checking
    "Report minimally",          # lead with the outcome, no fluff
]

_ORCHESTRATION = [
    # Live task list, exactly one in-progress, complete immediately.
    "task_update(task_id, status='in_progress')",
    "Never have ≥2 tasks",
    # The eight-pattern playbook pillars that define HOW work is approached.
    "Pre-probe over assume",
    "Verify-after-modify",
    "Stop-trigger awareness",
    "Parallel independent tool calls",
    # Autonomy with a genuine ask-when-uncertain brake.
    "Autonomy ≠ guessing",
    # Orient at session start the way a careful engineer does.
    "git status",
    "git log --oneline",
    # Always verify a code change by running it.
    "After every code edit",
]


@pytest.mark.parametrize("phrase", _WORKING_LOOP + _ORCHESTRATION)
def test_solo_prompt_keeps_orchestration_systematic(phrase):
    assert phrase in _text(), (
        f"solo_agent orchestration regressed — missing {phrase!r}. "
        "The agent must keep approaching tasks the Claude-Code way; restore it "
        "or update this contract consciously."
    )


def test_solo_prompt_exists_and_is_substantial():
    # A truncated prompt would pass individual phrase checks against a cache;
    # assert the real file is the comprehensive playbook.
    body = _text()
    assert len(body) > 5000
    assert body.lstrip().startswith("# Solo Agent")
