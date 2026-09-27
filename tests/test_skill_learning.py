"""Tests for delfin/agent/skill_learning.py (package 2).

The message forms below mirror real DELFIN transcripts saved by
session_store.save_session: roles are user/system/thinking/tool/assistant;
tool messages embed an HTML ``<span class="tool-name">`` chip plus a
``<details><summary> -> {json}`` result block. The samples here are
anonymised (no paths, accounts or hostnames) but keep those shapes.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin.agent.skill_learning import qualifies  # noqa: E402


def _tool(name: str, param: str, result: str = '{"status": "ok"}') -> dict:
    return {"role": "tool", "content":
            f'<span class="tool-name">{name}</span>  '
            f'<span class="tool-param">{param}</span>'
            f'<details><summary> &rarr; {result}</details>'}


def _bash(cmd: str, exit_code: int = 0, out: str = "") -> dict:
    res = f'{{"exit_code": {exit_code}, "stdout": "{out}"}}'
    return {"role": "tool", "content":
            f'<span class="tool-name">$</span> {cmd}'
            f'<details><summary> &rarr; {res}</details>'}


TRIVIAL = [
    {"role": "user", "content": "What does X mean?"},
    {"role": "assistant", "content": "X does Y."},
]

# Five tool calls in one working stretch, then a green gate run.
DEEP_WORK = [
    {"role": "user", "content": "Fix the failing test in module X."},
    {"role": "thinking", "content": "Grep first."},
] + [
    _bash("grep -n KEYWORD delfin/x.py"),
    _tool("Read", "delfin/x.py", "file body"),
    _tool("Edit", "delfin/x.py old new", '{"status": "ok"}'),
    _bash("ruff check delfin/x.py"),
    {"role": "tool", "content":
     '<span class="tool-name">TestRunner</span>  <span class="tool-param">'
     'tests/test_x.py</span><details><summary> &rarr; '
     '{"summary": {"passed": 5, "failed": 0}}</details>'},
]

# A tool call failed, then the same step went green.
ERROR_RECOVERY = [
    {"role": "user", "content": "Run the tests for module X."},
    _bash("gate tests/test_x.py", exit_code=1, out="1 failed"),
    _bash("gate tests/test_x.py", exit_code=0, out="all passed"),
]

# The user corrected the agent's course.
USER_CORRECTION = [
    {"role": "user", "content": "Do it in file A."},
    {"role": "assistant", "content": "I edited B."},
    {"role": "user", "content": "No, that is wrong — I meant A, please revert and use A."},
    {"role": "assistant", "content": "Reverted, edited A."},
]

# Green gate run carried by a tool message — evidence, but no qualification.
NO_QUALIFY_WITH_EVIDENCE = [
    {"role": "user", "content": "Run tests/test_x.py once."},
    {"role": "tool", "content":
     '<span class="tool-name">TestRunner</span>  <span class="tool-param">'
     'tests/test_x.py</span><details><summary> &rarr; '
     '{"summary": {"passed": 3, "failed": 0}}</details>'},
]


class TestQualifies:
    def test_trivial_session_qualifies_for_nothing(self):
        assert qualifies(TRIVIAL) == []

    def test_five_tool_calls_in_one_stretch_qualify(self):
        reasons = qualifies(DEEP_WORK)
        assert any("tool" in r.lower() for r in reasons), reasons
        assert len(reasons) >= 1

    def test_error_then_recovery_qualifies(self):
        reasons = qualifies(ERROR_RECOVERY)
        assert any("recover" in r.lower() for r in reasons), reasons

    def test_user_correction_qualifies(self):
        reasons = qualifies(USER_CORRECTION)
        assert any("correct" in r.lower() for r in reasons), reasons

    def test_green_run_alone_is_not_a_reason(self):
        # Evidence without depth: a single green test run is not, by
        # itself, a reason to propose a skill.
        assert qualifies(NO_QUALIFY_WITH_EVIDENCE) == []

    def test_empty_and_none_inputs(self):
        assert qualifies([]) == []
        assert qualifies(None) == []

    def test_reasons_are_named_english_strings(self):
        for r in qualifies(DEEP_WORK):
            assert isinstance(r, str) and 4 < len(r) < 200, r
