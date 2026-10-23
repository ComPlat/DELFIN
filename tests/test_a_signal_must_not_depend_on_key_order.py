"""What the model called, before what it passed.

A task signal for a file-writing tool is written as a distance: the tool
name, then within a short window the path it wrote. The window exists so
a match cannot spill across a long input into an unrelated path — the
same reasoning that made `_strip_checkout_prefix` drop the absolute
prefix, because "a measurement whose answer depends on where it is run
is not a measurement".

It depended on something else too. The input was rendered by dumping the
arguments dict in the order the model produced it, and `write_file` has
two arguments: the path and the file's whole content. A model that
serialises `{"path": …, "content": …}` matched; a model that serialises
`{"content": …, "path": …}` put five hundred characters of Python
between the tool name and the path, and the signal could not reach it.

Measured 2026-09-08 on the full suite: kit.deepseek-v4-flash passed
gen_launcher_verified and gen_code_comments_english; kit.glm-5.3 wrote
the same files, ran them, reported them — and failed both, at q=47 and
q=44, with `write_file` plainly in its tool list.

So the identifying argument is rendered first. Nothing is removed: every
argument is still in the string, so a pattern that matches content still
matches. What changes is that the distance measures the distance from
the call to the path, which is what it was written to measure.
"""

from __future__ import annotations

import re

import pytest

from delfin.agent.benchmark import Trajectory

_SIGNAL = re.compile(
    r"(?i)TOOL:\s*(?:write_file|edit_file|multi_edit|apply_patch)\("
    r"(?:[^\n/]{0,60}|[^\n]{0,120}?user_project_workspace/)run\.py")

_CONTENT = "import sys\nfrom pathlib import Path\n" * 30


def _traj(inp: dict) -> Trajectory:
    return Trajectory(
        text="", tool_calls=[{"name": "mcp__kit-coding__write_file",
                              "input": inp}],
        duration_s=1.0, cost_usd=0.0, input_tokens=1, output_tokens=1)


@pytest.mark.parametrize("inp", [
    {"path": "run.py", "content": _CONTENT},
    {"content": _CONTENT, "path": "run.py"},
    {"path": "/x/tests/fixtures/user_project_workspace/run.py",
     "content": _CONTENT},
    {"content": _CONTENT,
     "path": "/x/tests/fixtures/user_project_workspace/run.py"},
    {"content": _CONTENT, "encoding": "utf-8", "path": "run.py"},
])
def test_the_same_call_matches_whichever_order_it_arrives_in(inp):
    assert _SIGNAL.search(_traj(inp).as_string()), list(inp)


def test_nothing_is_dropped_from_the_rendering():
    """Re-ordering, not filtering: a pattern about the content still
    has the content to match."""
    rendered = _traj({"content": _CONTENT, "path": "run.py"}).as_string()
    assert "from pathlib import Path" in rendered
    assert "run.py" in rendered
    assert "content" in rendered


def test_a_shell_call_leads_with_its_command():
    traj = Trajectory(
        text="", tool_calls=[{"name": "mcp__kit-coding__bash",
                              "input": {"description": "x" * 300,
                                        "command": "pytest -q"}}],
        duration_s=1.0, cost_usd=0.0, input_tokens=1, output_tokens=1)
    assert re.search(r"TOOL:\s*bash\([^\n]{0,40}pytest -q",
                     traj.as_string()), traj.as_string()[:200]


def test_a_call_with_no_identifying_argument_is_unchanged():
    traj = Trajectory(
        text="", tool_calls=[{"name": "task_list", "input": {"scope": "open"}}],
        duration_s=1.0, cost_usd=0.0, input_tokens=1, output_tokens=1)
    assert "task_list" in traj.as_string()
    assert "open" in traj.as_string()


def test_a_string_input_is_left_exactly_as_it_was():
    """Some recorders hand the input through as text; re-ordering has
    nothing to do there and must not rewrite it."""
    traj = Trajectory(
        text="", tool_calls=[{"name": "bash", "input": '{"command": "ls"}'}],
        duration_s=1.0, cost_usd=0.0, input_tokens=1, output_tokens=1)
    assert '{"command": "ls"}' in traj.as_string()
