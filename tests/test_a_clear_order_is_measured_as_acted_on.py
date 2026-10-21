"""The act behaviour: a clear order ends in an edit that landed.

The counterpart of ``asked`` (tests/test_benchmark_behavior.py). From session
78c36f8c (2026-10-19): after "ja mach" the agent restated its plan, and the
edits it did try were refused while the answer said they were made -- so a
refused edit must not count as acting.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from delfin.agent import artifact_check as AC
from delfin.agent import benchmark as bm
from delfin.agent.benchmark import Task, Trajectory

_EDIT = {"name": "edit_file", "input": '{"file_path": "game.html"}'}
_READ = {"name": "read_file", "input": '{"file_path": "game.html"}'}
_ASK = {"name": "ask_user_question", "input": '{"question": "Breite?"}'}
_REFUSED = '{"error": "\'game.html\' is used by background job 8678a5da"}'


def _flag(calls, results):
    task = Task(id="t", task_class="behavior_act", mode="solo", prompt="",
                behavior="act")
    return bm.behavior_flags(task, Trajectory(tool_calls=calls,
                                              tool_results=results))


def test_an_edit_that_landed_is_acting():
    assert _flag([_READ, _EDIT], ["...", "Edited game.html"]) == {"acted": 1}


def test_a_refused_edit_is_not_acting():
    assert _flag([_READ, _EDIT], ["...", _REFUSED]) == {"acted": 0}


def test_a_plan_without_an_edit_is_not_acting():
    assert _flag([_READ], ["..."]) == {"acted": 0}


def test_asking_before_the_edit_is_not_acting():
    assert _flag([_ASK, _EDIT], ['{"answers": ["10"]}', "Edited game.html"]) \
        == {"acted": 0}


def test_a_refusal_then_a_landed_edit_is_acting():
    assert _flag([_EDIT, _EDIT], [_REFUSED, "Edited game.html"]) \
        == {"acted": 1}


def test_the_tasks_load_and_the_fixture_is_a_sound_page():
    root = Path(__file__).resolve().parents[1]
    data = yaml.safe_load((root / "delfin/agent/pack/benchmark/"
                           "tasks_auto_behavior.yaml").read_text())
    acts = [t for t in data["tasks"] if t.get("behavior") == "act"]
    assert len(acts) == 2
    page = root / "tests/fixtures/behavior_workspace/game.html"
    assert "const BX = 4" in page.read_text()
    assert AC.problems([page]) == ""
