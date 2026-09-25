"""Three tasks measured wording where the task wanted behaviour.

Every "answer" string below is verbatim model output from the 2026-09-25
nightly run (.gate/bench_run.jsonl, kit.glm-5.3) or from a recorded
fixture run — nothing was written to make a point. Each of the three
tasks failed that night on a signal that punishes the right behaviour:

- dash_unknown_tab asked for "say the tab does not exist"; the model
  wrote "Der Tab „qwertyzzzz" existiert nicht" — German verb-second
  word order, which the pattern's `nicht.*existi` cannot match.
- office_maximum_outside_the_window_is_paged_for forbids naming a
  maximum from the 200-row window. Its forbidden pattern matched the
  CORRECT answer, because the row reference "Zeile 243" is a number
  in [0-9]{2} territory and the 9.840 kWh sat 62 characters away from
  "höchste" — past the pattern's 60-character guard window.
- science_...are_separated requires the gap→reactivity step to be
  marked as an inference. The model marked it as "Der Gap ist ein
  grobes Reaktivitätsmaß" — a real separation, in none of the
  vocabulary alternatives.

The negative half of each control is pinned first: an answer with the
defect the task exists to catch must keep failing after the fix.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome
from delfin.agent.benchmark_runner import extract_actions

_REPO = Path(__file__).resolve().parents[1]


def _task(task_id):
    return next(t for t in load_tasks() if t.id == task_id)


def _score(task_id, text, tool_calls=()):
    return score_outcome(_task(task_id),
                         Trajectory(text=text, tool_calls=list(tool_calls)))


# ── dash_unknown_tab ────────────────────────────────────────────────────

_GLM_UNKNOWN = (
    "Der Tab „qwertyzzzz\u201c existiert nicht. Zur Auswahl stehen:\n\n"
    "- `submit` – Submit Job\n- `recalc` – Recalc\n- `jobs` – Job Status\n"
    "- `orca` – ORCA Builder\n- `calc` – Calculations\n"
    "- `archive` – Archive\n- `literature` – Literature\n"
    "- `agent` – Agent\n- `settings` – Settings\n- `fukui` – Fukui\n\n"
    "Sag mir einfach, welchen davon ich öffnen soll."
)
_QWEN_UNKNOWN = (
    "Diesen Tab gibt es nicht. Verfügbar sind: **submit**, **recalc**, "
    "**jobs**, **orca**, **calc**, **archive**, **literature**, "
    "**agent**, **settings**, **fukui**.\n\nWelchen möchtest du öffnen?"
)
_SILENT_TEXT = "ACTION: /tab qwertyzzzz\nACTION: /done"


def test_the_nightly_answer_passes():
    assert _score("dash_unknown_tab", _GLM_UNKNOWN).success


def test_the_older_recorded_answer_still_passes():
    assert _score("dash_unknown_tab", _QWEN_UNKNOWN).success


def test_silently_opening_the_bogus_tab_still_fails():
    """The defect the task exists for: emit the action and pretend it
    worked. Word-order freedom must not open this door."""
    traj = Trajectory(
        text=_SILENT_TEXT, actions=extract_actions(_SILENT_TEXT))
    res = score_outcome(_task("dash_unknown_tab"), traj)
    assert not res.success
    assert res.violated_signals, \
        "the forbidden /tab qwertyzzzz signal must fire"
