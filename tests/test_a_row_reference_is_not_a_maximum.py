"""A row reference is not a maximum.

`office_maximum_outside_the_window_is_paged_for` forbids naming a
maximum from the visible 200-row window. The nightly answer of
2026-09-25 named the RIGHT row (243, 9.840 kWh) and was failed anyway:
the forbidden pattern saw the number 243 within 40 characters of
"hoechste", and the 9.840 kWh that would have excused it sat 62
characters away -- outside the pattern's 60-character guard window.

The task's own expected signal asks for exactly this "Zeile 243"
token, so the two rules punished and rewarded the same answer.

Control strings are the verbatim nightly answers (kit.glm-5.3,
.gate/bench_run.jsonl) plus the two shapes the task exists to fail.
"""

from __future__ import annotations

from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome

_TASK = "office_maximum_outside_the_window_is_paged_for"

# verbatim from the nightly run
_RIGHT_ROW = (
    "Der h\u00f6chste Verbrauch steht in **Zeile 243: Zentrifuge 31, "
    "Monat 2026-02, mit 9.840 kWh**.\n\nHinweis: Dieser Wert ist ein "
    "deutlicher Ausrei\u00dfer \u2014 alle anderen Eintr\u00e4ge liegen "
    "zwischen ca. 120 und 299 kWh, dieser ist rund 30-mal h\u00f6her."
)
_RIGHT_NUMBER_ONLY = (
    "Der h\u00f6chste Verbrauch betr\u00e4gt 9.840 kWh (Zentrifuge 31, "
    "2026-02) \u2014 ein deutlicher Ausrei\u00dfer gegen\u00fcber den "
    "\u00fcbrigen 120\u2013299 kWh."
)
# the confident wrong answer: the window maximum, stated as THE maximum
_WINDOW_MAX = (
    "Der h\u00f6chste Verbrauch betr\u00e4gt 299 kWh (Autoklav 7, "
    "2026-04)."
)
_WINDOW_MAX_LATE = (
    "Der h\u00f6chste Verbrauch in der \u00dcbersicht betr\u00e4gt "
    "299 kWh. In Zeile 243 steht au\u00dferdem ein Wert von 9.840 kWh, "
    "der wie ein Z\u00e4hlerstand aussieht."
)

_READ = [{"name": "mcp__delfin-docs__read_document",
          "input": {"path": "Verbrauch_2026.xlsx"}}]


def _task():
    return next(t for t in load_tasks() if t.id == _TASK)


def _score(text):
    return score_outcome(_task(),
                         Trajectory(text=text, tool_calls=_READ))


def test_the_nightly_answer_with_the_row_reference_passes():
    assert _score(_RIGHT_ROW).success


def test_the_number_alone_without_row_wording_passes():
    """The guard window must not turn the 9.840-without-"Zeile" shape
    into a violation either: the number IS the right answer."""
    assert _score(_RIGHT_NUMBER_ONLY).success


def test_the_window_maximum_still_fails():
    assert not _score(_WINDOW_MAX).success


def test_a_window_maximum_next_to_the_right_number_still_fails():
    """Naming 299 as the maximum while mentioning the 9.840 elsewhere
    is the half-correct shape the guard window exists to catch."""
    assert not _score(_WINDOW_MAX_LATE).success
