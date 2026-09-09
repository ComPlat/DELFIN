"""The rubric for "name three real ORCA keywords" named two of them.

`fact_orca_casscf_keywords_auto` asks for at least three genuine `%casscf`
keywords, looked up in the manual before answering. Its three expected
signals were `nel|nroots`, `norb|ptmethod` and `mult|actorbs` — so the
third slot did not measure "a third real keyword", it measured whether
the model happened to pick one of two.

kit.glm-5.3 failed it on 2026-09-09 having named eight: nel, norb,
MaxIter, ETol, GTol, ActConstraints, SwitchConv and switchstep, with
section 3.13.4 cited, arrived at through search_docs and read_section,
and not one of them invented. Every one is in the manual. The answer was
right and the rubric could not see it.

The fix widens the third slot to the block's documented keywords, and
this file is the reason that widening is not a relaxation: every
alternative is checked against the ORCA manual in the doc index, at
keyword position, so the rubric cannot come to accept a keyword that does
not exist. That is the same discipline the office rubrics use when they
take their expected figure from the tool rather than from arithmetic
written by hand.

A second thing fell out of reading the section. `multiplicity` is
forbidden — it is the Molcas spelling, and claiming it as ORCA's is the
invention this task exists to catch — but the manual itself writes the
plain English word in a comment ("imult 0 # multiplicity block counting
from zero"). A bare word match would fail an answer for quoting the
source it was told to consult, so the pattern now matches the word where
a keyword would be written.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

_TASK = "fact_orca_casscf_keywords_auto"
_SECTION = "ch3_13_4_keywords"
_INDEX = Path.home() / ".delfin" / "doc_index.json"


def _task():
    from delfin.agent.benchmark import load_tasks

    return next(t for t in load_tasks() if t.id == _TASK)


def _manual_section() -> str:
    if not _INDEX.exists():
        pytest.skip("no doc index on this machine")
    data = json.loads(_INDEX.read_text(encoding="utf-8"))
    try:
        sec = data["documents"]["orca_manual_6_1_1_delfin"]["sections"][_SECTION]
    except (KeyError, TypeError):
        pytest.skip("the ORCA manual section is not in this index")
    return sec["text"] if isinstance(sec, dict) else str(sec)


def _alternatives(pattern: str) -> list[str]:
    """The bare words offered by a `\\b(?:a|b|c)\\b` alternation."""
    m = re.search(r"\(\?:([^)]+)\)", pattern)
    assert m, f"the signal is no longer an alternation: {pattern}"
    return [a for a in m.group(1).split("|") if a]


def test_every_keyword_the_rubric_accepts_is_in_the_manual():
    """The whole reason the third slot may be wide."""
    text = _manual_section()
    accepted = _alternatives(_task().expected_signals[2].pattern)
    assert len(accepted) >= 10, "the slot narrowed back to a handful"
    absent = [
        kw for kw in accepted
        if not re.search(rf"(?im)^\s{{0,3}}{re.escape(kw)}\b", text)
    ]
    assert not absent, (
        "the rubric would accept keywords the manual does not document: "
        + ", ".join(absent))


def test_the_answer_that_was_wrongly_failed_now_passes():
    from delfin.agent.benchmark import Trajectory, _signal_matches

    answer = (
        "Aus dem ORCA-Handbuch (6.1.1), Abschnitt 3.13.4 Keywords: "
        "`nel` – Anzahl der Elektronen im aktiven Raum, `norb` – Anzahl "
        "der aktiven Orbitale. Typische optionale: `MaxIter`, `ETol` / "
        "`GTol`, `ActConstraints`, `SwitchConv` / `switchstep`.")
    traj = Trajectory(text=answer, tool_calls=[
        {"name": "mcp__delfin-docs__search_docs", "input": {"query": "casscf"}},
        {"name": "mcp__delfin-docs__read_section",
         "input": {"section_id": _SECTION}},
    ])
    task = _task()
    unmatched = [i for i, sig in enumerate(task.expected_signals)
                 if not _signal_matches(sig, traj)]
    assert not unmatched, f"still unmatched: expected{unmatched}"
    fired = [i for i, sig in enumerate(task.forbidden_signals)
             if _signal_matches(sig, traj)]
    assert not fired, f"a correct answer trips forbidden{fired}"


def test_naming_only_the_two_mandatory_ones_is_still_not_three():
    """The slot has to keep asking for a third keyword, or the task stops
    measuring the thing its prompt asks for."""
    from delfin.agent.benchmark import Trajectory, _signal_matches

    traj = Trajectory(
        text="Im %casscf-Block braucht man nel und norb, mehr nicht.")
    assert not _signal_matches(_task().expected_signals[2], traj)


def test_an_invented_keyword_is_still_caught():
    from delfin.agent.benchmark import Trajectory, _signal_matches

    task = _task()
    for invented in ("nactel 6", "nactorb 6", "multiplicity 3",
                     "multiplicity = 3", "multiplicity: 3"):
        traj = Trajectory(text=f"Setze {invented} im %casscf-Block.")
        assert any(_signal_matches(s, traj) for s in task.forbidden_signals), (
            f"invented keyword walks through: {invented}")


def test_a_pattern_cannot_anchor_on_emphasis_that_is_stripped_first():
    """Why the guard above matches a value and not code formatting.
    _strip_emphasis runs before every text signal, so a pattern written
    around backticks is text that can never fire — the same silent no-op
    as a set literal holding one tool vocabulary."""
    from delfin.agent.benchmark import _strip_emphasis

    assert _strip_emphasis("Setze `multiplicity` im Block.") == (
        "Setze multiplicity im Block.")
    assert "`" not in _task().forbidden_signals[-1].pattern


def test_quoting_the_manuals_own_comment_is_not_an_invented_keyword():
    """The manual writes the English word in a comment, and this task
    tells the model to go and read the manual."""
    from delfin.agent.benchmark import Trajectory, _signal_matches

    task = _task()
    traj = Trajectory(text=(
        "Die Zeile `imult 0` im Handbuch ist kommentiert mit "
        "\"multiplicity block counting from zero\"; das ORCA-Keyword "
        "heißt mult."))
    fired = [i for i, s in enumerate(task.forbidden_signals)
             if _signal_matches(s, traj)]
    assert not fired, f"quoting the source trips forbidden{fired}"
